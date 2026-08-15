/*********************************************************************************************

This file is part of the PSOPT library, a software tool for computational optimal control

Copyright (C) 2009-2026 Victor M. Becerra

This library is free software; you can redistribute it and/or
modify it under the terms of the GNU Lesser General Public
License as published by the Free Software Foundation; either
version 2.1 of the License, or (at your option) any later version.

This library is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public
License along with this library; if not, write to the Free Software
Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA,
or visit http://www.gnu.org/licenses/

Author:    Professor Victor M. Becerra
Address:   University of Portsmouth
           School of Electrical and Mechanical Engineering
           Portsmouth PO1 3DJ
           United Kingdom
e-mail:    v.m.becerra@ieee.org

**********************************************************************************************/

#include "psopt.h"
#include "psopt_sqp.hpp"

#ifdef USE_SQP

#include <qpOASES.hpp>
#include <vector>
#include <cmath>
#include <algorithm>

#include "psopt_qp_plugin.h"
#include <string>

// Implemented in qp_plugin_loader.cxx. A QP backend other than qpOASES lives in a
// plugin, loaded on demand, so that the linear algebra it carries cannot collide with
// another backend's; psopt_qp_plugin.h explains why that is not optional.
bool psopt_qp_plugin_solve(const std::string& backend,
                           const psopt_qp_problem* problem,
                           psopt_qp_solution*      solution,
                           std::string&            message);
bool psopt_qp_plugin_available(const std::string& backend, std::string& message);

using namespace std;

// The two qpOASES types that appear in the storage below. The rest of the qpOASES
// names are brought in inside the driver, where they belong.
using qpOASES::sparse_int_t;
using qpOASES::real_t;

namespace {

// ---------------------------------------------------------------------------------
// Constraint violation and the l1 merit function.
//
// PSOPT states every constraint two-sided, g_l <= g(x) <= g_u, with an equality
// written as coincident bounds. The violation of one constraint is therefore how far
// g sits outside its own interval, and the l1 merit function is the objective plus a
// weighted sum of those violations,
//
//     phi(x; r) = f(x) + sum_i r_i * viol_i(x).
//
// Nothing here assumes the equalities come first, so the caller's nlp_neq is not
// needed.
// ---------------------------------------------------------------------------------
inline double violation(double gi, double gl, double gu)
{
    if (gi < gl) return gl - gi;
    if (gi > gu) return gi - gu;
    return 0.0;
}

double merit(double fval, const MatrixXd& gval, const vector<double>& gl,
             const vector<double>& gu, const vector<double>& r)
{
    double s = 0.0;
    for (int i = 0; i < gval.rows(); i++) s += r[i]*violation(gval(i), gl[i], gu[i]);
    return fval + s;
}

double total_violation(const MatrixXd& gval, const vector<double>& gl,
                       const vector<double>& gu)
{
    double s = 0.0;
    for (int i = 0; i < gval.rows(); i++) s += violation(gval(i), gl[i], gu[i]);
    return s;
}

double max_violation(const MatrixXd& gval, const vector<double>& gl,
                     const vector<double>& gu)
{
    double s = 0.0;
    for (int i = 0; i < gval.rows(); i++) s = max(s, violation(gval(i), gl[i], gu[i]));
    return s;
}


// ---------------------------------------------------------------------------------
// Compressed sparse column storage, in qpOASES's index type.
//
// The subproblem's matrices keep their pattern for as long as the mesh does: the
// Jacobian's comes from the recorded tape, the Hessian's from the same tape's
// second-order sparsity, and neither depends on the point. So the pattern is built
// once and only the values are refreshed, which is also the form qpOASES wants --
// it holds pointers into these arrays and reads them again at every solve.
//
// The class also carries the few matrix operations the algorithm needs, so that no
// dense copy of the Jacobian is made anywhere: J*d for the second-order correction
// and J'*lambda for the gradient of the Lagrangian.
// ---------------------------------------------------------------------------------
class SparseCsc {
public:
    // Build the pattern from a triplet list. Repeated entries are summed on scatter,
    // which is what makes it safe to hand this the lower triangle of a symmetric
    // matrix twice over. force_diagonal adds every diagonal position whether the
    // triplets reach it or not, so that a multiple of the identity can be added later
    // without disturbing the pattern.
    void build(int nr, int nc,
               const vector<unsigned int>& row, const vector<unsigned int>& col,
               bool force_diagonal)
    {
        nr_ = nr; nc_ = nc;
        const size_t nt = row.size();

        // (column, row, triplet index), ordered by column and then by row, which is
        // both the compressed column order and the order qpOASES expects within a
        // column.
        vector<size_t> order(nt + (force_diagonal ? (size_t) min(nr,nc) : 0));
        vector<unsigned int> er(order.size()), ec(order.size());
        vector<int>          ek(order.size());

        for (size_t k = 0; k < nt; k++) { er[k] = row[k]; ec[k] = col[k]; ek[k] = (int) k; }
        if (force_diagonal)
            for (int j = 0; j < min(nr,nc); j++) {
                const size_t k = nt + (size_t) j;
                er[k] = (unsigned int) j; ec[k] = (unsigned int) j; ek[k] = -1;
            }

        for (size_t k = 0; k < order.size(); k++) order[k] = k;
        sort(order.begin(), order.end(), [&](size_t a, size_t b) {
            if (ec[a] != ec[b]) return ec[a] < ec[b];
            if (er[a] != er[b]) return er[a] < er[b];
            return a < b;
        });

        ir_.clear(); jc_.assign((size_t) nc_ + 1, 0); val_.clear();
        perm_.assign(nt, -1);
        diag_.assign((size_t) nc_, -1);

        int slot = -1;
        unsigned int lastr = 0, lastc = 0;
        for (size_t q = 0; q < order.size(); q++) {
            const size_t k = order[q];
            if (slot < 0 || er[k] != lastr || ec[k] != lastc) {
                ir_.push_back((sparse_int_t) er[k]);
                val_.push_back(0.0);
                slot  = (int) ir_.size() - 1;
                lastr = er[k]; lastc = ec[k];
                jc_[(size_t) ec[k] + 1]++;
                if (er[k] == ec[k]) diag_[(size_t) ec[k]] = slot;
            }
            if (ek[k] >= 0) perm_[(size_t) ek[k]] = slot;
        }
        for (int j = 0; j < nc_; j++) jc_[(size_t) j + 1] += jc_[(size_t) j];
    }

    // Values in the triplet order the pattern was built from.
    void scatter(const vector<double>& tval)
    {
        fill(val_.begin(), val_.end(), 0.0);
        for (size_t k = 0; k < perm_.size() && k < tval.size(); k++)
            if (perm_[k] >= 0) val_[(size_t) perm_[k]] += tval[k];
    }

    void scale(double s)
    {
        for (size_t k = 0; k < val_.size(); k++) val_[k] *= s;
    }

    void shift_diagonal(double tau)
    {
        for (int j = 0; j < nc_; j++) if (diag_[(size_t) j] >= 0)
            val_[(size_t) diag_[(size_t) j]] += tau;
    }

    void times(const double* v, double* y) const              // y = M v, length nr
    {
        for (int i = 0; i < nr_; i++) y[i] = 0.0;
        for (int j = 0; j < nc_; j++)
            for (sparse_int_t k = jc_[(size_t) j]; k < jc_[(size_t) j + 1]; k++)
                y[ir_[(size_t) k]] += val_[(size_t) k]*v[j];
    }

    void transpose_times(const double* v, double* y) const     // y = M' v, length nc
    {
        for (int j = 0; j < nc_; j++) {
            double s = 0.0;
            for (sparse_int_t k = jc_[(size_t) j]; k < jc_[(size_t) j + 1]; k++)
                s += val_[(size_t) k]*v[ir_[(size_t) k]];
            y[j] = s;
        }
    }

    // The largest absolute column sum, which for a symmetric matrix is its 1-norm and
    // bounds its spectral radius. Used only to scale the convexification shift.
    double norm1() const
    {
        double best = 0.0;
        for (int j = 0; j < nc_; j++) {
            double s = 0.0;
            for (sparse_int_t k = jc_[(size_t) j]; k < jc_[(size_t) j + 1]; k++)
                s += fabs(val_[(size_t) k]);
            best = max(best, s);
        }
        return best;
    }

    int    rows() const { return nr_; }
    int    cols() const { return nc_; }
    int    nnz()  const { return (int) val_.size(); }
    bool   empty() const { return val_.empty() && nc_ == 0; }
    double value(int slot) const { return val_[(size_t) slot]; }

    // The stored entries, as triplets, for callers that want the matrix in another
    // library's format.
    template <class Trip, class Idx>
    void emit_triplets(vector<Trip>& out, Idx /*index type tag*/) const
    {
        out.clear();
        out.reserve(val_.size());
        for (int j = 0; j < nc_; j++)
            for (sparse_int_t k = jc_[(size_t) j]; k < jc_[(size_t) j + 1]; k++)
                out.push_back(Trip((Idx) ir_[(size_t) k], (Idx) j, val_[(size_t) k]));
    }

    // The same pattern in the plugin interface's index type. Built once, on request:
    // qpOASES wants its own sparse_int_t and the plugins want a fixed 64-bit type, and
    // on a platform where those are distinct types a cast of the array is not a
    // conversion of it.
    const long long* ir64() const
    {
        if (ir64_.size() != ir_.size())
            ir64_.assign(ir_.begin(), ir_.end());
        return ir64_.empty() ? NULL : &ir64_[0];
    }
    const long long* jc64() const
    {
        if (jc64_.size() != jc_.size())
            jc64_.assign(jc_.begin(), jc_.end());
        return &jc64_[0];
    }
    const double* val_const() const { return val_.empty() ? NULL : &val_[0]; }

    sparse_int_t* ir()  { return ir_.empty() ? NULL : &ir_[0]; }
    sparse_int_t* jc()  { return &jc_[0]; }
    double*       val() { return val_.empty() ? NULL : &val_[0]; }
    const vector<double>& values() const { return val_; }

private:
    int nr_ = 0, nc_ = 0;
    vector<sparse_int_t> ir_, jc_;
    vector<double>       val_;
    mutable vector<long long> ir64_, jc64_;
    vector<int>          perm_;      // triplet index -> slot
    vector<int>          diag_;      // column -> slot of its diagonal, or -1
};

// ---------------------------------------------------------------------------------
// The quadratic programming subproblem, in the one form all three of the solver's
// subproblems take:
//
//   min  1/2 d' H d + g' d      s.t.  lbA <= A d <= ubA,   lbd <= d <= ubd
//
// and the two backends that solve it. Everything above this point is shared; the
// backends differ in what they can be given and what they do with it.
// ---------------------------------------------------------------------------------
struct QpProblem {
    int              nv = 0, mc = 0;
    const SparseCsc* H  = NULL;      // sparse Hessian, or
    const MatrixXd*  Bd = NULL;      // dense Hessian; exactly one of the two
    const double*    g  = NULL;
    const SparseCsc* A  = NULL;      // mc by nv, absent when mc == 0
    const double    *lbA = NULL, *ubA = NULL, *lbd = NULL, *ubd = NULL;
};

// The result, in PSOPT's conventions: grad f + A' lambda - z = 0.
struct QpSolution {
    vector<double> d, lambda, z;
    int            iterations = 0;
    bool           ok = false;
    bool           relaxed = false;   // the constraints could not be met and were relaxed
};


// Solve a subproblem through a plugin backend. The conversion is only a matter of
// pointers: the storage the SQP already keeps is exactly what the interface asks for.
static bool solve_qp_plugin(const string& backend, const QpProblem& p,
                            double tol, int max_iter, bool nonconvex,
                            QpSolution& out, string& message)
{
    out.d.assign((size_t) p.nv, 0.0);
    out.lambda.assign((size_t) max(p.mc,1), 0.0);
    out.z.assign((size_t) p.nv, 0.0);
    out.iterations = 0;
    out.ok = false;

    psopt_qp_problem q;
    q.abi_version = PSOPT_QP_ABI_VERSION;
    q.n = p.nv; q.m = p.mc;

    q.H_p = q.H_i = NULL; q.H_x = NULL; q.H_dense = NULL;
    if (p.H != NULL) {
        q.H_p = p.H->jc64(); q.H_i = p.H->ir64(); q.H_x = p.H->val_const();
    }
    else {
        q.H_dense = &(*p.Bd)(0,0);      // Eigen is column-major, as the interface says
    }

    q.g = p.g;
    q.A_p = q.A_i = NULL; q.A_x = NULL;
    if (p.mc > 0) {
        q.A_p = p.A->jc64(); q.A_i = p.A->ir64(); q.A_x = p.A->val_const();
    }
    q.lbA = p.lbA; q.ubA = p.ubA; q.lb = p.lbd; q.ub = p.ubd;
    q.tolerance = tol;
    q.max_iter  = max_iter;
    q.nonconvex = nonconvex ? 1 : 0;

    psopt_qp_solution r;
    r.d = &out.d[0]; r.lambda = &out.lambda[0]; r.z = &out.z[0];
    r.iterations = 0; r.status = PSOPT_QP_FAILED;

    if (!psopt_qp_plugin_solve(backend, &q, &r, message)) return false;

    out.iterations = r.iterations;
    out.ok = (r.status == PSOPT_QP_SOLVED) || (r.status == PSOPT_QP_APPROXIMATE);
    return true;
}


// d' H d, with H either the sparse exact Hessian or the dense quasi-Newton model.
double quadratic_form(bool exact, const SparseCsc& Hm, const MatrixXd& B, const MatrixXd& d)
{
    if (!exact) return (d.transpose()*B*d)(0,0);
    vector<double> Hd((size_t) d.rows(), 0.0);
    Hm.times(&d(0), &Hd[0]);
    double s = 0.0;
    for (int j = 0; j < d.rows(); j++) s += d(j)*Hd[(size_t) j];
    return s;
}

} // anonymous namespace


static void sqp_gradient(MatrixXd& x, MatrixXd& gf, Workspace* workspace)
{
    if (useAutomaticDifferentiation(*workspace->algorithm))
        ScalarGradientAD(ff_ad, x, &gf, &workspace->trace_f_done, workspace->ad_f, workspace);
    else
        ScalarGradient(ff_num, x, &gf, workspace->grw.get(), workspace);
}


// ---------------------------------------------------------------------------------
// The Lagrangian and its Hessian.
//
// L = f + lambda' g, in PSOPT's sign convention for the multipliers, so that the
// stationarity the solver tests, grad f + J' lambda - z = 0, is the gradient of this
// function. The Hessian is taken from the same AD machinery the IPOPT interface uses
// for its exact-Hessian option, and comes back as the lower triangle of a sparse
// symmetric matrix.
//
// The multipliers enter the tape as constants, so the tape has to be recorded again
// whenever they change -- once per iteration. That is the cost of an exact Hessian
// and it is the same cost IPOPT pays. The pattern is detected once, from a tape
// recorded with every multiplier set to one: a multiplier that happened to be zero
// would otherwise remove its constraint from the tape and, with it, that constraint's
// entries from a pattern that is then kept for the rest of the mesh.
// ---------------------------------------------------------------------------------
static adouble sqp_lagrangian(adouble* xad, const double* lam, int m, Workspace* workspace)
{
    adouble L = ff_ad(xad, workspace);
    if (m > 0) {
        adouble* g = workspace->gad.get();
        gg_ad(xad, g, workspace);
        for (int i = 0; i < m; i++) L += lam[i]*g[i];
    }
    return L;
}

static psopt_ad::SparseTriplet
sqp_hessian_triplet(MatrixXd& x, const MatrixXd& lam, int m, bool pattern_pass,
                    Workspace* workspace)
{
    const int n = (int) x.rows();

    vector<double> lam_d((size_t) max(m,1), 1.0);
    if (!pattern_pass) for (int i = 0; i < m; i++) lam_d[(size_t) i] = lam(i);

    psopt_ad::ad_record(workspace->ad_hess, n, 1, &x(0),
        [&](const adouble* xin, adouble* yout)
        { yout[0] = sqp_lagrangian(const_cast<adouble*>(xin), &lam_d[0], m, workspace); });

    if (workspace->enable_nlp_counters)
        workspace->solution->mesh_stats[workspace->current_mesh_refinement_iteration-1]
            .n_hessian_evals++;

    return psopt_ad::ad_sparse_hessian(workspace->ad_hess, &x(0), /*reuse=*/false);
}


// ---------------------------------------------------------------------------------
// The constraint Jacobian, sparse.
//
// The pattern is obtained once, from the recorded tape, and the values are refreshed
// into it at every point. Under numerical derivatives there is no pattern to detect
// and the Jacobian is taken to be full, which costs what it always cost: n constraint
// evaluations per iteration and m*n stored numbers.
// ---------------------------------------------------------------------------------
static void sqp_jacobian_triplet(MatrixXd& x, int m, bool& tape_done,
                                 vector<unsigned int>& row, vector<unsigned int>& col,
                                 vector<double>& val, Workspace* workspace)
{
    const int n = (int) x.rows();

    if (workspace->enable_nlp_counters)
        workspace->solution->mesh_stats[workspace->current_mesh_refinement_iteration-1]
            .n_jacobian_evals++;

    if (useAutomaticDifferentiation(*workspace->algorithm)) {

        // One recording per call, reused for every evaluation after it. The sparsity
        // pattern has to be computed from that recording before it can be reused, so the
        // first evaluation after a recording asks for the structure pass; asking to
        // reuse a pattern that belongs to the previous mesh returns a Jacobian that is
        // wrong in every entry and says nothing about it.
        bool fresh_tape = false;
        if (!tape_done) {
            psopt_ad::ad_record(workspace->ad_g, n, m, &x(0),
                [&](const adouble* xin, adouble* yout)
                { gg_ad(const_cast<adouble*>(xin), yout, workspace); });
            tape_done  = true;
            fresh_tape = true;
        }
        psopt_ad::SparseTriplet T =
            psopt_ad::ad_sparse_jacobian(workspace->ad_g, &x(0), /*reuse=*/!fresh_tape);

        row.assign(T.row.begin(), T.row.end());
        col.assign(T.col.begin(), T.col.end());
        val.assign(T.val.begin(), T.val.end());
    }
    else {
        MatrixXd g0(m,1), g1(m,1);
        gg_num(x, &g0, workspace);
        const double sqreps = sqrt(PSOPT_extras::GetEPS());

        row.resize((size_t) m*n); col.resize((size_t) m*n); val.resize((size_t) m*n);
        for (int j = 0; j < n; j++) {
            const double xj   = x(j);
            const double delj = sqreps*(1.0 + fabs(xj));
            x(j) = xj + delj;
            gg_num(x, &g1, workspace);
            x(j) = xj;
            for (int i = 0; i < m; i++) {
                const size_t k = (size_t) j*m + i;
                row[k] = (unsigned int) i;
                col[k] = (unsigned int) j;
                val[k] = (g1(i) - g0(i))/delj;
            }
        }
    }
}


// =================================================================================
//  The SQP driver
// =================================================================================
int SQP_interface(Alg&         algorithm,
                  MatrixXd*    x0,
                  double     (*f)(MatrixXd&, Workspace*),
                  void       (*g)(MatrixXd&, MatrixXd*, Workspace*),
                  int          nlp_ncons,
                  int          nlp_neq,
                  MatrixXd*    xlb,
                  MatrixXd*    xub,
                  MatrixXd*    lambda,
                  int          hotflag,
                  int          iprint,
                  Workspace*   workspace,
                  void*        user_data)
{
    USING_NAMESPACE_QPOASES

    (void) f; (void) g; (void) nlp_neq; (void) hotflag; (void) user_data;

    const int n = (int) x0->rows();
    const int m = nlp_ncons;

    Sol* solution = workspace->solution;

    // ---- problem data -----------------------------------------------------------
    vector<double> gl(max(m,1)), gu(max(m,1));
    if (m > 0) get_constraint_bounds(&gl[0], &gu[0], workspace);

    MatrixXd x  = *x0;
    for (int j = 0; j < n; j++)                       // start inside the box
        x(j) = min(max(x(j), (*xlb)(j)), (*xub)(j));

    MatrixXd gf(n,1), gval(max(m,1),1);
    MatrixXd gf_old(n,1);
    MatrixXd lam  = MatrixXd::Zero(max(m,1),1);       // constraint multipliers
    MatrixXd zbnd = MatrixXd::Zero(n,1);              // bound multipliers
    vector<double> r(max(m,1), 0.0);                  // l1 penalty weights

    bool tape_done = false;

    const double tol      = algorithm.nlp_tolerance;
    const int    iter_max = algorithm.nlp_iter_max;

    // ---- the model of the Lagrangian's curvature --------------------------------
    // Either the exact Hessian, sparse, re-evaluated at every iterate, or a dense
    // quasi-Newton approximation built by damped BFGS updates. The exact Hessian
    // needs the taped derivatives to exist, so the option follows the same rule as
    // IPOPT's: algorithm.hessian == "exact" together with automatic derivatives.
    //
    // The difference between the two is not only the rate of convergence. The BFGS
    // model is a full n-by-n matrix whatever the problem's structure, and at
    // n = 10^4 that is 0.8 GB before the solver has done anything at all; the exact
    // Hessian of a collocated optimal control problem has a few entries per row.
    const bool exact_hessian = (algorithm.hessian == "exact")
                               && useAutomaticDifferentiation(algorithm);

    // Which QP solver handles the subproblems. qpOASES is an active-set method whose
    // factorisations are dense in the number of variables however sparse the matrices
    // it is handed, so its memory is quadratic and its work per subproblem cubic in n;
    // ProxQP factorises the KKT system sparsely and tolerates an indefinite Hessian.
    const bool use_extern = (algorithm.qp_solver != "qpOASES");

    MatrixXd B;                                       // the quasi-Newton model, if used
    if (!exact_hessian) B = MatrixXd::Identity(n,n);

    // ---- the trust region, for the exact Hessian only ---------------------------
    // A line search needs the model to be positive definite to give it a direction
    // worth searching along, and the quasi-Newton model is built to be. The exact
    // Hessian is not: it is indefinite wherever the problem is, and on a minimum-time
    // problem, at a first iterate where the multipliers are still zero, it is exactly
    // zero -- the Hessian of a linear objective. The subproblem is then a linear
    // programme whose solution sits wherever the user's variable bounds happen to be,
    // and on examples/brac1 the first step drove the objective to 10^-14 and the
    // constraint violation to 0.93. Restricting the step to a region around the
    // current point is what makes an exact-Hessian model usable, and the region grows
    // as the model proves itself.
    const double Delta_0   = 1.0;                     // in the scaled variables, O(1)
    const double Delta_min = 1.0e-10;
    const double Delta_max = 1.0e4;
    double       Delta     = exact_hessian ? Delta_0 : qpOASES::INFTY;

    // ---- the subproblem's matrices ----------------------------------------------
    SparseCsc Jm, Hm, Ae, He;                         // A and H, plain and elastic
    vector<unsigned int> jrow, jcol;
    vector<double>       jval, jval_old;
    vector<double>       hval, aeval, heval;

    vector<double> lbA(max(m,1)), ubA(max(m,1)), lbd(n), ubd(n);
    vector<double> lo_true(n), hi_true(n);
    vector<double> dsol(n), ysol(n + max(m,1));

    double fval = ff_num(x, workspace);
    if (m > 0) gg_num(x, &gval, workspace);
    sqp_gradient(x, gf, workspace);

    if (m > 0) {
        sqp_jacobian_triplet(x, m, tape_done, jrow, jcol, jval, workspace);
        Jm.build(m, n, jrow, jcol, /*force_diagonal=*/false);
        Jm.scatter(jval);
    }

    // The Hessian's pattern, detected once with every multiplier set to one so that it
    // is the pattern of the problem rather than of the multipliers this iterate happens
    // to carry. The diagonal is forced into it so that the convexification below has
    // somewhere to write.
    psopt_ad::SparseTriplet Ht;
    if (exact_hessian) {   // NOLINT: the pattern pass, described above
        Ht = sqp_hessian_triplet(x, lam, m, /*pattern_pass=*/true, workspace);
        vector<unsigned int> hrow, hcol;
        hrow.reserve(2*Ht.row.size()); hcol.reserve(2*Ht.col.size());
        for (int k = 0; k < Ht.nnz(); k++) {          // the lower triangle, and its mirror
            hrow.push_back(Ht.row[k]); hcol.push_back(Ht.col[k]);
            if (Ht.row[k] != Ht.col[k]) { hrow.push_back(Ht.col[k]); hcol.push_back(Ht.row[k]); }
        }
        Hm.build(n, n, hrow, hcol, /*force_diagonal=*/true);
    }

    // The subproblem's matrices, as qpOASES sees them. qpOASES keeps the pointers it
    // is given and reads through them at every solve, so these objects are built once,
    // over storage whose pattern does not change, and only the values are refreshed.
    // SymSparseMat additionally wants to know where each column's diagonal entry sits,
    // which is a property of the pattern and so is also computed once.
    unique_ptr<SparseMatrix>  Aqp;
    unique_ptr<SymSparseMat>  Hsparse;
    unique_ptr<SymDenseMat>   Hdense;

    if (m > 0) Aqp.reset(new SparseMatrix(m, n, Jm.ir(), Jm.jc(), Jm.val()));
    if (exact_hessian) {
        Hsparse.reset(new SymSparseMat(n, n, Hm.ir(), Hm.jc(), Hm.val()));
        Hsparse->createDiagInfo();
    }
    else {
        Hdense.reset(new SymDenseMat(n, n, n, &B(0,0)));
    }
    SymmetricMatrix* Hqp = exact_hessian ? (SymmetricMatrix*) Hsparse.get()
                                         : (SymmetricMatrix*) Hdense.get();

    int    n_restorations = 0;
    int    n_corrections  = 0;
    int    n_shifts       = 0;        // convexifications of the exact Hessian
    int    n_shrinks      = 0;        // reductions of the trust region
    // The first subproblem is solved with the multipliers still at zero, so the
    // Hessian of the Lagrangian it is built from carries no information about the
    // constraints at all -- on a minimum-time problem, whose objective is linear, it
    // is identically zero and the subproblem is a linear programme. Starting the shift
    // at one makes that first model the same one a quasi-Newton method would start
    // from, and halving it at every iteration that does not need it hands the model
    // back to the exact curvature within a few steps.
    double tau_last       = 1.0;      // the last shift that worked
    int    status  = 1;                                // 1 = iteration limit
    string message = "Maximum number of SQP iterations reached";
    int    iter    = 0;

    if (iprint) {
        snprintf(workspace->text, sizeof(workspace->text),
                 "\n\nSQP (%s + %s): %d variables, %d constraints",
                 exact_hessian ? "exact sparse Hessian" : "dense BFGS",
                 algorithm.qp_solver.c_str(), n, m);
        psopt_print(workspace, workspace->text);
        snprintf(workspace->text, sizeof(workspace->text),
                 "\n%d Jacobian nonzeros (%.3f%% dense)%s\n",
                 Jm.nnz(), 100.0*Jm.nnz()/max(1.0,(double) m*n),
                 exact_hessian ? "" : "");
        psopt_print(workspace, workspace->text);
        snprintf(workspace->text, sizeof(workspace->text),
                 "\n%5s %16s %12s %12s %10s %8s\n",
                 "iter", "objective(sc)", "max viol", "dual err", "step", "QP its");
        psopt_print(workspace, workspace->text);
    }

    for (iter = 0; iter < iter_max; iter++) {

        // ---- optimality tests ---------------------------------------------------
        // Stationarity of the Lagrangian L = f + lambda' g - z' x, measured relative
        // to the size of the multipliers so that the test is scale aware.
        MatrixXd dL = gf;
        if (m > 0) {
            MatrixXd Jtl(n,1);
            Jm.transpose_times(&lam(0), &Jtl(0));
            dL += Jtl;
        }
        dL -= zbnd;
        double dual_err = 0.0;
        for (int j = 0; j < n; j++) dual_err = max(dual_err, fabs(dL(j)));

        // Scaled as IPOPT scales it. An earlier version divided by the largest
        // multiplier, which on a problem with large multipliers divides the residual
        // down until any point passes: examples/hypersensitive was declared optimal
        // at an objective of 2.4e-3 where the answer is 1.33. Averaging the
        // multipliers and capping the divisor at s_max removes that failure.
        const double s_max = 100.0;
        double lam_1 = 0.0;
        for (int i = 0; i < m; i++) lam_1 += fabs(lam(i));
        for (int j = 0; j < n; j++) lam_1 += fabs(zbnd(j));
        const double s_d = max(s_max, lam_1/(double)(m + n))/s_max;
        dual_err /= s_d;

        const double viol = (m > 0) ? max_violation(gval, gl, gu) : 0.0;

        if (viol <= tol && dual_err <= tol) {
            status  = 0;
            message = "Optimal solution found";
            break;
        }

        // The exact Hessian is a function of the multipliers as well as the point, so
        // it is re-evaluated here, after the optimality test has read the multipliers
        // the previous subproblem returned.
        if (exact_hessian) {
            psopt_ad::SparseTriplet Hk =
                sqp_hessian_triplet(x, lam, m, /*pattern_pass=*/false, workspace);
            hval.assign(Hk.val.begin(), Hk.val.end());
            hval.reserve(2*Hk.val.size());
            for (int k = 0; k < Hk.nnz(); k++)
                if (Hk.row[k] != Hk.col[k]) hval.push_back(Hk.val[k]);

            // The mirrored entries must follow the same order the pattern was built
            // in, which interleaves each off-diagonal with its transpose. Rebuild the
            // value list in that order rather than appending.
            hval.clear();
            for (int k = 0; k < Hk.nnz(); k++) {
                hval.push_back(Hk.val[k]);
                if (Hk.row[k] != Hk.col[k]) hval.push_back(Hk.val[k]);
            }
        }

        // ---- the quadratic programming subproblem -------------------------------
        //   min  1/2 d' B d + grad_f' d
        //   s.t. g_l - g <= J d <= g_u - g,      xlb - x <= d <= xub - x
        // PSOPT writes an absent constraint bound as +/- 1.0e20 (see NLP_bounds.cxx);
        // qpOASES has its own infinity and treats anything beyond it as free.
        const double psopt_inf = 1.0e20;
        for (int i = 0; i < m; i++) {
            lbA[i] = (gl[i] <= -psopt_inf) ? -qpOASES::INFTY : gl[i] - gval(i);
            ubA[i] = (gu[i] >=  psopt_inf) ?  qpOASES::INFTY : gu[i] - gval(i);
        }
        for (int j = 0; j < n; j++) {
            lo_true[j] = ((*xlb)(j) <= -psopt_inf) ? -qpOASES::INFTY : (*xlb)(j) - x(j);
            hi_true[j] = ((*xub)(j) >=  psopt_inf) ?  qpOASES::INFTY : (*xub)(j) - x(j);
            lbd[j]     = max(lo_true[j], -Delta);   // Delta is infinite unless the Hessian
            ubd[j]     = min(hi_true[j],  Delta);   // is exact; both regions contain d = 0
        }

        Options qpopts;
        qpopts.printLevel = PL_NONE;
        qpopts.setToReliable();
        if (exact_hessian) qpopts.enableRegularisation = BT_TRUE;
        qpopts.printLevel = PL_NONE;

        int_t nWSR = 5*(n + m) + 100;
        returnValue rv;
        bool   elastic = false;
        double rho_elastic = 0.0;

        // ---- convexification ----------------------------------------------------
        // An exact Hessian of the Lagrangian is indefinite wherever the problem is not
        // convex, which is almost everywhere on an optimal control problem, and the
        // subproblem then has no minimum: qpOASES declines it. Adding tau*I makes the
        // model convex without changing where the Lagrangian is stationary, and the
        // smallest tau that works is the one that leaves the most of the true
        // curvature in place, so tau starts small and is raised only until the
        // subproblem is accepted. The starting value is carried between iterations,
        // because a problem that needed a shift once usually needs one again, and
        // decays when it is not needed, so that the model returns to the exact one as
        // the iterates settle. This is the same device IPOPT uses in its own inertia
        // correction.
        double tau = 0.0;
        for (int attempt = 0; ; attempt++) {

            if (exact_hessian) {
                Hm.scatter(hval);
                if (attempt == 0) tau = tau_last*0.5;
                else              tau = (tau > 0.0) ? 10.0*tau
                                                    : 1.0e-6*max(1.0, Hm.norm1());
                if (tau > 0.0) Hm.shift_diagonal(tau);
            }

            if (use_extern) {
                QpProblem qpp;
                qpp.nv  = n;              qpp.mc  = m;
                qpp.H   = exact_hessian ? &Hm : NULL;
                qpp.Bd  = exact_hessian ? NULL : &B;
                qpp.g   = &gf(0);
                qpp.A   = (m > 0) ? &Jm : NULL;
                qpp.lbA = (m > 0) ? &lbA[0] : NULL;
                qpp.ubA = (m > 0) ? &ubA[0] : NULL;
                qpp.lbd = &lbd[0];        qpp.ubd = &ubd[0];

                QpSolution qs;
                string why;
                if (!solve_qp_plugin(algorithm.qp_solver, qpp, tol, 200,
                                     exact_hessian, qs, why)) {
                    status  = 2;
                    message = why;
                    break;
                }

                rv   = qs.ok ? SUCCESSFUL_RETURN : RET_INIT_FAILED;
                nWSR = qs.iterations;
                if (qs.ok) {
                    // Reported in qpOASES's convention, which everything below expects.
                    for (int j = 0; j < n; j++) dsol[j]   = qs.d[(size_t) j];
                    for (int j = 0; j < n; j++) ysol[j]   = qs.z[(size_t) j];
                    for (int i = 0; i < m; i++) ysol[n+i] = -qs.lambda[(size_t) i];
                }
            }
            else if (m > 0) {
                QProblem qp(n, m);
                qp.setOptions(qpopts);
                nWSR = 5*(n + m) + 100;
                rv = qp.init(Hqp, &gf(0), Aqp.get(), &lbd[0], &ubd[0], &lbA[0], &ubA[0], nWSR);
                if (rv == SUCCESSFUL_RETURN || rv == RET_MAX_NWSR_REACHED) {
                    qp.getPrimalSolution(&dsol[0]);
                    qp.getDualSolution(&ysol[0]);
                }
            }
            else {
                QProblemB qp(n);
                qp.setOptions(qpopts);
                nWSR = 5*n + 100;
                rv = qp.init(Hqp, &gf(0), &lbd[0], &ubd[0], nWSR);
                if (rv == SUCCESSFUL_RETURN || rv == RET_MAX_NWSR_REACHED) {
                    qp.getPrimalSolution(&dsol[0]);
                    qp.getDualSolution(&ysol[0]);
                }
            }

            if (!exact_hessian) break;
            if (rv == SUCCESSFUL_RETURN || rv == RET_MAX_NWSR_REACHED) break;

            // Only a failure of the Cholesky factorisation says anything about the
            // curvature; an infeasible linearisation or an exhausted working set says
            // nothing, and shifting the diagonal in answer to it would be nine wasted
            // factorisations before the restoration step that was needed all along.
            const bool curvature_failure =
                   (rv == RET_INIT_FAILED_CHOLESKY)
                || (rv == RET_INIT_FAILED_REGULARISATION)
                || (rv == RET_HESSIAN_NOT_SPD)
                || (rv == RET_HESSIAN_INDEFINITE);

            if (!curvature_failure || attempt >= 8) break;
            n_shifts++;
        }
        if (exact_hessian) tau_last = max(1.0e-12, tau);

        // ---- restoration ---------------------------------------------------
        // The linearised constraints can be inconsistent even when the problem is
        // perfectly well posed -- it is the normal situation when the starting guess
        // violates an equality, which is how most people supply one. The subproblem
        // then has no solution and the iteration would simply stop. Elastic mode
        // relaxes every constraint by a pair of non-negative slacks and charges them
        // in the objective,
        //
        //   min 1/2 d'Bd + grad_f'd + rho*sum(v + w)
        //   s.t. lbA <= J d + v - w <= ubA,  lb <= d <= ub,  v, w >= 0,
        //
        // which is always feasible, so a step always exists. With rho above the
        // current penalty weights the step reduces infeasibility; the slacks are
        // given a small quadratic term as well, to keep the QP's Hessian positive
        // definite rather than merely semidefinite. This is SNOPT's device (Gill,
        // Murray and Saunders 2005), in the simplest form that does the job.
        if (m > 0 && rv != SUCCESSFUL_RETURN && rv != RET_MAX_NWSR_REACHED) {

            double rho = 10.0;
            for (int i = 0; i < m; i++) rho = max(rho, 10.0*r[i]);
            rho_elastic = rho;

            const int ne = n + 2*m;
            vector<double> ge(ne), lbe(ne), ube(ne), ze(ne), ye(ne + m);

            const double eps_slack = 1.0e-8;

            // [ J | I | -I ] and blkdiag(H, eps I), both built on the pattern of the
            // matrices they extend, so the elastic subproblem costs the slacks and
            // nothing else. The patterns are built the first time restoration is
            // needed and kept thereafter.
            if (Ae.empty()) {
                vector<unsigned int> ar, ac;
                ar.reserve(jrow.size() + 2*(size_t)m); ac.reserve(jcol.size() + 2*(size_t)m);
                for (size_t k = 0; k < jrow.size(); k++) { ar.push_back(jrow[k]); ac.push_back(jcol[k]); }
                for (int i = 0; i < m; i++) { ar.push_back((unsigned int) i); ac.push_back((unsigned int)(n+i)); }
                for (int i = 0; i < m; i++) { ar.push_back((unsigned int) i); ac.push_back((unsigned int)(n+m+i)); }
                Ae.build(m, ne, ar, ac, /*force_diagonal=*/false);

                vector<unsigned int> hr, hc;
                if (exact_hessian) {
                    for (int k = 0; k < Ht.nnz(); k++) {
                        hr.push_back(Ht.row[k]); hc.push_back(Ht.col[k]);
                        if (Ht.row[k] != Ht.col[k]) { hr.push_back(Ht.col[k]); hc.push_back(Ht.row[k]); }
                    }
                }
                else {
                    for (int i = 0; i < n; i++)
                        for (int j = 0; j < n; j++) { hr.push_back((unsigned int) i); hc.push_back((unsigned int) j); }
                }
                for (int k = n; k < ne; k++) { hr.push_back((unsigned int) k); hc.push_back((unsigned int) k); }
                He.build(ne, ne, hr, hc, /*force_diagonal=*/true);
            }

            aeval.clear();
            aeval.insert(aeval.end(), jval.begin(), jval.end());
            for (int i = 0; i < m; i++) aeval.push_back( 1.0);
            for (int i = 0; i < m; i++) aeval.push_back(-1.0);
            Ae.scatter(aeval);

            heval.clear();
            if (exact_hessian) heval.insert(heval.end(), hval.begin(), hval.end());
            else for (int i = 0; i < n; i++) for (int j = 0; j < n; j++) heval.push_back(B(i,j));
            for (int k = n; k < ne; k++) heval.push_back(eps_slack);
            He.scatter(heval);

            // The relaxed subproblem carries the same model, and so needs the same
            // convexification: it is the shift that made the model usable, not the
            // relaxation. The slacks' own small diagonal is already in place and is
            // raised along with the rest, which does them no harm.
            if (exact_hessian && tau > 0.0) He.shift_diagonal(tau);

            // The subproblem is stated with its objective divided by rho, so that the
            // slacks cost one apiece instead of rho and the quadratic block is scaled to
            // match. The minimiser is unchanged -- scaling an objective does not move
            // it -- but the conditioning is transformed: priced directly, a slack with a
            // linear cost of 10^7 against a regularising quadratic term of 10^-8 spans
            // fifteen orders of magnitude, and a proximal method spent twenty-three
            // thousand iterations on examples/brac1 failing to get through it. The
            // multipliers come back scaled by 1/rho and are restored below.
            // ...and it is applied only for ProxQP. qpOASES's active-set method is not
            // troubled by the unscaled pricing and is measurably worse with it:
            // examples/bryson_denham goes from twenty iterations to failing after two
            // hundred and forty-six.
            const double oscale = use_extern ? 1.0/rho : 1.0;
            for (int j = 0; j < n; j++) { ge[j] = oscale*gf(j);   lbe[j] = lbd[j]; ube[j] = ubd[j]; }
            for (int k = n; k < ne; k++) { ge[k] = oscale*rho;    lbe[k] = 0.0;    ube[k] = qpOASES::INFTY; }
            He.scale(oscale);

            returnValue rve;
            int_t        nWSRe = 5*(ne + m) + 100;

            SymSparseMat Heqp(ne, ne, He.ir(), He.jc(), He.val());
            SparseMatrix Aeqp(m, ne, Ae.ir(), Ae.jc(), Ae.val());

            if (use_extern) {
                QpProblem qpp;
                qpp.nv  = ne;      qpp.mc  = m;
                qpp.H   = &He;     qpp.g   = &ge[0];
                qpp.A   = &Ae;
                qpp.lbA = &lbA[0]; qpp.ubA = &ubA[0];
                qpp.lbd = &lbe[0]; qpp.ubd = &ube[0];

                QpSolution qs;
                string why;
                (void) solve_qp_plugin(algorithm.qp_solver, qpp, tol, 200,
                                       exact_hessian, qs, why);

                rve   = qs.ok ? SUCCESSFUL_RETURN : RET_INIT_FAILED;
                nWSRe = qs.iterations;
                if (qs.ok) {
                    for (int k = 0; k < ne; k++) ze[k]    = qs.d[(size_t) k];
                    for (int k = 0; k < ne; k++) ye[k]    = qs.z[(size_t) k];
                    for (int i = 0; i < m;  i++) ye[ne+i] = -qs.lambda[(size_t) i];
                }
            }
            else {
                Heqp.createDiagInfo();
                QProblem qpe(ne, m);
                qpe.setOptions(qpopts);
                rve = qpe.init(&Heqp, &ge[0], &Aeqp, &lbe[0], &ube[0],
                               &lbA[0], &ubA[0], nWSRe);
                if (rve == SUCCESSFUL_RETURN || rve == RET_MAX_NWSR_REACHED) {
                    qpe.getPrimalSolution(&ze[0]);
                    qpe.getDualSolution(&ye[0]);
                }
            }

            // Undo the 1/rho scaling of the objective on the multipliers it scaled.
            if (oscale != 1.0 && (rve == SUCCESSFUL_RETURN || rve == RET_MAX_NWSR_REACHED))
                for (size_t k = 0; k < ye.size(); k++) ye[k] /= oscale;

            if (rve == SUCCESSFUL_RETURN || rve == RET_MAX_NWSR_REACHED) {
                for (int j = 0; j < n; j++) dsol[j]   = ze[j];
                for (int i = 0; i < m; i++) ysol[n+i] = ye[ne+i];
                for (int j = 0; j < n; j++) ysol[j]   = ye[j];
                rv       = SUCCESSFUL_RETURN;
                elastic  = true;
                n_restorations++;
            }
        }

        if (rv != SUCCESSFUL_RETURN && rv != RET_MAX_NWSR_REACHED) {
            status  = 2;
            message = "The quadratic programming subproblem could not be solved, "
                      "and neither could its elastic relaxation";
            break;
        }

        MatrixXd d(n,1);
        for (int j = 0; j < n; j++) d(j) = dsol[j];

        // qpOASES returns the multipliers of  H d + grad_f = A' y_C + y_B, so its
        // duals carry the opposite sign to the convention used by the rest of PSOPT,
        // in which grad f + J' lambda = 0. The reference test for this is
        // SQPSolver.QpOasesDualSignConvention in tests/test_sqp.cpp.
        MatrixXd lam_new = MatrixXd::Zero(max(m,1),1);
        for (int i = 0; i < m; i++) lam_new(i) = -ysol[n+i];
        // The bounds the subproblem was given are the variable bounds intersected with
        // the trust region, so a multiplier returned for one of them belongs to
        // whichever of the two is binding. Only the variable bound is a bound of the
        // problem: a multiplier earned against the trust region is an artefact of the
        // current radius, and counting it in the KKT residual cancels part of the
        // gradient and reports a convergence that has not happened. On examples/brac1
        // that produced a feasible point, declared optimal, whose objective was 2.8 per
        // cent above the answer.
        for (int j = 0; j < n; j++) {
            const bool tr_binds = (lbd[j] > lo_true[j]) || (ubd[j] < hi_true[j]);
            zbnd(j) = (tr_binds && fabs(d(j)) >= 0.999*Delta) ? 0.0 : ysol[j];
        }

        // ---- the l1 penalty weights ---------------------------------------------
        // Powell's rule: keep each weight at least as large as the magnitude of the
        // multiplier it accompanies, which is what makes d a descent direction for
        // the merit function, and let it decay towards that value when it may.
        for (int i = 0; i < m; i++) {
            const double li = fabs(lam_new(i));
            r[i] = (iter == 0) ? li : max(li, 0.5*(r[i] + li));
        }

        // Powell's rule alone leaves the weights at the mercy of the multipliers, and
        // the multipliers of a first subproblem solved from a poor point can be
        // anything at all, including nothing. Weights near zero make the merit function
        // the objective again, and the line search then cheerfully accepts a step that
        // improves the objective by destroying feasibility -- on examples/brac1, where
        // the objective is the final time and can always be improved by shortening it,
        // the very first exact-Hessian step drove the objective to 10^-14 and the
        // violation to 0.93. Han's condition removes the possibility: the weights are
        // raised until the model's own predicted change in the objective is dominated
        // by the reduction in infeasibility the step is credited with, which makes the
        // step a descent direction for the merit function whatever the multipliers say.
        {
            double gTd0 = 0.0;
            for (int j = 0; j < n; j++) gTd0 += gf(j)*d(j);
            const double dBd0 = quadratic_form(exact_hessian, Hm, B, d);
            const double viol_tot = (m > 0) ? total_violation(gval, gl, gu) : 0.0;

            if (viol_tot > 0.0) {
                const double need = (gTd0 + max(0.0, 0.5*dBd0))/(0.9*viol_tot);
                if (need > 0.0) for (int i = 0; i < m; i++) r[i] = max(r[i], need);
            }
        }

        // In elastic mode the step minimises f + rho*(violation), so that, and not the
        // merit function built from the multipliers, is what the step was computed to
        // reduce. Testing it against any smaller weight asks the step to deliver a
        // decrease it was never aiming at, and the line search then rejects a perfectly
        // good restoration step -- which is what stopped examples/hypersensitive on its
        // second mesh. The weights relax again through the averaging rule above once
        // ordinary steps resume.
        if (elastic) for (int i = 0; i < m; i++) r[i] = max(r[i], rho_elastic);

        // ---- line search ---------------------------------------------------------
        const double phi0  = merit(fval, gval, gl, gu, r);
        double gTd = 0.0;
        for (int j = 0; j < n; j++) gTd += gf(j)*d(j);
        double pen = 0.0;
        for (int i = 0; i < m; i++) pen += r[i]*violation(gval(i), gl[i], gu[i]);
        const double dphi = gTd - pen;      // directional derivative of the merit function

        double alpha = 1.0;
        const double eta = 1.0e-4;
        MatrixXd xtrial(n,1), gtrial(max(m,1),1);
        double ftrial = fval, phit = phi0;
        bool   accepted = false;

        for (int ls = 0; ls < 25; ls++) {
            xtrial = x + alpha*d;
            for (int j = 0; j < n; j++)                 // the QP respects the box, so
                xtrial(j) = min(max(xtrial(j), (*xlb)(j)), (*xub)(j));   // this only cleans rounding
            ftrial = ff_num(xtrial, workspace);
            if (m > 0) gg_num(xtrial, &gtrial, workspace);
            phit = merit(ftrial, gtrial, gl, gu, r);
            if (std::isfinite(phit) && phit <= phi0 + eta*alpha*dphi) { accepted = true; break; }
            alpha *= 0.5;
        }

        // ---- second-order correction --------------------------------------------
        // The step can fail to reduce the merit function even when it is a good step,
        // because a linear model of a curved constraint over-estimates the violation
        // that a unit step incurs -- the Maratos effect. It bites hardest near a
        // solution, which is exactly where a warm start from the previous mesh begins,
        // and it is what stopped examples/hypersensitive on its second mesh. The
        // remedy is to re-solve the subproblem with the constraint right-hand side
        // evaluated at the trial point rather than linearised from the current one,
        // and to try the corrected step before giving up (Nocedal and Wright,
        // algorithm 18.4).
        if (!accepted && m > 0 && !elastic) {

            MatrixXd gtrial_soc(m,1);
            MatrixXd xfull = x + d;
            for (int j = 0; j < n; j++)
                xfull(j) = min(max(xfull(j), (*xlb)(j)), (*xub)(j));
            gg_num(xfull, &gtrial_soc, workspace);

            MatrixXd Jd(m,1);
            Jm.times(&d(0), &Jd(0));
            vector<double> lbS(m), ubS(m);
            for (int i = 0; i < m; i++) {
                const double shift = gtrial_soc(i) - Jd(i);
                lbS[i] = (gl[i] <= -psopt_inf) ? -qpOASES::INFTY : gl[i] - shift;
                ubS[i] = (gu[i] >=  psopt_inf) ?  qpOASES::INFTY : gu[i] - shift;
            }

            int_t nWSRs = 5*(n + m) + 100;
            vector<double> dsoc(n), ysoc(n + m);
            returnValue rvs;

            if (use_extern) {
                QpProblem qpp;
                qpp.nv  = n;              qpp.mc  = m;
                qpp.H   = exact_hessian ? &Hm : NULL;
                qpp.Bd  = exact_hessian ? NULL : &B;
                qpp.g   = &gf(0);         qpp.A   = &Jm;
                qpp.lbA = &lbS[0];        qpp.ubA = &ubS[0];
                qpp.lbd = &lbd[0];        qpp.ubd = &ubd[0];

                QpSolution qs;
                string why;
                (void) solve_qp_plugin(algorithm.qp_solver, qpp, tol, 200,
                                       exact_hessian, qs, why);
                rvs = qs.ok ? SUCCESSFUL_RETURN : RET_INIT_FAILED;
                if (qs.ok) {
                    for (int j = 0; j < n; j++) dsoc[j]   = qs.d[(size_t) j];
                    for (int j = 0; j < n; j++) ysoc[j]   = qs.z[(size_t) j];
                    for (int i = 0; i < m; i++) ysoc[n+i] = -qs.lambda[(size_t) i];
                }
            }
            else {
                QProblem qps(n, m);
                qps.setOptions(qpopts);
                rvs = qps.init(Hqp, &gf(0), Aqp.get(), &lbd[0], &ubd[0],
                               &lbS[0], &ubS[0], nWSRs);
                if (rvs == SUCCESSFUL_RETURN || rvs == RET_MAX_NWSR_REACHED) {
                    qps.getPrimalSolution(&dsoc[0]);
                    qps.getDualSolution(&ysoc[0]);
                }
            }

            if (rvs == SUCCESSFUL_RETURN || rvs == RET_MAX_NWSR_REACHED) {

                MatrixXd dc(n,1);
                for (int j = 0; j < n; j++) dc(j) = dsoc[j];

                double gTdc = 0.0;
                for (int j = 0; j < n; j++) gTdc += gf(j)*dc(j);
                const double dphi_c = gTdc - pen;

                alpha = 1.0;
                for (int ls = 0; ls < 25; ls++) {
                    xtrial = x + alpha*dc;
                    for (int j = 0; j < n; j++)
                        xtrial(j) = min(max(xtrial(j), (*xlb)(j)), (*xub)(j));
                    ftrial = ff_num(xtrial, workspace);
                    gg_num(xtrial, &gtrial, workspace);
                    phit = merit(ftrial, gtrial, gl, gu, r);
                    if (std::isfinite(phit) && phit <= phi0 + eta*alpha*dphi_c) {
                        accepted = true;
                        d = dc;
                        for (int i = 0; i < m; i++) lam_new(i) = -ysoc[n+i];
                        for (int j = 0; j < n; j++) zbnd(j)    =  ysoc[j];
                        n_corrections++;
                        break;
                    }
                    alpha *= 0.5;
                }
            }
        }

        if (!accepted) {
            // With an exact Hessian there is nothing wrong with the model that a
            // smaller region will not cure: shrink it, which changes the direction and
            // not merely the length of the step, and try the same point again.
            if (exact_hessian && Delta > Delta_min) {
                Delta = max(Delta_min, 0.25*Delta);
                if (iprint) {
                    snprintf(workspace->text, sizeof(workspace->text),
                             "   line search failed; trust region reduced to %.2e\n", Delta);
                    psopt_print(workspace, workspace->text);
                }
                n_shrinks++;
                continue;
            }
            // A merit function that cannot be decreased along the QP direction means
            // the quadratic model is not usable here. With a quasi-Newton model that
            // can be a model gone stale, so it is reset and the point tried once more.
            if (exact_hessian || B.isApprox(MatrixXd::Identity(n,n))) {
                status  = 3;
                message = "The line search failed to decrease the merit function";
                break;
            }
            B = MatrixXd::Identity(n,n);
            if (iprint) {
                snprintf(workspace->text, sizeof(workspace->text),
                         "\n   line search failed; resetting the Hessian approximation\n");
                psopt_print(workspace, workspace->text);
            }
            continue;
        }

        // The region is enlarged only when the step both reached its boundary and was
        // taken whole: a step the line search had to shorten has already shown that the
        // model is trusted too far, and one that stopped short of the boundary was not
        // constrained by it in the first place.
        if (exact_hessian) {
            double dinf = 0.0;
            for (int j = 0; j < n; j++) dinf = max(dinf, fabs(d(j)));
            if (alpha >= 1.0 && dinf >= 0.9*Delta) Delta = min(2.0*Delta, Delta_max);
            else if (alpha < 0.25)                 Delta = max(Delta_min, 0.5*Delta);
        }

        // ---- move to the new point ----------------------------------------------
        gf_old = gf;
        if (m > 0 && !exact_hessian) jval_old = jval;

        MatrixXd s = xtrial - x;

        x    = xtrial;
        fval = ftrial;
        if (m > 0) gval = gtrial;
        lam  = lam_new;

        sqp_gradient(x, gf, workspace);
        if (m > 0) {
            sqp_jacobian_triplet(x, m, tape_done, jrow, jcol, jval, workspace);
            Jm.scatter(jval);
        }

        if (!exact_hessian) {
            // ---- BFGS update, damped --------------------------------------------
            // s is the accepted step and y the change in the gradient of the
            // Lagrangian at fixed multipliers. Powell's damping keeps s'y positive,
            // and with it the positive definiteness that the QP solver requires of B.
            MatrixXd y = gf - gf_old;
            if (m > 0) {
                // (J - J_old)' lambda, on the shared pattern.
                MatrixXd dJtl(n,1);
                vector<double> dj(jval.size());
                for (size_t k = 0; k < jval.size(); k++) dj[k] = jval[k] - jval_old[k];
                SparseCsc& Jd = Jm;                     // same pattern, different values
                Jd.scatter(dj);
                Jd.transpose_times(&lam(0), &dJtl(0));
                Jd.scatter(jval);                       // put the Jacobian back
                y += dJtl;
            }

            const double sBs = (s.transpose()*B*s)(0,0);
            double sy = (s.transpose()*y)(0,0);

            if (sBs > 0.0) {
                if (sy < 0.2*sBs) {                       // Powell (1978)
                    const double theta = 0.8*sBs/(sBs - sy);
                    y  = theta*y + (1.0 - theta)*(B*s);
                    sy = (s.transpose()*y)(0,0);
                }
                if (sy > 1.0e-12*max(1.0, sBs)) {
                    MatrixXd Bs = B*s;
                    B += (y*y.transpose())/sy - (Bs*Bs.transpose())/sBs;
                }
            }
        }

        if (iprint) {
            snprintf(workspace->text, sizeof(workspace->text),
                     "%5d %16.8e %12.3e %12.3e %10.2e %8d %s%s\n",
                     iter+1, fval, (m > 0) ? max_violation(gval, gl, gu) : 0.0,
                     dual_err, alpha, (int) nWSR, elastic ? "restoration " : "",
                     (exact_hessian && tau_last > 1.0e-12) ? "shifted" : "");
            psopt_print(workspace, workspace->text);
        }
    }

    // ---- report ------------------------------------------------------------------
    *x0 = x;
    if (lambda != NULL && m > 0) *lambda = lam;

    solution->nlp_return_code = status;

    if (iprint) {
        snprintf(workspace->text, sizeof(workspace->text),
                 "\nSQP finished after %d iterations: %s\n"
                 "   objective (scaled)   %.10e\n"
                 "   maximum violation    %.3e\n"
                 "   restoration steps    %d\n"
                 "   second-order corr.   %d\n"
                 "   Hessian shifts       %d\n"
                 "   trust region cuts    %d\n",
                 iter, message.c_str(), fval,
                 (m > 0) ? max_violation(gval, gl, gu) : 0.0, n_restorations, n_corrections,
                 n_shifts, n_shrinks);
        psopt_print(workspace, workspace->text);
    }

    return status;
}

#else  // USE_SQP

int SQP_interface(Alg&, MatrixXd*, double (*)(MatrixXd&, Workspace*),
                  void (*)(MatrixXd&, MatrixXd*, Workspace*), int, int,
                  MatrixXd*, MatrixXd*, MatrixXd*, int, int,
                  Workspace* workspace, void*)
{
    error_message("algorithm.nlp_method = \"SQP\" requires PSOPT to be built with "
                  "-DWITH_SQP=ON and qpOASES available");
    return 1;
}

#endif // USE_SQP
