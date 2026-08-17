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

#include <Eigen/SparseCore>
#include <qpOASES.hpp>

extern "C" {
#include <dmumps_c.h>
}
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

// Betts's merit function (3rd ed., section 2.4), after Gill, Murray, Saunders and
// Wright: an augmented Lagrangian in the constraints,
//
//     M(x, lambda, s) = f + lambda'(c - s) + 1/2 (c - s)' Theta (c - s),
//
// in which the slacks s are the values the constraints would take if the step were
// exact, so that (c - s) measures the deviation from linearity rather than the
// violation. Two properties matter and neither belongs to the l1 function it replaces.
// It is smooth, so a line search has a derivative to work with and the Maratos effect
// does not arise. And it carries the multipliers explicitly, so the penalty weights do
// not have to dominate them: they can be chosen as small as the descent condition
// allows and recomputed from nothing at each iteration, where the l1 weights had to
// exceed the largest multiplier and so could only ratchet upwards -- on
// examples/hypersensitive they reached 10^10 and stayed there.
//
// The bound terms of Betts's formulation are omitted. He carries them because his
// iterates may sit outside their bounds; ours never do, since the subproblem is given
// the bounds directly and the trial point is clipped into the box, so (x - t) is
// identically zero and the terms contribute nothing.
//
// The sign convention is PSOPT's throughout: grad f + J'lambda - z = 0, where Betts
// writes g - G'lambda - nu = 0, so his lambda is the negative of the one used here and
// the signs below differ from his in print.
double merit_al(double fval, const MatrixXd& gval, const MatrixXd& lamv,
                const MatrixXd& sv, const vector<double>& theta)
{
    double M = fval;
    for (int i = 0; i < gval.rows(); i++) {
        const double r = gval(i) - sv(i);
        M += lamv(i)*r + 0.5*theta[(size_t) i]*r*r;
    }
    return M;
}

// The slacks at the start of a step, from Betts's (2.28): the constraint value shifted
// by the multiplier and clipped to its own bounds, which is the choice that minimises
// the merit function over s for the current x and lambda. With a weight still at zero
// the shift is unbounded and the limit is the projection of c itself.
void merit_slacks(const MatrixXd& gval, const MatrixXd& lamv, const vector<double>& theta,
                  const vector<double>& gl, const vector<double>& gu, MatrixXd& sv)
{
    for (int i = 0; i < gval.rows(); i++) {
        double c = gval(i);
        if (theta[(size_t) i] > 1.0e-300) c += lamv(i)/theta[(size_t) i];
        sv(i) = min(max(c, gl[(size_t) i]), gu[(size_t) i]);
    }
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

    // Gerschgorin's bound on the most negative eigenvalue: every eigenvalue lies in a
    // disc centred on a diagonal entry with radius that row's off-diagonal absolute sum.
    // It normalises the convexification shift, turning it from a quantity in the units of
    // the problem into a number between zero and one. Betts, 3rd ed., eq. (2.40).
    double gerschgorin_lower_bound() const
    {
        vector<double> dg(max(nc_,1), 0.0), rad(max(nc_,1), 0.0);
        for (int j = 0; j < nc_; j++)
            for (sparse_int_t k = jc_[(size_t) j]; k < jc_[(size_t) j+1]; k++) {
                const int i = (int) ir_[(size_t) k];
                if (i == j) dg[(size_t) j]  = val_[(size_t) k];
                else        rad[(size_t) j] += fabs(val_[(size_t) k]);
            }
        double sigma = 0.0;
        for (int j = 0; j < nc_; j++) sigma = min(sigma, dg[(size_t) j] - rad[(size_t) j]);
        return sigma;
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
// A vector with a NaN or an infinity in it is not an answer, whatever the solver's
// return code says. See the note in solve_qp_plugin.
static bool all_finite(const double* v, int k)
{
    for (int i = 0; i < k; i++) if (!std::isfinite(v[i])) return false;
    return true;
}

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

    // A backend that stops early can return a vector that is not a number. On
    // examples/zpm, GALAHAD's QPA came back from a relaxed subproblem with multipliers
    // of 2.1e+50 and then, at the next iteration, with exactly zero -- an overflow and
    // the NaN cascade after it. Those multipliers go into the Hessian of the
    // Lagrangian, whose Gerschgorin bound then reads -2.8e+50, and the Levenberg shift
    // that is scaled by it is carried into every iteration that follows. A single
    // non-finite entry poisons the next fifty steps.
    //
    // So a solution with a non-finite entry in it is not a solution, and is reported as
    // a failure: the restoration that follows is a step the solver knows how to take.
    if (out.ok) {
        for (int j = 0; j < p.nv && out.ok; j++)
            if (!std::isfinite(out.d[(size_t) j]) || !std::isfinite(out.z[(size_t) j]))
                out.ok = false;
        for (int i = 0; i < p.mc && out.ok; i++)
            if (!std::isfinite(out.lambda[(size_t) i])) out.ok = false;
    }

    // An approximate answer is welcome; an inadmissible one is not. A backend that
    // stops at its iteration limit returns wherever it had got to, and for an
    // interior-point or proximal method that is a point which need not satisfy the
    // simple bounds at all. GALAHAD's QPA, given the elastic subproblem of
    // examples/launch -- 1778 variables and an iteration limit it does not reach the
    // end of -- came back with a step twice the length of the trust region and a set
    // of slacks summing to minus twenty-two, the slacks being the variables bounded
    // below by zero. Nothing downstream can make sense of that: the line search
    // evaluates a point outside the box, the merit function judges it, the trust
    // region is not a trust region any more, and the diagnosis of whatever happens
    // next starts from a false premise.
    //
    // The bounds are a box, so the honest repair is to project onto it. It costs
    // nothing when the backend was right, which is nearly always, and it turns a
    // meaningless step into the nearest admissible one when it was not.
    if (out.ok)
        for (int j = 0; j < p.nv; j++)
            out.d[(size_t) j] = min(max(out.d[(size_t) j], p.lbd[j]), p.ubd[j]);

    return true;
}


// ---------------------------------------------------------------------------------
// Betts's constraint violation merit function, (2.51):
//
//   Mv(x) = sum_i chi^2(cL_i, c_i, cU_i) + sum_j chi^2(xL_j, x_j, xU_j),
//   chi(l, y, u) = max[0, l - y, y - u].
//
// A sum of squares, not the sum of absolute violations the optimality phase uses for
// reporting. The square is what makes it differentiable where it matters and what makes
// a least-distance step a descent direction for it.
// ---------------------------------------------------------------------------------
double violation_merit(const MatrixXd& gval, const vector<double>& gl,
                       const vector<double>& gu)
{
    double s = 0.0;
    for (int i = 0; i < gval.rows(); i++) {
        const double v = violation(gval(i), gl[(size_t) i], gu[(size_t) i]);
        s += v*v;
    }
    return s;
}

// The inertia of the KKT matrix
//
//     K = [ H   J' ]
//         [ J   0  ]
//
// which is what Betts's Hessian modification is actually steered by, and what none of
// the QP backends will tell us. The reduced Hessian of the Lagrangian is positive
// definite on the null space of the constraints exactly when K has n positive and m
// negative eigenvalues and none zero (Gould 1985; Betts, 3rd ed., eq. (2.41)). That is a
// sharp test, unlike a Cholesky of H alone: the full Hessian of a collocated optimal
// control problem is essentially never positive definite, so requiring it to be drives
// the Levenberg parameter up until the model is a scaled identity and the second-order
// information is gone -- measured, on examples/lts, as a run that reached the iteration
// limit still creeping in the eighth decimal place.
//
// MUMPS provides it, and PSOPT already links MUMPS: IPOPT's default linear solver is the
// same libdmumps_seq, so this adds no library, no symbol and none of the collision risk
// that the QP backends had to be separated to avoid. INFOG(12) counts the negative
// pivots and INFOG(28), with ICNTL(24) on, the null ones.
//
// GALAHAD's SLS also advertises this, through inform.negative_eigenvalues and
// inform.rank. On the same matrices, with inertia known by hand, MUMPS was right three
// times out of three and SLS returned values like (-23, 26, 0) -- which is either a
// misuse of its interface or a defect in it, and either way not something to build on
// without understanding it first. MUMPS is used here for that reason.
bool kkt_inertia(const SparseCsc& H, const SparseCsc& J, int n, int m,
                 int& npos, int& nneg, int& nzero)
{
    vector<MUMPS_INT> irn, jcn;
    vector<double>    val;
    irn.reserve((size_t) H.nnz()/2 + J.nnz());
    jcn.reserve(irn.capacity());
    val.reserve(irn.capacity());

    // The lower triangle of K, in coordinate form, one-based as MUMPS wants it. H is
    // stored with both triangles, so half of it is dropped here.
    vector<Eigen::Triplet<double, int> > th, tj;
    H.emit_triplets(th, int(0));
    for (size_t k = 0; k < th.size(); k++)
        if (th[k].row() >= th[k].col()) {
            irn.push_back(th[k].row() + 1); jcn.push_back(th[k].col() + 1);
            val.push_back(th[k].value());
        }
    if (m > 0) {
        J.emit_triplets(tj, int(0));
        for (size_t k = 0; k < tj.size(); k++) {          // J sits below the diagonal
            irn.push_back(n + tj[k].row() + 1); jcn.push_back(tj[k].col() + 1);
            val.push_back(tj[k].value());
        }
    }
    if (irn.empty()) return false;

    DMUMPS_STRUC_C id;
    id.comm_fortran = -987654;              // MPI_COMM_WORLD, in the sequential build
    id.par = 1;
    id.sym = 2;                             // general symmetric
    id.job = -1;
    dmumps_c(&id);

    id.icntl[0] = -1; id.icntl[1] = -1; id.icntl[2] = -1; id.icntl[3] = 0;  // silent
    id.icntl[23] = 1;                       // ICNTL(24): detect null pivots

    id.n   = (MUMPS_INT) (n + m);
    id.nnz = (MUMPS_INT8) irn.size();
    id.irn = &irn[0]; id.jcn = &jcn[0]; id.a = &val[0];
    id.job = 4;                             // analyse and factorize
    dmumps_c(&id);

    const bool ok = (id.infog[0] >= 0);
    nneg  = (int) id.infog[11];
    nzero = (int) id.infog[27];
    npos  = (n + m) - nneg - nzero;

    id.job = -2;
    dmumps_c(&id);
    return ok;
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
//  Betts's feasible point strategy, section 2.8
// =================================================================================
//
//  The first phase of his FM strategy, and the thing the optimality phase falls back on
//  when it concludes that the trouble is in the constraints rather than in the model.
//  Its whole point is what it leaves out: no objective, no multipliers, no Hessian of
//  the Lagrangian, no merit function that is trying to do two things at once. Section
//  2.7 states the premise plainly -- segregate difficulties caused by the constraints
//  from difficulties caused by the objective -- and this is the half that deals with the
//  constraints.
//
//  Two subproblems, in Betts's order of preference.
//
//  The primary one is a least distance program, (2.46)-(2.47): the shortest step that
//  satisfies the linearised constraints,
//
//      min 1/2 p'p   s.t.  bL <= (Gp; p) <= bU,
//
//  which for a problem of equalities alone is p = -G# b, the pseudoinverse solution, and
//  a generalisation of the Newton step. Its Hessian is the identity, so it is convex
//  whatever the problem is doing, and its reduced Hessian is positive definite whatever
//  the rank of G -- which is why no inertia test appears below. What it cannot survive
//  is an inconsistent linearisation: then it has no solution at all.
//
//  The fallback is the relaxation, (2.49)-(2.50): give the constraints a residual vector
//  and charge it quadratically,
//
//      min 1/2 p'p + rho/2 u'u   s.t.  bL <= (Gp + u; p) <= bU,
//
//  which always has a solution because u can absorb anything. Betts uses rho = 10^6 and
//  notes that the KKT matrix is then conditioned like rho^2, which he treats with
//  iterative refinement; here the value is smaller, because the backends have no such
//  refinement and because measurement did not reward the larger one.
//
//  The step is judged by the constraint violation and by nothing else -- Mv from (2.51),
//  a sum of squares. This is the part that has to be got right and that an earlier
//  attempt at this phase got wrong: with the objective removed from the subproblem but
//  the augmented Lagrangian still judging the step, the line search cut a perfectly good
//  feasibility step to a thirty-second because the objective had gone the wrong way, and
//  the violation fell by one per cent an iteration. Half of this phase is worse than
//  none of it.
//
//  The switching between the two is Betts's, section 2.8.2: prefer the least distance
//  step; when it fails, use the relaxation and keep using it; when the relaxation takes
//  a full step, which is approximately a Newton step, try the least distance step again.
//  He calls the logic somewhat ad hoc and reports that it works, which is also what is
//  found here.
//
//  Returns 0 if a feasible point was found, 1 on the iteration limit, 2 if neither
//  subproblem could be solved, 3 if no step reduced the violation.

static int feasible_point_phase(MatrixXd& x, int m,
                                const vector<double>& gl, const vector<double>& gu,
                                MatrixXd* xlb, MatrixXd* xub,
                                bool& tape_done,
                                vector<unsigned int>& jrow, vector<unsigned int>& jcol,
                                vector<double>& jval, SparseCsc& Jm,
                                MatrixXd& gval,
                                const string& backend, bool use_extern,
                                double tol, int max_iter, int qp_iter_max, int iprint,
                                int& n_iters, int& n_relaxed, Workspace* workspace)
{
    const int    n         = (int) x.rows();
    const int    nr        = n + m;                 // the relaxed subproblem's size
    const double psopt_inf = 1.0e20;
    const double rho       = 1.0e4;                 // Betts's 10^6, tempered; see above

    if (m <= 0) return 0;

    // The identity, for the least distance program, and [J | I] with the quadratic
    // penalty block, for the relaxation. Both patterns are fixed for the mesh.
    SparseCsc Ip, Ar, Hr;
    {
        vector<unsigned int> ir, ic;
        for (int j = 0; j < n; j++) { ir.push_back((unsigned int) j); ic.push_back((unsigned int) j); }
        Ip.build(n, n, ir, ic, /*force_diagonal=*/true);
        Ip.scatter(vector<double>((size_t) n, 1.0));

        vector<unsigned int> ar, ac;
        ar.reserve(jrow.size() + (size_t) m); ac.reserve(jcol.size() + (size_t) m);
        for (size_t k = 0; k < jrow.size(); k++) { ar.push_back(jrow[k]); ac.push_back(jcol[k]); }
        for (int i = 0; i < m; i++) { ar.push_back((unsigned int) i); ac.push_back((unsigned int)(n+i)); }
        Ar.build(m, nr, ar, ac, /*force_diagonal=*/false);

        vector<unsigned int> hr, hc;
        for (int k = 0; k < nr; k++) { hr.push_back((unsigned int) k); hc.push_back((unsigned int) k); }
        Hr.build(nr, nr, hr, hc, /*force_diagonal=*/true);
        vector<double> hv((size_t) nr);
        for (int k = 0; k < n;  k++) hv[(size_t) k] = 1.0;
        for (int k = n; k < nr; k++) hv[(size_t) k] = rho;
        Hr.scatter(hv);
    }

    vector<double> lbA((size_t) m), ubA((size_t) m);
    vector<double> lbd((size_t) nr), ubd((size_t) nr);
    vector<double> gzero((size_t) nr, 0.0);
    vector<double> dsol((size_t) nr), ysol((size_t) (nr + m));
    vector<double> arval;

    MatrixXd d(n,1), xtrial(n,1), gtrial(m,1);

    bool   primary    = true;                       // the least distance step, first
    int    status     = 1;
    double alpha_prev = 1.0;                        // reported, not used

    for (int iter = 0; iter < max_iter; iter++) {

        gg_num(x, &gval, workspace);
        const double viol = max_violation(gval, gl, gu);
        const double mv0  = violation_merit(gval, gl, gu);

        if (iprint) {
            snprintf(workspace->text, sizeof(workspace->text),
                     "%5d  %14s  %12.3e  %9.2e\n", iter + 1,
                     primary ? "least distance" : "relaxation", viol, alpha_prev);
            psopt_print(workspace, workspace->text);
        }

        if (viol <= tol) { status = 0; break; }

        sqp_jacobian_triplet(x, m, tape_done, jrow, jcol, jval, workspace);
        Jm.scatter(jval);

        for (int i = 0; i < m; i++) {
            lbA[(size_t) i] = (gl[(size_t) i] <= -psopt_inf) ? -qpOASES::INFTY
                                                             :  gl[(size_t) i] - gval(i);
            ubA[(size_t) i] = (gu[(size_t) i] >=  psopt_inf) ?  qpOASES::INFTY
                                                             :  gu[(size_t) i] - gval(i);
        }
        for (int j = 0; j < n; j++) {
            lbd[(size_t) j] = ((*xlb)(j) <= -psopt_inf) ? -qpOASES::INFTY : (*xlb)(j) - x(j);
            ubd[(size_t) j] = ((*xub)(j) >=  psopt_inf) ?  qpOASES::INFTY : (*xub)(j) - x(j);
        }
        for (int k = n; k < nr; k++) {              // the residuals are free
            lbd[(size_t) k] = -qpOASES::INFTY;
            ubd[(size_t) k] =  qpOASES::INFTY;
        }

        arval.clear();
        arval.insert(arval.end(), jval.begin(), jval.end());
        for (int i = 0; i < m; i++) arval.push_back(1.0);
        Ar.scatter(arval);

        // ---- the subproblem, primary or relaxed -----------------------------
        bool solved = false;
        for (int attempt = 0; attempt < 2 && !solved; attempt++) {

            const int  nv = primary ? n : nr;
            SparseCsc& H  = primary ? Ip : Hr;
            SparseCsc& A  = primary ? Jm : Ar;

            if (use_extern) {
                QpProblem qpp;
                qpp.nv  = nv;             qpp.mc  = m;
                qpp.H   = &H;             qpp.Bd  = NULL;
                qpp.g   = &gzero[0];      qpp.A   = &A;
                qpp.lbA = &lbA[0];        qpp.ubA = &ubA[0];
                qpp.lbd = &lbd[0];        qpp.ubd = &ubd[0];

                QpSolution qs;
                string why;
                if (solve_qp_plugin(backend, qpp, tol, qp_iter_max, /*nonconvex=*/false, qs, why)
                    && qs.ok) {
                    for (int j = 0; j < n; j++) dsol[(size_t) j] = qs.d[(size_t) j];
                    solved = true;
                }
            }
            else {
                qpOASES::Options qpopts;
                qpopts.printLevel = qpOASES::PL_NONE;
                qpopts.setToReliable();
                qpopts.printLevel = qpOASES::PL_NONE;

                qpOASES::SymSparseMat Hq(nv, nv, H.ir(), H.jc(), H.val());
                qpOASES::SparseMatrix Aq(m, nv, A.ir(), A.jc(), A.val());
                Hq.createDiagInfo();

                qpOASES::QProblem qp(nv, m);
                qp.setOptions(qpopts);
                qpOASES::int_t nWSR = 5*(nv + m) + 100;
                const qpOASES::returnValue rv =
                    qp.init(&Hq, &gzero[0], &Aq, &lbd[0], &ubd[0], &lbA[0], &ubA[0], nWSR);
                if (rv == qpOASES::SUCCESSFUL_RETURN || rv == qpOASES::RET_MAX_NWSR_REACHED) {
                    qp.getPrimalSolution(&dsol[0]);
                    solved = true;
                }
            }

            // Betts, 2.8.2 step 2(a)ii: the least distance program failing is the
            // signal to change strategy, not to give up.
            if (!solved && primary) { primary = false; n_relaxed++; }
            else if (!solved)       break;
        }

        if (!solved) { status = 2; break; }

        for (int j = 0; j < n; j++) d(j) = dsol[(size_t) j];

        // ---- line search on the violation alone ------------------------------
        // Sufficient decrease in Mv, and nothing else consulted. The relaxation step
        // gets a quadratic interpolation as well -- Betts's "accurate line search" --
        // because it is the step being taken when the linearisation is poor and the
        // one whose best steplength is least likely to be a power of one half.
        double alpha = 1.0, mvt = mv0;
        bool   taken = false;
        const double kappa = 1.0e-4;

        for (int ls = 0; ls < 30; ls++) {
            xtrial = x + alpha*d;
            for (int j = 0; j < n; j++)
                xtrial(j) = min(max(xtrial(j), (*xlb)(j)), (*xub)(j));
            gg_num(xtrial, &gtrial, workspace);
            mvt = violation_merit(gtrial, gl, gu);
            if (std::isfinite(mvt) && mvt <= (1.0 - kappa*alpha)*mv0) { taken = true; break; }
            alpha *= 0.5;
        }

        if (taken && !primary && alpha < 1.0) {
            // One quadratic through Mv(0), Mv'(0) approximated by -2*mv0, and Mv(alpha).
            const double denom = 2.0*(mvt - mv0 + 2.0*mv0*alpha);
            if (fabs(denom) > 0.0) {
                const double a_q = 2.0*mv0*alpha*alpha/denom;
                if (a_q > 0.1*alpha && a_q < 1.0) {
                    MatrixXd xq(n,1), gq(m,1);
                    xq = x + a_q*d;
                    for (int j = 0; j < n; j++)
                        xq(j) = min(max(xq(j), (*xlb)(j)), (*xub)(j));
                    gg_num(xq, &gq, workspace);
                    const double mq = violation_merit(gq, gl, gu);
                    if (std::isfinite(mq) && mq < mvt) { alpha = a_q; xtrial = xq; mvt = mq; }
                }
            }
        }

        if (!taken) { status = 3; break; }

        x = xtrial;
        n_iters++;
        alpha_prev = alpha;

        // Betts, 2.8.2 step 4(b): a full relaxation step is approximately a Newton
        // step, so it is worth trying the least distance program again.
        if (!primary && alpha >= 1.0) primary = true;
    }

    return status;
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
    vector<double> r(max(m,1), 0.0);                  // penalty weights, Betts's Theta
    MatrixXd smerit = MatrixXd::Zero(max(m,1),1);     // the merit function's slacks

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
    const bool bound_rho_by_multipliers = (algorithm.elastic_penalty == "multipliers");

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
    // The initial radius. It was 1.0, on the reasoning that the variables are scaled
    // and so an O(1) region is the natural first guess. Measured rather than reasoned,
    // 0.3 is better, and the reason is visible in what the QP backend is being asked to
    // do: on examples/zpm, counting GALAHAD's return codes over the first mesh, a radius
    // of 1.0 has it solving 44 subproblems, declaring 38 primal infeasible and reaching
    // its iteration limit on 22 -- a 58 per cent failure rate -- while a radius of 0.1
    // has it solving 386, declaring 19 infeasible and 6 limited, which is 6 per cent,
    // and gets through four times as many subproblems in the same wall clock. A region
    // the linearisation cannot be satisfied inside is a subproblem the backend cannot
    // solve, and a subproblem it cannot solve is a restoration the solver did not need.
    //
    // 0.3 rather than 0.1 because the small dense examples want the larger of the two:
    // at 0.1, brac1 takes 23 iterations against 17 and lts 28 against 20. At 0.3 those
    // are 17 and 18, and what the reduction buys elsewhere is kept -- launch goes from
    // 42 iterations and 195 seconds to 26 and 126, and manutec under qpOASES from 103
    // iterations and 106 seconds to 13 and 15.
    const double Delta_0   = 0.3;
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

    // ---- the feasibility phase, if the strategy asks for one --------------------
    // Betts's FM: locate a point feasible with respect to the constraints before any
    // attempt is made to optimise, so that the optimality phase begins somewhere its
    // model is worth trusting. Everything computed above is recomputed afterwards,
    // because the point has moved.
    int n_feas_iters = 0, n_feas_relaxed = 0, feas_status = 0;
    if (m > 0 && (algorithm.sqp_strategy == "FM" || algorithm.sqp_strategy == "F")) {

        if (iprint) {
            snprintf(workspace->text, sizeof(workspace->text),
                     "\n Locating a feasible point before optimising (strategy %s)\n"
                     "\n iter          method     max viol       step\n",
                     algorithm.sqp_strategy.c_str());
            psopt_print(workspace, workspace->text);
        }

        feas_status = feasible_point_phase(x, m, gl, gu, xlb, xub, tape_done,
                                           jrow, jcol, jval, Jm, gval,
                                           algorithm.qp_solver, use_extern,
                                           tol, algorithm.nlp_iter_max,
                                           algorithm.qp_iter_max, iprint,
                                           n_feas_iters, n_feas_relaxed, workspace);

        if (iprint) {
            static const char* why[4] = { "a feasible point was found",
                                          "the iteration limit was reached",
                                          "neither subproblem could be solved",
                                          "no step reduced the violation" };
            snprintf(workspace->text, sizeof(workspace->text),
                     "\n Feasibility phase finished after %d iterations: %s\n"
                     "   relaxation steps  %d\n\n",
                     n_feas_iters, why[feas_status <= 3 ? feas_status : 1], n_feas_relaxed);
            psopt_print(workspace, workspace->text);
        }

        // The point has moved, so everything evaluated at the old one is stale.
        fval = ff_num(x, workspace);
        gg_num(x, &gval, workspace);
        sqp_gradient(x, gf, workspace);
        sqp_jacobian_triplet(x, m, tape_done, jrow, jcol, jval, workspace);
        Jm.scatter(jval);

        if (algorithm.sqp_strategy == "F") {
            // Strategy F stops here. Reporting it as optimal would be a lie; reporting
            // it as a failure would be one too, when finding this point is what was
            // asked for.
            *x0 = x;
            if (lambda != NULL && m > 0) *lambda = MatrixXd::Zero(m,1);
            solution->nlp_return_code = (feas_status == 0) ? 0 : 1;
            if (iprint) {
                snprintf(workspace->text, sizeof(workspace->text),
                         " Strategy F: stopping at the feasible point, as requested.\n");
                psopt_print(workspace, workspace->text);
            }
            return (feas_status == 0) ? 0 : 1;
        }
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
    // Betts's Levenberg parameter: H = H_L + tau*(|sigma| + 1) I, tau in [0,1], sigma the
    // Gerschgorin bound on the most negative eigenvalue. It starts at zero -- the exact
    // Hessian -- because the inertia of the KKT matrix says at once when that is not
    // usable, which is the signal the previous two attempts at this lacked.
    double tau          = 0.0;
    int    n_inertia    = 0;          // KKT factorisations spent on the inertia test

    // The first subproblem is solved with H = I for its multipliers alone: with the
    // multipliers at zero the constraints contribute nothing to the Hessian of the
    // Lagrangian, and on a minimum-time problem the model is then identically zero.
    bool multiplier_pass = exact_hessian;
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

        // Betts drives the Levenberg parameter from the agreement between the predicted
        // and the actual reduction, in the manner of a trust region, and an earlier
        // version of this file did the same. It was dead code: the parameter it wrote
        // was shadowed by a declaration of the same name inside the subproblem loop
        // below, so the value computed here was never read and the shift began at
        // nothing at every iteration. Removing the shadow rather than the rule turned
        // out to matter more than the rule -- see the inertia loop below -- and the rule
        // itself is not reinstated: the inertia of the KKT matrix answers the same
        // question directly, and answers it before a step has to be taken and judged.

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
        for (int attempt = 0; ; attempt++) {

            if (exact_hessian && multiplier_pass) {
                Hm.scatter(vector<double>(hval.size(), 0.0));
                Hm.shift_diagonal(1.0);                   // H = I, for multipliers only
            }
            else if (exact_hessian) {
                // Raise the Levenberg parameter until the KKT matrix has the inertia
                // that says the reduced Hessian is positive definite. This is Betts's
                // inertia control, with the inertia obtained from MUMPS.
                //
                // The size of the shift matters as much as the fact of it, and how it is
                // arrived at matters more than either. Betts's parameter multiplies
                // |sigma| + 1, where sigma is the Gerschgorin bound on the most negative
                // eigenvalue of the Hessian, so a shift is always a fraction of a
                // worst-case bound rather than of anything the model actually needs; and
                // sigma moves with the multipliers, which on examples/bryson_denham took
                // it from -8.8e+02 to -5.4e+08 within three iterations. Started afresh
                // from a fixed fraction of that at every iteration, the shift is either
                // far too large or far too small, and neither is recoverable within the
                // iteration: too large flattens the curvature the model does have along
                // with the direction that has none, and the Newton step becomes a scaled
                // gradient step -- 150 iterations on bryson_denham, converging linearly,
                // against 19 or fewer on every other example; too small leaves the model
                // near-singular, and on brac1 the first full step took the scaled
                // objective of a minimum-time problem to zero.
                //
                // What removes the choice is not the starting value but the memory. The
                // shift that worked at the last iteration is a far better guess than any
                // fraction of a bound, so it is carried forward and tried at a fifth of
                // its size: if the model has improved the shift decays towards zero and
                // the exact Hessian returns, and if it has not the escalation begins
                // from somewhere useful instead of from nothing. This is IPOPT's inertia
                // correction, with Betts's bound setting the scale of the first attempt
                // and the ceiling.
                //
                // Sweeping the starting fraction over 1.0e-10 to 1.0e-4 and the
                // escalation over 10 and 100, with and without the carry-over, the
                // carry-over decides every case: with it all eight combinations solve
                // all five small examples, without it six of the eight fail or stall on
                // one or another. The decay factor is a real choice rather than a free
                // one: at 0.1 the shift decays faster than the model earns its curvature
                // back and bryson_denham fails, while 0.3333 to 0.6 all solve everything
                // under qpOASES. Measured across three backends rather than one, 0.2 is
                // better than 0.5 nearly everywhere -- brac1 in 14, 14 and 20 iterations
                // against 27, 58 and 40, manutec in 103, 12 and 12 against 193, 12 and
                // 12 -- and the one exception, bryson_denham under qpOASES, is 91
                // iterations against 19. Counts move by a factor of two on changes of a
                // few per cent, so they are not a quantity to tune against; the ends of
                // the working range are.
                Hm.scatter(hval);
                const double sigma     = Hm.gerschgorin_lower_bound();
                const double delta_max = fabs(sigma) + 1.0;
                // The shift that worked last time is the first thing to try, but it
                // has to be a shift for *this* matrix. delta_max moves with the
                // multipliers, and on examples/zpm it went from 1.5e+01 to 2.8e+50 and
                // back within three iterations when a subproblem returned multipliers
                // of 1e+50. Carried forward uncapped, a shift of 5.6e+49 then sat on a
                // matrix whose own bound was 15, decaying by a fifth each iteration and
                // needing seventy of them to come back -- seventy iterations of a model
                // that is a multiple of the identity, a step of zero, an objective
                // static to nine figures and a dual error that does not move. That is
                // the plateau seen on zpm and on low_thrust. Capping the carried value
                // at the current ceiling costs nothing when the scale is steady and
                // ends the plateau when it is not.
                double delta = (tau > 0.0) ? 0.2*tau : 0.0;
                for (;;) {
                    Hm.scatter(hval);
                    if (delta > 0.0) Hm.shift_diagonal(delta);

                    int npos = 0, nneg = 0, nzero = 0;
                    const bool got = kkt_inertia(Hm, Jm, n, m, npos, nneg, nzero);
                    n_inertia++;
                    // Gould's condition is usually written In(K) = (n, m, 0), and that
                    // is what was tested here. It is the right condition only when the
                    // Jacobian has full row rank. In general, with r = rank(J), the
                    // inertia of K is (r, r, m-r) plus that of the reduced Hessian on
                    // the null space of J, so a rank deficiency of m - r shows up as
                    // m - r zero eigenvalues and m - r fewer negative ones however the
                    // Hessian is shifted -- the deficiency is in J, and no change to H
                    // reaches it. Demanding zero zeros therefore asks for something a
                    // shift cannot deliver, and the loop answers by escalating to its
                    // ceiling, at every iteration, for ever.
                    //
                    // That is what examples/glider had been doing since inertia control
                    // went in: a Jacobian one row short of full rank, a shift pinned at
                    // |sigma| + 1, a model that was therefore a multiple of the identity,
                    // and a steepest-descent step of fixed length taken 1000 times per
                    // mesh. It solved in 117 iterations before, and stopped solving at
                    // all.
                    //
                    // The rank-tolerant form of the same condition is npos == n: the
                    // reduced Hessian is positive definite exactly when the positive
                    // eigenvalues number n, whatever the rank of J. It reduces to
                    // In(K) = (n, m, 0) when J has full rank, so nothing is given up.
                    if (!got || npos == n) break;
                    if (delta >= delta_max) break;
                    delta = (delta == 0.0) ? 1.0e-10*delta_max
                                           : min(100.0*delta, delta_max);
                    n_shifts++;
                }
                tau = delta;                     // carried into the next iteration
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
                if (!solve_qp_plugin(algorithm.qp_solver, qpp, tol, algorithm.qp_iter_max,
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
                    if (!all_finite(&dsol[0], n) || !all_finite(&ysol[0], n + m))
                        rv = RET_INIT_FAILED;
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
                    if (!all_finite(&dsol[0], n)) rv = RET_INIT_FAILED;
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

            // The inertia was corrected before the subproblem saw the model, so a
            // refusal now is about the constraints, and restoration is the answer.
            (void) curvature_failure;
            break;
        }

        // ---- restoration ---------------------------------------------------
        // The linearised constraints can be inconsistent even when the problem is
        // perfectly well posed -- it is the normal situation when the starting guess
        // violates an equality, which is how most people supply one. The subproblem
        // then has no solution and the iteration would simply stop.
        //
        // Two relaxations are offered, and psopt.h records what is known about when
        // each is the right one.
        //
        // Elastic mode, which is the default, gives every constraint its own pair of
        // non-negative slacks and charges them in the objective (Gill, Murray and
        // Saunders):
        //
        //   min 1/2 d'Bd + grad_f'd + rho*sum(v + w)
        //   s.t. lbA <= J d + v - w <= ubA,  lb <= d <= ub,  v, w >= 0.
        //
        // Betts's relaxation uses one variable instead of 2m (3rd ed., section 2.7).
        // For each row let a_i be how far d = 0 falls short of its lower bound and b_i
        // how far it exceeds its upper -- at most one can be positive, since the two
        // together would mean an empty interval -- and put
        //
        //   min 1/2 d'Bd + grad_f'd + rho*xi
        //   s.t. lbA <= J d + c xi <= ubA,  lb <= d <= ub,  0 <= xi <= 1,
        //
        // with c_i = a_i where the lower bound is the one violated and -b_i where it is
        // the upper. At xi = 1 the point d = 0 satisfies every row by construction, so
        // the subproblem is feasible for the same reason elastic mode is; at xi = 0 it
        // is the original subproblem exactly. Minimising rho*xi asks for the smallest
        // fraction of the current infeasibility that has to be tolerated for the
        // linearisation to be consistent.
        //
        // Both are always feasible, so a step always exists; both give the slacks a
        // small quadratic term as well, to keep the subproblem's Hessian positive
        // definite rather than merely semidefinite.
        if (m > 0 && rv != SUCCESSFUL_RETURN && rv != RET_MAX_NWSR_REACHED) {

            const bool relax = (algorithm.qp_restoration == "relaxation");

            // The price of a unit of infeasibility. See psopt.h for why the default is
            // taken from the merit weights alone although the multipliers are the
            // quantity the theory names, and what has to move with it when they are.
            double rho = 10.0;
            for (int i = 0; i < m; i++) rho = max(rho, 10.0*r[i]);
            if (bound_rho_by_multipliers) {
                for (int i = 0; i < m; i++) rho = max(rho, 10.0*fabs(lam(i)));
                for (int j = 0; j < n; j++) rho = max(rho, 10.0*fabs(zbnd(j)));
            }
            rho_elastic = rho;

            const int ne = relax ? (n + 1) : (n + 2*m);
            vector<double> ge(ne), lbe(ne), ube(ne), ze(ne), ye(ne + m);

            // Small relative to what the slacks cost, not small in absolute terms: once
            // rho is bounded below by the multipliers it is no longer near ten, and a
            // linear cost of 10^5 against a fixed quadratic term of 10^-8 is thirteen
            // orders of magnitude in one subproblem. qpOASES declined the relaxation on
            // the first mesh of examples/lts for that reason and no other.
            const double eps_slack = bound_rho_by_multipliers ? 1.0e-8*rho : 1.0e-8;

            // [ J | I | -I ] or [ J | c ], and blkdiag(H, eps I), built on the pattern
            // of the matrices they extend so that the relaxed subproblem costs the
            // extra columns and nothing else. The patterns are built the first time
            // restoration is needed and kept thereafter; the relaxation column carries
            // an entry in every row whether that row needs one at this iterate or not,
            // because the rows that need relaxing change from step to step and the
            // pattern does not.
            if (Ae.empty()) {
                vector<unsigned int> ar, ac;
                const size_t extra = relax ? (size_t) m : 2*(size_t) m;
                ar.reserve(jrow.size() + extra); ac.reserve(jcol.size() + extra);
                for (size_t k = 0; k < jrow.size(); k++) { ar.push_back(jrow[k]); ac.push_back(jcol[k]); }
                if (relax) {
                    for (int i = 0; i < m; i++) { ar.push_back((unsigned int) i); ac.push_back((unsigned int) n); }
                }
                else {
                    for (int i = 0; i < m; i++) { ar.push_back((unsigned int) i); ac.push_back((unsigned int)(n+i)); }
                    for (int i = 0; i < m; i++) { ar.push_back((unsigned int) i); ac.push_back((unsigned int)(n+m+i)); }
                }
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
            if (relax) {
                for (int i = 0; i < m; i++) {
                    const double lo = lbA[i], hi = ubA[i];
                    double c = 0.0;
                    if      (lo > 0.0 && lo <  qpOASES::INFTY) c = lo;   // short of the lower bound
                    else if (hi < 0.0 && hi > -qpOASES::INFTY) c = hi;   // past the upper bound
                    aeval.push_back(c);
                }
            }
            else {
                for (int i = 0; i < m; i++) aeval.push_back( 1.0);
                for (int i = 0; i < m; i++) aeval.push_back(-1.0);
            }
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

            // The subproblem is stated with its objective divided by rho, so that a
            // slack costs one apiece instead of rho and the quadratic block is scaled to
            // match. The minimiser is unchanged -- scaling an objective does not move
            // it -- but the conditioning is transformed: priced directly, a slack with a
            // linear cost of 10^7 against a regularising quadratic term of 10^-8 spans
            // fifteen orders of magnitude, and a proximal method spent twenty-three
            // thousand iterations on examples/brac1 failing to get through it.
            // ...and it is applied only for the plugin backends. qpOASES's active-set
            // method is not troubled by the unscaled pricing and is measurably worse
            // with it: examples/bryson_denham goes from twenty iterations to failing
            // after two hundred and forty-six.
            const double oscale = use_extern ? 1.0/rho : 1.0;
            for (int j = 0; j < n; j++) { ge[j] = oscale*gf(j);   lbe[j] = lbd[j]; ube[j] = ubd[j]; }
            for (int k = n; k < ne; k++) { ge[k] = oscale*rho;    lbe[k] = 0.0;
                                          ube[k] = relax ? 1.0 : qpOASES::INFTY; }
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
                (void) solve_qp_plugin(algorithm.qp_solver, qpp, tol, algorithm.qp_iter_max,
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
                    if (!all_finite(&ze[0], ne) || !all_finite(&ye[0], ne + m))
                        rve = RET_INIT_FAILED;
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

        if (multiplier_pass) {
            lam             = lam_new;
            multiplier_pass = false;
            iter--;
            continue;
        }

        // ---- the penalty weights -------------------------------------------------
        // The slacks first, from the weights the previous iteration settled on, then the
        // weights themselves. Betts's condition is that the directional derivative of
        // the merit function be at least as negative as -1/2 p'Hp, which is what makes
        // the step a descent direction; one condition does not determine m weights, so
        // the smallest set that satisfies it is taken, in the sense of least norm.
        merit_slacks(gval, lam, r, gl, gu, smerit);

        MatrixXd dlam(max(m,1),1), Jd(max(m,1),1), ds(max(m,1),1), rdev(max(m,1),1);
        if (m > 0) {
            Jm.times(&d(0), &Jd(0));
            for (int i = 0; i < m; i++) {
                dlam(i) = lam_new(i) - lam(i);
                rdev(i) = gval(i) - smerit(i);
                ds(i)   = Jd(i) + rdev(i);          // the predicted slack step
            }
        }

        {
            double gTd0 = 0.0;
            for (int j = 0; j < n; j++) gTd0 += gf(j)*d(j);
            const double dHd = quadratic_form(exact_hessian, Hm, B, d);

            // M'(0) = g'p + sum[ (dlam_i - lam_i) rdev_i - theta_i rdev_i^2 ], so the
            // condition M'(0) <= -1/2 p'Hp reads sum theta_i rdev_i^2 >= varsigma.
            double varsigma = gTd0 + 0.5*max(0.0, dHd);
            for (int i = 0; i < m; i++) varsigma += (dlam(i) - lam(i))*rdev(i);

            double aTa = 0.0;
            for (int i = 0; i < m; i++) aTa += rdev(i)*rdev(i)*rdev(i)*rdev(i);

            const double psi0 = PSOPT_extras::GetEPS();
            for (int i = 0; i < m; i++) r[i] = psi0;
            if (varsigma > 0.0 && aTa > 0.0)
                for (int i = 0; i < m; i++)
                    r[i] = psi0 + rdev(i)*rdev(i)*varsigma/aTa;

            // The l1 merit function is exact only while its weights dominate the
            // multipliers, so after a restoration they were raised to the price of the
            // relaxation. Betts's augmented Lagrangian is not a penalty function in
            // that sense -- it carries the multipliers explicitly and its weights are
            // the least-norm set that gives descent -- so with the price bounded below
            // by the multipliers this rule does not make the merit function exact, it
            // makes it stiff: weights of 10^5 on every constraint of
            // examples/bryson_denham, and a line search that cannot find a step. The
            // two rules were reading the same quantity for opposite purposes, and only
            // worked together while it was the l1 weights, which were the right size
            // for neither.
            if (elastic && !bound_rho_by_multipliers)
                for (int i = 0; i < m; i++) r[i] = max(r[i], rho_elastic);

            // The slacks were computed against the previous weights; recompute them
            // against these, so that the function the line search sees is the one the
            // weights were chosen for.
            merit_slacks(gval, lam, r, gl, gu, smerit);
            for (int i = 0; i < m; i++) {
                rdev(i) = gval(i) - smerit(i);
                ds(i)   = Jd(i) + rdev(i);
            }
        }

        // ---- line search ---------------------------------------------------------
        // The step moves the multipliers and the slacks alongside x, which is what makes
        // the merit function's value at alpha = 1 the value at the point the subproblem
        // actually proposed. Betts, (2.30).
        const double phi0 = merit_al(fval, gval, lam, smerit, r);

        double gTd = 0.0;
        for (int j = 0; j < n; j++) gTd += gf(j)*d(j);
        double dphi = gTd;
        for (int i = 0; i < m; i++)
            dphi += (dlam(i) - lam(i))*rdev(i) - r[i]*rdev(i)*rdev(i);

        // Betts fits a quadratic and a cubic to the merit function rather than halving,
        // and imposes the Wolfe condition to stop steplengths becoming too small
        // (section 2.6.1). Halving is the crudest thing that works, and the traces show
        // what it costs: steplengths of 1.95e-03 and 1.19e-07 on the harder examples,
        // each one reached by throwing away a full evaluation of the objective and every
        // constraint at every power of a half on the way down.
        //
        // What is fitted here is the standard safeguarded interpolation (Nocedal and
        // Wright, section 3.5): a quadratic through phi(0), phi'(0) and the first
        // rejected point, a cubic through those and the second, each new trial confined
        // to [0.1, 0.5] of the last so that the search cannot stall or leap. Betts's
        // Wolfe condition needs phi'(alpha) at the trial point, which for this merit
        // function means the gradient and the whole Jacobian there -- one full
        // derivative evaluation per trial, which is not worth it when the purpose the
        // condition serves is to stop the steplength collapsing and a floor does that
        // for nothing. The floor is alpha_min.
        double alpha = 1.0;
        const double eta       = 1.0e-4;      // Betts's kappa_1
        const double alpha_min = 1.0e-8;
        MatrixXd xtrial(n,1), gtrial(max(m,1),1);
        MatrixXd lam_a(max(m,1),1), s_a(max(m,1),1);
        double ftrial = fval, phit = phi0;
        bool   accepted = false;

        double a_prev = 0.0, phi_prev = 0.0;   // the last rejected trial, for the cubic

        for (int ls = 0; ls < 25; ls++) {
            xtrial = x + alpha*d;
            for (int j = 0; j < n; j++)
                xtrial(j) = min(max(xtrial(j), (*xlb)(j)), (*xub)(j));
            ftrial = ff_num(xtrial, workspace);
            if (m > 0) gg_num(xtrial, &gtrial, workspace);
            for (int i = 0; i < m; i++) {
                lam_a(i) = lam(i) + alpha*dlam(i);
                s_a(i)   = smerit(i) + alpha*ds(i);
            }
            phit = merit_al(ftrial, gtrial, lam_a, s_a, r);
            if (std::isfinite(phit) && phit <= phi0 + eta*alpha*dphi) { accepted = true; break; }
            if (alpha <= alpha_min) break;

            // The next trial, by interpolation where the numbers allow it and by
            // halving where they do not. dphi is the slope at zero and is negative by
            // construction of the penalty weights, so a finite minimiser exists unless
            // the model is degenerate.
            double a_new = 0.5*alpha;
            if (std::isfinite(phit) && dphi < 0.0) {
                if (a_prev == 0.0) {
                    // quadratic through phi(0), phi'(0), phi(alpha)
                    const double den = 2.0*(phit - phi0 - dphi*alpha);
                    if (den > 0.0) a_new = -dphi*alpha*alpha/den;
                }
                else {
                    // cubic through phi(0), phi'(0), phi(alpha), phi(a_prev)
                    const double d1 = alpha - a_prev;
                    if (fabs(d1) > 0.0 && std::isfinite(phi_prev)) {
                        const double u = phit     - phi0 - dphi*alpha;
                        const double v = phi_prev - phi0 - dphi*a_prev;
                        const double aa = ( u/(alpha*alpha) - v/(a_prev*a_prev))/d1;
                        const double bb = (-u*a_prev/(alpha*alpha)
                                           + v*alpha/(a_prev*a_prev))/d1;
                        const double disc = bb*bb - 3.0*aa*dphi;
                        if (disc >= 0.0) {
                            if (fabs(aa) < 1.0e-300) { if (bb > 0.0) a_new = -dphi/(2.0*bb); }
                            else                     a_new = (-bb + sqrt(disc))/(3.0*aa);
                        }
                    }
                }
            }
            if (!std::isfinite(a_new)) a_new = 0.5*alpha;
            // Nocedal and Wright's band. It is not a free parameter: at [0.25, 0.5]
            // the search cuts too gently and bryson_denham stops converging, while
            // examples/launch, which wants long steps, takes 30 iterations instead of
            // 42. The wider band is kept because what it buys on bryson_denham is an
            // answer rather than a count.
            a_new = min(max(a_new, 0.1*alpha), 0.5*alpha);

            a_prev   = alpha;
            phi_prev = phit;
            alpha    = a_new;
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
                (void) solve_qp_plugin(algorithm.qp_solver, qpp, tol, algorithm.qp_iter_max,
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
                    if (!all_finite(&dsoc[0], n) || !all_finite(&ysoc[0], n + m))
                        rvs = RET_INIT_FAILED;
                }
            }

            if (rvs == SUCCESSFUL_RETURN || rvs == RET_MAX_NWSR_REACHED) {

                MatrixXd dc(n,1);
                for (int j = 0; j < n; j++) dc(j) = dsoc[j];

                // The corrected step's slope, in the same merit function. The slack
                // step follows the corrected direction, so the deviation term is the
                // same as before; only the objective part changes.
                double gTdc = 0.0;
                for (int j = 0; j < n; j++) gTdc += gf(j)*dc(j);
                double dphi_c = gTdc;
                for (int i = 0; i < m; i++)
                    dphi_c += (dlam(i) - lam(i))*rdev(i) - r[i]*rdev(i)*rdev(i);

                alpha = 1.0;
                for (int ls = 0; ls < 25; ls++) {
                    xtrial = x + alpha*dc;
                    for (int j = 0; j < n; j++)
                        xtrial(j) = min(max(xtrial(j), (*xlb)(j)), (*xub)(j));
                    ftrial = ff_num(xtrial, workspace);
                    gg_num(xtrial, &gtrial, workspace);
                    for (int i = 0; i < m; i++) {
                        lam_a(i) = lam(i) + alpha*dlam(i);
                        s_a(i)   = smerit(i) + alpha*ds(i);
                    }
                    phit = merit_al(ftrial, gtrial, lam_a, s_a, r);
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
                     (exact_hessian && tau > 0.0) ? "shifted" : "");
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
