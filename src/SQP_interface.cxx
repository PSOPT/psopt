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

using namespace std;

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

} // anonymous namespace


// ---------------------------------------------------------------------------------
// The constraint Jacobian, dense.
//
// The values come from PSOPT's own machinery, so the SQP sees exactly the derivatives
// the rest of the library sees: a sparse automatic-differentiation evaluation when
// algorithm.derivatives is "automatic", and one-sided finite differences otherwise.
// Only the storage is dense -- the sparse stage keeps the evaluation and drops the
// scatter.
// ---------------------------------------------------------------------------------
static void sqp_jacobian(MatrixXd& x, MatrixXd& J, Workspace* workspace, bool& tape_done)
{
    const int n = (int) x.rows();
    const int m = (int) J.rows();

    J.setZero();

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
        for (int k = 0; k < T.nnz(); k++)
            J((int) T.row[k], (int) T.col[k]) = T.val[k];
    }
    else {
        // One-sided differences, one column per variable. This is the expensive path
        // and it is why the dense stage is a stage: it costs n constraint evaluations
        // per SQP iteration.
        MatrixXd g0(m,1), g1(m,1);
        gg_num(x, &g0, workspace);
        const double sqreps = sqrt(PSOPT_extras::GetEPS());
        for (int j = 0; j < n; j++) {
            const double xj   = x(j);
            const double delj = sqreps*(1.0 + fabs(xj));
            x(j) = xj + delj;
            gg_num(x, &g1, workspace);
            x(j) = xj;
            for (int i = 0; i < m; i++) J(i,j) = (g1(i) - g0(i))/delj;
        }
    }
}


static void sqp_gradient(MatrixXd& x, MatrixXd& gf, Workspace* workspace)
{
    if (useAutomaticDifferentiation(*workspace->algorithm))
        ScalarGradientAD(ff_ad, x, &gf, &workspace->trace_f_done, workspace->ad_f, workspace);
    else
        ScalarGradient(ff_num, x, &gf, workspace->grw.get(), workspace);
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

    MatrixXd gf(n,1), gval(max(m,1),1), J(max(m,1), n);
    MatrixXd gf_old(n,1), J_old(max(m,1), n);
    MatrixXd lam  = MatrixXd::Zero(max(m,1),1);       // constraint multipliers
    MatrixXd zbnd = MatrixXd::Zero(n,1);              // bound multipliers
    MatrixXd B    = MatrixXd::Identity(n,n);          // the quasi-Newton model
    vector<double> r(max(m,1), 0.0);                  // l1 penalty weights

    bool tape_done = false;

    const double tol      = algorithm.nlp_tolerance;
    const int    iter_max = algorithm.nlp_iter_max;

    // qpOASES wants row-major storage; B is symmetric so its layout is immaterial,
    // but the Jacobian's is not.
    typedef Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> RowMajorXd;
    RowMajorXd Arow(max(m,1), n);

    vector<double> lbA(max(m,1)), ubA(max(m,1)), lbd(n), ubd(n);
    vector<double> dsol(n), ysol(n + max(m,1));

    double fval = ff_num(x, workspace);
    if (m > 0) gg_num(x, &gval, workspace);
    sqp_gradient(x, gf, workspace);
    if (m > 0) sqp_jacobian(x, J, workspace, tape_done);

    int    n_restorations = 0;
    int    n_corrections  = 0;
    int    status  = 1;                                // 1 = iteration limit
    string message = "Maximum number of SQP iterations reached";
    int    iter    = 0;

    if (iprint) {
        snprintf(workspace->text, sizeof(workspace->text),
                 "\n\nSQP (dense, BFGS + qpOASES): %d variables, %d constraints\n", n, m);
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
        if (m > 0) dL += J.transpose()*lam;
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

        // ---- the quadratic programming subproblem -------------------------------
        //   min  1/2 d' B d + grad_f' d
        //   s.t. g_l - g <= J d <= g_u - g,      xlb - x <= d <= xub - x
        // PSOPT writes an absent constraint bound as +/- 1.0e20 (see NLP_bounds.cxx);
        // qpOASES has its own infinity and treats anything beyond it as free.
        const double psopt_inf = 1.0e20;
        for (int i = 0; i < m; i++) {
            lbA[i] = (gl[i] <= -psopt_inf) ? -qpOASES::INFTY : gl[i] - gval(i);
            ubA[i] = (gu[i] >=  psopt_inf) ?  qpOASES::INFTY : gu[i] - gval(i);
            for (int j = 0; j < n; j++) Arow(i,j) = J(i,j);
        }
        for (int j = 0; j < n; j++) {
            lbd[j] = ((*xlb)(j) <= -psopt_inf) ? -qpOASES::INFTY : (*xlb)(j) - x(j);
            ubd[j] = ((*xub)(j) >=  psopt_inf) ?  qpOASES::INFTY : (*xub)(j) - x(j);
        }

        Options qpopts;
        qpopts.printLevel = PL_NONE;
        qpopts.setToReliable();
        qpopts.printLevel = PL_NONE;

        int_t nWSR = 5*(n + m) + 100;
        returnValue rv;
        bool   elastic = false;
        double rho_elastic = 0.0;

        if (m > 0) {
            QProblem qp(n, m);
            qp.setOptions(qpopts);
            rv = qp.init(&B(0,0), &gf(0), &Arow(0,0), &lbd[0], &ubd[0], &lbA[0], &ubA[0], nWSR);
            if (rv == SUCCESSFUL_RETURN || rv == RET_MAX_NWSR_REACHED) {
                qp.getPrimalSolution(&dsol[0]);
                qp.getDualSolution(&ysol[0]);
            }
        }
        else {
            QProblemB qp(n);
            qp.setOptions(qpopts);
            rv = qp.init(&B(0,0), &gf(0), &lbd[0], &ubd[0], nWSR);
            if (rv == SUCCESSFUL_RETURN || rv == RET_MAX_NWSR_REACHED) {
                qp.getPrimalSolution(&dsol[0]);
                qp.getDualSolution(&ysol[0]);
            }
        }

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
            RowMajorXd Ae = RowMajorXd::Zero(m, ne);
            MatrixXd   He = MatrixXd::Zero(ne, ne);
            vector<double> ge(ne), lbe(ne), ube(ne), ze(ne), ye(ne + m);

            const double eps_slack = 1.0e-8;
            He.topLeftCorner(n,n) = B;
            for (int k = n; k < ne; k++) He(k,k) = eps_slack;

            for (int i = 0; i < m; i++) {
                for (int j = 0; j < n; j++) Ae(i,j) = J(i,j);
                Ae(i, n + i)     =  1.0;      // v
                Ae(i, n + m + i) = -1.0;      // w
            }
            for (int j = 0; j < n; j++) { ge[j] = gf(j); lbe[j] = lbd[j]; ube[j] = ubd[j]; }
            for (int k = n; k < ne; k++) { ge[k] = rho;  lbe[k] = 0.0;    ube[k] = qpOASES::INFTY; }

            QProblem qpe(ne, m);
            qpe.setOptions(qpopts);
            int_t nWSRe = 5*(ne + m) + 100;
            returnValue rve = qpe.init(&He(0,0), &ge[0], &Ae(0,0), &lbe[0], &ube[0],
                                       &lbA[0], &ubA[0], nWSRe);

            if (rve == SUCCESSFUL_RETURN || rve == RET_MAX_NWSR_REACHED) {
                qpe.getPrimalSolution(&ze[0]);
                qpe.getDualSolution(&ye[0]);
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
        for (int j = 0; j < n; j++) zbnd(j)    =  ysol[j];

        // ---- the l1 penalty weights ---------------------------------------------
        // Powell's rule: keep each weight at least as large as the magnitude of the
        // multiplier it accompanies, which is what makes d a descent direction for
        // the merit function, and let it decay towards that value when it may.
        for (int i = 0; i < m; i++) {
            const double li = fabs(lam_new(i));
            r[i] = (iter == 0) ? li : max(li, 0.5*(r[i] + li));
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

            MatrixXd Jd = J*d;
            vector<double> lbS(m), ubS(m);
            for (int i = 0; i < m; i++) {
                const double shift = gtrial_soc(i) - Jd(i);
                lbS[i] = (gl[i] <= -psopt_inf) ? -qpOASES::INFTY : gl[i] - shift;
                ubS[i] = (gu[i] >=  psopt_inf) ?  qpOASES::INFTY : gu[i] - shift;
            }

            QProblem qps(n, m);
            qps.setOptions(qpopts);
            int_t nWSRs = 5*(n + m) + 100;
            vector<double> dsoc(n), ysoc(n + m);
            returnValue rvs = qps.init(&B(0,0), &gf(0), &Arow(0,0), &lbd[0], &ubd[0],
                                       &lbS[0], &ubS[0], nWSRs);

            if (rvs == SUCCESSFUL_RETURN || rvs == RET_MAX_NWSR_REACHED) {
                qps.getPrimalSolution(&dsoc[0]);
                qps.getDualSolution(&ysoc[0]);

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
            // A merit function that cannot be decreased along the QP direction means
            // the quadratic model is not usable here. Reset it and try once more from
            // the same point before giving up.
            if (B.isApprox(MatrixXd::Identity(n,n))) {
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

        // ---- BFGS update, damped ------------------------------------------------
        // s is the accepted step and y the change in the gradient of the Lagrangian
        // at fixed multipliers. Powell's damping keeps s'y positive, and with it the
        // positive definiteness that the QP solver requires of B.
        gf_old = gf;
        if (m > 0) J_old = J;

        MatrixXd s = xtrial - x;

        x    = xtrial;
        fval = ftrial;
        if (m > 0) gval = gtrial;
        lam  = lam_new;

        sqp_gradient(x, gf, workspace);
        if (m > 0) sqp_jacobian(x, J, workspace, tape_done);

        MatrixXd y = gf - gf_old;
        if (m > 0) y += (J - J_old).transpose()*lam;

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

        if (iprint) {
            snprintf(workspace->text, sizeof(workspace->text),
                     "%5d %16.8e %12.3e %12.3e %10.2e %8d %s\n",
                     iter+1, fval, (m > 0) ? max_violation(gval, gl, gu) : 0.0,
                     dual_err, alpha, (int) nWSR, elastic ? "restoration" : "");
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
                 "   second-order corr.   %d\n",
                 iter, message.c_str(), fval,
                 (m > 0) ? max_violation(gval, gl, gu) : 0.0, n_restorations, n_corrections);
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
