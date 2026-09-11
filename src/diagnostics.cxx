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
e-mail:    vmbecerra@vmb1.com

**********************************************************************************************/

// diagnostics.cxx -- optional post-solve solution diagnostics for PSOPT.
//
// This module is entirely gated behind Alg::diagnostic_level (0 = off, the default, in which
// case solution_diagnostics() is never called and the solve is bit-identical to before). It
// reads only quantities that PSOPT has already computed for the converged solution (per-interval
// relative ODE error, states, costates and the stored Hamiltonian) and reports an interpretation;
// it never alters the solution.
//
//   level 1 : per-interval discretization-error localisation + a per-state spectral smoothness
//             indicator (this increment, D1); costate-structure summary (D2a).
//   level 2 : additionally the Hamiltonian-constancy and stationarity (dH/du) residuals (D2b),
//             and the rank and conditioning of the constraint Jacobian at the final iterate (D3).
//
// A note on the smoothness indicator: hp_refine.cxx contains its own legendre_decay_rate used to
// drive mesh refinement. We deliberately keep a private copy here (diag_legendre_decay_rate)
// rather than share a symbol, so that the refinement path is left completely untouched by this
// optional diagnostic feature. The two are similar in spirit but serve different purposes: the
// refinement version fits each interval exactly (low, well-conditioned degree), whereas the
// diagnostic version fits the whole phase in a least-squares sense with a capped degree, giving a
// single robust "is this state's trajectory globally smooth?" read per state.

#include "psopt.h"
#include <Eigen/SVD>
#include <cmath>
#include <vector>
#include <algorithm>
using namespace Eigen;
using namespace std;

// Exponential decay rate sigma of the Legendre-coefficient envelope of the samples (xi, vals),
// with xi in [-1,1]. The envelope behaves like |a_k| ~ exp(-sigma * k); a larger sigma means
// faster spectral decay, i.e. a smoother function. A least-squares fit to a capped Legendre degree
// keeps the estimate well-conditioned irrespective of how many nodes the phase has. Returns -1.0
// when there are too few points to judge, and a large value for an essentially constant trajectory.
static double diag_legendre_decay_rate(const VectorXd& xi, const VectorXd& vals, double sigma_min)
{
    const int npts = (int) xi.size();
    if (npts < 4) return -1.0;                       // too few points to estimate

    // Essentially constant trajectory -> treat as perfectly smooth.
    const double vspan = vals.maxCoeff() - vals.minCoeff();
    const double vscale = vals.cwiseAbs().maxCoeff() + 1.0;
    if (vspan < 1.0e-12 * vscale) return 50.0;

    const int m = std::min(npts - 1, 24);            // capped fit degree

    // Legendre Vandermonde V(npts, m+1) via the three-term recurrence.
    MatrixXd V(npts, m + 1);
    for (int r = 0; r < npts; r++) {
        const double x = xi(r);
        double Pkm1 = 1.0, Pk = x;
        V(r, 0) = 1.0;
        if (m >= 1) V(r, 1) = x;
        for (int k = 1; k < m; k++) {
            const double Pkp1 = ((2.0 * k + 1.0) * x * Pk - k * Pkm1) / (k + 1.0);
            V(r, k + 1) = Pkp1;
            Pkm1 = Pk;
            Pk = Pkp1;
        }
    }

    // Least-squares Legendre coefficients a_0..a_m.
    const VectorXd a = V.colPivHouseholderQr().solve(vals);

    // A well-resolved smooth state's coefficients decay until they reach the round-off plateau
    // (~1e-14 relative) and then flatten; fitting that flat noise tail would spuriously report a
    // slow decay. So fit only the genuine decay region: from k=1 up to the last coefficient that
    // sits above the round-off floor.
    const double amax = a.cwiseAbs().maxCoeff() + 1.0e-300;
    const double floor_level = 1.0e-11 * amax;
    int kmax = 0;
    for (int k = 0; k <= m; k++)
        if (std::fabs(a(k)) > floor_level) kmax = k;
    if (kmax < 2) return 50.0;          // spectrum collapses within <2 modes => very smooth

    // Slope of ln|a_k| over [1, kmax]; -slope is the exponential decay rate sigma. Exponential
    // (smooth) decay gives a steep negative slope (large sigma); algebraic (non-smooth) decay is
    // shallow on this linear-in-k axis (small sigma).
    int cnt = 0;
    double sx = 0.0, sy = 0.0, sxx = 0.0, sxy = 0.0;
    for (int k = 1; k <= kmax; k++) {
        const double y = std::log(std::max(std::fabs(a(k)), floor_level) / amax);
        sx += k; sy += y; sxx += (double) k * k; sxy += (double) k * y;
        cnt++;
    }
    if (cnt < 2) return -1.0;
    const double denom = cnt * sxx - sx * sx;
    if (std::fabs(denom) < 1.0e-300) return sigma_min;
    const double slope = (cnt * sxy - sx * sy) / denom;
    double sigma = -slope;
    if (!(sigma > sigma_min)) sigma = sigma_min;     // also catches NaN
    return sigma;
}

// ---------------------------------------------------------------------------------------
// D3: the rank and conditioning of the constraint Jacobian at the final iterate.
//
// An NLP whose constraint Jacobian is rank-deficient has no unique set of multipliers, so
// the dual residual grad f + J'lambda - z has no well-defined minimum and what a solver
// reports for it is decided by its regularisation and by rounding. The symptom is a run
// whose objective and constraint violation settle to machine precision while the dual
// error wanders over orders of magnitude and finishes wherever it happens to be -- and,
// because nothing else looks wrong, a user has no way to find the cause.
//
// examples/dae_i3 is the case this was written for. It is the pendulum in index-3 form:
// two positions, two velocities, the holonomic constraint L^2 - x1^2 - x2^2 = 0 imposed at
// every node, and its multiplier carried as a control. Collocated directly, its constraint
// Jacobian is square and rank-deficient by three at every mesh from 20 nodes to 80 -- the
// three smallest singular values at 30 nodes are 6.2e-15, 4.1e-15 and 1.8e-27 against a
// largest of 4.0e+02 -- which is the textbook consequence of collocating an index-3 DAE
// without reducing it: the constraint and the two hidden constraints obtained by
// differentiating it are not independent of the dynamics. The same example under the
// integrated-residual transcription has full column rank and a condition number of 9.1e+07.
//
// The check is a dense SVD and is therefore opt-in and size-limited. It is the right price
// only for someone who is already asking why a run behaves the way it does.
// ---------------------------------------------------------------------------------------
static void constraint_jacobian_conditioning(Prob& problem, Alg& algorithm,
                                             Workspace* workspace)
{
    char* t = workspace->text;
    const size_t tn = sizeof(workspace->text);

    const int nv = get_number_nlp_vars(problem, workspace);
    const int nc = get_number_nlp_constraints(problem, workspace);
    if (nv <= 0 || nc <= 0) return;

    psopt_print(workspace, "\n  Constraint Jacobian at the final iterate");

    // A dense SVD costs O(m n min(m,n)). Past a few million entries that is minutes rather
    // than seconds, and a diagnostic that hangs is worse than one that declines.
    const double entries = (double) nc * (double) nv;
    if (entries > 4.0e6) {
        snprintf(t, tn, "\n    %d rows by %d columns: too large for the dense factorisation this"
                        "\n    check uses, and it is skipped. Reduce the mesh to ask the question.",
                 nc, nv);
        psopt_print(workspace, t);
        return;
    }

    MatrixXd& x = *workspace->x0;          // the final iterate: both solvers write it back

    psopt_ad::ad_record(workspace->ad_gc, nv, nc, &x(0),
        [&](const adouble* xin, adouble* yout){ gg_ad(const_cast<adouble*>(xin), yout, workspace); });
    psopt_ad::SparseTriplet Jt = psopt_ad::ad_sparse_jacobian(workspace->ad_gc, &x(0), false);

    // Only the columns Ipopt keeps. A variable pinned by coincident bounds is removed from
    // the problem (fixed_variable_treatment defaults to make_parameter), so the rank that
    // decides whether the multipliers are unique is the rank over the rest.
    MatrixXd& blb = *(workspace->xlb);
    MatrixXd& bub = *(workspace->xub);
    std::vector<int> keep;
    for (int j = 0; j < nv; j++) if (blb(j) != bub(j)) keep.push_back(j);
    const int nf = (int) keep.size();
    if (nf <= 0) return;

    std::vector<int> col_of((size_t) nv, -1);
    for (int j = 0; j < nf; j++) col_of[(size_t) keep[(size_t) j]] = j;

    MatrixXd J = MatrixXd::Zero(nc, nf);
    for (int k = 0; k < Jt.nnz(); k++) {
        const int c = col_of[(size_t) Jt.col[k]];
        if (c >= 0) J(Jt.row[k], c) += Jt.val[k];
    }

    Eigen::BDCSVD<MatrixXd> svd(J);
    const MatrixXd sv = svd.singularValues();
    const int    ns   = (int) sv.rows();
    if (ns <= 0) return;

    const double smax = sv(0);
    const double smin = sv(ns - 1);
    // The standard numerical rank: singular values below smax * max(dim) * eps are zero.
    const double rtol = smax * (double) std::max(nc, nf) * 2.220446049250313e-16;
    int rank = 0;
    for (int k = 0; k < ns; k++) if (sv(k) > rtol) rank++;
    const int full = std::min(nc, nf);
    const int def  = full - rank;

    snprintf(t, tn, "\n    %d rows, %d free columns (%d of %d variables fixed by coincident bounds)",
             nc, nf, nv - nf, nv);
    psopt_print(workspace, t);
    if (def > 0) {
        // The gap is the evidence, so show it: the smallest singular value that counts as
        // rank against the largest that does not. A clean separation says the deficiency is
        // structural rather than an artefact of the tolerance.
        const double kept    = (rank > 0)  ? sv(rank - 1) : 0.0;
        const double dropped = (rank < ns) ? sv(rank)     : 0.0;
        snprintf(t, tn, "\n    rank %d of a possible %d;  sigma_max %10.3e"
                        "\n    RANK DEFICIENT BY %d, at a tolerance of %10.3e:"
                        "\n      smallest singular value counted as rank  %10.3e"
                        "\n      largest  singular value counted as zero  %10.3e",
                 rank, full, smax, def, rtol, kept, dropped);
        psopt_print(workspace, t);
        psopt_print(workspace,
            "\n    Read: the constraints are not independent, so the multipliers are not unique"
            "\n          and the dual residual has no well-defined minimum. Expect the objective"
            "\n          and the constraint violation to settle while the dual error wanders, and"
            "\n          expect the verdict to move under perturbations that do not change the"
            "\n          problem. The usual causes are a constraint stated twice, a boundary"
            "\n          condition already implied by the dynamics, and a differential-algebraic"
            "\n          system of index two or higher collocated without being index-reduced.");
    }
    else {
        snprintf(t, tn, "\n    rank %d of a possible %d, which is full;  sigma_max %10.3e,"
                        " sigma_min %10.3e, ratio %10.3e",
                 rank, full, smax, smin, (smin > 0.0) ? smax/smin : PSOPT::inf);
        psopt_print(workspace, t);
        if (smin > 0.0 && smax/smin > 1.0e12)
            psopt_print(workspace,
                "\n    Read: full rank, but the ratio says the constraints are close to being"
                "\n          dependent. The dual error may be limited by that rather than by"
                "\n          anything the solver is doing.");
    }
}


void solution_diagnostics(Prob& problem, Alg& algorithm, Sol& solution, Workspace* workspace)
{
    if (algorithm.diagnostic_level <= 0) return;

    const bool ps = use_global_collocation(algorithm);
    const double sigma_min = 0.05;
    const double sigma_lo  = 0.5;     // below this an interval/state is flagged non-smooth
    char* t = workspace->text;
    const size_t tn = sizeof(workspace->text);

    psopt_print(workspace, "\n*******************************************************************************");
    psopt_print(workspace, "\n                              Solution diagnostics                             ");
    psopt_print(workspace, "\n*******************************************************************************\n");

    for (int iphase = 1; iphase <= problem.nphases; iphase++) {
        MatrixXd& err    = solution.relative_errors[iphase - 1];   // 1 x norder
        MatrixXd& time   = solution.get_time_in_phase(iphase);     // 1 x (norder+1)
        MatrixXd& states = solution.get_states_in_phase(iphase);   // nstates x (norder+1)
        MatrixXd& sig    = solution.smoothness[iphase - 1];        // 1 x nstates
        const int nstates = problem.phase[iphase - 1].nstates;
        const int norder  = (int) err.cols();
        const int N       = (int) states.cols();

        snprintf(t, tn, "\nPhase %i", iphase);
        psopt_print(workspace, t);

        // ---- per-state spectral smoothness (global pseudospectral methods only) ----
        if (ps && N >= 4) {
            const double t0 = time(0), tf = time(N - 1);
            VectorXd xi(N);
            for (int k = 0; k < N; k++)
                xi(k) = (tf > t0) ? (2.0 * (time(k) - t0) / (tf - t0) - 1.0) : 0.0;
            for (int s = 0; s < nstates; s++) {
                VectorXd v = states.row(s).transpose();
                sig(0, s) = diag_legendre_decay_rate(xi, v, sigma_min);
            }
        } else {
            sig.setConstant(-1.0);     // not applicable
        }

        // ---- worst discretization-error intervals ----
        const int kshow = std::min(3, norder);
        if (norder > 0) {
            std::vector<std::pair<double,int> > ev(norder);
            for (int i = 0; i < norder; i++) ev[i] = std::make_pair(err(0, i), i);
            std::partial_sort(ev.begin(), ev.begin() + kshow, ev.end(),
                              [](const std::pair<double,int>& a, const std::pair<double,int>& b)
                              { return a.first > b.first; });
            snprintf(t, tn, "\n  Largest relative ODE error over %i intervals:", norder);
            psopt_print(workspace, t);
            for (int j = 0; j < kshow; j++) {
                const int i = ev[j].second;
                const double tm = 0.5 * (time(i) + time(i + 1));
                snprintf(t, tn, "\n    interval %4i   t=%12.5e   error=%12.5e", i + 1, tm, ev[j].first);
                psopt_print(workspace, t);
            }
        }

        // ---- per-state smoothness read ----
        if (ps) {
            psopt_print(workspace, "\n  State smoothness (Legendre-decay sigma; sigma < 0.50 => non-smooth):");
            for (int s = 0; s < nstates; s++) {
                const double sg = sig(0, s);
                const char* lab = (sg < 0.0) ? "undetermined"
                                : (sg < sigma_lo) ? "NON-SMOOTH" : "smooth";
                if (sg < 0.0)
                    snprintf(t, tn, "\n    state %3i:  sigma=  n/a        (%s)", s + 1, lab);
                else
                    snprintf(t, tn, "\n    state %3i:  sigma=%10.3e   (%s)", s + 1, sg, lab);
                psopt_print(workspace, t);
            }
            psopt_print(workspace,
                "\n  Read: large error with high sigma => smooth but under-resolved (raise the order);"
                "\n        large error with low sigma  => genuinely non-smooth (subdivide or regularise).");
        } else {
            psopt_print(workspace,
                "\n  (Spectral smoothness indicator is reported for the global pseudospectral methods only.)");
        }

        // ---- costate-structure summary (D2a) ----
        {
            MatrixXd& lam = solution.dual.costates[iphase - 1];   // nstates x ncols
            const int lc = (int) lam.cols();
            if (lc > 0 && (int) lam.rows() >= nstates) {
                const double tau    = 1.0e-3;     // near-null threshold, relative to the dominant peak
                const double afloor = 1.0e-12;    // absolute floor (sign tolerance / "all ~ 0" guard)
                double maxpeak = 0.0;
                for (int j = 0; j < nstates; j++) {
                    const double pk = lam.row(j).cwiseAbs().maxCoeff();
                    if (pk > maxpeak) maxpeak = pk;
                }
                psopt_print(workspace, "\n  Costate structure (lambda per state):");
                for (int j = 0; j < nstates; j++) {
                    const VectorXd r = lam.row(j).transpose();
                    const double pk   = r.cwiseAbs().maxCoeff();
                    const double rms  = std::sqrt(r.squaredNorm() / (double) lc);
                    const double rmin = r.minCoeff(), rmax = r.maxCoeff();
                    const char* sgn = (rmin >= -afloor) ? "+  "
                                    : (rmax <=  afloor) ? "-  " : "+/-";
                    const bool nearnull = (maxpeak > afloor) && (pk < tau * maxpeak);
                    snprintf(t, tn, "\n    state %3i:  peak=%10.3e  rms=%10.3e  sign=%s  %s",
                             j + 1, pk, rms, sgn, nearnull ? "NEAR-NULL" : "active");
                    psopt_print(workspace, t);
                }
                const int npath = problem.phase[iphase - 1].npath;
                if (npath > 0 && solution.dual.path[iphase - 1].cols() > 0) {
                    MatrixXd& mu = solution.dual.path[iphase - 1];
                    psopt_print(workspace, "\n  Path-constraint multipliers (peak |mu| per constraint):");
                    for (int j = 0; j < npath && j < (int) mu.rows(); j++) {
                        const double pk = mu.row(j).cwiseAbs().maxCoeff();
                        snprintf(t, tn, "\n    path %3i:   peak=%10.3e   %s",
                                 j + 1, pk, (pk > afloor) ? "active" : "inactive");
                        psopt_print(workspace, t);
                    }
                }
                psopt_print(workspace,
                    "\n  Read: a costate negligible against the dominant one (NEAR-NULL) flags a state"
                    "\n        whose adjoint is structurally inactive; a constant-sign costate is a clean"
                    "\n        shadow price, a sign-changing one is not.");
            }
        }

        // ---- consistency residuals (D2b, diagnostic_level >= 2) ----
        if (algorithm.diagnostic_level >= 2) {
            // (a) Hamiltonian constancy from the stored H = L + lambda.f
            MatrixXd& H = solution.dual.Hamiltonian[iphase - 1];   // 1 x ncols
            if (H.cols() > 0) {
                const double Hmax = H.maxCoeff(), Hmin = H.minCoeff(), Hmean = H.mean();
                const double dH = (Hmax - Hmin) / (std::fabs(Hmean) + 1.0e-10);
                snprintf(t, tn, "\n  Hamiltonian constancy:  dH=(max-min)/(|mean|+eps)=%10.3e", dH);
                psopt_print(workspace, t);
                psopt_print(workspace,
                    "\n        (an optimality indicator for autonomous problems; for time-varying"
                    "\n         problems H need not be constant, so a large dH is not necessarily a defect)");
            }

            // (b) stationarity residual dH/du by central finite differences on the user dae/integrand.
            // H = L + lambda.f; dH/du vanishes only at interior controls. At an active bound it equals
            // the bound multiplier (non-zero by design), so it is summarised separately.
            const int nc   = problem.phase[iphase - 1].ncontrols;
            const int ns   = nstates;
            const int npar  = problem.phase[iphase - 1].nparameters;
            const int np    = problem.phase[iphase - 1].npath;
            MatrixXd& Rst = solution.stationarity_residual[iphase - 1];   // nc x ncols
            MatrixXd& lam = solution.dual.costates[iphase - 1];
            MatrixXd& Xs  = solution.states[iphase - 1];
            MatrixXd& Us  = solution.controls[iphase - 1];
            MatrixXd& Tn  = solution.nodes[iphase - 1];
            const int ncols = (int) Us.cols();

            std::vector<adouble> st(std::max(ns,1)),  ct(std::max(nc,1)),
                                 pa(std::max(npar,1)), de(std::max(ns,1)), pth(std::max(np,1));
            for (int l = 0; l < npar; l++) pa[l] = solution.parameters[iphase-1](l);

            double max_interior = 0.0;
            int n_bound_active = 0;
            int n_path_active = 0;
            bool any_interior = false;
            const double pfloor = 1.0e-9;

            for (int k = 0; k < ncols; k++) {
                adouble tm = Tn(0, k);
                for (int l = 0; l < ns; l++) st[l] = Xs(l, k);
                for (int c = 0; c < nc; c++) ct[c] = Us(c, k);

                // is a path constraint active at this node?
                bool path_active_k = false;
                if (np > 0 && solution.dual.path[iphase-1].cols() == ncols) {
                    for (int p = 0; p < np && p < (int) solution.dual.path[iphase-1].rows(); p++)
                        if (std::fabs(solution.dual.path[iphase-1](p, k)) > pfloor) path_active_k = true;
                }

                for (int c = 0; c < nc; c++) {
                    const double u0 = Us(c, k);
                    const double h  = 1.0e-6 * std::max(1.0, std::fabs(u0));

                    ct[c] = u0 + h;
                    double Hp = (problem.integrand_cost)
                        ? problem.integrand_cost(&st[0], &ct[0], &pa[0], tm, solution.xad, iphase, workspace).value() : 0.0;
                    problem.dae(&de[0], &pth[0], &st[0], &ct[0], &pa[0], tm, solution.xad, iphase, workspace);
                    for (int j = 0; j < ns; j++) Hp += lam(j, k) * de[j].value();

                    ct[c] = u0 - h;
                    double Hm = (problem.integrand_cost)
                        ? problem.integrand_cost(&st[0], &ct[0], &pa[0], tm, solution.xad, iphase, workspace).value() : 0.0;
                    problem.dae(&de[0], &pth[0], &st[0], &ct[0], &pa[0], tm, solution.xad, iphase, workspace);
                    for (int j = 0; j < ns; j++) Hm += lam(j, k) * de[j].value();

                    ct[c] = u0;                              // restore
                    const double dHduc = (Hp - Hm) / (2.0 * h);
                    Rst(c, k) = dHduc;

                    const double lo = problem.phase[iphase-1].bounds.lower.controls(c);
                    const double up = problem.phase[iphase-1].bounds.upper.controls(c);
                    const double btol = 1.0e-6 * (1.0 + std::max(std::fabs(lo), std::fabs(up)));
                    const bool at_bound = (std::fabs(u0 - lo) < btol) || (std::fabs(up - u0) < btol);
                    if (at_bound) {
                        n_bound_active++;
                    } else if (path_active_k) {
                        n_path_active++;
                    } else {
                        any_interior = true;
                        if (std::fabs(dHduc) > max_interior) max_interior = std::fabs(dHduc);
                    }
                }
            }

            if (nc > 0) {
                if (any_interior)
                    snprintf(t, tn, "\n  Stationarity |dH/du| at free interior controls:  max=%10.3e", max_interior);
                else
                    snprintf(t, tn, "\n  Stationarity |dH/du| at free interior controls:  none (all at a bound or active-path node)");
                psopt_print(workspace, t);
                snprintf(t, tn, "\n    control nodes excluded: %d at an active bound, %d at an active-path node (of %d)",
                         n_bound_active, n_path_active, nc * ncols);
                psopt_print(workspace, t);
                psopt_print(workspace,
                    "\n  Read: dH/du = 0 is the stationarity condition only for a free interior control with no"
                    "\n        active path constraint. The max above is taken over exactly those nodes, so ~0"
                    "\n        confirms stationarity; at bound/active-path nodes dH/du is balanced by the"
                    "\n        corresponding multiplier and is expected nonzero.");
            }
        }
        psopt_print(workspace, "\n");
    }

    // Not per phase: the Jacobian is of the whole NLP.
    if (algorithm.diagnostic_level >= 2) {
        constraint_jacobian_conditioning(problem, algorithm, workspace);
        psopt_print(workspace, "\n");
    }
}
