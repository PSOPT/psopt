//  Clarabel as a PSOPT QP backend plugin.
//
//  Built as a shared object with hidden visibility; see psopt_qp_plugin.h for why the
//  backends are kept apart. Clarabel is a Rust library reached through its C interface,
//  which is exactly the kind of dependency the plugin boundary was built for.
//
//  Clarabel states the problem in conic form,
//
//    min 1/2 x'Px + q'x   s.t.  A x + s = b,  s in K
//
//  so a two-sided row has to be split. With K the zero cone for the equalities followed
//  by the nonnegative orthant for everything else, a row l <= a'x <= u contributes
//  a'x + s = u and -a'x + s = -l, a simple bound contributes the same two rows with a
//  unit coefficient, and a row whose bounds coincide contributes one row in the zero
//  cone. Rows free on both sides contribute nothing.
//
//  Its stationarity condition is P x + q + A'z = 0, PSOPT's sign convention, so with the
//  splitting above the multiplier of an original row is z_eq for an equality and
//  z_upper - z_lower for an inequality, and the bound multiplier PSOPT wants -- which
//  enters its own condition with a minus -- is z_lower - z_upper. This is pinned by
//  SQPSolver.ClarabelPluginSolvesAndUsesTheRightDualSign, on the same QP as every other
//  backend.
//
//  Clarabel is a primal-dual interior-point method with the same requirement as the
//  other interior-point and splitting backends here, and unlike GALAHAD's QPA: P must be
//  positive semidefinite, and the SQP raises its shift until the model is accepted. Like
//  PIQP, and for the reason set out in that plugin, an exhausted iteration limit is
//  reported as a failure rather than as an approximate step; unlike PIQP it also has a
//  reduced-tolerance verdict of its own, which is a genuine near-solution rather than an
//  arbitrary point on the central path, and that one is passed on as approximate. That
//  distinction was measured rather than assumed: refusing the almost-solved verdict too
//  costs catmix outright and takes moon from 33 SQP iterations to 57.
//
//  Clarabel is also the one backend built here that solves conic programmes rather than
//  only quadratic ones, and so the one that can take PSOPT's Euclidean trust region --
//  algorithm.trust_region = "l2" -- which is a second-order cone and is not expressible in
//  a quadratic programme at all. Every other plugin refuses such a subproblem. The cone
//  does not relieve Clarabel of needing P positive semidefinite, so the SQP's shift ladder
//  is unchanged by it: what the region changes is the shape of the restriction on the
//  step, not the convexity of the model.

#include "psopt_qp_plugin.h"

// Clarabel's C header carries no extern "C" guard of its own, so the declarations would
// otherwise be mangled as C++ and fail to find the Rust library's symbols.
extern "C" {
#include <clarabel.h>
}

#include <algorithm>
#include <vector>

namespace {

struct Csc {
    std::vector<uintptr_t> p, i;
    std::vector<double>    x;
};

// Column-compress a set of triplets that are gathered column by column already, so that
// only the rows within each column need ordering.
void compress(int nrows, int ncols,
              std::vector<std::vector<std::pair<int,double> > >& cols, Csc& out)
{
    (void) nrows;
    out.p.assign((size_t) ncols + 1, 0);
    for (int j = 0; j < ncols; j++) {
        std::sort(cols[(size_t) j].begin(), cols[(size_t) j].end());
        for (size_t k = 0; k < cols[(size_t) j].size(); k++) {
            out.i.push_back((uintptr_t) cols[(size_t) j][k].first);
            out.x.push_back(cols[(size_t) j][k].second);
        }
        out.p[(size_t) j + 1] = (uintptr_t) out.i.size();
    }
}

} // anonymous namespace


extern "C" {

__attribute__((visibility("default")))
int psopt_qp_abi_version(void) { return PSOPT_QP_ABI_VERSION; }

__attribute__((visibility("default")))
const char* psopt_qp_name(void) { return "Clarabel"; }

__attribute__((visibility("default")))
int psopt_qp_solve(const psopt_qp_problem* p, psopt_qp_solution* s)
{
    s->iterations = 0;
    s->status     = PSOPT_QP_FAILED;

    if (p == NULL || s == NULL || p->abi_version != PSOPT_QP_ABI_VERSION)
        return s->status;

    const int n = p->n, m = p->m;
    if (n <= 0) return s->status;

    // ---- the Euclidean trust region, if one was asked for -------------------------------
    // ||d[0..k)||_2 <= Delta is a second-order cone of dimension k+1, and Clarabel is the
    // one backend built here that has such cones. In Clarabel's form the cone variable is
    // s = b - A x, so a block whose rows are b = (Delta, 0, ..., 0) and A = (0'; -I_k)
    // gives s = (Delta, d_0, ..., d_{k-1}), and s in the second-order cone says exactly
    // Delta >= ||d[0..k)||_2. The block goes last, after the orthant, because Clarabel
    // reads the cones in the order their rows appear.
    const bool   has_tr = (p->trust_radius < PSOPT_QP_INFINITY);
    const int    ktr    = has_tr ? ((p->trust_dim > 0 && p->trust_dim <= n) ? p->trust_dim : n) : 0;
    const double Rtr    = has_tr ? p->trust_radius : 0.0;
    // A region of zero or negative radius admits only d = 0, which is not a step; the SQP
    // never asks for one, and a backend should not be the place that discovers it did.
    if (has_tr && !(Rtr > 0.0)) return s->status;

    // ---- the row layout ---------------------------------------------------------------
    // The zero cone first, then the nonnegative orthant, because Clarabel takes the cones
    // in the order the rows appear.
    std::vector<int> eq_of(m > 0 ? (size_t) m : 1, -1);
    std::vector<int> up_of(m > 0 ? (size_t) m : 1, -1);
    std::vector<int> lo_of(m > 0 ? (size_t) m : 1, -1);
    std::vector<int> bup_of((size_t) n, -1), blo_of((size_t) n, -1);

    int neq = 0;
    for (int i = 0; i < m; i++)
        if (p->lbA[i] == p->ubA[i] && p->lbA[i] > -PSOPT_QP_INFINITY
                                   && p->lbA[i] <  PSOPT_QP_INFINITY)
            eq_of[(size_t) i] = neq++;

    int row = neq;
    for (int i = 0; i < m; i++) {
        if (eq_of[(size_t) i] >= 0) continue;
        if (p->ubA[i] <  PSOPT_QP_INFINITY) up_of[(size_t) i] = row++;
        if (p->lbA[i] > -PSOPT_QP_INFINITY) lo_of[(size_t) i] = row++;
    }
    for (int j = 0; j < n; j++) {
        if (p->ub[j] <  PSOPT_QP_INFINITY) bup_of[(size_t) j] = row++;
        if (p->lb[j] > -PSOPT_QP_INFINITY) blo_of[(size_t) j] = row++;
    }
    const int nnneg   = row - neq;
    const int tr_row0 = row;                 // first row of the cone block, if any
    if (has_tr) row += ktr + 1;
    const int nrows = row;

    // ---- P, its upper triangle ---------------------------------------------------------
    // Clarabel takes the whole symmetric matrix and reduces it itself; handing over the
    // triangle saves it that copy.
    Csc P;
    {
        std::vector<std::vector<std::pair<int,double> > > cols((size_t) n);
        for (int j = 0; j < n; j++) {
            if (p->H_p != NULL) {
                for (long long k = p->H_p[j]; k < p->H_p[j+1]; k++)
                    if (p->H_i[k] <= j)
                        cols[(size_t) j].push_back(std::make_pair((int) p->H_i[k], p->H_x[k]));
            }
            else {
                for (int i = 0; i <= j; i++) {
                    const double v = p->H_dense[(size_t) j*n + i];   // column-major
                    if (v != 0.0) cols[(size_t) j].push_back(std::make_pair(i, v));
                }
            }
        }
        compress(n, n, cols, P);
    }

    // ---- A and b -------------------------------------------------------------------------
    Csc A;
    std::vector<double> b((size_t) (nrows > 0 ? nrows : 1), 0.0);
    {
        std::vector<std::vector<std::pair<int,double> > > cols((size_t) n);
        if (m > 0) {
            for (int j = 0; j < n; j++)
                for (long long k = p->A_p[j]; k < p->A_p[j+1]; k++) {
                    const int i = (int) p->A_i[k];
                    const double v = p->A_x[k];
                    if (eq_of[(size_t) i] >= 0)
                        cols[(size_t) j].push_back(std::make_pair(eq_of[(size_t) i],  v));
                    else {
                        if (up_of[(size_t) i] >= 0)
                            cols[(size_t) j].push_back(std::make_pair(up_of[(size_t) i],  v));
                        if (lo_of[(size_t) i] >= 0)
                            cols[(size_t) j].push_back(std::make_pair(lo_of[(size_t) i], -v));
                    }
                }
        }
        for (int j = 0; j < n; j++) {
            if (bup_of[(size_t) j] >= 0)
                cols[(size_t) j].push_back(std::make_pair(bup_of[(size_t) j],  1.0));
            if (blo_of[(size_t) j] >= 0)
                cols[(size_t) j].push_back(std::make_pair(blo_of[(size_t) j], -1.0));
        }
        // The cone block: its leading row has no coefficients at all, since s_0 is the
        // constant radius, and the k that follow carry -I.
        if (has_tr)
            for (int j = 0; j < ktr; j++)
                cols[(size_t) j].push_back(std::make_pair(tr_row0 + 1 + j, -1.0));
        compress(nrows, n, cols, A);
    }

    for (int i = 0; i < m; i++) {
        if (eq_of[(size_t) i] >= 0) b[(size_t) eq_of[(size_t) i]] = p->lbA[i];
        else {
            if (up_of[(size_t) i] >= 0) b[(size_t) up_of[(size_t) i]] =  p->ubA[i];
            if (lo_of[(size_t) i] >= 0) b[(size_t) lo_of[(size_t) i]] = -p->lbA[i];
        }
    }
    for (int j = 0; j < n; j++) {
        if (bup_of[(size_t) j] >= 0) b[(size_t) bup_of[(size_t) j]] =  p->ub[j];
        if (blo_of[(size_t) j] >= 0) b[(size_t) blo_of[(size_t) j]] = -p->lb[j];
    }
    if (has_tr) b[(size_t) tr_row0] = Rtr;      // the rest of the block is zero already

    std::vector<double> q((size_t) n);
    for (int j = 0; j < n; j++) q[(size_t) j] = p->g[j];

    // ---- solve ---------------------------------------------------------------------------
    ClarabelCscMatrix Pm, Am;
    clarabel_CscMatrix_init(&Pm, (uintptr_t) n, (uintptr_t) n, &P.p[0],
                            P.i.empty() ? NULL : &P.i[0], P.x.empty() ? NULL : &P.x[0]);
    clarabel_CscMatrix_init(&Am, (uintptr_t) nrows, (uintptr_t) n, &A.p[0],
                            A.i.empty() ? NULL : &A.i[0], A.x.empty() ? NULL : &A.x[0]);

    std::vector<ClarabelSupportedConeT> cones;
    if (neq > 0) {
        ClarabelSupportedConeT c;
        c.tag = ClarabelZeroConeT_Tag;
        c.zero_cone_t = (uintptr_t) neq;
        cones.push_back(c);
    }
    if (nnneg > 0) {
        ClarabelSupportedConeT c;
        c.tag = ClarabelNonnegativeConeT_Tag;
        c.nonnegative_cone_t = (uintptr_t) nnneg;
        cones.push_back(c);
    }
    if (has_tr) {
        ClarabelSupportedConeT c;
        c.tag = ClarabelSecondOrderConeT_Tag;
        c.second_order_cone_t = (uintptr_t) (ktr + 1);
        cones.push_back(c);
    }
    if (cones.empty()) return s->status;     // an unconstrained subproblem is not one

    ClarabelDefaultSettings settings = clarabel_DefaultSettings_default();
    settings.verbose  = false;
    // The SQP asks for a subproblem two orders tighter than the NLP's own tolerance.
    // Clarabel's reduced tolerances stay at their defaults, and are what its
    // almost-solved verdict is measured against.
    settings.tol_gap_abs = std::max(1.0e-10, 1.0e-2*p->tolerance);
    settings.tol_feas    = std::max(1.0e-10, 1.0e-2*p->tolerance);
    // An interior-point method that has not converged in a few hundred iterations is not
    // going to. The caller's budget is an active-set budget and far too generous here.
    settings.max_iter = (uint32_t) std::min(std::max(200, p->max_iter), 500);

    ClarabelDefaultSolver* solver =
        clarabel_DefaultSolver_new(&Pm, &q[0], &Am, &b[0],
                                   (uintptr_t) cones.size(), &cones[0], &settings);
    if (solver == NULL) return s->status;

    clarabel_DefaultSolver_solve(solver);
    ClarabelDefaultSolution sol = clarabel_DefaultSolver_solution(solver);

    s->iterations = (int) sol.iterations;

    const bool solved      = (sol.status == ClarabelSolved);
    const bool approximate = (sol.status == ClarabelAlmostSolved);

    if ((solved || approximate) && sol.x != NULL && sol.z != NULL) {
        for (int j = 0; j < n; j++) s->d[j] = sol.x[j];

        for (int i = 0; i < m; i++) {
            if (eq_of[(size_t) i] >= 0) s->lambda[i] = sol.z[eq_of[(size_t) i]];
            else {
                const double zu = (up_of[(size_t) i] >= 0) ? sol.z[up_of[(size_t) i]] : 0.0;
                const double zl = (lo_of[(size_t) i] >= 0) ? sol.z[lo_of[(size_t) i]] : 0.0;
                s->lambda[i] = zu - zl;
            }
        }
        // The bound multipliers, and only those. The cone block's own multiplier is not
        // reported and must not be: it is the price of the trust region, an artefact of
        // the current radius rather than of the problem, and adding it here would cancel
        // part of the gradient in the SQP's dual residual and report a convergence that
        // has not happened. Leaving it out is what makes a step held back by the region
        // show up as a step that is not yet optimal, which is what it is.
        for (int j = 0; j < n; j++) {
            const double zu = (bup_of[(size_t) j] >= 0) ? sol.z[bup_of[(size_t) j]] : 0.0;
            const double zl = (blo_of[(size_t) j] >= 0) ? sol.z[blo_of[(size_t) j]] : 0.0;
            s->z[j] = zl - zu;
        }

        s->status = solved ? PSOPT_QP_SOLVED : PSOPT_QP_APPROXIMATE;
    }

    clarabel_DefaultSolver_free(solver);
    return s->status;
}

} // extern "C"
