//  PIQP as a PSOPT QP backend plugin.
//
//  Built as a shared object with hidden visibility; see psopt_qp_plugin.h for why the
//  backends are kept apart.
//
//  PIQP states the problem as
//
//    min 1/2 x'Px + c'x   s.t.  A x = b,  h_l <= G x <= h_u,  x_l <= x <= x_u
//
//  which is the first backend here to take equalities as equalities. A collocation
//  problem is mostly equalities -- every defect constraint is one -- so the rows are
//  split rather than all handed over as two-sided inequalities: an interior-point
//  method has nothing to eliminate in a row whose two bounds coincide, and PIQP's KKT
//  system keeps the equality block separate.
//
//  The simple bounds are appended to G as identity rows rather than passed as x_l and
//  x_u. PIQP stores the bound multipliers compressed to the variables that have a
//  finite bound, and the map from that back to the variable it belongs to is internal;
//  as identity rows the multipliers come back by row and the mapping is the identity.
//  The trust region makes most bounds finite most of the time, so little is lost.
//
//  Its stationarity condition is
//
//    P x + c + A'y + G'(z_u - z_l) - z_bl + z_bu = 0
//
//  which is PSOPT's sign convention, so the multipliers pass through unnegated: the
//  equality rows carry y and the inequality rows z_u - z_l. This is pinned by
//  SQPSolver.PiqpPluginSolvesAndUsesTheRightDualSign, on the same QP as every other
//  backend.
//
//  PIQP is a proximal interior-point method: the proximal terms mean it does not need
//  the constraint gradients to be linearly independent, which on a collocation mesh
//  with a rank-deficient Jacobian is the property that matters. Like every
//  interior-point and splitting method here, and unlike GALAHAD's QPA, it requires P to
//  be positive semidefinite; the SQP raises its shift until the model is accepted.

#include "psopt_qp_plugin.h"

#include <piqp/piqp.hpp>

#include <algorithm>
#include <limits>
#include <vector>

namespace {

const double PIQP_PLUGIN_INF = std::numeric_limits<double>::infinity();

inline double lower_of(double v) { return (v <= -PSOPT_QP_INFINITY) ? -PIQP_PLUGIN_INF : v; }
inline double upper_of(double v) { return (v >=  PSOPT_QP_INFINITY) ?  PIQP_PLUGIN_INF : v; }

} // anonymous namespace


extern "C" {

__attribute__((visibility("default")))
int psopt_qp_abi_version(void) { return PSOPT_QP_ABI_VERSION; }

__attribute__((visibility("default")))
const char* psopt_qp_name(void) { return "PIQP"; }

__attribute__((visibility("default")))
int psopt_qp_solve(const psopt_qp_problem* p, psopt_qp_solution* s)
{
    s->iterations = 0;
    s->status     = PSOPT_QP_FAILED;

    if (p == NULL || s == NULL || p->abi_version != PSOPT_QP_ABI_VERSION)
        return s->status;

    const int n = p->n, m = p->m;
    if (n <= 0) return s->status;

    // ---- which rows are equalities -------------------------------------------------
    // A row whose bounds coincide is an equality. Everything else, including a row
    // bounded on one side only, is an inequality.
    // A row free on both sides is left out altogether rather than handed over as a row
    // of G that constrains nothing: PIQP zeroes such a row and says so on every solve,
    // and the SQP's bound vectors are full of infinities whenever the trust region is.
    std::vector<int> eq_row(m > 0 ? (size_t) m : 1, -1);
    std::vector<int> in_row(m > 0 ? (size_t) m : 1, -1);
    int neq = 0, nin = 0;
    for (int i = 0; i < m; i++) {
        const bool lo = (p->lbA[i] > -PSOPT_QP_INFINITY);
        const bool hi = (p->ubA[i] <  PSOPT_QP_INFINITY);
        if (lo && hi && p->lbA[i] == p->ubA[i]) eq_row[(size_t) i] = neq++;
        else if (lo || hi)                      in_row[(size_t) i] = nin++;
    }

    // The bounds ride along as identity rows below the general inequalities, and on the
    // same terms: a variable free on both sides contributes no row.
    std::vector<int> bnd_row((size_t) n, -1);
    int nbnd = 0;
    for (int j = 0; j < n; j++)
        if (p->lb[j] > -PSOPT_QP_INFINITY || p->ub[j] < PSOPT_QP_INFINITY)
            bnd_row[(size_t) j] = nin + nbnd++;

    const int ng = nin + nbnd;

    // ---- the matrices ---------------------------------------------------------------
    typedef Eigen::SparseMatrix<double> Sparse;
    typedef Eigen::Triplet<double>      Tri;

    Sparse P((Eigen::Index) n, (Eigen::Index) n);
    {
        std::vector<Tri> t;
        if (p->H_p != NULL) {
            t.reserve((size_t) p->H_p[n]);
            for (int j = 0; j < n; j++)
                for (long long k = p->H_p[j]; k < p->H_p[j+1]; k++)
                    t.push_back(Tri((int) p->H_i[k], j, p->H_x[k]));
        }
        else {
            t.reserve((size_t) n*(size_t) n);
            for (int j = 0; j < n; j++)
                for (int i = 0; i < n; i++) {
                    const double v = p->H_dense[(size_t) j*n + i];   // column-major
                    if (v != 0.0) t.push_back(Tri(i, j, v));
                }
        }
        // PIQP takes the upper triangle of whatever it is given, so the whole symmetric
        // matrix can go over as it stands.
        P.setFromTriplets(t.begin(), t.end());
    }
    P.makeCompressed();

    Sparse A((Eigen::Index) neq, (Eigen::Index) n);
    Sparse G((Eigen::Index) ng,  (Eigen::Index) n);
    {
        std::vector<Tri> ta, tg;
        if (m > 0) {
            for (int j = 0; j < n; j++)
                for (long long k = p->A_p[j]; k < p->A_p[j+1]; k++) {
                    const int i = (int) p->A_i[k];
                    if      (eq_row[(size_t) i] >= 0) ta.push_back(Tri(eq_row[(size_t) i], j, p->A_x[k]));
                    else if (in_row[(size_t) i] >= 0) tg.push_back(Tri(in_row[(size_t) i], j, p->A_x[k]));
                    // a row free on both sides has no place in either block
                }
        }
        for (int j = 0; j < n; j++)
            if (bnd_row[(size_t) j] >= 0) tg.push_back(Tri(bnd_row[(size_t) j], j, 1.0));
        A.setFromTriplets(ta.begin(), ta.end());
        G.setFromTriplets(tg.begin(), tg.end());
    }
    A.makeCompressed();
    G.makeCompressed();

    Eigen::VectorXd c((Eigen::Index) n), b((Eigen::Index) neq),
                    h_l((Eigen::Index) ng), h_u((Eigen::Index) ng);
    for (int j = 0; j < n; j++) c((Eigen::Index) j) = p->g[j];
    for (int i = 0; i < m; i++) {
        if (eq_row[(size_t) i] >= 0) b((Eigen::Index) eq_row[(size_t) i]) = p->lbA[i];
        else if (in_row[(size_t) i] >= 0) {
            const Eigen::Index r = (Eigen::Index) in_row[(size_t) i];
            h_l(r) = lower_of(p->lbA[i]);
            h_u(r) = upper_of(p->ubA[i]);
        }
    }
    for (int j = 0; j < n; j++) {
        if (bnd_row[(size_t) j] < 0) continue;
        const Eigen::Index r = (Eigen::Index) bnd_row[(size_t) j];
        h_l(r) = lower_of(p->lb[j]);
        h_u(r) = upper_of(p->ub[j]);
    }

    // ---- solve ------------------------------------------------------------------------
    piqp::SparseSolver<double> solver;
    solver.settings().verbose         = false;
    solver.settings().compute_timings = false;
    // The SQP asks for a subproblem two orders tighter than the NLP's own tolerance,
    // which is what eps_abs carries. eps_rel is left at PIQP's own default rather than
    // set to zero: an absolute-only test on a badly scaled subproblem is one an
    // interior-point method may never pass, and the consequence is not a slow solve but
    // a wrong answer -- see the status mapping below.
    solver.settings().eps_abs         = std::max(1.0e-10, 1.0e-2*p->tolerance);
    solver.settings().eps_rel         = 1.0e-9;
    // An interior-point method that has not converged in a few hundred iterations is
    // not going to. The caller's budget is an active-set budget and is far too generous
    // to be worth spending here.
    solver.settings().max_iter        = std::min(std::max(200, p->max_iter), 500);

    piqp::Status status;
    try {
        solver.setup(P, c, A, b, G, h_l, h_u, piqp::nullopt, piqp::nullopt);
        status = solver.solve();
    }
    catch (...) {
        // The ABI says a subproblem a backend cannot take is reported, not escalated.
        return s->status;
    }

    s->iterations = (int) solver.result().info.iter;

    // An exhausted iteration limit is reported as a failure and not as an approximate
    // step, which is where this differs from the OSQP plugin and where it has to. A
    // first-order method stopped early is somewhere near the solution and its step is
    // worth taking. An interior-point method stopped early is on its central path, and
    // the step it returns can be arbitrarily short while the multipliers that come with
    // it still satisfy the regularised stationarity condition to machine precision.
    // The SQP measures its dual error with exactly those multipliers, so a capped solve
    // offered as usable reads as a dual error of 1e-12 next to a constraint violation
    // of 19, and on examples/chance_covariance that was accepted as convergence: the
    // run stopped at iteration zero and reported an objective half again too large.
    // Wrong answers are worse than refused subproblems, and the SQP has an answer to a
    // refused subproblem.
    const bool solved = (status == piqp::Status::PIQP_SOLVED);

    if (!solved) return s->status;

    for (int j = 0; j < n; j++) s->d[j] = solver.result().x((Eigen::Index) j);

    // The inequality multiplier of a row is z_u - z_l, by the stationarity condition
    // quoted above; the identity rows carry the bound multipliers, and PSOPT's
    // convention has H d + g + A'lambda - z = 0, so those enter with the opposite sign.
    for (int i = 0; i < m; i++) {
        if (eq_row[(size_t) i] >= 0)
            s->lambda[i] = solver.result().y((Eigen::Index) eq_row[(size_t) i]);
        else if (in_row[(size_t) i] >= 0) {
            const Eigen::Index r = (Eigen::Index) in_row[(size_t) i];
            s->lambda[i] = solver.result().z_u(r) - solver.result().z_l(r);
        }
        else s->lambda[i] = 0.0;                       // the row constrained nothing
    }
    for (int j = 0; j < n; j++) {
        if (bnd_row[(size_t) j] < 0) { s->z[j] = 0.0; continue; }
        const Eigen::Index r = (Eigen::Index) bnd_row[(size_t) j];
        s->z[j] = solver.result().z_l(r) - solver.result().z_u(r);
    }

    s->status = PSOPT_QP_SOLVED;
    return s->status;
}

} // extern "C"
