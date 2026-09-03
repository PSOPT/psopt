//  ProxQP as a PSOPT QP backend plugin.
//
//  Built as a shared object with hidden visibility, so that the only symbols it exports
//  are the three of psopt_qp_plugin.h and everything ProxQP brings with it -- its
//  linear algebra, its ordering -- stays inside. See psopt_qp_plugin.h for why that
//  matters.
//
//  ProxQP states the problem as
//
//    min 1/2 x'Hx + g'x   s.t.  A_eq x = b,   l <= C x <= u
//
//  with no separate provision for simple bounds in its sparse interface, so the bounds
//  go into C as identity rows. Its stationarity condition is
//  H x + g + A_eq' y + C' z = 0, which is already PSOPT's sign convention, so the
//  multipliers pass through unnegated.
//
//  Equality rows are separated from inequality rows rather than passed as inequalities
//  with coincident bounds: a proximal method treats the two blocks differently, and
//  there is no reason to hide from it that most of a collocated problem's constraints
//  are equalities.

#include "psopt_qp_plugin.h"

#include <proxsuite/proxqp/sparse/sparse.hpp>

#include <algorithm>
#include <vector>

namespace {

typedef long long Idx;
typedef Eigen::SparseMatrix<double, Eigen::ColMajor, Idx> SpMat;
typedef Eigen::Triplet<double, Idx>                       Trip;

void csc_to_triplets(const long long* p, const long long* i, const double* x,
                     int nc, int row_offset, int col_offset, std::vector<Trip>& out)
{
    for (int j = 0; j < nc; j++)
        for (long long k = p[j]; k < p[j+1]; k++)
            out.push_back(Trip((Idx)(i[k] + row_offset), (Idx)(j + col_offset), x[k]));
}

} // anonymous namespace


extern "C" {

__attribute__((visibility("default")))
int psopt_qp_abi_version(void) { return PSOPT_QP_ABI_VERSION; }

__attribute__((visibility("default")))
const char* psopt_qp_name(void) { return "ProxQP"; }

__attribute__((visibility("default")))
int psopt_qp_solve(const psopt_qp_problem* p, psopt_qp_solution* s)
{
    s->iterations = 0;
    s->status     = PSOPT_QP_FAILED;

    if (p == NULL || s == NULL || p->abi_version != PSOPT_QP_ABI_VERSION)
        return s->status;

    // A quadratic programming solver cannot impose a Euclidean trust region: the ball is
    // a second-order cone and this backend has no cones. Refusing is required of it --
    // solving the subproblem without the region would return an unrestricted step as
    // though it were a restricted one. See psopt_qp_plugin.h.
    if (p->trust_radius < PSOPT_QP_INFINITY) return s->status;

    try {
        const int n = p->n, m = p->m;

        // Which rows are equalities. The partition follows the bounds, so it is a
        // property of the subproblem and is recomputed here rather than carried across
        // the interface.
        std::vector<int> eq_of_row(std::max(m,1), -1), in_of_row(std::max(m,1), -1);
        int n_eq = 0, n_in = 0;
        for (int i = 0; i < m; i++) {
            if (p->ubA[i] - p->lbA[i] <= 0.0) eq_of_row[i] = n_eq++;
            else                              in_of_row[i] = n_in++;
        }
        const int n_in_total = n_in + n;      // the bounds occupy the last n rows

        std::vector<Trip> th, ta, tc;

        if (p->H_p != NULL) csc_to_triplets(p->H_p, p->H_i, p->H_x, n, 0, 0, th);
        else {
            th.reserve((size_t) n*n);
            for (int j = 0; j < n; j++)
                for (int i = 0; i < n; i++) {
                    const double v = p->H_dense[(size_t) j*n + i];   // column-major
                    if (v != 0.0) th.push_back(Trip(i, j, v));
                }
        }

        if (m > 0) {
            std::vector<Trip> tj;
            csc_to_triplets(p->A_p, p->A_i, p->A_x, n, 0, 0, tj);
            ta.reserve(tj.size()); tc.reserve(tj.size() + (size_t) n);
            for (size_t k = 0; k < tj.size(); k++) {
                const int i = (int) tj[k].row();
                if (eq_of_row[i] >= 0)
                    ta.push_back(Trip(eq_of_row[i], tj[k].col(), tj[k].value()));
                else
                    tc.push_back(Trip(in_of_row[i], tj[k].col(), tj[k].value()));
            }
        }
        for (int j = 0; j < n; j++) tc.push_back(Trip(n_in + j, j, 1.0));

        SpMat H(n, n), A(std::max(n_eq,0), n), C(n_in_total, n);
        H.setFromTriplets(th.begin(), th.end());
        if (!ta.empty()) A.setFromTriplets(ta.begin(), ta.end());
        C.setFromTriplets(tc.begin(), tc.end());

        Eigen::VectorXd g(n), b(std::max(n_eq,0)), l(n_in_total), u(n_in_total);
        for (int j = 0; j < n; j++) g(j) = p->g[j];
        for (int i = 0; i < m; i++) {
            if (eq_of_row[i] >= 0) b(eq_of_row[i]) = p->lbA[i];
            else { l(in_of_row[i]) = p->lbA[i]; u(in_of_row[i]) = p->ubA[i]; }
        }
        for (int j = 0; j < n; j++) { l(n_in + j) = p->lb[j]; u(n_in + j) = p->ub[j]; }

        proxsuite::proxqp::sparse::QP<double, Idx> qp(n, n_eq, n_in_total);
        qp.settings.verbose  = false;
        qp.settings.eps_abs  = std::max(1.0e-12, 1.0e-2*p->tolerance);
        qp.settings.eps_rel  = 1.0e-10;
        qp.settings.max_iter = (proxsuite::proxqp::isize) std::max(50, p->max_iter);
        qp.init(H, g, A, b, C, l, u);
        qp.solve();

        s->iterations = (int) qp.results.info.iter;

        const bool solved =
            qp.results.info.status == proxsuite::proxqp::QPSolverOutput::PROXQP_SOLVED;
        const bool approximate =
            qp.results.info.status ==
                proxsuite::proxqp::QPSolverOutput::PROXQP_MAX_ITER_REACHED;

        if (!solved && !approximate) return s->status;   // already PSOPT_QP_FAILED

        for (int j = 0; j < n; j++) s->d[j] = qp.results.x(j);
        for (int i = 0; i < m; i++)
            s->lambda[i] = (eq_of_row[i] >= 0) ? qp.results.y(eq_of_row[i])
                                               : qp.results.z(in_of_row[i]);
        // The bound rows sit in C, so their multipliers arrive with the sign of a
        // general constraint and are turned back into bound multipliers.
        for (int j = 0; j < n; j++) s->z[j] = -qp.results.z(n_in + j);

        s->status = solved ? PSOPT_QP_SOLVED : PSOPT_QP_APPROXIMATE;
    }
    catch (...) {
        // A subproblem the backend cannot take is reported, never escalated: the
        // caller has a restoration step for exactly this.
        s->status = PSOPT_QP_FAILED;
    }

    return s->status;
}

} // extern "C"
