//  QPALM as a PSOPT QP backend plugin.
//
//  Built as a shared object with hidden visibility. That is not a tidiness measure
//  here: QPALM factorises through LADEL, which links SuiteSparse's AMD and exports
//  amd_order and amd_l_order built for 64-bit indices, and linking those into a program
//  that already has an AMD of another width is what made PSOPT's ProxQP test crash
//  inside ProxQP. Kept inside a plugin loaded with RTLD_LOCAL, they are nobody else's
//  business. See psopt_qp_plugin.h.
//
//  QPALM states the problem as
//
//    min 1/2 x'Qx + q'x   s.t.  bmin <= A x <= bmax
//
//  with one two-sided constraint block and no separate provision for simple bounds, so
//  the bounds are appended to A as identity rows. Its stationarity condition is
//  Q x + q + A' y = 0, PSOPT's convention, so the multipliers pass through unnegated.

#include "psopt_qp_plugin.h"

#include <qpalm.hpp>
#include <qpalm/constants.h>

#include <algorithm>
#include <vector>

extern "C" {

__attribute__((visibility("default")))
int psopt_qp_abi_version(void) { return PSOPT_QP_ABI_VERSION; }

__attribute__((visibility("default")))
const char* psopt_qp_name(void) { return "QPALM"; }

__attribute__((visibility("default")))
int psopt_qp_solve(const psopt_qp_problem* p, psopt_qp_solution* s)
{
    s->iterations = 0;
    s->status     = PSOPT_QP_FAILED;

    if (p == NULL || s == NULL || p->abi_version != PSOPT_QP_ABI_VERSION)
        return s->status;

    try {
        const int n = p->n, m = p->m;
        const qpalm::index_t nrows = (qpalm::index_t)(m + n);

        std::vector<qpalm::triplet_t> tq, ta;

        if (p->H_p != NULL) {
            for (int j = 0; j < n; j++)
                for (long long k = p->H_p[j]; k < p->H_p[j+1]; k++)
                    tq.push_back(qpalm::triplet_t((qpalm::sp_index_t) p->H_i[k],
                                                  (qpalm::sp_index_t) j, p->H_x[k]));
        }
        else {
            tq.reserve((size_t) n*n);
            for (int j = 0; j < n; j++)
                for (int i = 0; i < n; i++) {
                    const double v = p->H_dense[(size_t) j*n + i];   // column-major
                    if (v != 0.0) tq.push_back(qpalm::triplet_t(i, j, v));
                }
        }

        if (m > 0)
            for (int j = 0; j < n; j++)
                for (long long k = p->A_p[j]; k < p->A_p[j+1]; k++)
                    ta.push_back(qpalm::triplet_t((qpalm::sp_index_t) p->A_i[k],
                                                  (qpalm::sp_index_t) j, p->A_x[k]));
        for (int j = 0; j < n; j++)
            ta.push_back(qpalm::triplet_t((qpalm::sp_index_t)(m + j),
                                          (qpalm::sp_index_t) j, 1.0));

        qpalm::sparse_mat_t Q(n, n), A(nrows, n);
        Q.setFromTriplets(tq.begin(), tq.end());
        A.setFromTriplets(ta.begin(), ta.end());

        qpalm::Data data((qpalm::index_t) n, nrows);
        data.set_Q(Q);
        data.set_A(A);
        data.q = qpalm::vec_t::Zero(n);
        for (int j = 0; j < n; j++) data.q(j) = p->g[j];
        data.c = 0.0;

        data.bmin = qpalm::vec_t::Zero(nrows);
        data.bmax = qpalm::vec_t::Zero(nrows);
        for (int i = 0; i < m; i++) {
            data.bmin(i) = std::max(p->lbA[i], -QPALM_INFTY);
            data.bmax(i) = std::min(p->ubA[i],  QPALM_INFTY);
        }
        for (int j = 0; j < n; j++) {
            data.bmin(m + j) = std::max(p->lb[j], -QPALM_INFTY);
            data.bmax(m + j) = std::min(p->ub[j],  QPALM_INFTY);
        }

        qpalm::Settings settings;
        settings.verbose   = 0;
        settings.eps_abs   = std::max(1.0e-12, 1.0e-2*p->tolerance);
        settings.eps_rel   = 1.0e-10;
        settings.max_iter  = std::max(50, p->max_iter);
        settings.nonconvex = p->nonconvex ? 1 : 0;

        // QPALM's C++ interface reports a refused problem by throwing, and its setup
        // refuses subproblems the other backends accept.
        qpalm::Solver solver(data, settings);
        solver.solve();

        s->iterations = (int) solver.get_info().iter;

        const bool solved      = solver.get_info().status_val == QPALM_SOLVED;
        const bool approximate = solver.get_info().status_val == QPALM_MAX_ITER_REACHED;
        if (!solved && !approximate) return s->status;

        const qpalm::SolutionView sol = solver.get_solution();
        for (int j = 0; j < n; j++) s->d[j]      = sol.x(j);
        for (int i = 0; i < m; i++) s->lambda[i] = sol.y(i);
        for (int j = 0; j < n; j++) s->z[j]      = -sol.y(m + j);

        s->status = solved ? PSOPT_QP_SOLVED : PSOPT_QP_APPROXIMATE;
    }
    catch (...) {
        s->status = PSOPT_QP_FAILED;
    }

    return s->status;
}

} // extern "C"
