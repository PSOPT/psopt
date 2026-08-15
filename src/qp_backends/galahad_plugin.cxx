//  GALAHAD's QPA as a PSOPT QP backend plugin.
//
//  Built as a shared object with hidden visibility; see psopt_qp_plugin.h for why the
//  backends are kept apart. GALAHAD brings a Fortran runtime and its own sparse
//  symmetric solvers, so it is exactly the kind of dependency the plugin boundary was
//  built for.
//
//  QPA states the problem as
//
//    min 1/2 x'Hx + g'x + f   s.t.  c_l <= A x <= c_u,   x_l <= x <= x_u
//
//  which is PSOPT's own statement of the subproblem, bounds included, so unlike the
//  other three backends nothing has to be smuggled in as identity rows. H is given by
//  its lower triangle in coordinate form.
//
//  Its multipliers carry the opposite sign to PSOPT's: QPA's stationarity is
//  H x + g - A' y - z = 0, so lambda = -y and the bound multipliers pass through
//  unchanged. That is qpOASES's convention rather than ProxQP's and QPALM's, and it is
//  pinned by the same test as all the others.
//
//  QPA is an active-set method for a possibly indefinite Hessian, working on an l1
//  penalty formulation of the constraints. It is the one backend here that needs no
//  convexification at all, which makes it the natural comparison for what the shift of
//  the diagonal costs the others.

#include "psopt_qp_plugin.h"

extern "C" {
#include "galahad_qpa.h"
}

#include <algorithm>
#include <cstring>
#include <vector>

extern "C" {

__attribute__((visibility("default")))
int psopt_qp_abi_version(void) { return PSOPT_QP_ABI_VERSION; }

__attribute__((visibility("default")))
const char* psopt_qp_name(void) { return "GALAHAD"; }

__attribute__((visibility("default")))
int psopt_qp_solve(const psopt_qp_problem* p, psopt_qp_solution* s)
{
    s->iterations = 0;
    s->status     = PSOPT_QP_FAILED;

    if (p == NULL || s == NULL || p->abi_version != PSOPT_QP_ABI_VERSION)
        return s->status;

    const int n = p->n, m = p->m;

    // GALAHAD requires m > 0; a subproblem with no general constraints is a bound
    // constrained QP, which is not what QPA is for and is not what the SQP produces
    // for any problem with dynamics.
    if (n <= 0 || m <= 0) return s->status;

    // ---- H, lower triangle, coordinate form ------------------------------------
    std::vector<ipc_> H_row, H_col, A_row, A_col;
    std::vector<rpc_> H_val, A_val;

    if (p->H_p != NULL) {
        for (int j = 0; j < n; j++)
            for (long long k = p->H_p[j]; k < p->H_p[j+1]; k++)
                if (p->H_i[k] >= j) {                       // lower triangle
                    H_row.push_back((ipc_) p->H_i[k]);
                    H_col.push_back((ipc_) j);
                    H_val.push_back((rpc_) p->H_x[k]);
                }
    }
    else {
        for (int j = 0; j < n; j++)
            for (int i = j; i < n; i++) {
                const double v = p->H_dense[(size_t) j*n + i];   // column-major
                if (v != 0.0) {
                    H_row.push_back(i); H_col.push_back(j); H_val.push_back(v);
                }
            }
    }

    for (int j = 0; j < n; j++)
        for (long long k = p->A_p[j]; k < p->A_p[j+1]; k++) {
            A_row.push_back((ipc_) p->A_i[k]);
            A_col.push_back((ipc_) j);
            A_val.push_back((rpc_) p->A_x[k]);
        }

    if (H_val.empty()) { H_row.push_back(0); H_col.push_back(0); H_val.push_back(0.0); }
    if (A_val.empty()) return s->status;

    std::vector<rpc_> g(n), c_l(m), c_u(m), x_l(n), x_u(n);
    std::vector<rpc_> x(n, 0.0), c(m, 0.0), y(m, 0.0), z(n, 0.0);
    std::vector<ipc_> x_stat(n, 0), c_stat(m, 0);

    // GALAHAD regards a bound as absent when it is *larger* than control.infinity in
    // modulus, so a bound of exactly the infinity value is a bound, not an absence.
    // PSOPT writes an absent bound as exactly +/- 1.0e20, which would therefore be
    // taken literally: every free variable would acquire a pair of enormous but finite
    // bounds, and the subproblem would be a different problem. They are pushed past the
    // threshold here.
    // GALAHAD's own default for control.infinity is smaller than PSOPT's 1.0e20, so an
    // absent bound written as +/- 1.0e20 is already past it and is read as absent. The
    // default is left in place: setting control.infinity to 1.0e20 and pushing the
    // bounds beyond it, which looks equivalent, made QPA stall at its iteration limit
    // on a two-variable problem it otherwise solves in four.
    const rpc_ beyond = (rpc_) PSOPT_QP_INFINITY;

    for (int j = 0; j < n; j++) {
        g[j]   = p->g[j];
        x_l[j] = (p->lb[j] <= -PSOPT_QP_INFINITY) ? -beyond : (rpc_) p->lb[j];
        x_u[j] = (p->ub[j] >=  PSOPT_QP_INFINITY) ?  beyond : (rpc_) p->ub[j];
    }
    for (int i = 0; i < m; i++) {
        c_l[i] = (p->lbA[i] <= -PSOPT_QP_INFINITY) ? -beyond : (rpc_) p->lbA[i];
        c_u[i] = (p->ubA[i] >=  PSOPT_QP_INFINITY) ?  beyond : (rpc_) p->ubA[i];
    }

    struct qpa_control_type control;
    struct qpa_inform_type  inform;
    void* data   = NULL;
    ipc_  status = 0;

    qpa_initialize(&data, &control, &status);

    control.f_indexing  = 0;        // C indexing, which is what is assembled above
    control.print_level = 0;
    control.out         = 0;
    control.error       = 0;

    // The default symmetric solver prints a warning per factorisation unless
    // OMP_CANCELLATION and OMP_PROC_BIND are set in the environment, which a library
    // cannot arrange for its caller. Naming a different one is tempting and wrong:
    // asking for a solver GALAHAD was not built with makes it refuse the problem, and
    // the refusal arrives as an ordinary solve failure with nothing to distinguish it
    // from a hard subproblem. The default is left alone.

    qpa_import(&control, &data, &status, (ipc_) n, (ipc_) m,
               "coordinate", (ipc_) H_val.size(), &H_row[0], &H_col[0], NULL,
               "coordinate", (ipc_) A_val.size(), &A_row[0], &A_col[0], NULL);

    if (status == 0) {
        status = 1;                                  // cold start
        qpa_solve_qp(&data, &status, (ipc_) n, (ipc_) m,
                     (ipc_) H_val.size(), &H_val[0], &g[0], (rpc_) 0.0,
                     (ipc_) A_val.size(), &A_val[0], &c_l[0], &c_u[0],
                     &x_l[0], &x_u[0],
                     &x[0], &c[0], &y[0], &z[0], &x_stat[0], &c_stat[0]);

        ipc_ info_status = 0;
        qpa_information(&data, &inform, &info_status);
        s->iterations = (int) inform.iter;

        // 0 is a solution; -18 is QPA's iteration limit, which yields an approximate
        // step, and whether a step is any good is the caller's business rather than
        // something to escalate.
        if (status == 0 || status == -18) {
            for (int j = 0; j < n; j++) s->d[j] = (double) x[j];
            // QPA's multipliers carry the opposite sign to PSOPT's.
            for (int i = 0; i < m; i++) s->lambda[i] = -(double) y[i];
            for (int j = 0; j < n; j++) s->z[j]      =  (double) z[j];

            s->status = (status == 0) ? PSOPT_QP_SOLVED : PSOPT_QP_APPROXIMATE;
        }
    }

    ipc_ term_status = 0;
    qpa_terminate(&data, &control, &inform);
    (void) term_status;

    return s->status;
}

} // extern "C"
