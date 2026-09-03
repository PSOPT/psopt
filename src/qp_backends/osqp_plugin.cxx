//  OSQP as a PSOPT QP backend plugin.
//
//  Built as a shared object with hidden visibility; see psopt_qp_plugin.h for why the
//  backends are kept apart.
//
//  OSQP states the problem as
//
//    min 1/2 x'Px + q'x   s.t.  l <= A x <= u
//
//  with one two-sided constraint block and no separate provision for simple bounds, so
//  the bounds are appended to A as identity rows. Its stationarity condition is
//  P x + q + A' y = 0, PSOPT's sign convention, so the multipliers pass through
//  unnegated. P must be given as its upper triangle, and must be positive semidefinite:
//  OSQP is a first-order operator-splitting method and has no way to cope with
//  indefinite curvature, so the SQP's shift of the Hessian's diagonal is doing real
//  work here rather than merely satisfying a factorisation.
//
//  What OSQP is good at is warm starts and very large sparse problems; what it is not
//  good at is the last few digits, being first order. The SQP asks for a tolerance two
//  orders below the NLP's own, which OSQP reaches quickly and then approaches slowly;
//  an exhausted iteration limit is reported as an approximate step rather than a
//  failure, which is what the interface is for.

#include "psopt_qp_plugin.h"

// OSQP's exported CMake target sets the include directory to <prefix>/include/osqp,
// so the header it means is <osqp.h>. Writing <osqp/osqp.h> works only where
// <prefix>/include also happens to be on the default search path -- true of
// /usr/local on a Linux box, false of Homebrew's /opt/homebrew, where this failed to
// build. Older layouts put the header a directory up, so both spellings are tried
// rather than picking one and being wrong somewhere else.
#if defined(__has_include)
#  if __has_include(<osqp.h>)
#    include <osqp.h>
#  elif __has_include(<osqp/osqp.h>)
#    include <osqp/osqp.h>
#  else
#    error "OSQP headers were not found. Check that osqp_DIR points at the directory holding osqpConfig.cmake."
#  endif
#else
#  include <osqp.h>
#endif

#include <algorithm>
#include <vector>

namespace {

// Compressed columns in OSQP's index type, accumulated from triplets that arrive in
// column order already, so no sort is needed -- only the upper triangle of the Hessian
// and, for the constraint matrix, every entry.
struct Csc {
    std::vector<OSQPInt>   p, i;
    std::vector<OSQPFloat> x;
};

} // anonymous namespace


extern "C" {

__attribute__((visibility("default")))
int psopt_qp_abi_version(void) { return PSOPT_QP_ABI_VERSION; }

__attribute__((visibility("default")))
const char* psopt_qp_name(void) { return "OSQP"; }

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

    const int n = p->n, m = p->m;
    const OSQPInt nrows = (OSQPInt)(m + n);

    // ---- P, upper triangle -----------------------------------------------------
    Csc P;
    P.p.resize((size_t) n + 1, 0);
    for (int j = 0; j < n; j++) {
        if (p->H_p != NULL) {
            for (long long k = p->H_p[j]; k < p->H_p[j+1]; k++)
                if (p->H_i[k] <= j) {                      // upper triangle only
                    P.i.push_back((OSQPInt) p->H_i[k]);
                    P.x.push_back((OSQPFloat) p->H_x[k]);
                }
        }
        else {
            for (int i = 0; i <= j; i++) {
                const double v = p->H_dense[(size_t) j*n + i];   // column-major
                if (v != 0.0) { P.i.push_back(i); P.x.push_back(v); }
            }
        }
        P.p[(size_t) j + 1] = (OSQPInt) P.i.size();
    }

    // ---- A, with the bounds as identity rows below the constraints --------------
    Csc A;
    A.p.resize((size_t) n + 1, 0);
    for (int j = 0; j < n; j++) {
        if (m > 0)
            for (long long k = p->A_p[j]; k < p->A_p[j+1]; k++) {
                A.i.push_back((OSQPInt) p->A_i[k]);
                A.x.push_back((OSQPFloat) p->A_x[k]);
            }
        A.i.push_back((OSQPInt)(m + j));                   // the bound on d_j
        A.x.push_back(1.0);
        A.p[(size_t) j + 1] = (OSQPInt) A.i.size();
    }

    std::vector<OSQPFloat> q((size_t) n), l((size_t) nrows), u((size_t) nrows);
    for (int j = 0; j < n; j++) q[(size_t) j] = p->g[j];
    for (int i = 0; i < m; i++) {
        l[(size_t) i] = std::max(p->lbA[i], -OSQP_INFTY);
        u[(size_t) i] = std::min(p->ubA[i],  OSQP_INFTY);
    }
    for (int j = 0; j < n; j++) {
        l[(size_t)(m + j)] = std::max(p->lb[j], -OSQP_INFTY);
        u[(size_t)(m + j)] = std::min(p->ub[j],  OSQP_INFTY);
    }

    OSQPCscMatrix* Pm = OSQPCscMatrix_new((OSQPInt) n, (OSQPInt) n,
                                          (OSQPInt) P.x.size(),
                                          P.x.empty() ? NULL : &P.x[0],
                                          P.i.empty() ? NULL : &P.i[0], &P.p[0]);
    OSQPCscMatrix* Am = OSQPCscMatrix_new(nrows, (OSQPInt) n,
                                          (OSQPInt) A.x.size(),
                                          A.x.empty() ? NULL : &A.x[0],
                                          A.i.empty() ? NULL : &A.i[0], &A.p[0]);

    OSQPSettings* settings = OSQPSettings_new();
    if (Pm == NULL || Am == NULL || settings == NULL) {
        OSQPCscMatrix_free(Pm); OSQPCscMatrix_free(Am); OSQPSettings_free(settings);
        return s->status;
    }

    osqp_set_default_settings(settings);
    settings->verbose        = 0;
    settings->eps_abs        = std::max(1.0e-10, 1.0e-2*p->tolerance);
    settings->eps_rel        = 0.0;
    settings->max_iter       = std::max(200, 20*p->max_iter);
    settings->polishing      = 1;    // a cheap Newton step on the identified active set,
                                     // which is what lifts a first-order answer to
                                     // something an SQP can use
    settings->scaled_termination = 0;

    OSQPSolver* solver = NULL;
    const OSQPInt rc = osqp_setup(&solver, Pm, q.empty() ? NULL : &q[0], Am,
                                  &l[0], &u[0], nrows, (OSQPInt) n, settings);

    if (rc == 0 && solver != NULL) {
        osqp_solve(solver);

        s->iterations = (int) solver->info->iter;

        const bool solved = (solver->info->status_val == OSQP_SOLVED)
                         || (solver->info->status_val == OSQP_SOLVED_INACCURATE);
        const bool approximate = (solver->info->status_val == OSQP_MAX_ITER_REACHED);

        if (solved || approximate) {
            for (int j = 0; j < n; j++) s->d[j]      = solver->solution->x[j];
            for (int i = 0; i < m; i++) s->lambda[i] = solver->solution->y[i];
            // The bound rows sit in A, so their multipliers arrive with the sign of a
            // general constraint and are turned back into bound multipliers.
            for (int j = 0; j < n; j++) s->z[j]      = -solver->solution->y[m + j];

            s->status = (solver->info->status_val == OSQP_SOLVED)
                            ? PSOPT_QP_SOLVED : PSOPT_QP_APPROXIMATE;
        }
    }

    osqp_cleanup(solver);
    OSQPCscMatrix_free(Pm);
    OSQPCscMatrix_free(Am);
    OSQPSettings_free(settings);

    return s->status;
}

} // extern "C"
