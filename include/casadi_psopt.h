// casadi_psopt.h - CasADi NLP backend for PSOPT (guarded by USE_CASADI).
//
// Routes PSOPT's already-transcribed sparse NLP through CasADi's IPOPT plugin,
// so the solve can use HSL (ma57/ma97) as well as MUMPS. Reuses PSOPT's
// objective/constraint callbacks and its variable/constraint bounds; CasADi
// supplies derivatives by finite differences over those callbacks (exact
// PSOPT-AD derivatives could be wired in later via Callback::get_jacobian).
//
// Enabled by CMake option PSOPT_WITH_CASADI (which links CasADi and defines
// USE_CASADI). Requires a CasADi build with the IPOPT plugin; the PSOPT
// superbuild (-DPSOPT_SUPERBUILD=ON -DPSOPT_WITH_CASADI=ON) provides it.

#ifndef CASADI_PSOPT_H
#define CASADI_PSOPT_H

#ifdef USE_CASADI

#include "psopt.h"

// Solve the current NLP held in `workspace` via CasADi -> IPOPT.
// Mirrors the argument list of the IPOPT branch in NLP_interface().
// Returns 0 on success. On return, *x0 holds the primal solution and
// *lambda the constraint multipliers.
int psopt_casadi_solve(
    Alg&       algorithm,
    MatrixXd*  x0,
    double   (*f)(MatrixXd&, Workspace*),
    void     (*g)(MatrixXd&, MatrixXd*, Workspace*),
    int        nlp_ncons,
    int        nlp_neq,
    MatrixXd*  xlb,
    MatrixXd*  xub,
    MatrixXd*  lambda,
    Workspace* workspace,
    void*      user_data);

#endif // USE_CASADI
#endif // CASADI_PSOPT_H
