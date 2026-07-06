/*********************************************************************************************
 SCIP_interface.cxx  --  mixed-integer NLP backend for PSOPT via SCIP.

 Routes PSOPT's transcribed problem to SCIP so that declared-integer controls and
 parameters are solved to GLOBAL optimality (branch-and-bound), which the
 continuous IPOPT/CasADi backends cannot do.

 SCIP cannot branch on PSOPT's opaque ADOL-C nonlinear callbacks, so this backend
 targets the LINEAR-TRANSCRIPTION subclass: linear dynamics/constraints and a
 linear objective (use auxiliary variables for |u|, etc., exactly as in a MILP).
 For such problems the constraint map is g(x)=J x + b with a CONSTANT Jacobian J,
 which is recovered by unit perturbations and handed to SCIP as explicit linear
 constraints; integrality comes from problem.phase[i].integer_controls /
 integer_parameters. A linearity self-check warns if the transcription is nonlinear
 (then the SCIP model is only a linearization at 0).

 Compiled into libpsopt only when PSOPT_WITH_SCIP=ON (guarded by USE_SCIP).
 Author: added on branch ecbrown, 2026-07.
**********************************************************************************************/

#include "psopt.h"

#ifdef USE_SCIP

#include <scip/scip.h>
#include <scip/scipdefplugins.h>
#include <vector>
#include <cmath>

using Eigen::MatrixXd;

// From the PSOPT core:
void get_constraint_bounds(double* g_l, double* g_u, Workspace* workspace);
int  get_nvars_phase_i(Prob& problem, int i, Workspace* workspace);
void psopt_print(Workspace* workspace, const char* msg);

int psopt_scip_solve(
    Alg&       algorithm,
    MatrixXd*  x0,
    double   (*f)(MatrixXd&, Workspace*),
    void     (*g)(MatrixXd&, MatrixXd*, Workspace*),
    int        nlp_ncons,
    int        /*nlp_neq*/,
    MatrixXd*  xlb,
    MatrixXd*  xub,
    MatrixXd*  lambda,
    Workspace* workspace,
    void*      /*user_data*/)
{
    Prob& problem = *workspace->problem;
    const int nx = workspace->nvars;
    const int ng = nlp_ncons;
    const double tol = 1e-10;

    // ---- constraint bounds (scaled NLP space, same as IPOPT/CasADi see) --------
    std::vector<double> gl(ng > 0 ? ng : 1), gu(ng > 0 ? ng : 1);
    if (ng > 0) get_constraint_bounds(gl.data(), gu.data(), workspace);

    // ---- linear model about the FEASIBLE base point x0 (not 0): the start/end
    //      times t0,tf are decision variables that multiply the dynamics/cost, so
    //      x=0 (=> tf=0) degenerately zeroes every integral. Linearising at x0
    //      (correct t0,tf) is exact for a transcription that is linear in the free
    //      variables.  g(x) ~= b + J (x - x0),  obj ~= f0 + c'(x - x0).
    MatrixXd base = *x0;
    MatrixXd gb(ng > 0 ? ng : 1, 1);
    if (ng > 0) g(base, &gb, workspace);
    double f0 = f(base, workspace);

    std::vector<double> c(nx, 0.0);
    std::vector<std::vector<std::pair<int,double> > > rows(ng > 0 ? ng : 1); // per-constraint (var,coeff)
    MatrixXd xj(nx, 1), gj(ng > 0 ? ng : 1, 1);
    for (int j = 0; j < nx; ++j) {
        xj = base; xj(j, 0) += 1.0;
        c[j] = f(xj, workspace) - f0;
        if (ng > 0) {
            g(xj, &gj, workspace);
            for (int i = 0; i < ng; ++i) {
                double v = gj(i, 0) - gb(i, 0);
                if (v > tol || v < -tol) rows[i].push_back(std::make_pair(j, v));
            }
        }
    }
    // J*x0 per row (to convert the shifted model back to absolute variables x)
    std::vector<double> Jx0(ng > 0 ? ng : 1, 0.0);
    for (int i = 0; i < ng; ++i)
        for (size_t t = 0; t < rows[i].size(); ++t) Jx0[i] += rows[i][t].second * base(rows[i][t].first, 0);

    // ---- linearity self-check: perturb the free (non-time) variables and verify
    //      g reproduces the linear prediction (time vars held fixed to avoid the
    //      benign t0,tf bilinearity). Warns if the dynamics/cost are nonlinear.
    {
        std::vector<char> istime(nx, 0);
        { int xo = 0; for (int p = 0; p < problem.nphases; ++p) { int nv = get_nvars_phase_i(problem, p, workspace); istime[xo+nv-2]=1; istime[xo+nv-1]=1; xo += nv; } }
        MatrixXd xd = base, gd(ng > 0 ? ng : 1, 1);
        for (int j = 0; j < nx; ++j) if (!istime[j]) xd(j, 0) += 0.5;
        if (ng > 0) g(xd, &gd, workspace);
        double worst = 0.0;
        for (int i = 0; i < ng; ++i) {
            double lin = gb(i, 0);
            for (size_t t = 0; t < rows[i].size(); ++t) if (!istime[rows[i][t].first]) lin += 0.5 * rows[i][t].second;
            worst = std::max(worst, std::fabs(gd(i, 0) - lin));
        }
        if (worst > 1e-6) {
            snprintf(workspace->text, sizeof(workspace->text),
                "\n*** SCIP backend: transcription appears NONLINEAR (residual %.2e); the MILP\n"
                "*** is only a linearization. Use linear dynamics + a linear objective for\n"
                "*** exact mixed-integer results.\n", worst);
            psopt_print(workspace, workspace->text);
        }
    }

    // ---- integrality: map declared integer controls/params to NLP indices ------
    std::vector<char> isInt(nx, 0);
    int xoff = 0;
    for (int p = 0; p < problem.nphases; ++p) {
        int norder = problem.phase[p].current_number_of_intervals;
        int nc = problem.phase[p].ncontrols, ns = problem.phase[p].nstates;
        MatrixXd& ic = problem.phase[p].integer_controls;
        for (int r = 0; r < ic.size(); ++r) {
            int cj = (int) ic(r);
            for (int k = 0; k <= norder; ++k) { int idx = xoff + k*nc + cj; if (idx >= 0 && idx < nx) isInt[idx] = 1; }
        }
        MatrixXd& ip = problem.phase[p].integer_parameters;
        int pbase = xoff + (nc + ns)*(norder + 1);
        for (int r = 0; r < ip.size(); ++r) { int pj = (int) ip(r); int idx = pbase + pj; if (idx >= 0 && idx < nx) isInt[idx] = 1; }
        xoff += get_nvars_phase_i(problem, p, workspace);
    }

    // ---- build and solve the SCIP MILP -----------------------------------------
    SCIP* scip = NULL;
    SCIP_CALL( SCIPcreate(&scip) );
    SCIP_CALL( SCIPincludeDefaultPlugins(scip) );
    SCIP_CALL( SCIPcreateProbBasic(scip, "psopt_miop") );
    SCIP_CALL( SCIPsetObjsense(scip, SCIP_OBJSENSE_MINIMIZE) );
    SCIP_CALL( SCIPsetRealParam(scip, "limits/gap", 1e-9) );

    std::vector<SCIP_VAR*> var(nx);
    char nm[32];
    for (int j = 0; j < nx; ++j) {
        snprintf(nm, sizeof(nm), "x%d", j);
        SCIP_VARTYPE t = isInt[j] ? SCIP_VARTYPE_INTEGER : SCIP_VARTYPE_CONTINUOUS;
        double lo = (*xlb)(j, 0), hi = (*xub)(j, 0);
        if (lo <= -1e19) lo = -SCIPinfinity(scip);
        if (hi >=  1e19) hi =  SCIPinfinity(scip);
        SCIP_CALL( SCIPcreateVarBasic(scip, &var[j], nm, lo, hi, c[j], t) );
        SCIP_CALL( SCIPaddVar(scip, var[j]) );
    }
    for (int i = 0; i < ng; ++i) {
        SCIP_CONS* cons; snprintf(nm, sizeof(nm), "g%d", i);
        // g_l <= b + J(x - x0) <= g_u  =>  g_l - b + J*x0 <= J*x <= g_u - b + J*x0
        double lhs = gl[i] - gb(i, 0) + Jx0[i], rhs = gu[i] - gb(i, 0) + Jx0[i];
        if (lhs <= -1e19) lhs = -SCIPinfinity(scip);
        if (rhs >=  1e19) rhs =  SCIPinfinity(scip);
        std::vector<SCIP_VAR*> vs; std::vector<double> cs;
        for (size_t t = 0; t < rows[i].size(); ++t) { vs.push_back(var[rows[i][t].first]); cs.push_back(rows[i][t].second); }
        SCIP_CALL( SCIPcreateConsBasicLinear(scip, &cons, nm, (int) vs.size(),
                                             vs.empty() ? NULL : vs.data(),
                                             cs.empty() ? NULL : cs.data(), lhs, rhs) );
        SCIP_CALL( SCIPaddCons(scip, cons) );
        SCIP_CALL( SCIPreleaseCons(scip, &cons) );
    }

    SCIP_CALL( SCIPsolve(scip) );

    SCIP_SOL* sol = SCIPgetBestSol(scip);
    if (sol != NULL) {
        for (int j = 0; j < nx; ++j) (*x0)(j, 0) = SCIPgetSolVal(scip, sol, var[j]);
        snprintf(workspace->text, sizeof(workspace->text),
                 "\n>>> SCIP mixed-integer solve: %s, objective (scaled) = %.8g\n",
                 SCIPgetStatus(scip) == SCIP_STATUS_OPTIMAL ? "optimal" : "feasible",
                 SCIPgetSolOrigObj(scip, sol));
        psopt_print(workspace, workspace->text);
    } else {
        snprintf(workspace->text, sizeof(workspace->text),
                 "\n*** SCIP: no feasible mixed-integer solution found\n");
        psopt_print(workspace, workspace->text);
    }

    for (int j = 0; j < nx; ++j) SCIP_CALL( SCIPreleaseVar(scip, &var[j]) );
    SCIP_CALL( SCIPfree(&scip) );

    if (lambda) lambda->setZero();   // integer OC has no classical costates; duals not mapped
    return 0;
}

#endif  // USE_SCIP
