/*********************************************************************************************
 SCIP_interface.cxx  --  mixed-integer NLP backend for PSOPT via SCIP.

 Routes PSOPT's transcribed problem to SCIP so that declared-integer controls and
 parameters are solved to GLOBAL optimality (branch-and-bound), which the
 continuous IPOPT/CasADi backends cannot do.

 SCIP cannot branch on PSOPT's opaque ADOL-C nonlinear callbacks, so the backend
 works on LINEARISATIONS of the transcription:

  - LINEAR transcription (linear dynamics/constraints + linear objective -- use
    auxiliary controls for |u| etc., as in a MILP): the constraint map g(x)=Jx+b
    has a CONSTANT Jacobian, recovered by unit perturbations about the feasible
    guess x0 and handed to SCIP as explicit linear constraints. ONE MILP solve is
    globally optimal and exact. Integrality comes from integer_controls /
    integer_parameters.

  - NONLINEAR transcription: a Sequential-Linear-Programming loop -- re-linearise
    about the incumbent, let SCIP re-optimise the MILP within a trust region on the
    continuous variables, iterate; a second phase pins the integer schedule and
    polishes the continuous states. This is a first-order HEURISTIC: it reaches a
    near-feasible integer solution but not guaranteed tight/global feasibility.
    A linearity self-check picks the exact one-shot path vs. the SLP loop.

 Linearising about x0 (not 0) is essential: the transcription is bilinear in the
 free times t0,tf, so x=0 (=> tf=0) degenerately zeroes every integral.

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
    void*      user_data)
{
    Prob& problem = *workspace->problem;
    const int nx = workspace->nvars;
    const int ng = nlp_ncons;
    const double tol = 1e-10;

    // ---- constraint bounds (scaled NLP space, same as IPOPT/CasADi see) --------
    std::vector<double> gl(ng > 0 ? ng : 1), gu(ng > 0 ? ng : 1);
    if (ng > 0) get_constraint_bounds(gl.data(), gu.data(), workspace);

    // ---- time-variable mask (t0,tf: last two of each phase) --------------------
    std::vector<char> istime(nx, 0);
    { int xo = 0; for (int p = 0; p < problem.nphases; ++p) { int nv = get_nvars_phase_i(problem, p, workspace); istime[xo+nv-2]=1; istime[xo+nv-1]=1; xo += nv; } }

    // ---- integrality: map declared integer controls/params to NLP indices ------
    std::vector<char> isInt(nx, 0);
    { int xoff = 0;
      for (int p = 0; p < problem.nphases; ++p) {
        int norder = problem.phase[p].current_number_of_intervals;
        int nc = problem.phase[p].ncontrols, ns = problem.phase[p].nstates;
        MatrixXd& ic = problem.phase[p].integer_controls;
        for (int r = 0; r < ic.size(); ++r) { int cj = (int) ic(r);
            for (int k = 0; k <= norder; ++k) { int idx = xoff + k*nc + cj; if (idx >= 0 && idx < nx) isInt[idx] = 1; } }
        MatrixXd& ip = problem.phase[p].integer_parameters;
        int pbase = xoff + (nc + ns)*(norder + 1);
        for (int r = 0; r < ip.size(); ++r) { int pj = (int) ip(r); int idx = pbase + pj; if (idx >= 0 && idx < nx) isInt[idx] = 1; }
        xoff += get_nvars_phase_i(problem, p, workspace);
      } }

    // Nonlinear-constraint violation of a point (scaled space).
    MatrixXd gtmp(ng > 0 ? ng : 1, 1);
    #define VIOL_OF(X, OUT) do { OUT = 0.0; if (ng>0){ g((X), &gtmp, workspace); \
        for (int _i=0;_i<ng;_i++){ double lo=gl[_i]-gtmp(_i,0), hi=gtmp(_i,0)-gu[_i]; \
            if (lo>OUT) OUT=lo; if (hi>OUT) OUT=hi; } } } while(0)

    // ---- Sequential-Linear-Programming with a trust region ---------------------
    // Each iteration linearises the transcription about the current point xc,
    // g(x) ~= b + J(x-xc), obj ~= f0 + c'(x-xc), and lets SCIP re-optimise the
    // MILP (integers free, continuous non-time vars restricted to a trust region).
    // For a LINEAR transcription the first solve is already exact (one iteration).
    // Base point x0: essential because the raw transcription is bilinear in the
    // free times t0,tf, so xc must carry a feasible (nonzero) tf.
    MatrixXd xc = *x0, best = *x0;
    double bestViol = 1e30, bestObj = 1e30, bestMerit = 1e30;
    const double MERIT_RHO = 1e4, FEASTOL = 1e-6;
    const int    MAXITER = 40;
    double Delta = 5.0;                    // trust-region radius (scaled units)
    bool linear = false;
    int  status_optimal = 0;

    // Two-phase SLP: iterations < PIN_AT search the integer schedule (integers
    // free); afterwards, if still infeasible, the incumbent integer schedule is
    // PINNED and the remaining iterations polish the continuous states (LP masters
    // -> Newton-like, drives the nonlinear defects to feasibility).
    const int PIN_AT = 12;
    bool pinned = false;
    std::vector<double> pinval(nx, 0.0);

    for (int iter = 0; iter < MAXITER; ++iter)
    {
        if (!linear && !pinned && iter >= PIN_AT && bestViol > FEASTOL) {
            pinned = true; xc = best; Delta = 3.0; bestMerit = 1e30;
            for (int j = 0; j < nx; ++j) if (isInt[j]) pinval[j] = std::floor(best(j,0)+0.5);
            snprintf(workspace->text, sizeof(workspace->text),
                     ">>> SLP+SCIP: pinning integer schedule, polishing continuous states\n");
            psopt_print(workspace, workspace->text);
        }
        // ---- linearise at xc -------------------------------------------------
        MatrixXd gb(ng > 0 ? ng : 1, 1);
        if (ng > 0) g(xc, &gb, workspace);
        double f0 = f(xc, workspace);
        std::vector<double> c(nx, 0.0);
        std::vector<std::vector<std::pair<int,double> > > rows(ng > 0 ? ng : 1);
        MatrixXd xj(nx, 1), gj(ng > 0 ? ng : 1, 1);
        for (int j = 0; j < nx; ++j) {
            xj = xc; xj(j, 0) += 1.0;
            c[j] = f(xj, workspace) - f0;
            if (ng > 0) { g(xj, &gj, workspace);
                for (int i = 0; i < ng; ++i) { double v = gj(i,0)-gb(i,0); if (v>tol||v<-tol) rows[i].push_back(std::make_pair(j,v)); } }
        }
        std::vector<double> Jxc(ng > 0 ? ng : 1, 0.0);
        for (int i = 0; i < ng; ++i)
            for (size_t t = 0; t < rows[i].size(); ++t) Jxc[i] += rows[i][t].second * xc(rows[i][t].first, 0);

        // linearity detection on the first pass (hold time vars fixed)
        if (iter == 0) {
            MatrixXd xd = xc, gd(ng > 0 ? ng : 1, 1);
            for (int j = 0; j < nx; ++j) if (!istime[j]) xd(j,0) += 0.5;
            if (ng > 0) g(xd, &gd, workspace);
            double worst = 0.0;
            for (int i = 0; i < ng; ++i) { double lin = gb(i,0);
                for (size_t t = 0; t < rows[i].size(); ++t) if (!istime[rows[i][t].first]) lin += 0.5*rows[i][t].second;
                worst = std::max(worst, std::fabs(gd(i,0)-lin)); }
            linear = (worst < 1e-7);
        }

        // ---- master build+solve, optionally with an SOC constraint shift -----
        // `shift[i]` moves constraint i's bounds by -shift[i]; used for the
        // second-order correction below. (No SCIP_CALL inside a bool lambda.)
        auto solveMaster = [&](const std::vector<double>& shift, MatrixXd& xout, SCIP_STATUS& sout) -> bool {
            SCIP* scip = NULL;
            SCIPcreate(&scip); SCIPincludeDefaultPlugins(scip);
            SCIPcreateProbBasic(scip, "psopt_miop"); SCIPsetObjsense(scip, SCIP_OBJSENSE_MINIMIZE);
            SCIPsetIntParam(scip, "display/verblevel", 0);
            if (linear) { SCIPsetRealParam(scip, "limits/gap", 1e-9); }
            else { SCIPsetRealParam(scip, "limits/gap", 1e-4); SCIPsetRealParam(scip, "limits/time", 30.0); }
            std::vector<SCIP_VAR*> var(nx);
            char nm[32];
            for (int j = 0; j < nx; ++j) {
                snprintf(nm, sizeof(nm), "x%d", j);
                SCIP_VARTYPE ty = isInt[j] ? SCIP_VARTYPE_INTEGER : SCIP_VARTYPE_CONTINUOUS;
                double lo = (*xlb)(j, 0), hi = (*xub)(j, 0);
                if (pinned && isInt[j]) { lo = hi = pinval[j]; }
                else if (!linear && !isInt[j] && !istime[j]) { lo = std::max(lo, xc(j,0)-Delta); hi = std::min(hi, xc(j,0)+Delta); }
                if (lo <= -1e19) lo = -SCIPinfinity(scip);
                if (hi >=  1e19) hi =  SCIPinfinity(scip);
                if (lo > hi) { double m=0.5*(lo+hi); lo=hi=m; }
                SCIPcreateVarBasic(scip, &var[j], nm, lo, hi, c[j], ty); SCIPaddVar(scip, var[j]);
            }
            for (int i = 0; i < ng; ++i) {
                SCIP_CONS* cons; snprintf(nm, sizeof(nm), "g%d", i);
                double sh = shift.empty() ? 0.0 : shift[i];
                double lhs = gl[i] - gb(i, 0) + Jxc[i] - sh, rhs = gu[i] - gb(i, 0) + Jxc[i] - sh;
                if (lhs <= -1e19) lhs = -SCIPinfinity(scip);
                if (rhs >=  1e19) rhs =  SCIPinfinity(scip);
                std::vector<SCIP_VAR*> vs; std::vector<double> cs;
                for (size_t t = 0; t < rows[i].size(); ++t) { vs.push_back(var[rows[i][t].first]); cs.push_back(rows[i][t].second); }
                SCIPcreateConsBasicLinear(scip, &cons, nm, (int) vs.size(),
                                          vs.empty()?NULL:vs.data(), cs.empty()?NULL:cs.data(), lhs, rhs);
                SCIPaddCons(scip, cons); SCIPreleaseCons(scip, &cons);
            }
            SCIPsolve(scip);
            SCIP_SOL* sol = SCIPgetBestSol(scip);
            sout = SCIPgetStatus(scip);
            bool hs = (sol != NULL);
            xout = xc;
            if (hs) for (int j = 0; j < nx; ++j) xout(j,0) = SCIPgetSolVal(scip, sol, var[j]);
            for (int j = 0; j < nx; ++j) SCIPreleaseVar(scip, &var[j]);
            SCIPfree(&scip);
            return hs;
        };

        std::vector<double> noshift;
        SCIP_STATUS sstat; MatrixXd xstar;
        bool haveSol = solveMaster(noshift, xstar, sstat);

        // ---- SECOND-ORDER CORRECTION (SQP-style): the linear step ignores
        //      constraint curvature, which is what stalls first-order SLP. Evaluate
        //      the TRUE constraints at the trial point, shift the model by that
        //      residual, and re-solve -- one extra solve that captures curvature and
        //      drives feasibility much faster. Kept only if it lowers the violation.
        if (haveSol && !linear && ng > 0) {
            MatrixXd gx(ng, 1); g(xstar, &gx, workspace);
            std::vector<double> r(ng);
            for (int i = 0; i < ng; ++i) {
                double gmodel = gb(i, 0);
                for (size_t t = 0; t < rows[i].size(); ++t)
                    gmodel += rows[i][t].second * (xstar(rows[i][t].first,0) - xc(rows[i][t].first,0));
                r[i] = gx(i, 0) - gmodel;         // nonlinear (curvature) residual at xstar
            }
            MatrixXd xsoc; SCIP_STATUS s2;
            if (solveMaster(r, xsoc, s2)) {
                double v1; VIOL_OF(xstar, v1);
                double v2; VIOL_OF(xsoc,  v2);
                if (v2 < v1) xstar = xsoc;        // accept the second-order-corrected step
            }
        }
        if (haveSol) status_optimal = (sstat == SCIP_STATUS_OPTIMAL);

        if (!haveSol) {
            // An INFEASIBLE master means the trust region is too tight to satisfy the
            // linearised constraints -> GROW Delta; any other no-solution -> shrink.
            bool infeas = (sstat == SCIP_STATUS_INFEASIBLE || sstat == SCIP_STATUS_INFORUNBD);
            snprintf(workspace->text, sizeof(workspace->text),
                     ">>> SLP+SCIP iter %2d: master %s -> %s Delta\n", iter,
                     infeas ? "infeasible" : "no solution", infeas ? "growing" : "shrinking");
            psopt_print(workspace, workspace->text);
            if (infeas) { Delta *= 2.5; if (Delta > 200.0) break; }
            else        { Delta *= 0.4; if (Delta < 1e-6) break; }
            continue;
        }

        // ---- evaluate the TRUE nonlinear model at xstar ----------------------
        double viol; VIOL_OF(xstar, viol);
        double objx = f(xstar, workspace);
        double merit = objx + MERIT_RHO * viol;

        if (linear) { best = xstar; bestObj = objx; bestViol = viol; break; }  // exact single shot

        double step = 0.0; for (int j = 0; j < nx; ++j) step = std::max(step, std::fabs(xstar(j,0)-xc(j,0)));

        bool accepted = (merit < bestMerit - 1e-12);
        // Always retain the least-infeasible incumbent (feasibility first).
        if (viol < bestViol - 1e-12 || (viol < FEASTOL && objx < bestObj)) { best = xstar; bestViol = viol; bestObj = objx; }
        snprintf(workspace->text, sizeof(workspace->text),
                 ">>> SLP+SCIP iter %2d: viol=%.2e obj=%.5g step=%.2e Delta=%.2g %s\n",
                 iter, viol, objx, step, Delta, accepted ? "accept" : "reject");
        psopt_print(workspace, workspace->text);
        if (accepted) {                           // accepted step
            bestMerit = merit; xc = xstar;
            Delta = std::min(Delta * 1.5, 50.0);
            if (bestViol < FEASTOL && step < 1e-6) break;   // converged
        } else {                                  // rejected: shrink (gently) and retry from xc
            Delta *= 0.6; if (Delta < 1e-4) break;
        }
    }

    // ---- Outer-Approximation driver (tight nonlinear MINLP) --------------------
    // Alternate: (master) a SCIP MILP over the integers with accumulated OA
    // linearisation cuts g(p)+J_p(x-p) and an objective epigraph, giving an integer
    // schedule + lower bound; (subproblem) pin that schedule and solve the continuous
    // NLP EXACTLY with IPOPT (true nonlinear dynamics) for an upper bound; add its
    // linearisation as a new cut. This finds a FEASIBLE integer schedule even when the
    // SLP incumbent was infeasible. Nested IPOPT is safe now that the SCIP branch
    // allocates the Jacobian-sparsity workspace (see workspace.cxx).
    if (!linear && bestViol > FEASTOL) {
        extern int psopt_ipopt_resolve(Workspace*, void*);
        double v_before = bestViol;
        // linearise f and g at a point p -> gradf (in cpt) and rows (Jp) and gp
        auto linearize = [&](const MatrixXd& p, std::vector<double>& cpt,
                             std::vector<std::vector<std::pair<int,double> > >& rws, MatrixXd& gp, double& f0){
            gp.resize(ng>0?ng:1,1); if (ng>0) g(const_cast<MatrixXd&>(p),&gp,workspace);
            f0 = f(const_cast<MatrixXd&>(p),workspace);
            cpt.assign(nx,0.0); rws.assign(ng>0?ng:1, std::vector<std::pair<int,double> >());
            MatrixXd xj(nx,1), gj(ng>0?ng:1,1);
            for (int j=0;j<nx;++j){ xj=p; xj(j,0)+=1.0; cpt[j]=f(xj,workspace)-f0;
                if (ng>0){ g(xj,&gj,workspace); for(int i=0;i<ng;++i){ double v=gj(i,0)-gp(i,0); if(v>tol||v<-tol) rws[i].push_back(std::make_pair(j,v)); } } } };

        std::vector<MatrixXd> cutPts;                 // points we have linearised at
        std::vector<std::vector<double> > cutC;       // gradf at each
        std::vector<std::vector<std::vector<std::pair<int,double> > > > cutJ;
        std::vector<MatrixXd> cutG; std::vector<double> cutF0;
        { std::vector<double> c0; std::vector<std::vector<std::pair<int,double> > > r0; MatrixXd g0; double f0;
          linearize(best,c0,r0,g0,f0); cutPts.push_back(best); cutC.push_back(c0); cutJ.push_back(r0); cutG.push_back(g0); cutF0.push_back(f0); }

        const int OA_MAX = 8;
        for (int oa = 0; oa < OA_MAX && bestViol > FEASTOL; ++oa) {
            // ---- MASTER: min eta s.t. OA cuts + integrality --------------------
            SCIP* scip=NULL; SCIPcreate(&scip); SCIPincludeDefaultPlugins(scip);
            SCIPcreateProbBasic(scip,"oa_master"); SCIPsetObjsense(scip,SCIP_OBJSENSE_MINIMIZE);
            SCIPsetIntParam(scip,"display/verblevel",0); SCIPsetRealParam(scip,"limits/gap",1e-6);
            SCIPsetRealParam(scip,"limits/time",30.0);
            std::vector<SCIP_VAR*> var(nx); SCIP_VAR* eta; char nm[40];
            for (int j=0;j<nx;++j){ snprintf(nm,sizeof(nm),"x%d",j);
                SCIP_VARTYPE ty = isInt[j]?SCIP_VARTYPE_INTEGER:SCIP_VARTYPE_CONTINUOUS;
                double lo=(*xlb)(j,0), hi=(*xub)(j,0);
                if (lo<=-1e19) lo=-SCIPinfinity(scip); if (hi>=1e19) hi=SCIPinfinity(scip);
                SCIPcreateVarBasic(scip,&var[j],nm,lo,hi,0.0,ty); SCIPaddVar(scip,var[j]); }
            SCIPcreateVarBasic(scip,&eta,"eta",-SCIPinfinity(scip),SCIPinfinity(scip),1.0,SCIP_VARTYPE_CONTINUOUS);
            SCIPaddVar(scip,eta);
            for (size_t k=0;k<cutPts.size();++k){ const MatrixXd& p=cutPts[k];
                // objective epigraph:  eta - sum gradf_j x_j >= f0 - sum gradf_j p_j
                { std::vector<SCIP_VAR*> vs; std::vector<double> cs; double rhs=cutF0[k];
                  vs.push_back(eta); cs.push_back(1.0);
                  for (int j=0;j<nx;++j){ if (cutC[k][j]>tol||cutC[k][j]<-tol){ vs.push_back(var[j]); cs.push_back(-cutC[k][j]); rhs -= cutC[k][j]*p(j,0);} }
                  SCIP_CONS* c; snprintf(nm,sizeof(nm),"obj%zu",k);
                  SCIPcreateConsBasicLinear(scip,&c,nm,(int)vs.size(),vs.data(),cs.data(),rhs,SCIPinfinity(scip));
                  SCIPaddCons(scip,c); SCIPreleaseCons(scip,&c); }
                // constraint cuts:  gl <= g(p)+Jp(x-p) <= gu
                for (int i=0;i<ng;++i){ double Jp=0.0; for (size_t t=0;t<cutJ[k][i].size();++t) Jp+=cutJ[k][i][t].second*p(cutJ[k][i][t].first,0);
                    double lhs=gl[i]-cutG[k](i,0)+Jp, rhs=gu[i]-cutG[k](i,0)+Jp;
                    if (lhs<=-1e19) lhs=-SCIPinfinity(scip); if (rhs>=1e19) rhs=SCIPinfinity(scip);
                    std::vector<SCIP_VAR*> vs; std::vector<double> cs;
                    for (size_t t=0;t<cutJ[k][i].size();++t){ vs.push_back(var[cutJ[k][i][t].first]); cs.push_back(cutJ[k][i][t].second); }
                    if (vs.empty()) continue;
                    SCIP_CONS* c; snprintf(nm,sizeof(nm),"g%zu_%d",k,i);
                    SCIPcreateConsBasicLinear(scip,&c,nm,(int)vs.size(),vs.data(),cs.data(),lhs,rhs);
                    SCIPaddCons(scip,c); SCIPreleaseCons(scip,&c); } }
            SCIPsolve(scip); SCIP_SOL* sol=SCIPgetBestSol(scip);
            MatrixXd xmas = best; bool hm=(sol!=NULL);
            if (hm) for (int j=0;j<nx;++j) xmas(j,0)=SCIPgetSolVal(scip,sol,var[j]);
            SCIPreleaseVar(scip,&eta); for (int j=0;j<nx;++j) SCIPreleaseVar(scip,&var[j]); SCIPfree(&scip);
            if (!hm) break;                                    // master infeasible: no schedule left

            // ---- SUBPROBLEM: pin the schedule, solve continuous NLP with IPOPT --
            MatrixXd xlb_s=*xlb, xub_s=*xub;
            for (int j=0;j<nx;++j) if (isInt[j]){ double v=std::floor(xmas(j,0)+0.5); (*xlb)(j,0)=v; (*xub)(j,0)=v; }
            *x0 = xmas;                                        // warm start
            int st = psopt_ipopt_resolve(workspace, user_data);
            *xlb=xlb_s; *xub=xub_s;
            (void) st;
            MatrixXd xnlp = *x0; double vnlp; VIOL_OF(xnlp,vnlp); double onlp=f(xnlp,workspace);
            // The IPOPT subproblem returns the best continuous completion of the pinned
            // schedule; keep it whenever it is less infeasible (or feasible & cheaper).
            if (vnlp < bestViol - 1e-12 || (vnlp < FEASTOL && onlp < bestObj)) { best=xnlp; bestViol=vnlp; bestObj=onlp; }
            // add this NLP point as a new OA cut so the master avoids it next round
            { std::vector<double> ck; std::vector<std::vector<std::pair<int,double> > > rk; MatrixXd gk; double fk;
              linearize(xnlp,ck,rk,gk,fk); cutPts.push_back(xnlp); cutC.push_back(ck); cutJ.push_back(rk); cutG.push_back(gk); cutF0.push_back(fk); }
            snprintf(workspace->text,sizeof(workspace->text),
                     ">>> OA iter %d: NLP subproblem (IPOPT, integers pinned) viol=%.2e obj=%.5g -> best viol %.2e\n",
                     oa, vnlp, onlp, bestViol);
            psopt_print(workspace, workspace->text);
        }
        snprintf(workspace->text, sizeof(workspace->text),
                 ">>> OA driver (SCIP master <-> IPOPT NLP subproblem): violation %.2e -> %.2e\n", v_before, bestViol);
        psopt_print(workspace, workspace->text);
    }

    for (int j = 0; j < nx; ++j) (*x0)(j, 0) = best(j, 0);
    if (linear || bestViol < FEASTOL) {
        snprintf(workspace->text, sizeof(workspace->text),
                 "\n>>> SCIP mixed-integer solve (%s): objective (scaled) = %.8g, max constraint violation = %.2e\n",
                 linear ? (status_optimal ? "MILP, global optimum" : "MILP, feasible")
                        : "SLP+SCIP converged", bestObj, bestViol);
    } else if (bestViol < 2.0e-2) {
        snprintf(workspace->text, sizeof(workspace->text),
                 "\n>>> SCIP mixed-integer solve (SLP+SCIP, NEAR-feasible heuristic): objective (scaled) = %.8g,\n"
                 ">>> max constraint violation = %.2e. First-order SLP; refine with a finer mesh / better guess,\n"
                 ">>> or use a linear transcription for exact results.\n", bestObj, bestViol);
    } else {
        snprintf(workspace->text, sizeof(workspace->text),
                 "\n*** SCIP mixed-integer solve did not reach feasibility (best violation %.2e).\n"
                 "*** For nonlinear dynamics the SLP+SCIP loop is a heuristic; try a better guess,\n"
                 "*** a finer mesh, or a linear transcription.\n", bestViol);
    }
    psopt_print(workspace, workspace->text);

    #undef VIOL_OF
    if (lambda) lambda->setZero();   // integer OC has no classical costates; duals not mapped
    return 0;
}

#endif  // USE_SCIP
