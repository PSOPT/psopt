//
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



#include "psopt.h"

// Bring std names into this translation unit (formerly leaked via psopt.h).
using namespace std;

// ===========================================================================================
// Nie-Kerrigan flexible-order local basis. On each mesh element the state (and control) is a
// degree-d Lagrange polynomial through d+1 local LGL nodes on [0,1] (endpoints shared between
// elements, so the state is C0). This routine precomputes, on the reference element, the local
// LGL nodes and the Lagrange basis l_r and its derivative l'_r evaluated at the residual GL
// quadrature points s_q, so the residual sweep is a pair of small matrix-vector products per
// element (the basis is element-independent; only the physical width h_e scales the derivative).
// ===========================================================================================
static void build_ir_local_basis(int d, const MatrixXd& gl, MatrixXd& lgl01, MatrixXd& lgl_w,
                                  MatrixXd& Bval, MatrixXd& Bder, Workspace* workspace)
{
    MatrixXd x, w, P, D;
    lglnodes(d, x, w, P, D, workspace);          // d+1 LGL nodes/weights on [-1,1]
    int np = d + 1;
    lgl01.resize(np,1); lgl_w.resize(np,1);
    for (int r=0;r<np;r++) { lgl01(r) = 0.5*(x(r) + 1.0); lgl_w(r) = w(r); }  // map nodes [-1,1]->[0,1]
    for (int a=0;a<np;a++) for (int b=a+1;b<np;b++)               // sort nodes ascending,
        if (lgl01(b) < lgl01(a)) { std::swap(lgl01(a), lgl01(b)); // carrying the weights along
                                   std::swap(lgl_w(a), lgl_w(b)); }

    int m = gl.rows();
    Bval.resize(m,np); Bder.resize(m,np);
    for (int q=0;q<m;q++) {
        double s = gl(q);
        for (int r=0;r<np;r++) {
            double tr  = lgl01(r);
            double val = 1.0;
            for (int j=0;j<np;j++) if (j!=r) val *= (s - lgl01(j))/(tr - lgl01(j));
            Bval(q,r) = val;
            double der = 0.0;                                     // l'_r(s) = sum_k 1/(tr-tk) prod_{j!=r,k}(...)
            for (int k=0;k<np;k++) if (k!=r) {
                double term = 1.0/(tr - lgl01(k));
                for (int j=0;j<np;j++) if (j!=r && j!=k) term *= (s - lgl01(j))/(tr - lgl01(j));
                der += term;
            }
            Bder(q,r) = der;
        }
    }
}

// ===========================================================================================
// Robust-DAIR (Option B) adjoint costate recovery.
//
// For the residual-box solve the dynamics are imposed by the box |xdot-f| <= delta rather than
// by defect equalities, so the defect multipliers are ~0 and carry no costate information (the
// box inequality multipliers are sensitivities to delta, not the discrete adjoint, and were
// found to be structurally unusable). The costates are instead recovered by post-processing:
// integrate the adjoint
//        lambda_dot = -( dL/dx + (df/dx)^T lambda )
// backward along the recovered primal from the transversality boundary value
//        lambda(tf) = dphi/dxf + sum_e nu_e de_e/dxf,
// where phi is the Mayer cost, e the event vector and nu_e its multipliers (already recovered
// for the box solve, since the box does not touch the event constraints). This single boundary
// formula covers the free-terminal (dphi/dxf = 0, no terminal event -> lambda(tf)=0), Mayer
// (dphi/dxf) and fixed-terminal (event xf - x_f = 0 -> lambda(tf) = dphi/dxf + nu) cases without
// classifying them, and mixed events through de_e/dxf. All Jacobians (df/dx, dL/dx, dphi/dxf,
// de_e/dxf) are formed by central finite differences on the user functions (self-contained, no
// dependence on the AD constraint tape); the backward sweep is RK4 on a sub-grid per interval
// with the primal x(t), u(t) linearly interpolated between nodes.
// ===========================================================================================
// ===========================================================================================
// The solved element widths, written back into the stored mesh.
//
// snodes is the source of truth for everything that does not need a derivative: the reported
// node times, the local error estimator, the plots, and the interpolation that hot-starts the
// next mesh. All of it reads doubles, and all of it keeps working unchanged -- provided the
// doubles describe the mesh that was actually solved on. Left at the uniform mesh the phase
// started from, every node of a moved element is reported at a time it does not have, and the
// error estimator measures the residual of the solved trajectory against a mesh nobody solved
// on. That is the second half of the flexible mesh: the widths drive snodes, and snodes drives
// everything else.
//
// Called once per NLP solve, on the returned primal, before anything reads the solution.
// ===========================================================================================
static void ir_write_back_snodes(MatrixXd& x, Prob& problem, Alg& algorithm, Workspace* workspace)
{
    for (int i = 0; i < problem.nphases; i++) {

        const int norder = problem.phase[i].current_number_of_intervals;
        const int nflex  = ir_flex_mesh_vars(norder, algorithm);
        if ( nflex == 0 ) continue;

        const int d = algorithm.ir_local_order;
        const int M = norder/d;

        const int iphase_offset = get_iphase_offset(problem, i+1, workspace);
        const int nvars_phase_i = get_nvars_phase_i(problem, i, workspace);
        const int base          = iphase_offset + nvars_phase_i - 2 - nflex;

        MatrixXd& sn    = workspace->snodes[i];
        MatrixXd& lgl01 = workspace->ir_lgl01;

        double a = -1.0;
        for (int e = 0; e < M; e++) {
            const double h = x(base+e);        // the widths carry no scale factor, by design
            for (int r = 0; r < d; r++) sn(e*d + r) = a + lgl01(r)*h;
            a += h;
        }
        // The sum equality holds to the NLP's tolerance, not exactly, so the last node is
        // pinned rather than accumulated to. The phase ends at tf; a reported final time that
        // missed it by the constraint violation would be a worse answer than the solver gave.
        sn(norder) = 1.0;
    }
}

static void recover_costates_adjoint(Prob& problem, Alg& algorithm, Sol& solution, Workspace* workspace)
{
    const double fd   = 1.0e-6;   // central-difference step
    const int    Nsub = 20;       // RK4 sub-steps per mesh interval

    for (int i=0; i<problem.nphases; i++) {
        int iphase    = i+1;
        int iph       = (problem.multi_segment_flag || workspace->auto_linked_flag) ? 1 : iphase;
        int nstates   = problem.phase[i].nstates;
        int ncontrols = problem.phase[i].ncontrols;
        int nparam    = problem.phase[i].nparameters;
        int nevents   = problem.phase[i].nevents;
        int norder    = problem.phase[i].current_number_of_intervals;

        adouble* st  = workspace->states[i].get();
        adouble* ct  = workspace->controls[i].get();
        adouble* pr  = workspace->parameters[iph-1].get();
        adouble* dv  = workspace->derivatives[i].get();
        adouble* pth = workspace->path[i].get();
        for (int l=0;l<nparam;l++) pr[l] = (solution.parameters[i])(l);

        double t0 = (solution.nodes[i])(0);
        double tf = (solution.nodes[i])(norder);

        // ---- point evaluators on the user functions (value extraction) ----
        auto eval_f = [&](const double* x,const double* u,double t,double* fout){
            for(int l=0;l<nstates;l++)   st[l]=x[l];
            for(int l=0;l<ncontrols;l++) ct[l]=u[l];
            adouble tt=t;
            problem.dae(dv,pth,st,ct,pr,tt,solution.xad,iphase,workspace);
            for(int j=0;j<nstates;j++) fout[j]=dv[j].value();
        };
        auto eval_L = [&](const double* x,const double* u,double t)->double{
            // A Mayer-only performance index leaves integrand_cost null (see the
            // zero_cost_integrand flag set in psopt_main); the running cost is then zero.
            if ( problem.integrand_cost == NULL ) return 0.0;
            for(int l=0;l<nstates;l++)   st[l]=x[l];
            for(int l=0;l<ncontrols;l++) ct[l]=u[l];
            adouble tt=t;
            return problem.integrand_cost(st,ct,pr,tt,solution.xad,iphase,workspace).value();
        };

        std::vector<adouble> as0(nstates), asf(nstates), aev(nevents>0?nevents:1);
        auto eval_phi = [&](const double* X0,const double* XF)->double{
            for(int l=0;l<nstates;l++){ as0[l]=X0[l]; asf[l]=XF[l]; }
            adouble at0=t0, atf=tf;
            return problem.endpoint_cost(as0.data(),asf.data(),pr,at0,atf,solution.xad,iphase,workspace).value();
        };
        auto eval_ev = [&](const double* X0,const double* XF,double* eout){
            for(int l=0;l<nstates;l++){ as0[l]=X0[l]; asf[l]=XF[l]; }
            adouble at0=t0, atf=tf;
            problem.events(aev.data(),as0.data(),asf.data(),pr,at0,atf,solution.xad,iphase,workspace);
            for(int e=0;e<nevents;e++) eout[e]=aev[e].value();
        };

        // ---- transversality: lambda(tf) ----
        std::vector<double> x0v(nstates), xf(nstates);
        for(int l=0;l<nstates;l++){ x0v[l]=(solution.states[i])(l,0); xf[l]=(solution.states[i])(l,norder); }

        MatrixXd lam(nstates, norder+1);
        std::vector<double> lamf(nstates,0.0), xfp(nstates), xfm(nstates);
        for(int l=0;l<nstates;l++){
            xfp=xf; xfm=xf; xfp[l]+=fd; xfm[l]-=fd;
            lamf[l] += ( eval_phi(x0v.data(),xfp.data()) - eval_phi(x0v.data(),xfm.data()) )/(2*fd);
        }
        if (nevents>0){
            std::vector<double> ep(nevents), em(nevents);
            for(int l=0;l<nstates;l++){
                xfp=xf; xfm=xf; xfp[l]+=fd; xfm[l]-=fd;
                eval_ev(x0v.data(),xfp.data(),ep.data());
                eval_ev(x0v.data(),xfm.data(),em.data());
                for(int e=0;e<nevents;e++)
                    lamf[l] += (solution.dual.events[i])(e) * ( ep[e]-em[e] )/(2*fd);
            }
        }
        for(int l=0;l<nstates;l++) lam(l,norder)=lamf[l];

        // ---- adjoint RHS: lambda_dot = -( dL/dx + (df/dx)^T lambda ) ----
        auto rhs = [&](double t,const double* x,const double* u,const double* L_,double* ld){
            std::vector<double> xp(nstates),xm(nstates),fp(nstates),fm(nstates),g(nstates);
            for(int l=0;l<nstates;l++){
                for(int q=0;q<nstates;q++){xp[q]=x[q];xm[q]=x[q];}
                xp[l]+=fd; xm[l]-=fd;
                g[l]=( eval_L(xp.data(),u,t)-eval_L(xm.data(),u,t) )/(2*fd);   // dL/dx_l
            }
            std::vector<double> JTl(nstates,0.0);
            for(int l=0;l<nstates;l++){
                for(int q=0;q<nstates;q++){xp[q]=x[q];xm[q]=x[q];}
                xp[l]+=fd; xm[l]-=fd;
                eval_f(xp.data(),u,t,fp.data());
                eval_f(xm.data(),u,t,fm.data());
                // column l of df/dx is d f / d x_l; accumulate (J^T lambda)_l = sum_j (df_j/dx_l) lambda_j
                for(int j=0;j<nstates;j++)
                    JTl[l] += ( (fp[j]-fm[j])/(2*fd) ) * L_[j];
            }
            for(int l=0;l<nstates;l++) ld[l] = -( g[l] + JTl[l] );
        };

        // ---- backward RK4, interval by interval ----
        std::vector<double> L0(nstates),xa(nstates),ua(ncontrols),
                            k1(nstates),k2(nstates),k3(nstates),k4(nstates),tmp(nstates);
        for(int k=norder-1;k>=0;k--){
            double ta=(solution.nodes[i])(k), tb=(solution.nodes[i])(k+1);
            for(int l=0;l<nstates;l++) L0[l]=lam(l,k+1);
            auto interp=[&](double t){
                double w=(tb>ta)?(t-ta)/(tb-ta):0.0;
                for(int l=0;l<nstates;l++)   xa[l]=(solution.states[i])(l,k)+w*((solution.states[i])(l,k+1)-(solution.states[i])(l,k));
                for(int l=0;l<ncontrols;l++) ua[l]=(solution.controls[i])(l,k)+w*((solution.controls[i])(l,k+1)-(solution.controls[i])(l,k));
            };
            double dt=(tb-ta)/Nsub, h=-dt;
            for(int s=0;s<Nsub;s++){
                double t1=tb - s*dt;
                interp(t1);              rhs(t1,      xa.data(),ua.data(),L0.data(),k1.data());
                for(int l=0;l<nstates;l++) tmp[l]=L0[l]+0.5*h*k1[l];
                interp(t1+0.5*h);        rhs(t1+0.5*h,xa.data(),ua.data(),tmp.data(),k2.data());
                for(int l=0;l<nstates;l++) tmp[l]=L0[l]+0.5*h*k2[l];
                                         rhs(t1+0.5*h,xa.data(),ua.data(),tmp.data(),k3.data());
                for(int l=0;l<nstates;l++) tmp[l]=L0[l]+h*k3[l];
                interp(t1+h);            rhs(t1+h,    xa.data(),ua.data(),tmp.data(),k4.data());
                for(int l=0;l<nstates;l++) L0[l]+=(h/6.0)*(k1[l]+2*k2[l]+2*k3[l]+k4[l]);
            }
            for(int l=0;l<nstates;l++) lam(l,k)=L0[l];
        }

        solution.dual.costates[i] = lam;

        // recompute the Hamiltonian with the recovered costates (Xdot[i] holds f)
        MatrixXd Temp1 = solution.dual.costates[i].cwiseProduct(workspace->Xdot[i]);
        solution.dual.Hamiltonian[i] = solution.integrand_cost[i] + sum_columns(Temp1);
    }
}




// ---------------------------------------------------------------------------------
// Gauss: put the terminal point into the reported solution.
//
// The Gauss (Legendre-Gauss) scheme collocates strictly interior points. PSOPT stores
// norder+1 nodes per phase -- the initial breakpoint plus the norder Gauss points --
// and x(+1) is an appended NLP variable, which is where the event constraints are
// correctly imposed. It was not, however, in anything the solution accessors returned,
// so get_states_in_phase, get_controls_in_phase and get_time_in_phase handed back a
// trajectory that stopped at the last Gauss node.
//
// That is a long way short. The largest Legendre-Gauss node on 40 points is
// tau = 0.99814738, so on the linear tangent steering problem of the book the reported
// trajectory ended 0.39 s and 3 km early and its last state read y = 407.999044,
// vx = 7.653981, vy = 4.919e-3 against required terminal values of 408, 7.66 and 0.
// That looks exactly like a converged-to-the-wrong-answer failure and is nothing of the
// kind: solution.cost agreed with the Radau run to eight decimals, and integrating the
// returned control on to tf reproduced the Radau endpoint. On a three-interval hp mesh
// it is worse -- the last stored node sits at tau = 0.98809 -- because the final
// interval carries fewer Gauss points.
//
// Radau has the same non-collocated terminal point and does report it, so the two
// siblings disagreed inside one library. This is also the same species as the Gauss
// breakpoint controls, which used to be reported as the barrier's artefact rather than
// as the control the dynamics saw.
//
// The terminal point is appended here, after the mesh loop, so that nothing which
// drives the solve can see the wider arrays: the hot start, the mesh refinement, the
// error estimate and the costate recovery have all finished with them, and
// workspace->prev_states and its siblings keep the unaugmented copies. Every array
// that solution_diagnostics pairs with the trajectory is widened together with it, or
// its column loops would run off the end of the ones left behind.
//
// The values are exact where an exact value exists. The state is the NLP variable; the
// costate is lambda(+1), already recovered as dual.terminal_costates. The control at
// tau = +1 is not a variable -- it enters no defect -- so it is the Lagrange interpolant
// of the last interval's own collocation controls, which is what the dynamics saw and
// the same construction used for the breakpoints. The path multiplier is extrapolated
// linearly, as the Lobatto endpoints already are, since no path constraint is imposed
// there either.
// ---------------------------------------------------------------------------------
void append_gauss_terminal_point(Prob& problem, Alg& algorithm, Sol& solution,
                                 Workspace* workspace)
{
    if ( algorithm.collocation_method != "Gauss" ) return;
    if ( solution.terminal_states == NULL )        return;

    for (int i = 0; i < problem.nphases; i++) {

        const int  nstates   = problem.phase[i].nstates;
        const int  ncontrols = problem.phase[i].ncontrols;
        const int  npath     = problem.phase[i].npath;
        const int  norder    = problem.phase[i].current_number_of_intervals;
        const long M         = solution.nodes[i].cols();

        if ( M < 2 ) continue;
        if ( solution.terminal_states[i].rows() != nstates ) continue;   // never captured

        const double tf = (*workspace->prev_tf)(i);
        const double t_last = solution.nodes[i](0, M-1);
        // Nothing to do if the terminal point is already the last stored node, which is
        // how a second call -- or a scheme that stores it -- is recognised.
        if ( t_last >= tf - 1.0e-13*(1.0 + std::fabs(tf)) ) continue;

        // ---- the last interval, whose interpolant carries tau = +1 ----
        const MatrixXd& sn = workspace->snodes[i];
        const int K  = hp_mesh_active(problem.phase[i]) ? hp_num_intervals(problem.phase[i]) : 1;
        int s = 0, nj = norder;
        for (int j = 0; j < K; j++) {
            nj = hp_mesh_active(problem.phase[i]) ? hp_interval_order(problem.phase[i], j) : norder;
            if ( j == K-1 ) break;
            s += nj + 1;
        }
        const bool interp_ok = ( nj >= 1 && s + nj <= norder && sn.size() >= norder+1 );

        // ---- nodes ----
        { MatrixXd tmp(1, M+1); tmp.leftCols(M) = solution.nodes[i]; tmp(0,M) = tf;
          solution.nodes[i] = tmp; }

        // ---- states ----
        { MatrixXd tmp(nstates, M+1); tmp.leftCols(M) = solution.states[i];
          tmp.col(M) = solution.terminal_states[i];
          solution.states[i] = tmp; }

        // ---- controls: the last interval's Lagrange interpolant at tau = +1 ----
        if ( ncontrols > 0 ) {
            MatrixXd tmp(ncontrols, M+1); tmp.leftCols(M) = solution.controls[i];
            for (int l = 0; l < ncontrols; l++) {
                double val = (solution.controls[i])(l, M-1);      // fallback: hold
                if ( interp_ok ) {
                    val = 0.0;
                    for (int m = s+1; m <= s+nj; m++) {
                        double wL = 1.0;
                        for (int q = s+1; q <= s+nj; q++) if (q != m) wL *= (1.0 - sn(q))/(sn(m) - sn(q));
                        val += wL*(solution.controls[i])(l, m);
                    }
                }
                tmp(l, M) = val;
            }
            solution.controls[i] = tmp;
        }

        // ---- costates: lambda(+1), already recovered ----
        if ( solution.dual.costates != NULL && solution.dual.costates[i].cols() == M ) {
            MatrixXd tmp(nstates, M+1); tmp.leftCols(M) = solution.dual.costates[i];
            if ( solution.dual.terminal_costates != NULL
                 && solution.dual.terminal_costates[i].rows() == nstates )
                tmp.col(M) = solution.dual.terminal_costates[i];
            else
                tmp.col(M) = solution.dual.costates[i].col(M-1);
            solution.dual.costates[i] = tmp;
        }

        // ---- the running cost and the Hamiltonian at the terminal point ----
        {
            std::vector<adouble> st(std::max(nstates,1)), ct(std::max(ncontrols,1)),
                                 pa(std::max(problem.phase[i].nparameters,1)),
                                 de(std::max(nstates,1)), pth(std::max(npath,1));
            for (int l = 0; l < nstates;   l++) st[l] = (solution.states[i])(l, M);
            for (int c = 0; c < ncontrols; c++) ct[c] = (solution.controls[i])(c, M);
            for (int l = 0; l < problem.phase[i].nparameters; l++)
                pa[l] = (solution.parameters[i])(l);
            adouble tm = tf;
            double L = (problem.integrand_cost)
                ? problem.integrand_cost(&st[0], &ct[0], &pa[0], tm, solution.xad, i+1, workspace).value()
                : 0.0;
            problem.dae(&de[0], &pth[0], &st[0], &ct[0], &pa[0], tm, solution.xad, i+1, workspace);

            if ( solution.integrand_cost != NULL && solution.integrand_cost[i].cols() == M ) {
                MatrixXd tmp(1, M+1); tmp.leftCols(M) = solution.integrand_cost[i];
                tmp(0,M) = L; solution.integrand_cost[i] = tmp;
            }
            if ( solution.dual.Hamiltonian != NULL && solution.dual.Hamiltonian[i].cols() == M ) {
                double H = L;
                for (int j = 0; j < nstates; j++)
                    H += (solution.dual.costates[i])(j, M) * de[j].value();
                MatrixXd tmp(1, M+1); tmp.leftCols(M) = solution.dual.Hamiltonian[i];
                tmp(0,M) = H; solution.dual.Hamiltonian[i] = tmp;
            }
        }

        // ---- path multipliers: no constraint is imposed at tau = +1, so extrapolate ----
        if ( npath > 0 && solution.dual.path != NULL && solution.dual.path[i].cols() == M ) {
            MatrixXd tmp(npath, M+1); tmp.leftCols(M) = solution.dual.path[i];
            const double d = solution.nodes[i](0,M-1) - solution.nodes[i](0,M-2);
            const double w = ( d > 0.0 ) ? (tf - solution.nodes[i](0,M-1))/d : 0.0;
            tmp.col(M) = solution.dual.path[i].col(M-1)
                       + w*( solution.dual.path[i].col(M-1) - solution.dual.path[i].col(M-2) );
            solution.dual.path[i] = tmp;
        }

        // ---- stationarity residual: written by solution_diagnostics over every column ----
        if ( solution.stationarity_residual != NULL && ncontrols > 0 )
            (solution.stationarity_residual[i]).resize(ncontrols, M+1);
    }
}





int psopt(Sol& solution, Prob& problem, Alg& algorithm)
{
    // psopt() is total: no exception may escape it. The Workspace construction and
    // initialize_solution() are inside the try as well, so a failure in set-up or in
    // validate_user_input (which throws ErrorHandler) is reported through
    // solution.error_flag rather than propagating out to the caller's main().
    solution.error_flag = 0;
    solution.error_msg  = "";
    // Safe default for the post-failure accessor/utility policy, in case the failure
    // occurs before initialize_solution() resolves algorithm.on_error.
    solution.on_error_fast = true;

    try {
        // Integer-control outer convexification (no-op unless a phase declares
        // an integer control). Expands the problem here, before the Workspace is
        // sized from it, and restores the user layout on scope exit.
        IntegerControlExpansionGuard psopt_ic_guard(problem);

#ifdef PSOPT_ALLOW_ENV_OVERRIDES
        // Before the Workspace, not merely before the mesh loop. get_max_nodes reads the
        // node schedule differently under "manual" than under "automatic" -- manual takes
        // the last entry the user listed, automatic takes an a-priori growth ceiling -- and
        // every array in the Workspace is sized from the answer. Applying the override after
        // the allocation, as this did, leaves the arrays sized for the schedule that was not
        // used, and a prescribed mesh larger than the automatic ceiling then writes past
        // them. examples/manutec is such a case: it prescribes 20, 30, 40, 60, 80 nodes
        // against an automatic ceiling of 76, and its fifth mesh wrote off the end of
        // workspace->states_traj. glibc absorbs that silently; macOS traps on it, which is
        // how it was found -- IPOPT exiting 133 (SIGTRAP) on the fifth mesh with no message.
        psopt_apply_pre_workspace_environment_overrides(algorithm);
#endif

        unique_ptr<Workspace> workspace_up{ new Workspace{problem, algorithm, solution} };

        initialize_solution(solution, problem, algorithm, workspace_up.get());

        psopt_main(solution, problem, algorithm, workspace_up);
    }
    catch (ErrorHandler& handler)
    {
        solution.error_msg  = handler.error_message;
        solution.error_flag = 1;
    }

    return solution.error_flag;
}


void psopt_main(Sol& solution, Prob& problem, Alg& algorithm,  unique_ptr<Workspace>& workspace_up)
{
// PSOPT:  main algorithm


Workspace* workspace = workspace_up.get();

string startup_message= "\n *******************************************************************************\n * This is PSOPT, an optimal control solver based on pseudospectral and local  *\n * collocation methods, together with large scale nonlinear programming        *";

snprintf(workspace->text,sizeof(workspace->text), "%s %s %s", "\n *******************************************************************************\n * PSOPT release number: ", PSOPT_RELEASE_STRING, "                              *");
string release_message= workspace->text;
snprintf(workspace->text,sizeof(workspace->text), "%s %s %s", "\n * PSOPT build date: ", PSOPT_BUILD_DATE, "                                     *");

string build_date= workspace->text;

string license_notice=  "\n * Copyright (C) 2010-2025  Victor M. Becerra.                                 *\n *                                                                             *\n * This library is free software; you can redistribute it and/or               *\n * modify it under the terms of the GNU Lesser General Public License          *\n * as published by the Free Software Foundation;  version 2.1.                 *\n * This library is distributed in the hope that it will be useful,             *\n * but WITHOUT ANY WARRANTY; without even the implied warranty of              *\n * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU           *\n * Lesser General Public License for more details.                             *\n * You should have received a copy of the GNU Lesser General Public            *\n * License along with this library;                                            *\n * If not please visit http://www.gnu.org/licenses                             *\n *                                                                             *";

string contact_notice=  "\n * The author can be contacted at his email address:    vmbecerra@vmb1.com   *\n *                                                                             *\n *******************************************************************************\n\n";


  PSOPT_extras::tic();

  get_local_time( solution.start_date_and_time );


  int nlp_ncons;
  int nlp_neq;
  int offset = 0;   // accumulated in the per-phase loop before use; init to 0
                    // silences a -Wmaybe-uninitialized false positive (the
                    // compiler can't prove the phase loop always runs).
  int nphases = problem.nphases;

  int number_of_mesh_refinement_iterations = get_number_of_mesh_refinement_iterations(problem,algorithm);


  int i, k;
  int iter_nodes;
  int hotflag = 0;
  int x_phase_offset = 0;
  int lam_phase_offset = 0;



  snprintf(workspace->text,sizeof(workspace->text),"%s",startup_message.c_str());
  psopt_print(workspace,workspace->text);
  snprintf(workspace->text,sizeof(workspace->text),"%s",release_message.c_str());
  psopt_print(workspace,workspace->text);
  snprintf(workspace->text,sizeof(workspace->text),"%s",build_date.c_str());
  psopt_print(workspace,workspace->text);
  snprintf(workspace->text,sizeof(workspace->text),"%s",license_notice.c_str());
  psopt_print(workspace,workspace->text);
  snprintf(workspace->text,sizeof(workspace->text),"%s",contact_notice.c_str());
  psopt_print(workspace,workspace->text);

  validate_user_input(problem,algorithm, workspace);

  PSOPT_extras::SetPrintLevel( algorithm.print_level );


  if (problem.integrand_cost == NULL )  {
      for(i=1;i<=problem.nphases;i++) {
	  problem.phase[i-1].zero_cost_integrand = true;
      }
  }
  
  for (iter_nodes=1; iter_nodes<= number_of_mesh_refinement_iterations; iter_nodes++) {


    workspace->current_mesh_refinement_iteration = iter_nodes;

    MatrixXd& x0     = *workspace->x0;
    MatrixXd& lambda = *workspace->lambda;
    MatrixXd& xlb    = *workspace->xlb;
    MatrixXd& xub    = *workspace->xub;



    if (algorithm.collocation_method=="trapezoidal") {
             workspace->differential_defects = "trapezoidal";
    }

    else if (algorithm.collocation_method == "Hermite-Simpson") {
             workspace->differential_defects = "Hermite-Simpson";
    }

    else if (algorithm.collocation_method == "Radau") {
             workspace->differential_defects = "Radau";
    }

    else if (algorithm.collocation_method == "Gauss") {
             workspace->differential_defects = "Gauss";
    }

    else if (use_global_collocation(algorithm)) {
          if (algorithm.diff_matrix=="standard") {
             workspace->differential_defects = "standard";
	  }

	  else if (algorithm.diff_matrix=="reduced-roundoff") {
             workspace->differential_defects = "reduced-roundoff";
	  }
    }



    // Integrated-residual machinery: record the mode and build the per-interval
    // Gauss-Legendre residual grid once (shared across phases and intervals). The grid is
    // needed both by the integrated-residual transcription (increment 1) and by
    // integrated-residual regularisation (increment 2, ir_regularization>0 under collocation).
    workspace->transcription_method = algorithm.transcription_method;
    if ( algorithm.transcription_method == "integrated-residual"
         || algorithm.ir_regularization > 0.0 ) {
        workspace->ir_m = algorithm.ir_residual_nodes;
        gauss_legendre_unit(workspace->ir_m, workspace->ir_nodes, workspace->ir_weights);
        if ( algorithm.ir_local_order >= 2 )   // Nie-Kerrigan flexible-order local basis (p-refinement)
            build_ir_local_basis(algorithm.ir_local_order, workspace->ir_nodes,
                                  workspace->ir_lgl01, workspace->ir_lgl_w, workspace->ir_Bval, workspace->ir_Bder, workspace);
    }
    else {
        workspace->ir_m = 0;
    }

    if (algorithm.mesh_refinement == "manual" )
    {
			for (i=0; i<nphases; i++)
			{
	  			// The loop bound comes from the length of this vector, so the index is in
	  			// range unless something has changed the option after that bound was
	  			// taken. Eigen does not check the index in a release build, and what
	  			// follows an out-of-range read here is a negative interval count and a
	  			// bad_alloc several functions away, so it is checked.
	  			if ( iter_nodes > (int) problem.phase[i].nodes.cols() ) {
	  			    error_message("manual mesh refinement asked for more iterations than "
	  			                  "there are entries in problem.phases(i).nodes");
	  			}
	  			problem.phase[i].current_number_of_intervals    = ( (int) problem.phase[i].nodes(iter_nodes-1)) -1;
			}
    }

    else if ( hp_auto_active(algorithm) ) {
        // hp-adaptive automatic (Route B, Liu-Hager-Rao ph): the driver owns the mesh
        // schedule. At iteration 1 seed K0=1 interval at order N0 = nodes(0)-1 (the
        // discretisation the user already requests), unless an initial hp mesh was supplied
        // manually. From iteration 2 on, hp_refine_driver (called after the previous solve)
        // has already written hp_breakpoints/hp_orders, so nothing is set here; the N_eff
        // update happens through the hp_mesh_active branch a few lines below.
        if ( iter_nodes == 1 ) {
            for (i=0; i<nphases; i++) {
                if ( !hp_mesh_active(problem.phase[i]) ) {
                    problem.phase[i].hp_breakpoints.resize(0);
                    problem.phase[i].hp_orders.resize(1);
                    problem.phase[i].hp_orders(0) = ((int) problem.phase[i].nodes(0)) - 1;
                }
            }
        }
    }

    else  if (algorithm.mesh_refinement=="automatic" && use_local_collocation(algorithm) ) {
          // Local mesh refinement algorithm by Betts (2001)

	  // Step 3: Estimate primary order for the new mesh (either trapezoidal or Hermite-Simpson)
         if ( iter_nodes == 1 ) {
	    			for (i=0; i<nphases; i++)
					{
		    			problem.phase[i].current_number_of_intervals    = ( (int) problem.phase[i].nodes(iter_nodes-1)) -1;
					}
					if ( workspace->differential_defects != "trapezoidal" && algorithm.switch_order > 0) {
		   	  		workspace->differential_defects = "trapezoidal";
					}
	  		}	

	 		else if (iter_nodes == 2) {

	      	bool equi_error = check_for_equidistributed_error(problem,algorithm,solution);

	      	if (equi_error && algorithm.switch_order>0) {
				 	  workspace->differential_defects = "Hermite-Simpson";
	      	}
	 		}


	 		else if (workspace->differential_defects == "trapezoidal" && iter_nodes> algorithm.switch_order && algorithm.switch_order>0) {

	           workspace->differential_defects = "Hermite-Simpson";

    		}

	  		// Step 4: Estimate order reduction for each interval.

	 		if (iter_nodes>2) {
	      		estimate_order_reduction(problem,algorithm,solution, workspace);
	 		}

	 		else {
	      		zero_order_reduction(problem,algorithm,solution, workspace);
	 		}		

	 		// Step 5: construct new mesh

	 		if (iter_nodes>=2) {

	      		construct_new_mesh(problem,algorithm,solution, workspace);

	 		}



    }

    workspace->trace_f_done = false;

    // hp-adaptive fixed mesh (Route B, increment 1): a phase carrying an explicit
    // multi-interval Radau mesh (hp_orders populated) is discretised as a single
    // effective LGR block of order N_eff = sum(hp_orders) with shared breakpoints.
    // Set current_number_of_intervals = N_eff here, before the NLP var/constraint
    // counts and offsets below, all of which key off it; the composite node/weight/
    // differentiation assembler (lgr_nodes_multi) reproduces the single-block
    // storage layout for this N_eff, so the entire downstream pipeline is reused
    // unchanged. Phases without an hp mesh are untouched (bit-identical).
    for (i=0; i<nphases; i++) {
        if ( hp_mesh_active(problem.phase[i]) ) {
            int Nc = problem.phase[i].hp_orders.sum();              // total collocation points
            int K  = (int) problem.phase[i].hp_orders.size();       // number of intervals
            // Radau (and Legendre/Chebyshev): breakpoints are collocated/shared, storage is
            // Nc collocation + 1 terminal, so current_number_of_intervals = Nc. Gauss collocates
            // strictly interior points, so the K stored breakpoints tau_0..tau_{K-1} are extra
            // non-collocated storage nodes: storage = Nc + K, hence norder = Nc + K - 1.
            if ( algorithm.collocation_method == "Gauss" )
                problem.phase[i].current_number_of_intervals = Nc + K - 1;
            else
                problem.phase[i].current_number_of_intervals = Nc;
        }
    }

    // ---- Robust-DAIR (Option B): at each mesh-refinement level the first solve is the
    // feasibility step (min integral(||xdot-f||^2)); the optimality step (min J s.t. the
    // residual box) follows just after the NLP solve below. Under ir_dair the residual-box
    // block is always present (so the NLP dimension is constant across the two sub-solves and
    // no mid-iteration resize is needed); the feasibility step leaves the box non-binding.
    if ( algorithm.ir_dair && algorithm.transcription_method == "integrated-residual" ) {
        algorithm.ir_objective      = "residual";
        algorithm.ir_residual_bound = 1.0e20;     // box present but non-binding
        workspace->ir_delta_phase.clear();        // feasibility uses the scalar bound (all phases)
    }

    workspace->nvars     = get_number_nlp_vars(problem, workspace);

    nlp_ncons           = get_number_nlp_constraints(problem, workspace);

    workspace->ncons = nlp_ncons;

    resize_workspace_vars(problem,algorithm,solution, workspace);

    resize_solution(solution,problem,algorithm);

    // Compute the nodes for each phase

    if (algorithm.collocation_method == "Legendre") {
    	for(i=0; i<nphases; i++)
    	{
                 if ( hp_mesh_active(problem.phase[i]) ) {
                    // hp multi-interval LGL (shared, collocated breakpoints; Option A): the
                    // composite assembler returns ascending nodes and an M x (M+K-1) D whose
                    // leading M x M block is the positive-convention differentiation matrix
                    // (consumed by mtrx_mul_trans) and whose K-1 trailing columns carry the
                    // interface defects. Ascending output => NO sort_vector/rearrange_vector.
                    lgl_nodes_multi( problem.phase[i].hp_breakpoints, problem.phase[i].hp_orders,
                                     workspace->snodes[i], workspace->w[i], workspace->D[i] );
                 }
                 else {
        	 lglnodes( problem.phase[i].current_number_of_intervals, workspace->snodes[i], workspace->w[i], workspace->P[i], workspace->D[i], workspace);

                sort_vector(workspace->snodes[i],workspace->sindex[i]);

                rearrange_vector(workspace->w[i], workspace->sindex[i] );
                 }

    	}
    }

    else if ( algorithm.collocation_method == "Chebyshev" ) {

	    for(i=0; i<nphases; i++)
    	    {
                 if ( hp_mesh_active(problem.phase[i]) ) {
                    // hp multi-interval Chebyshev (shared collocated breakpoints; Option A), same
                    // layout as hp-LGL. cgl_nodes_multi returns ascending nodes, the M x (M+K-1)
                    // composite D, and the composite Clenshaw-Curtis weights (sum 2). The cost is
                    // then a plain weighted sum - no sqrt(1-x^2) factor. Ascending => no sort.
                    cgl_nodes_multi( problem.phase[i].hp_breakpoints, problem.phase[i].hp_orders,
                                     workspace->snodes[i], workspace->w[i], workspace->D[i] );
                 }
                 else {
         	cglnodes( problem.phase[i].current_number_of_intervals, workspace->snodes[i], workspace->w[i], workspace->D[i], workspace );

                sort_vector(workspace->snodes[i],workspace->sindex[i]);

                rearrange_vector(workspace->w[i], workspace->sindex[i] );

                // single-block Chebyshev cost uses Clenshaw-Curtis weights (integrate f directly,
                // spectral) instead of the legacy CGL pi-weights + sqrt(1-x^2) compensation. The
                // nodes/D from cglnodes are kept; only the quadrature weights are replaced. C-C
                // weights are symmetric, so the sorted (ascending) order needs no rearrangement.
                clenshaw_curtis_weights( (int)workspace->snodes[i].size() - 1, workspace->w[i] );
                 }
            }

    }

    else if ( algorithm.collocation_method == "Radau" ) {

	    for(i=0; i<nphases; i++)
    	    {
                // Radau nodes/weights and rectangular differentiation matrix.
                // snodes are returned ascending (N collocation pts incl. -1, then the
                // terminal +1), so NO sort_vector / rearrange_vector is applied here:
                // unlike lglnodes (descending output), the ordering is already correct
                // and D's rows are aligned to it.
                if ( hp_mesh_active(problem.phase[i]) )
                    // hp multi-interval mesh (shared breakpoints): the composite
                    // assembler returns the same (N_eff+1) storage layout as a single
                    // block of order N_eff = sum(hp_orders), block-structured D.
                    lgr_nodes_multi( problem.phase[i].hp_breakpoints, problem.phase[i].hp_orders,
                                     workspace->snodes[i], workspace->w[i], workspace->D[i] );
                else
        	    lgr_nodes( problem.phase[i].current_number_of_intervals, workspace->snodes[i], workspace->w[i], workspace->D[i] );
            }

    }

    else if ( algorithm.collocation_method == "Gauss" ) {

	    for(i=0; i<nphases; i++)
    	    {
                // Legendre-Gauss nodes/weights and rectangular differentiation matrix.
                // snodes are returned ascending: { -1 (initial, non-collocated), N Gauss
                // points }. Collocation is at rows 1..norder; row 0 (initial) is zero.
                // No sort_vector / rearrange_vector (ordering already correct).
                if ( hp_mesh_active(problem.phase[i]) )
                    // hp multi-interval Gauss (non-collocated breakpoints): composite storage
                    // [ -1, Gauss_1, tau_1, Gauss_2, ..., tau_{K-1}, Gauss_K ], block-diagonal D.
                    lg_nodes_multi( problem.phase[i].hp_breakpoints, problem.phase[i].hp_orders,
                                    workspace->snodes[i], workspace->w[i], workspace->D[i] );
                else
        	    lg_nodes( problem.phase[i].current_number_of_intervals, workspace->snodes[i], workspace->w[i], workspace->D[i] );
            }

    }

    else if ( ( use_local_collocation(algorithm) && (iter_nodes==1)) || (use_local_collocation(algorithm) && (iter_nodes>1) && (algorithm.mesh_refinement=="manual") )  ) {

	    for(i=0; i<nphases; i++)
    	    {
	        if ( algorithm.ir_local_order >= 2 ) {
	            // Nie-Kerrigan: the snodes are concatenated per-element LGL sub-nodes. With
	            // norder = current_number_of_intervals total sub-intervals and local degree d,
	            // there are M = norder/d elements of d+1 LGL nodes each (endpoints shared).
	            int d = algorithm.ir_local_order;
	            int norder = problem.phase[i].current_number_of_intervals;
	            if ( norder % d != 0 )
	                error_message("ir_local_order: (number of nodes - 1) must be divisible by ir_local_order ");
	            int M = norder / d;
	            MatrixXd& sn = workspace->snodes[i];
	            // A ROW, like every other producer of snodes. This branch resized it to a
	            // column, which nothing noticed because every reader indexes sn(k) and Eigen
	            // does not care -- until construct_new_mesh, which grows the array with
	            // block(0,0,1,cols) and then sorts it, and died in sort_vector() on "argument
	            // must be a column or row vector". A shape that only one consumer can see is
	            // a shape nobody maintains.
	            sn.resize(1,norder+1);
	            MatrixXd& lgl01 = workspace->ir_lgl01;       // d+1 reference LGL nodes on [0,1]
	            double H = 2.0 / (double) M;                 // element width on [-1,1]
	            for (int e=0; e<M; e++) {
	                double a_e = -1.0 + e*H;
	                for (int r=0; r<=d; r++) sn(e*d + r) = a_e + lgl01(r)*H;
	            }
	        }
	        else {
	            workspace->snodes[i] = linspace(-1.0, 1.0, problem.phase[i].current_number_of_intervals+1);
	        }
            }

    }


    // Define initial NLP guess
    if (iter_nodes==1) {
       hotflag = 0;
       define_initial_nlp_guess(x0, lambda, solution, problem, algorithm, workspace);
    }
    else {
        	hotflag = 1;
			hot_start_nlp_guess(x0, lambda, solution,problem,algorithm, workspace->prev_states.get(), workspace->prev_controls.get(), workspace->prev_costates.get(), workspace->prev_path.get(), workspace->prev_nodes.get(), workspace->prev_param.get(), *workspace->prev_t0, *workspace->prev_tf, workspace);
    }

  	// Define NLP bounds on the decision vector

    define_nlp_bounds(*workspace->xlb, *workspace->xub, problem, algorithm, workspace);

    nlp_neq = 0; //   equality constraints

    // Note: the decision variables are:
    // x = [vec(controls)' vec(states)'  parameters' t0 tf ... (this is repeated for each phase)]';
    // where "controls" is a real matrix with dimensions [ncontrols x norder+1]
    //       "states" is a real matrix with dimensions [nstates x norder+1]
    //       "parameters" is a real vector with dimensions [nparameters x 1]
    //       "t0" and "tf" are real numbers
    //       and the vec(.) operator stacks the columns of a matrix one below the next


    snprintf(workspace->text,sizeof(workspace->text),"\nProblem:\t\t\t\t\t\t%s", problem.name.c_str());
    psopt_print(workspace,workspace->text);


    snprintf(workspace->text,sizeof(workspace->text), "\nThis is mesh refinement iteration:\t\t\t%i", iter_nodes);
    psopt_print(workspace,workspace->text);
    if ( use_global_collocation(algorithm) ) {
      snprintf(workspace->text,sizeof(workspace->text), "\nCollocation method:\t\t\t\t\t%s", algorithm.collocation_method.c_str());
      psopt_print(workspace,workspace->text);
    }
    else {
      snprintf(workspace->text,sizeof(workspace->text), "\nCollocation method:\t\t\t\t\t%s", workspace->differential_defects.c_str());
      psopt_print(workspace,workspace->text);
    }
    if ( use_global_collocation(algorithm) ) {
	  snprintf(workspace->text,sizeof(workspace->text), "\nDifferentiation matrix:\t\t\t\t\t%s", algorithm.diff_matrix.c_str());
	  psopt_print(workspace,workspace->text);
    }

    snprintf(workspace->text,sizeof(workspace->text), "\nNumber of NLP variables\t\t\t\t\t%i", workspace->nvars );
    psopt_print(workspace,workspace->text);
    snprintf(workspace->text,sizeof(workspace->text), "\nNumber of NLP nonlinear constraints:\t\t\t%i", nlp_ncons);
    psopt_print(workspace,workspace->text);
    for(i=1;i<=problem.nphases;i++) {
      snprintf(workspace->text,sizeof(workspace->text), "\nNumber of nodes phase %i:\t\t\t\t%i", i,problem.phases(i).current_number_of_intervals+1);
      psopt_print(workspace,workspace->text);
    }
    snprintf(workspace->text,sizeof(workspace->text),"\n");
    psopt_print(workspace,workspace->text);

    // Once per mesh iteration, before the NLP sees the problem: does anything in gg_ad
    // write every row that get_ncons_phase_i counted? Three defects have been exactly this
    // disagreement, and none of them announced itself, because the constraint buffer is
    // zero-filled and zero is a plausible constraint value. One extra constraint evaluation
    // per mesh; see the note above check_constraint_coverage in NLP_constraints.cxx.
    check_constraint_coverage( x0, nlp_ncons, workspace );

    workspace->enable_nlp_counters = true;

    chronometer_tic(workspace);

    NLP_interface( algorithm, &x0,  ff_num, gg_num, nlp_ncons,  nlp_neq , &xlb, &xub, &lambda, hotflag, 1, workspace, problem.user_data   );

    // Where the elements ended up. Nothing downstream of here may read snodes before this.
    ir_write_back_snodes( x0, problem, algorithm, workspace );

    solution.mesh_stats[workspace->current_mesh_refinement_iteration-1].CPU_time = chronometer_toc(workspace);

    workspace->enable_nlp_counters = false;

    // ---- Robust-DAIR (Option B) optimality sub-solve ----
    // The solve above was the feasibility step (box non-binding). Read the mesh scale h from
    // the feasibility primal, set the box tolerance from it, and minimise the user cost J subject
    // to |(xdot-f)_{g,q,j}| <= delta, warm-started from the feasibility primal (x0) and duals
    // (lambda). The NLP dimension is unchanged (the box block is already present), so only the
    // objective and the box bounds change. As the mesh-refinement loop refines, h and delta fall
    // together and J converges to the optimum without penalty-continuation divergence.
    //
    // The schedule is order-aware. The residual of a degree-p local representation behaves like
    // h^p (the derivative error of a degree-p interpolant), so the box tolerance is matched to
    // that order: delta_p = ir_dair_delta_factor * h_node^expo, where h_node = Tlen/norder is the
    // mesh node spacing and expo is the representation order:
    //   cubic-Hermite (ir_local_order==0): expo = 2  (the legacy K*h^2 rule, unchanged);
    //   Nie-Kerrigan  (ir_local_order>=2): expo = d = ir_local_order.
    // The node spacing (not the element width Tlen/M = d*h_node) is the right scale empirically and
    // theoretically: the degree-d interpolation constant ~1/(d+1)! folds the d^d from the element
    // width back to O(1), so K*h_node^d tracks the achievable residual whereas K*h_element^d is far
    // too loose at coarse meshes. Each phase gets its own delta_p (ir_delta_phase); the scalar
    // ir_residual_bound is set to the largest delta_p for reporting and gating.
    if ( algorithm.ir_dair && algorithm.transcription_method == "integrated-residual" ) {
        copy_decision_variables(solution, x0, problem, algorithm, workspace);   // feasibility nodes
        workspace->ir_delta_phase.assign(nphases, 0.0);
        int    lo   = algorithm.ir_local_order;
        int    expo = ( lo >= 2 ) ? lo : 2;
        double delta_max = 0.0;
        for (i=0; i<nphases; i++) {
            int norder  = problem.phase[i].current_number_of_intervals;
            double Tlen = (solution.nodes[i])(norder) - (solution.nodes[i])(0);
            double h_rel = 1.0/norder;                       // dimensionless node spacing h_node/Tlen
            double d_p = algorithm.ir_dair_delta_factor;     // K*(h_node/Tlen)^expo, scale-invariant
            for (int e=0; e<expo; e++) d_p *= h_rel;          // horizon-normalised: no T^expo blow-up
            (void) Tlen;
            workspace->ir_delta_phase[i] = d_p;
            delta_max = std::max(delta_max, d_p);
        }
        algorithm.ir_objective      = "cost";
        algorithm.ir_residual_bound = delta_max;
        define_nlp_bounds(*workspace->xlb, *workspace->xub, problem, algorithm, workspace);
        workspace->trace_f_done = false;   // re-tape the objective gradient: residual -> cost J
        if (nphases == 1)
            snprintf(workspace->text,sizeof(workspace->text),
                     "\nDAIR optimality sub-solve: min J s.t. |r| <= %e  (delta = K*(h/T)^%d)", delta_max, expo);
        else
            snprintf(workspace->text,sizeof(workspace->text),
                     "\nDAIR optimality sub-solve: min J s.t. per-phase |r| <= delta_p (max %e, K*(h/T)^%d)", delta_max, expo);
        psopt_print(workspace,workspace->text);
        workspace->enable_nlp_counters = true;
        NLP_interface( algorithm, &x0, ff_num, gg_num, nlp_ncons, nlp_neq, &xlb, &xub, &lambda, 1, 1, workspace, problem.user_data );
        workspace->enable_nlp_counters = false;
        ir_write_back_snodes( x0, problem, algorithm, workspace );
    }

    // Copy the resultant decision vector into the relevant solution variables.

    copy_decision_variables(solution, x0, problem, algorithm, workspace);

    // Radau: the terminal node is non-collocated and carries no optimised control.
    // Pin the reported terminal control to the Lagrange interpolant of the collocation
    // controls evaluated at tau = +1 (the terminal node). hp multi-interval mesh: tau=+1
    // belongs to the LAST interval, so the interpolant must be built on that interval's
    // collocation nodes sn(m0..norder-1) with m0 = N_eff - n_last. A global interpolant
    // over all clustered multi-interval nodes is severely ill-conditioned (Sum|L| ~ 1e9)
    // and amplifies per-node noise into the terminal value. For a single interval m0 = 0,
    // so this reduces exactly to the legacy global interpolant (bit-identical). This mirrors
    // the interval-aware terminal-control pin in NLP_constraints.cxx.
    if ( algorithm.collocation_method == "Radau" ) {
        for(int ip=0; ip<nphases; ip++) {
            int ncontrols = problem.phase[ip].ncontrols;
            int norder    = problem.phase[ip].current_number_of_intervals;
            MatrixXd& sn  = workspace->snodes[ip];
            double xe = sn(norder);                       // terminal node, tau = +1
            int m0 = 0;                                   // first node of the interpolating interval
            if ( hp_mesh_active(problem.phase[ip]) )
                m0 = norder - hp_interval_order(problem.phase[ip], hp_num_intervals(problem.phase[ip]) - 1);
            MatrixXd Lw(norder,1); Lw.setZero();
            for(int m=m0;m<norder;m++) {                  // Lagrange weights at xe over the
                double wL = 1.0;                          // last interval's nodes sn(m0..norder-1)
                for(int j=m0;j<norder;j++) if(j!=m) wL *= (xe - sn(j))/(sn(m)-sn(j));
                Lw(m) = wL;
            }
            for(int l=0;l<ncontrols;l++) {
                double val = 0.0;
                for(int m=m0;m<norder;m++) val += Lw(m)*(solution.controls[ip])(l,m);
                (solution.controls[ip])(l,norder) = val;
            }
        }
    }

    // Gauss: the same problem as Radau's terminal node, at the other end and once per mesh
    // interval. Gauss collocates strictly interior points, so each interval's left breakpoint
    // is a stored node carrying a control that appears in no defect constraint. The variable
    // exists in the NLP but is unconstrained, and the barrier term alone decides it: it comes
    // back at the midpoint of the control bounds, which for a control bounded in [0,1] is a
    // reported value of exactly 0.5 at the start of every interval. Those are artefacts, and
    // they plot as spikes in a control trajectory that is otherwise correct.
    //
    // Report instead the control the dynamics actually saw: the Lagrange interpolant of that
    // interval's own collocation controls, evaluated at the breakpoint. Built per interval,
    // for the reason given above for Radau -- an interpolant over the clustered nodes of a
    // whole multi-interval mesh is severely ill-conditioned. This changes only what is
    // reported (and what the next mesh is hot-started from); it does not change the NLP.
    if ( algorithm.collocation_method == "Gauss" ) {
        for(int ip=0; ip<nphases; ip++) {
            int ncontrols = problem.phase[ip].ncontrols;
            if ( ncontrols < 1 ) continue;
            int norder    = problem.phase[ip].current_number_of_intervals;
            MatrixXd& sn  = workspace->snodes[ip];
            int K = hp_mesh_active(problem.phase[ip]) ? hp_num_intervals(problem.phase[ip]) : 1;
            int s = 0;                                   // storage index of interval's breakpoint
            for (int j=0; j<K; j++) {
                int nj = hp_mesh_active(problem.phase[ip])
                         ? hp_interval_order(problem.phase[ip], j) : norder;
                if ( nj < 1 || s + nj > norder ) break;  // defensive: leave the node as it is
                double xe = sn(s);                       // the non-collocated breakpoint
                for(int l=0;l<ncontrols;l++) {
                    double val = 0.0;
                    for(int m=s+1;m<=s+nj;m++) {         // interval's collocation nodes
                        double wL = 1.0;
                        for(int q=s+1;q<=s+nj;q++) if(q!=m) wL *= (xe - sn(q))/(sn(m)-sn(q));
                        val += wL*(solution.controls[ip])(l,m);
                    }
                    (solution.controls[ip])(l,s) = val;
                }
                s += nj + 1;
            }
        }
    }

    solution.cost = ff_num(x0, workspace)/problem.scale.objective;

    // Gauss: the running cost at each interval's left breakpoint has to be recomputed,
    // for the same reason the control there had to be. The loop above replaced the
    // barrier's artefact with the interval's own interpolant, but ff_num immediately
    // re-evaluates the objective from the decision vector and puts the artefact's value
    // straight back into solution.integrand_cost -- which is the array the Hamiltonian is
    // built from, so H at every Gauss breakpoint was L(u_artefact) + lambda.f rather than
    // L(u) + lambda.f. On examples/mineng_di, whose control is bounded in [-50,50] so the
    // unconstrained variable comes back at 0 and the running cost u^2/2 with it, that made
    // the reported H(t_0) = -36 against a true -18: a Hamiltonian that is supposed to be
    // constant, and is, reading as though it were not.
    //
    // Only the breakpoints are touched, and only in the reported array. The objective is
    // solution.integrated_cost, which NLP_objective forms from the quadrature sum and not
    // from this array, and the Gauss weight at a breakpoint is zero in any case, so no cost
    // and no derivative moves.
    if ( algorithm.collocation_method == "Gauss" ) {
        for(int ip=0; ip<nphases; ip++) {
            if ( problem.phase[ip].zero_cost_integrand || problem.integrand_cost == NULL ) continue;
            const int nstates   = problem.phase[ip].nstates;
            const int ncontrols = problem.phase[ip].ncontrols;
            const int nparam    = problem.phase[ip].nparameters;
            const int norder    = problem.phase[ip].current_number_of_intervals;
            if ( solution.integrand_cost[ip].cols() < norder+1 ) continue;
            const int K = hp_mesh_active(problem.phase[ip]) ? hp_num_intervals(problem.phase[ip]) : 1;
            std::vector<adouble> st(std::max(nstates,1)), ct(std::max(ncontrols,1)),
                                 pa(std::max(nparam,1));
            for (int l=0; l<nparam; l++) pa[l] = (solution.parameters[ip])(l);
            int sidx = 0;
            for (int j=0; j<K; j++) {
                const int nj = hp_mesh_active(problem.phase[ip])
                               ? hp_interval_order(problem.phase[ip], j) : norder;
                if ( nj < 1 || sidx > norder ) break;
                for (int l=0; l<nstates;   l++) st[l] = (solution.states[ip])(l, sidx);
                for (int c=0; c<ncontrols; c++) ct[c] = (solution.controls[ip])(c, sidx);
                adouble tb = (solution.nodes[ip])(0, sidx);
                (solution.integrand_cost[ip])(0, sidx) =
                    problem.integrand_cost(&st[0], &ct[0], &pa[0], tb, solution.xad, ip+1, workspace).value();
                sidx += nj + 1;
            }
        }
    }

    snprintf(workspace->text,sizeof(workspace->text),"\nReturned (unscaled) cost function value: %e", solution.cost);
    psopt_print(workspace,workspace->text);
    x_phase_offset   = 0;
    lam_phase_offset = 0;

    for(i=0; i<nphases; i++)
    {
        int nstates   = problem.phase[i].nstates;
        int norder    = problem.phase[i].current_number_of_intervals;
        int nevents   = problem.phase[i].nevents;
        int npath     = problem.phase[i].npath;
        int iphase = i+1;
        MatrixXd& deriv_scaling = problem.phase[i].scale.defects;
        MatrixXd pz(nstates,1);
        MatrixXd pint;
        MatrixXd tint(1,norder);
        MatrixXd ts;
        MatrixXd pl;
        MatrixXd pextra(1,norder+1);
        double hk, t0, tf;

        int nvars_phase_i = get_nvars_phase_i(problem, i, workspace);

        int ncons_phase_i =  get_ncons_phase_i(problem,i, workspace);

        solution.dual.costates[i] = lambda.block(lam_phase_offset,0,nstates*(norder+1),1);

		  solution.dual.costates[i] = reshape(solution.dual.costates[i], nstates, norder+1).eval();

		  workspace->prev_costates[i]      = solution.dual.costates[i];

    	  // The phase's own start and end. t0 read the *second* node rather than the
	     // first, which made (tf-t0) short by exactly one interval and scaled every
	     // local-collocation costate below by (M-1)/M on a mesh of M intervals. The
	     // error vanishes as the mesh refines, which is why it survived: on the
	     // 39-interval mesh of tests/test_costates.cpp it is 2.6 per cent, and by 300
	     // intervals it is a third of one per cent. The pseudospectral branches divide
	     // by the quadrature weights instead and were never affected.
    	  t0 = (solution.nodes[i])(0);

	     tf = (solution.nodes[i])(0, solution.nodes[i].cols()-1); 

	     if ( algorithm.collocation_method=="Legendre" ) {
	    		for(k=0;k<norder+1;k++) { // EIGEN_UPDATE: Index k shifted to start at 0

                   (solution.dual.costates[i]).block(0,k,nstates,1) = (solution.dual.costates[i]).block(0,k,nstates,1)/(workspace->w[i])(k);  // See PhD thesis by Huntington (2006).
	    		}
	     }	

	     if ( algorithm.collocation_method=="Radau" ) {
	    		// Radau covector mapping. Sign confirmed in-situ against an analytic adjoint
	    		// (scalar LQR, lambda(t)=sinh(T-t)/cosh(T)): in PSOPT's stored-dual convention
	    		// the mapping is  lambda_k = + lambda_tilde_k / w_k  at collocation points
	    		// k=0..norder-1 (same convention as Legendre). The terminal node (k=norder)
	    		// is the non-collocated boundary; its costate is carried by the boundary
	    		// multiplier (0 by transversality for a free terminal). w[i](norder)=0, so the
	    		// terminal column is deliberately not divided.
	    		for(k=0;k<norder;k++) {
                   (solution.dual.costates[i]).block(0,k,nstates,1) = (solution.dual.costates[i]).block(0,k,nstates,1)/(workspace->w[i])(k);
	    		}
	     }

	     if ( algorithm.collocation_method=="Gauss" ) {
	    		// Legendre-Gauss covector mapping (Benson canonical / Option 1).
	    		// Collocation is at the interior Gauss points; their costates are
	    		// lambda_k = + lambda_tilde_k / w_k (same convention as Legendre/Radau).
	    		// Non-collocated nodes (single-block: the initial node; multi-interval: every
	    		// stored breakpoint) carry w=0 and no defect, so they are NOT divided here -
	    		// they are recovered from the defining-constraint multipliers below.
	    		for(k=1;k<=norder;k++) {
                   if ( (workspace->w[i])(k) != 0.0 )
                       (solution.dual.costates[i]).block(0,k,nstates,1) = (solution.dual.costates[i]).block(0,k,nstates,1)/(workspace->w[i])(k);
	    		}
	     }

		  if (use_local_collocation(algorithm)) {
        		for(k=0;k<norder;k++) { // EIGEN_UPDATE: Index k shifted to start at 0
	         // For local collocation, this is an estimate of the costate.
	         	hk = (solution.nodes[i])(k+1)-(solution.nodes[i])(k);

	         	(solution.dual.costates[i]).block(0,k,nstates,1) = (solution.dual.costates[i]).block(0,k,nstates,1)/(2.0*hk); // This gives the costates at the midpoints between the collocation nodes.
               (solution.dual.costates[i]).block(0,k,nstates,1)= (solution.dual.costates[i]).block(0,k,nstates,1)*(tf-t0);
 	   		}
		  }		

    	  if ( algorithm.collocation_method == "Chebyshev" ) {
        		// Covector map lambda_k = nu_k / w_k with w_k the Clenshaw-Curtis cost weight.
        		// CGL collocates both endpoints and the C-C endpoint weights are nonzero, so the
        		// map runs over ALL nodes exactly as for Legendre (no endpoint extrapolation).
        		for(k=0;k<norder+1;k++) {
                    (solution.dual.costates[i]).block(0,k,nstates,1) = (solution.dual.costates[i]).block(0,k,nstates,1)/(workspace->w[i])(k);
        		}
    	  }

    	  if ( algorithm.scaling=="user" ) {
	     		for (k=0;k<norder+1;k++) {  // EIGEN_UPDATE: Index k shifted by -1
 //                  (solution.dual.costates[i])(colon(),k)=(solution.dual.costates[i])(colon(),k) & deriv_scaling;
                     (solution.dual.costates[i]).block(0,k,nstates,1) =(solution.dual.costates[i]).block(0,k,nstates,1).cwiseProduct(deriv_scaling);
            }
		  }

   	  solution.dual.costates[i] = solution.dual.costates[i]/problem.scale.objective;

   	  if ( algorithm.scaling=="automatic" ) {
	  			solution.dual.costates[i] = reshape(solution.dual.costates[i], nstates*(norder+1), 1 );
          	solution.dual.costates[i] = solution.dual.costates[i].cwiseProduct((*workspace->constraint_scaling).block(lam_phase_offset,0, nstates*(norder+1),1) );
	  			solution.dual.costates[i] = reshape(solution.dual.costates[i], nstates, norder+1);
   	  }

	  // Legendre/Chebyshev hp: correct the interior-breakpoint costates with the interface-
	  // defect multipliers. An interior Lobatto breakpoint is collocated from both sides; the
	  // primary (left-interval terminal) multiplier divided by the accumulated breakpoint
	  // weight, already stored in costates(bp) by the lambda_k = lambda_tilde_k/w_k loop above,
	  // is only the left-side contribution. The interface (right-interval initial) defect
	  // multiplier mu completes the Lobatto covector mapping at the shared node:
	  //     lambda(bp) = ( nu_primary + mu_interface ) / w_bp .
	  // mu is pushed through the same user/objective/automatic scaling the costate block
	  // received (using the interface constraint's own automatic-scaling entry), then divided
	  // by the accumulated weight. K=1 has no interfaces and is unchanged (bit-identical to
	  // single-block). LGL and CGL share this Option-A interface-defect structure exactly.
	  if ( ( algorithm.collocation_method=="Legendre" || algorithm.collocation_method=="Chebyshev" )
	       && hp_mesh_active(problem.phase[i]) ) {
	      int iface = lam_phase_offset + nstates*(norder+1) + nevents + npath*(norder+1);
	      int Kl = hp_num_intervals(problem.phase[i]);
	      int cbp = 0;
	      for (int e=0; e<Kl-1; e++) {
	          cbp += hp_interval_order(problem.phase[i], e);          // breakpoint storage index c_{e+1}
	          double wbp = (workspace->w[i])(cbp);
	          for (int j=0;j<nstates;j++) {
	              double v = lambda(iface + e*nstates + j);
	              if ( algorithm.scaling=="user" )      v *= deriv_scaling(j);
	              v /= problem.scale.objective;
	              if ( algorithm.scaling=="automatic" ) v *= (*workspace->constraint_scaling)(iface + e*nstates + j);
	              (solution.dual.costates[i])(j,cbp) += v / wbp;
	          }
	      }
	  }

	  // Radau: recover the terminal (non-collocated) costate by Lagrange extrapolation of
	  // the physical collocation costates to tau = +1. Scaling-robust (a pure function of the
	  // already-scaled collocation costates) and reuses the terminal-control-pin weights.
	  // hp multi-interval mesh: tau=+1 belongs to the last interval, so extrapolate only on
	  // that interval's nodes sn(m0..norder-1), m0 = N_eff - n_last (a global interpolant over
	  // clustered multi-interval nodes is severely ill-conditioned). Single interval m0 = 0,
	  // bit-identical to the legacy extrapolation.
	  if ( algorithm.collocation_method=="Radau" ) {
	      MatrixXd& sn = workspace->snodes[i];
	      double xe = sn(norder);
	      int m0 = 0;
	      if ( hp_mesh_active(problem.phase[i]) )
	          m0 = norder - hp_interval_order(problem.phase[i], hp_num_intervals(problem.phase[i]) - 1);
	      for (int r=0;r<nstates;r++) (solution.dual.costates[i])(r,norder) = 0.0;
	      for (int m=m0;m<norder;m++) {
	          double Lwm = 1.0;
	          for (int jj=m0;jj<norder;jj++) if (jj!=m) Lwm *= (xe - sn(jj))/(sn(m)-sn(jj));
	          for (int r=0;r<nstates;r++)
	              (solution.dual.costates[i])(r,norder) += Lwm * (solution.dual.costates[i])(r,m);
	      }
	  }

   	  // (Chebyshev endpoint costates are now mapped directly by the covector loop above -
   	  // the legacy near-endpoint linear-extrapolation workaround, needed only when the old
   	  // pi-weight + 1/sqrt(1-t^2) map blew up at the ends, has been removed.)

    	  if ( algorithm.collocation_method == "Legendre"
    	       || algorithm.collocation_method == "Chebyshev" ) {
                // ------------------------------------------------------------------
                // Smoothing the Lobatto costate, LGL and CGL alike.
                //
                // The covector map lambda_k = nu_k/w_k leaves the costate of a Lobatto
                // scheme with a node-to-node alternating component -- a genuine and well
                // known defect of those schemes, not of this implementation -- and a filter
                // is the usual remedy (Fahroo and Ross, "Costate estimation by a Legendre
                // pseudospectral method", J. Guidance, Control and Dynamics 24(2), 2001).
                // On the linear tangent steering problem at 40 nodes the raw costate
                // oscillates about the true value with an amplitude of 1.7 at the centre of
                // the mesh, growing to 5.5 near the ends, and the two endpoint values are
                // out by 13.7 on a costate whose true value is the constant -38.70.
                //
                // Chebyshev was left out of this for no reason anyone recorded, and it has
                // the same defect in the same shape: on examples/mineng_di, whose adjoint
                // lambda_1 is the constant -12 and which CGL solves exactly, the error at
                // 40 nodes alternates in sign at every node, is smallest at the centre of
                // the mesh and largest at the two ends, reaching 4.0e-6. Applying this
                // filter to it is a uniform improvement wherever it has been measured --
                // a factor of ten on both costates of mineng_di at 10, 20, 40 and 80 nodes,
                // three to ten on a three-interval hp mesh, and 10.3 on the linear tangent
                // steering costates at each of six initial guesses, never once worse.
                // examples/launch is the only shipped example that uses CGL and its
                // objective and mesh history are unchanged, as they must be: the costates
                // are computed after the solve and are not fed back.
                //
                // The filter used to be the fixed stencil (1/4, 1/2, 1/4) at the interior
                // nodes and a plain average of the last two values at each end. Both parts
                // were wrong, in different ways, and the two had to be repaired together.
                //
                // The interior stencil reproduces a linear function only on a UNIFORM grid,
                // and the LGL nodes are not uniform. So a linearly-varying costate came back
                // multiplied by a factor short of one -- 0.9729 on 10 nodes, 0.9935 on 20,
                // the same factor at every interior node -- while a constant costate came
                // back exact. On the minimum-energy double integrator of examples/mineng_di,
                // whose adjoint lambda_2 = 12t - 6 the raw map returns to nine figures and
                // which LGL solves exactly, that put the reported costate out by 1.4e-2 at
                // 40 nodes and drove PSOPT's own stationarity residual dH/du to the same
                // 1.4e-2 on a problem it had solved to machine precision. The error falls
                // like 1/N^2, so it read as discretization error and survived.
                //
                // The endpoint average is the natural way to kill an alternating mode at a
                // boundary -- averaging two successive extremes of the oscillation does
                // exactly that -- but it returns the value at the MIDPOINT of the first two
                // nodes rather than at the endpoint, which for a non-constant costate is
                // simply the wrong place. It also feeds the worst value in the vector, the
                // endpoint itself, straight back into its neighbour.
                //
                // What replaces them keeps both intentions and drops both errors. Every
                // stencil below has weights that sum to one, have their centroid at the node
                // being computed (so a linear costate passes through untouched) and satisfy
                // w_left - w_centre + w_right = 0 (so a constant-amplitude alternating mode
                // is annihilated exactly, as the old interior stencil also did).
                //
                //   interior, 2 <= k <= M-2: the three-point filter on the node's own
                //     neighbours. Writing h- and h+ for the intervals either side,
                //         a = (1/2) h+/(h- + h+),  c = (1/2) h-/(h- + h+),  centre = 1/2.
                //     Only the ratio of the two intervals enters, so this is the same filter
                //     read in physical time or in tau, and on a uniform grid a = c = 1/4 and
                //     it is the old stencil exactly.
                //
                //   k = 1 and k = M-1: the same three conditions on a one-sided stencil
                //     drawn from interior nodes only, which for target t* on (t_a,t_b,t_c) is
                //         w_b = 1/2,  w_a = [t* - (t_b+t_c)/2]/(t_a - t_c),  w_c = 1/2 - w_a.
                //
                //   k = 0 and k = M: linear extrapolation from the two repaired neighbours.
                //
                // The endpoint value therefore never enters any stencil, its own included.
                // That is the point: it is the least trustworthy entry in the vector and the
                // old filter spread it inward instead of replacing it.
                //
                // Measured against the closed form, 40 nodes, old filter -> this one:
                //
                //   double integrator   lambda_2  1.4e-02 -> 1.4e-09    dH/du 1.4e-02 -> 1.8e-07
                //   linear tangent      lambda_3  4.09    -> 1.45       lambda_4 3.43 -> 1.23
                //                       lambda_2  0.62    -> 0.22       lambda_1 5.5e-5 -> 1.9e-5
                //
                // Better on every costate of both problems. It does not make the Legendre
                // costate competitive with Radau or Gauss, and nothing applied after the fact
                // could; see the discussion in the book's direct collocation chapter.
                // ------------------------------------------------------------------

                if (norder >= 4) {

                pint.resize(nstates,norder+1);

                // interior: the node's own two neighbours
                for (k=1;k<norder;k++) {
                    const double hm = (solution.nodes[i])(k)   - (solution.nodes[i])(k-1);
                    const double hp = (solution.nodes[i])(k+1) - (solution.nodes[i])(k);
                    if (hm + hp <= 0.0) {   // coincident nodes: pass the value through
                        pint.block(0,k,nstates,1) = solution.dual.costates[i].block(0,k,nstates,1);
                        continue;
                    }
                    const double a = 0.5*hp/(hm+hp);
                    const double c = 0.5*hm/(hm+hp);
                    pint.block(0,k,nstates,1) =
                          a*solution.dual.costates[i].block(0,k-1,nstates,1)
                        + (1.0-a-c)*solution.dual.costates[i].block(0,k,nstates,1)
                        + c*solution.dual.costates[i].block(0,k+1,nstates,1);
                }

                // k = 1 and k = M-1, recomputed one-sided so that the endpoint value is not
                // used. w_a = [t* - (t_b+t_c)/2]/(t_a - t_c) with t* the node itself.
                {
                    const double ta=(solution.nodes[i])(1), tb=(solution.nodes[i])(2), tc=(solution.nodes[i])(3);
                    if (ta != tc) {
                        const double wa = (ta - 0.5*(tb+tc))/(ta-tc), wc = 0.5-wa;
                        pint.block(0,1,nstates,1) =
                              wa*solution.dual.costates[i].block(0,1,nstates,1)
                            + 0.5*solution.dual.costates[i].block(0,2,nstates,1)
                            + wc*solution.dual.costates[i].block(0,3,nstates,1);
                    }
                }
                {
                    const double ta=(solution.nodes[i])(norder-1), tb=(solution.nodes[i])(norder-2), tc=(solution.nodes[i])(norder-3);
                    if (ta != tc) {
                        const double wa = (ta - 0.5*(tb+tc))/(ta-tc), wc = 0.5-wa;
                        pint.block(0,norder-1,nstates,1) =
                              wa*solution.dual.costates[i].block(0,norder-1,nstates,1)
                            + 0.5*solution.dual.costates[i].block(0,norder-2,nstates,1)
                            + wc*solution.dual.costates[i].block(0,norder-3,nstates,1);
                    }
                }

                // the two endpoints, from the repaired neighbours
                {
                    const double t0e=(solution.nodes[i])(0), t1=(solution.nodes[i])(1), t2=(solution.nodes[i])(2);
                    if (t2 != t1)
                        pint.block(0,0,nstates,1) = pint.block(0,1,nstates,1)
                            + ((t0e-t1)/(t2-t1))*(pint.block(0,2,nstates,1) - pint.block(0,1,nstates,1));
                    else
                        pint.block(0,0,nstates,1) = solution.dual.costates[i].block(0,0,nstates,1);
                }
                {
                    const double tNe=(solution.nodes[i])(norder), u1=(solution.nodes[i])(norder-1), u2=(solution.nodes[i])(norder-2);
                    if (u2 != u1)
                        pint.block(0,norder,nstates,1) = pint.block(0,norder-1,nstates,1)
                            + ((tNe-u1)/(u2-u1))*(pint.block(0,norder-2,nstates,1) - pint.block(0,norder-1,nstates,1));
                    else
                        pint.block(0,norder,nstates,1) = solution.dual.costates[i].block(0,norder,nstates,1);
                }

                solution.dual.costates[i] = pint;

                }   // norder >= 4; below that there is no oscillation to filter

         }


			if ( use_local_collocation(algorithm) ) {
        // use linear extrapolation to approximate the costate values at the collocation nodes...

         	pint = solution.dual.costates[i].block(0,0,nstates,norder);
				for(int k=0;k<norder;k++) {  // EIGEN_UPDATE: k index shifted by -1.
		   		tint(k) = ( (workspace->snodes[i])(k)+(workspace->snodes[i])(k+1) )/2.0;
				}
        		for (int l=0;l<nstates;l++) { // EIGEN_UPDATE: l index shifted by -1.

                   ts = workspace->snodes[i].transpose();

                   tint = tint.transpose().eval();

                   pl = pint.row(l);
                   linear_interpolation(pextra, ts, tint, pl,norder); // EIGEN_UPDATE - Check function

                   solution.dual.costates[i].row(l) = pextra;
       		}		

        // use linear extrapolation to approximate costate values at end point (not perfect but better than nothing)

                pint = solution.dual.costates[i].block(0  ,norder-2  ,nstates  , 2 );

                tint = workspace->snodes[i].block(0  ,norder-2  ,1  , 2 );
        		for (int l=0;l<nstates;l++) {  // EIGEN_UPDATE: l index shifted by -1.

                   double tss = (workspace->snodes[i])(norder);

                   tint = tint.transpose().eval();

                   long ncols = pint.cols();
                   pl = pint.block(l, 0, 1, ncols);
                   linear_interpolation(pextra, tss, tint, pl,2);

                   solution.dual.costates[i](l,norder) = pextra(0);
        		}	

         }


		// Event (and, below, path) multipliers are gg_ad-aligned: the events block
		// begins immediately after the nstates*(norder+1) differential-defect multipliers,
		// at lam_phase_offset + nstates*(norder+1). (A previous +1 here skipped event 0 and
		// over-read into the constraint following the events block -- the t0-tf constraint
		// for Legendre/Chebyshev, the terminal-control pin for Radau, or the terminal-state
		// quadrature constraint for Gauss.)
		offset = lam_phase_offset+nstates*(norder+1);
	
		if (   algorithm.nlp_method == "IPOPT" || algorithm.nlp_method == "SQP"   ) {
	                solution.dual.costates[i] = -solution.dual.costates[i];
	    }
	

		// Gauss (Legendre-Gauss / Benson Option 1) boundary costate recovery, generalised
		// to the multi-interval mesh via the Darby-Garg-Rao covector mapping. In each
		// interval the raw interior mapping lambda_k = lambda_tilde_k/w_k recovers the
		// costate RELATIVE to that interval's right-endpoint costate, which is exactly the
		// multiplier of that interval's Gauss-quadrature defining constraint, lamf_j. So per
		// interval j we add lamf_j back to its interior costates; the interface breakpoint
		// tau_j then carries lamf_j (the costate at the right end of interval j = left end of
		// interval j+1, continuous for a smooth costate); the terminal costate is lamf_K; and
		// the non-collocated initial costate is a LOCAL Lagrange extrapolation of interval 1's
		// corrected interior costates to tau=-1 (a global extrapolation over the clustered
		// multi-interval nodes would be ill-conditioned). Single interval (K=1) reduces to the
		// original Benson recovery bit-identically.
		if ( algorithm.collocation_method=="Gauss" ) {
		    int quad = lam_phase_offset + nstates*(norder+1) + nevents + npath*(norder+1);
		    MatrixXd& sn = workspace->snodes[i];

		    int Kg = hp_mesh_active(problem.phase[i]) ? hp_num_intervals(problem.phase[i]) : 1;
		    std::vector<int> bp(Kg), gn(Kg);
		    int gidx = 0;
		    for (int jj=0; jj<Kg; jj++) {
		        gn[jj] = hp_mesh_active(problem.phase[i]) ? hp_interval_order(problem.phase[i], jj) : norder;
		        bp[jj] = gidx;
		        gidx  += 1 + gn[jj];
		    }

		    MatrixXd lamfK(nstates,1); lamfK.setZero();
		    for (int jint=0; jint<Kg; jint++) {
		        // interval jint's defining-constraint multiplier, pushed through the same
		        // scaling/sign pipeline as the defect costates (so it lands in physical units).
		        MatrixXd lamfj(nstates,1);
		        for (int j=0;j<nstates;j++) {
		            double v = lambda(quad + jint*nstates + j);
		            if ( algorithm.scaling=="user" )      v *= deriv_scaling(j);
		            v /= problem.scale.objective;
		            if ( algorithm.scaling=="automatic" ) v *= (*workspace->constraint_scaling)(quad + jint*nstates + j);
		            if ( algorithm.nlp_method=="IPOPT" || algorithm.nlp_method=="SQP" )  v = -v;
		            lamfj(j) = v;
		        }
		        // (1) undo the per-interval constant shift on this interval's interior costates
		        int kf = bp[jint]+1, kl = bp[jint]+gn[jint];
		        for (k=kf; k<=kl; k++)
		            (solution.dual.costates[i]).block(0,k,nstates,1) += lamfj;
		        // (2) interface breakpoint tau_{jint+1} carries this interval's right-end costate
		        if (jint < Kg-1) {
		            int bcol = bp[jint+1];
		            for (int r=0;r<nstates;r++) (solution.dual.costates[i])(r,bcol) = lamfj(r);
		        }
		        if (jint == Kg-1) lamfK = lamfj;     // terminal costate lambda(+1)
		    }

		    // (3) initial costate (column 0, tau=-1) by LOCAL extrapolation over interval 1's
		    //     corrected interior Gauss costates.
		    {
		        double xe = sn(0);
		        int kf = bp[0]+1, kl = bp[0]+gn[0];
		        for (int r=0;r<nstates;r++) (solution.dual.costates[i])(r,0) = 0.0;
		        for (int m=kf; m<=kl; m++) {
		            double Lwm = 1.0;
		            for (int jj=kf; jj<=kl; jj++) if (jj!=m) Lwm *= (xe - sn(jj))/(sn(m)-sn(jj));
		            for (int r=0;r<nstates;r++)
		                (solution.dual.costates[i])(r,0) += Lwm * (solution.dual.costates[i])(r,m);
		        }
		    }
		    workspace->prev_costates[i] = solution.dual.costates[i]; // keep hot-start consistent
		    solution.dual.terminal_costates[i] = lamfK;              // lambda(+1) for reporting
		}

		solution.dual.events[i]  = -lambda.block(offset,0,nevents,1);
		workspace->dual_events[i] = -(solution.dual.events[i]);
	
		if (algorithm.scaling=="user") {

			   solution.dual.events[i] = solution.dual.events[i].cwiseProduct( problem.phase[i].scale.events );
		}
	
	    if (algorithm.scaling=="automatic" && nevents>0) {

		   solution.dual.events[i] = solution.dual.events[i].cwiseProduct(  (*workspace->constraint_scaling).block(offset,0,nevents,1) );
	    }
		solution.dual.events[i] /= problem.scale.objective;
	
		if (   algorithm.nlp_method == "IPOPT" || algorithm.nlp_method == "SQP"   ) {
	                solution.dual.events[i] = -solution.dual.events[i];
	    }

		offset = offset+nevents;
		if (npath>0) {

	             solution.dual.path[i]      = -lambda.block(offset,0, npath*(norder+1), 1);
		     solution.dual.path[i]      = reshape(solution.dual.path[i], npath, norder+1);
		     workspace->prev_path[i] = -(solution.dual.path[i]);
		     if (use_local_collocation(algorithm)) {
	        	for (k=0;k<norder;k++) {  // EIGEN_UPDATE: k index shifted by -1
	        	    // Below are the path constraint adjoint estimates for trapezoidal discretization
	        	    // See Betts (2010), p. 176.
	        	    if ( algorithm.collocation_method == "trapezoidal") {
	                    double hk = (solution.nodes[i])(k+1)-(solution.nodes[i])(k);
	                    if (k==0) {

	                          (solution.dual.path[i]).block(0,k,npath,1) = -(solution.dual.path[i]).block(0,k,npath,1)*(2.0/hk);
	                    }
	                    else if (k==norder-1) {
	                          double hk_1 = (solution.nodes[i])(k)-(solution.nodes[i])(k-1);

	                          (solution.dual.path[i]).block(0,k,npath,1) = -(solution.dual.path[i]).block(0,k,npath,1)*(2.0/hk_1);
	                    }
	                    else {
	                        double hk_1 = (solution.nodes[i])(k)-(solution.nodes[i])(k-1);

	                          (solution.dual.path[i]).block(0,k,npath,1) = -(solution.dual.path[i]).block(0,k,npath,1)*(2.0/(hk+hk_1));
	        	    }
	                    }
	           	    if ( algorithm.collocation_method == "Hermite-Simpson") {
	           	        // These are the path constraint adjoint estimates for Hermite-Simpson discretization
	           	        // See Betts (2010), p. 177.
	           	        double hk = (solution.nodes[i])(k+1)-(solution.nodes[i])(k);
	                    if (k==0) {

	                          (solution.dual.path[i]).block(0,k,npath,1) = -(solution.dual.path[i]).block(0,k,npath,1)*(6.0/hk);
	                    }
	                    else if (k==norder-1) {
	                        double hk_1 = (solution.nodes[i])(k)-(solution.nodes[i])(k-1);

	                        (solution.dual.path[i]).block(0,k,npath,1) = -(solution.dual.path[i]).block(0,k,npath,1)*(6.0/hk_1);
	                    }
	                    else {
	                        double hk_1 = (solution.nodes[i])(k)-(solution.nodes[i])(k-1);

	                        (solution.dual.path[i]).block(0,k,npath,1) = -(solution.dual.path[i]).block(0,k,npath,1)*(6.0/(hk+hk_1));
	                    }
	           	    }
	
	
	        	}
             // use linear extrapolation to approximate multiplier values at end point (not perfect but better than nothing)

          	pint = solution.dual.path[i].block(0,norder-2,npath,2);

          	tint = workspace->snodes[i].block(0,norder-2,1,2);
             for (int l=0;l<npath;l++) {  //EIGEN_UPDATE: l index shifted by -1.

                double tss = (workspace->snodes[i])(norder);

                tint = tint.transpose().eval();

                pl   = pint.block(l,0,1,2);
                linear_interpolation(pextra, tss, tint, pl,2);

                  solution.dual.path[i](l,norder) = pextra(0);

             }

	     }

	     if ( algorithm.collocation_method == "Legendre") {
             for (k=0;k<norder+1;k++) {  // EIGEN_UPDATE: k index shifted by -1

		     		(solution.dual.path[i]).block(0,k,npath,1) = -(solution.dual.path[i]).block(0,k,npath,1)/(workspace->w[i])(k);  // See PhD thesis by Huntington (2006).

		     		(solution.dual.path[i]).block(0,k,npath,1) = (solution.dual.path[i]).block(0,k,npath,1)*(2.0/(tf-t0));
             }
	     }

	     if ( algorithm.collocation_method == "Chebyshev" ) {
             for (k=1; k< norder;k++) { // EIGEN_UPDATE: k index shifted by -1
                  double tk = (workspace->snodes[i])(k);

                    (solution.dual.path[i]).block(0,k,npath,1) = (solution.dual.path[i]).block(0,k,npath,1)/(workspace->w[i])(k);

                    (solution.dual.path[i]).block(0,k,npath,1) = -(1.0/sqrt(1.0 - tk*tk))*(solution.dual.path[i]).block(0,k,npath,1); // See PhD thesis by Pietz (2003)

                    (solution.dual.path[i]).block(0,k,npath,1) = (solution.dual.path[i]).block(0,k,npath,1)*(2.0/(tf-t0));

             }
            // use linear extrapolation to approximate multiplier values at both ends (not perfect but better than nothing)

                pint = solution.dual.path[i].block(0,2,npath,norder-3);


             tint.resize(1,norder-4);
             for (int ii=0; ii<(norder-4); ii++) {
                 tint(ii) = workspace->snodes[i](2+ii);
             }
            for (int l=0;l<npath;l++) {  //EIGEN_UPDATE: l index shifted by -1

                   ts = workspace->snodes[i].transpose();

                   tint = tint.transpose().eval();

                     long ncols = pint.cols();
                     pl = pint.block(l,0,1,ncols);
                   linear_interpolation(pextra, ts, tint, pl,norder-4);

                   for(k=0; k<length(pextra);k++) {
                	   solution.dual.path[i](l,k) = pextra(k);
                   }
            }
	     }

	     if ( algorithm.scaling == "user") {
	         for(k=0;k<norder+1;k++) { //EIGEN_UPDATE: k index shifted by -1

                      (solution.dual.path[i]).block(0,k,npath,1) = (solution.dual.path[i]).block(0,k,npath,1).cwiseProduct(problem.phase[i].scale.path);
	         }
         }

         if (algorithm.scaling=="automatic") {
             for (k=0;k<norder+1;k++) { //EIGEN_UPDATE: k index shifted by -1

                    solution.dual.path[i].block(0,k,npath,1)  =solution.dual.path[i].block(0,k,npath,1).cwiseProduct( (*workspace->constraint_scaling).block(offset+(k)*npath,0,npath,1) );
  	         }
         }
         solution.dual.path[i] /= problem.scale.objective;

	}

    offset = offset + npath*(norder+1);


    compute_derivatives_trajectory( workspace->Xdot[i], problem, solution, iphase, workspace );


       solution.dual.Hamiltonian[i] = (solution.integrand_cost[i]);
       MatrixXd Temp1 =  solution.dual.costates[i].cwiseProduct(workspace->Xdot[i]);
       MatrixXd Temp2 = sum_columns(Temp1);
       solution.dual.Hamiltonian[i]+= Temp2;



		workspace->prev_states[i]   =   solution.states[i];
		workspace->prev_controls[i] =   solution.controls[i];
		workspace->prev_nodes[i]    =   solution.nodes[i];
        if (problem.phase[i].nparameters) workspace->prev_param[i]    =   solution.parameters[i];

        (*workspace->prev_t0)(i) = x0(x_phase_offset+nvars_phase_i-2)/problem.phase[i].scale.time;

        (*workspace->prev_tf)(i) = x0(x_phase_offset+nvars_phase_i-1)/problem.phase[i].scale.time;

        x_phase_offset   += nvars_phase_i;
        lam_phase_offset += ncons_phase_i;

	// store scaled nodes
	workspace->old_snodes[i] = workspace->snodes[i];


    }

    if (problem.nlinkages) {

         *solution.dual.linkages = -lambda.block(lam_phase_offset, 0, problem.nlinkages, 1);
         if (algorithm.scaling=="user") {

             *solution.dual.linkages = (*solution.dual.linkages).cwiseProduct(problem.scale.linkages);

         }
         if (algorithm.scaling=="automatic") {

           *solution.dual.linkages =  (*solution.dual.linkages).cwiseProduct(   (*workspace->constraint_scaling).block(offset,0,problem.nlinkages,1)   );
         }
         *solution.dual.linkages /= problem.scale.objective;

         if (algorithm.nlp_method == "IPOPT" || algorithm.nlp_method == "SQP") {
             *solution.dual.linkages = -(*solution.dual.linkages);
         }
    }

    if (!useAutomaticDifferentiation(algorithm) && algorithm.nlp_method=="IPOPT")  {
//          deleteIndexGroups( workspace->igroup, workspace->nvars );
    }

    // ---- Robust-DAIR (Option B): adjoint costate recovery ----
    // For the residual-box solve the defect multipliers carry no costate information; recover
    // the costates by backward adjoint integration along the primal (see recover_costates_adjoint).
    if ( algorithm.transcription_method == "integrated-residual"
         && ( algorithm.ir_dair
              || ( algorithm.ir_objective == "cost" && algorithm.ir_residual_bound >= 0.0 ) ) ) {
        recover_costates_adjoint(problem, algorithm, solution, workspace);
    }

    evaluate_solution(problem, algorithm, solution, workspace);

    if ( algorithm.mesh_refinement == "automatic" ) {
       // Check satisfaction of mesh refinement tolerance
       int mr_phase_convergence_count = 0;
       for ( i=0; i< problem.nphases; i++ ) {
	    MatrixXd& emax_history = workspace->emax_history[i];

            if ( emax_history( iter_nodes-1, 1 ) <= algorithm.ode_tolerance )
	        mr_phase_convergence_count++;
       }



       if (mr_phase_convergence_count == problem.nphases ) {
	    psopt_print(workspace,"\n>>> PSOPT: automatic mesh refinement iterations converged as the maximum");
	    psopt_print(workspace,"\n>>> relative error in all phases is lower than algorithm.ode_tolerance\n");
	    break; // break the iterations.
       }

    }

    // hp-adaptive automatic (Route B, Liu-Hager-Rao ph): rewrite each phase's hp mesh from
    // the per-interval error, ready for the next iteration. This is the only automatic
    // refinement for the global pseudospectral methods; the local methods (Betts) build their
    // next mesh in the use_local_collocation pre-solve branch (construct_new_mesh).
    if ( hp_auto_active(algorithm) && iter_nodes < number_of_mesh_refinement_iterations )
    {
        hp_refine_driver( problem, algorithm, solution, workspace );
    }


  } // End of mesh refinement iterations loop

  // Before anything reports the solution, and after everything that drives the solve has
  // finished reading the unaugmented arrays.
  append_gauss_terminal_point(problem, algorithm, solution, workspace);

  solution.cpu_time = PSOPT_extras::toc();

  get_local_time( solution.end_date_and_time );

  if (algorithm.diagnostic_level > 0) {
     solution_diagnostics(problem, algorithm, solution, workspace);
  }



  // The parameter statistics are a result of the solve, not a side effect of printing:
  // they are computed here so that a caller reaches them through the solution object
  // whatever print_level is set to. The report below then formats what is already there.

  store_parameter_statistics(problem, algorithm, solution, workspace);

  if (algorithm.print_level>0) {

    print_algorithm_summary(problem, algorithm, solution, workspace);

    print_solution_summary(problem, algorithm, solution, workspace);

    print_constraint_summary(problem, solution, workspace);

    print_iterations_summary(problem,algorithm,solution, workspace);

    print_iterations_summary_tex(problem,algorithm,solution, workspace);

    print_psopt_summary(problem, algorithm, solution, workspace);

  }


//  if (workspace) delete  workspace;
  
  return;

}


