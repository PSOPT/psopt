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



adouble integrate( adouble (*integrand)(adouble*,adouble*,adouble*,adouble&,adouble*,int, Workspace* workspace), adouble* xad, int iphase, Workspace* workspace )
{
// Evaluates the integral of a user supplied function over the time
// span of phase iphase
    	   Prob& problem = *workspace->problem;
	     Alg&  algorithm = *workspace->algorithm;
        int i = iphase-1;

        int k;

        adouble* controls;
        adouble* states;
	     adouble* states_next;
        adouble* parameters;
        adouble time;
        adouble t0;
        adouble tf;
        adouble ieval;
        adouble retval = 0.0;

	     MatrixXd& w = workspace->w[i];

        int norder    = problem.phase[i].current_number_of_intervals;

	     states        = workspace->states[i].get();
	     states_next   = workspace->states_next[i].get();
        controls      = workspace->controls[i].get();
        parameters    = workspace->parameters[i].get();

        get_parameters(parameters, xad, iphase, workspace);

        get_times(&t0, &tf, xad, iphase, workspace);

	     if ( !use_local_collocation(algorithm) ) {

	        for(k=0; k<norder+1; k++)  // EIGEN_UPDATE: k index shifted by -1
	        {

		         get_controls(controls, xad, iphase, k, workspace);

		         get_states(states, xad, iphase, k, workspace);

		         time = convert_to_original_time_ad( (workspace->snodes[i])(k), t0, tf );

		         ieval = (*integrand)(states,controls,parameters,time,xad,iphase, workspace);

		         retval += ((tf-t0)/2.0)*ieval*w(k);
	       }

	}

	else if ( workspace->algorithm->ir_local_order >= 2
	          && norder >= workspace->algorithm->ir_local_order
	          && (norder % workspace->algorithm->ir_local_order) == 0 ) {

		  // Nie-Kerrigan local representation: the same per-element Lobatto quadrature
		  // that phase_running_cost uses for the objective, and for the same reason. The
		  // state and control on an element are degree-d polynomials through that
		  // element's own d+1 nodes, and the Hermite-Simpson midpoint of the branch below
		  // -- the cubic Hermite state and the midpoint control VARIABLE -- belongs to a
		  // different representation. The midpoint control is not read by anything in this
		  // transcription, so integrating against it returned the integral of a trajectory
		  // the solver never chose: on the test problem min int u^2 with xdot = u over
		  // [0,1] from 0 to 1, whose answer is 1, this routine returned about a half,
		  // because the midpoint control variables sit wherever the initial guess left
		  // them. The snodes within element e are that element's d+1 LGL nodes; the rule
		  // is exact to degree 2d-1 and so matched to the representation.
		  int d = workspace->algorithm->ir_local_order;
		  int M = norder / d;
		  MatrixXd& wl = workspace->ir_lgl_w;         // d+1 LGL weights on [-1,1], sum 2
		  // A flexible mesh changes both the element width the rule is scaled by and the
		  // times the integrand is sampled at; empty, and every time below is the stored
		  // node position it always was.
		  std::vector<adouble> ir_tau;
		  ir_node_taus(ir_tau, xad, iphase, workspace);
		  for (int e=0; e<M; e++) {
		      int base = e*d;
		      adouble te0 = ir_node_time( ir_tau, base,   t0, tf, workspace->snodes[i] );
		      adouble te1 = ir_node_time( ir_tau, base+d, t0, tf, workspace->snodes[i] );
		      adouble he  = te1 - te0;
		      for (int r=0; r<=d; r++) {
		          int gk = base + r;
		          get_element_controls(controls, xad, iphase, e, r, workspace);
		          get_states(states,     xad, iphase, gk, workspace);
		          adouble tnode = ir_node_time( ir_tau, gk, t0, tf, workspace->snodes[i] );
		          retval += (he/2.0) * wl(r)
		                    * (*integrand)(states,controls,parameters,tnode,xad,iphase, workspace);
		      }
		  }

	}

	else {


		  // Under Hermite-Simpson the midpoint state of Simpson's rule is the cubic
		  // Hermite value the transcription defines, not the arithmetic mean; see the
		  // note in phase_running_cost (NLP_objective.cxx), which had the same defect.
		  const bool     hs_int    = need_midpoint_controls(algorithm, workspace);
		  const int      ns_int    = problem.phase[i].nstates;
		  adouble* const derivs    = workspace->derivatives[i].get();
		  adouble* const derivs_nx = workspace->derivatives_next[i].get();
		  adouble* const path_scr  = workspace->path[i].get();
		  adouble* const path_scr2 = workspace->path_next[i].get();
		  adouble* const states_bar= workspace->states_bar[i].get();

		  // Non-empty only when the legacy integrated-residual form runs on a flexible
		  // mesh, where the interval ends are variables; see the same pairing in
		  // phase_running_cost.
		  std::vector<adouble> ir_tau_hs;
		  ir_node_taus(ir_tau_hs, xad, iphase, workspace);

		  for (k=0; k<norder;k++) {  // EIGEN_UPDATE: k index shifted by -1
		      int l;

		      adouble interval_value = 0.0;

		      get_controls(controls, xad, iphase, k, workspace);
		      get_states(states, xad, iphase, k, workspace);

		      adouble tk = ir_node_time( ir_tau_hs, k,   t0, tf, workspace->snodes[i] );
		      adouble tk1= ir_node_time( ir_tau_hs, k+1, t0, tf, workspace->snodes[i] );

		      adouble h = tk1-tk;

		      interval_value = (*integrand)(states,controls,parameters,tk,xad,iphase, workspace);

		      if (hs_int) {
		          if (k == 0)
		              problem.dae(derivs, path_scr, states, controls, parameters, tk, xad, iphase, workspace);
		          else
		              for (l=0; l<ns_int; l++) derivs[l] = derivs_nx[l];
		      }

		      get_controls(controls, xad, iphase,k+1 , workspace);
		      get_states(states_next, xad, iphase, k+1, workspace);


		      interval_value += (*integrand)(states_next,controls,parameters,tk1,xad,iphase, workspace);

		      if (hs_int) {

			         problem.dae(derivs_nx, path_scr2, states_next, controls, parameters, tk1, xad, iphase, workspace);

			         adouble tmiddle = (tk+tk1)/2.0;

			         for( l =0; l< ns_int; l++ ) {

			             states_bar[l] = 0.5*(states[l]+states_next[l])
			                             + h*(derivs[l]-derivs_nx[l])/8.0;

			        }

			         get_controls_bar(controls,xad,iphase,k, workspace);

			        interval_value += 4.0*(*integrand)(states_bar,controls,parameters,tmiddle,xad,iphase, workspace);


			        interval_value *= h/6.0;

		     }

		     else {
		           interval_value *= h/2.0;
		      }

		      retval+= interval_value;

		  }

       }

       return retval;

}



// ===========================================================================================
// One multiple-shooting segment, propagated.
//
// The scheme is classical RK4 at a fixed step, unrolled onto the same tape as the rest of the
// problem. That choice is not a convenience. A shooting method needs the derivative of a
// segment's end state with respect to its start state, its control, the parameters and its
// duration, and there are only two honest ways to get one: tape a fixed-step scheme, as here,
// or propagate sensitivities with the same steps and factorisations as the nominal trajectory,
// which is Bock's internal numerical differentiation. What must NOT be done is to difference
// an adaptive integrator, because that differentiates its step controller along with its
// dynamics and returns noise. A fixed-step scheme has no controller to differentiate, so
// taping it gives the exact derivative of the discrete map the transcription actually uses --
// which is the derivative the NLP wants, not the derivative of the ODE.
//
// The running cost is accumulated by the same RK4 stages rather than by a separate quadrature
// on the segment ends. That is state augmentation without the extra variable: the integral of
// L along the segment obeys dJ/dt = L, so the same scheme applied to the augmented system
// integrates it to the same order. A trapezoidal rule on the two segment ends would be second
// order against the trajectory's fourth, and the objective the NLP minimised would not be the
// objective the user wrote -- which is the defect that had to be fixed in the Hermite-Simpson
// running cost, and there is no reason to reintroduce it here.
//
// The control is piecewise constant over the segment, so it is read once and held.
// ===========================================================================================
// The two schemes, as Butcher tableaux. Writing the propagation against a tableau rather than
// unrolling one scheme is what makes a second scheme cost a table instead of a second loop.
//
// Classical fourth-order Runge-Kutta.
static const double ms_rk4_c[4] = { 0.0, 0.5, 0.5, 1.0 };
static const double ms_rk4_b[4] = { 1.0/6.0, 1.0/3.0, 1.0/3.0, 1.0/6.0 };
static const double ms_rk4_A[16] = {
    0.0, 0.0, 0.0, 0.0,
    0.5, 0.0, 0.0, 0.0,
    0.0, 0.5, 0.0, 0.0,
    0.0, 0.0, 1.0, 0.0 };

// Cooper and Verner's eleven-stage eighth-order formula (SIAM J. Numer. Anal. 9 (3), 1972).
// Eleven stages is the minimum for order eight, and every coefficient lies in Q(sqrt(21)),
// so the table is exact rather than a decimal truncation of one.
//
// The table was checked before it was used, which for a Runge-Kutta tableau is the whole of
// the verification problem: a mistyped coefficient gives a scheme of lower order that still
// converges to the right answer, so it does not fail, it just stops being what it claims to
// be. Four checks, none of which depends on PSOPT: the row sums equal c to 1.3e-15; sum b is
// 1 exactly; the linear conditions sum b_i c_i^(k-1) = 1/k hold to round-off for k = 1..8 and
// fail at k = 9; the stability polynomial matches the exponential series through z^8 and
// departs at z^9; and the measured order on a nonlinear non-autonomous system, in fifty-digit
// arithmetic, is 7.914, 7.975, 7.991 for successive halvings.
static const double MS_S21 = 4.58257569495584000658804719372800848898445657676797190260724212;

static double ms_rk8_c[11];
static double ms_rk8_b[11];
static double ms_rk8_A[121];
static bool   ms_rk8_built = false;

static void ms_build_rk8(void)
{
    if (ms_rk8_built) return;
    const double s = MS_S21;
    for (int q = 0; q < 121; q++) ms_rk8_A[q] = 0.0;
    for (int q = 0; q < 11;  q++) { ms_rk8_b[q] = 0.0; ms_rk8_c[q] = 0.0; }

#define AA(i,j) ms_rk8_A[(i)*11 + (j)]
    ms_rk8_c[1] = 0.5;                 AA(1,0) = 0.5;
    ms_rk8_c[2] = 0.5;                 AA(2,0) = 0.25;            AA(2,1) = 0.25;
    ms_rk8_c[3] = (7.0+s)/14.0;        AA(3,0) = 1.0/7.0;         AA(3,1) = (-7.0-3.0*s)/98.0;
                                       AA(3,2) = (21.0+5.0*s)/49.0;
    ms_rk8_c[4] = (7.0+s)/14.0;        AA(4,0) = (11.0+s)/84.0;   AA(4,2) = (18.0+4.0*s)/63.0;
                                       AA(4,3) = (21.0-s)/252.0;
    ms_rk8_c[5] = 0.5;                 AA(5,0) = (5.0+s)/48.0;    AA(5,2) = (9.0+s)/36.0;
                                       AA(5,3) = (-231.0+14.0*s)/360.0;
                                       AA(5,4) = (63.0-7.0*s)/80.0;
    ms_rk8_c[6] = (7.0-s)/14.0;        AA(6,0) = (10.0-s)/42.0;   AA(6,2) = (-432.0+92.0*s)/315.0;
                                       AA(6,3) = (633.0-145.0*s)/90.0;
                                       AA(6,4) = (-504.0+115.0*s)/70.0;
                                       AA(6,5) = (63.0-13.0*s)/35.0;
    ms_rk8_c[7] = (7.0-s)/14.0;        AA(7,0) = 1.0/14.0;        AA(7,4) = (14.0-3.0*s)/126.0;
                                       AA(7,5) = (13.0-3.0*s)/63.0;
                                       AA(7,6) = 1.0/9.0;
    ms_rk8_c[8] = 0.5;                 AA(8,0) = 1.0/32.0;        AA(8,4) = (91.0-21.0*s)/576.0;
                                       AA(8,5) = 11.0/72.0;
                                       AA(8,6) = (-385.0-75.0*s)/1152.0;
                                       AA(8,7) = (63.0+13.0*s)/128.0;
    ms_rk8_c[9] = (7.0+s)/14.0;        AA(9,0) = 1.0/14.0;        AA(9,4) = 1.0/9.0;
                                       AA(9,5) = (-733.0-147.0*s)/2205.0;
                                       AA(9,6) = (515.0+111.0*s)/504.0;
                                       AA(9,7) = (-51.0-11.0*s)/56.0;
                                       AA(9,8) = (132.0+28.0*s)/245.0;
    ms_rk8_c[10] = 1.0;                AA(10,4) = (-42.0+7.0*s)/18.0;
                                       AA(10,5) = (-18.0+28.0*s)/45.0;
                                       AA(10,6) = (-273.0-53.0*s)/72.0;
                                       AA(10,7) = (301.0+53.0*s)/72.0;
                                       AA(10,8) = (28.0-28.0*s)/45.0;
                                       AA(10,9) = (49.0-7.0*s)/18.0;
#undef AA

    ms_rk8_b[0]  = 1.0/20.0;
    ms_rk8_b[7]  = 49.0/180.0;
    ms_rk8_b[8]  = 16.0/45.0;
    ms_rk8_b[9]  = 49.0/180.0;
    ms_rk8_b[10] = 1.0/20.0;

    ms_rk8_built = true;
}

void ms_propagate_segment(adouble* xend, adouble* Lint, int k, adouble* xad, int iphase,
                          adouble& t0, adouble& tf, adouble* parameters, Workspace* workspace,
                          int nsteps_override,
                          adouble* xsamp, adouble* usamp, adouble* tsamp,
                          std::vector<adouble>* tau)
{
    Prob& problem   = *workspace->problem;
    Alg&  algorithm = *workspace->algorithm;
    const int i         = iphase-1;
    const int nstates   = problem.phase[i].nstates;
    const int ncontrols = problem.phase[i].ncontrols;

    int nsteps = ( nsteps_override > 0 ) ? nsteps_override : algorithm.ms_steps_per_segment;
    if ( nsteps < 1 ) nsteps = 1;

    // Local scratch, deliberately. Every one of the workspace's per-phase buffers is live in
    // the caller -- this is called from inside gg_ad's node loop, which is holding states,
    // controls, derivatives and path across the call -- and borrowing one of them here would
    // corrupt the row being written rather than fail.
    // The scheme. Everything below is written against the tableau, so RK4 and RK8 differ in
    // the table and in nothing else.
    const int     nstg = ms_integrator_stages(algorithm);
    const double* Atab;
    const double* btab;
    const double* ctab;
    if ( ms_rk8(algorithm) ) { ms_build_rk8(); Atab = ms_rk8_A; btab = ms_rk8_b; ctab = ms_rk8_c; }
    else                     { Atab = ms_rk4_A; btab = ms_rk4_b; ctab = ms_rk4_c; }

    const int npath = problem.phase[i].npath;
    std::vector<adouble> u_( (ncontrols>0) ? ncontrols : 1 );
    std::vector<adouble> xw_(nstates), xstg_(nstates);
    std::vector<adouble> K_(nstg*nstates);
    std::vector<adouble> Lstg_(nstg);
    std::vector<adouble> pscr_( (npath>0) ? npath : 1 );
    adouble* const u     = u_.data();
    adouble* const xw    = xw_.data();
    adouble* const xstg  = xstg_.data();
    adouble* const K     = K_.data();
    adouble* const Lstg  = Lstg_.data();
    adouble* const pscr  = pscr_.data();

    // The segment ends. Under a fixed partition these are stored node positions and constants
    // to the tape; under a flexible one they are expressions in the segment widths, and the
    // whole of the segment -- its duration, its step length, every stage time inside it, and
    // therefore the end state and the running cost -- becomes differentiable in them. That is
    // the entire cost of a moving partition here, and it is why the accessor is shared with
    // the integrated residual rather than reimplemented: it is the same parameterisation.
    std::vector<adouble> ms_tau_local;
    if ( tau == NULL ) { ir_node_taus(ms_tau_local, xad, iphase, workspace); tau = &ms_tau_local; }
    adouble tk  = ir_node_time( *tau, k,   t0, tf, workspace->snodes[i] );
    adouble tk1 = ir_node_time( *tau, k+1, t0, tf, workspace->snodes[i] );
    adouble dt  = (tk1 - tk)/((double) nsteps);

    get_states(xw, xad, iphase, k, workspace);

    // The control across the segment. Piecewise constant reads one value and holds it;
    // piecewise linear reads both ends and ramps between them; piecewise quadratic reads both
    // ends and the segment's own midpoint variable and carries the parabola through the three.
    // Every coefficient is an ordinary double -- the local coordinate of a stage is a fixed
    // fraction of the segment, whatever the segment's physical length turns out to be, and
    // that stays true when the partition itself is a decision variable -- so the higher forms
    // cost two or three multiplications per stage and nothing on the tape's structure.
    const bool linear_u = ms_linear_controls(algorithm);
    const bool quad_u   = ms_quadratic_controls(algorithm);
    const bool varying_u = ( linear_u || quad_u );
    std::vector<adouble> u0_( (ncontrols>0) ? ncontrols : 1 );
    std::vector<adouble> u1_( (ncontrols>0) ? ncontrols : 1 );
    std::vector<adouble> um_( (ncontrols>0) ? ncontrols : 1 );
    if (ncontrols > 0) {
        get_controls(u0_.data(), xad, iphase, k, workspace);
        if (varying_u) get_controls(u1_.data(), xad, iphase, k+1, workspace);
        if (quad_u)    get_controls_bar(um_.data(), xad, iphase, k, workspace);
    }

    // The control at local coordinate s in [0,1] across this segment, written once so that
    // the four RK4 stages, the interior path samples and the reader of the reported control
    // cannot drift apart. The quadratic weights are the same three Lagrange factors
    // get_interpolated_control uses, which is what makes the reported control the control the
    // integrator actually saw.
    auto eval_u = [&](double s) {
        if (ncontrols <= 0) return;
        if (quad_u) {
            const double w0 = (2.0*s-1.0)*(s-1.0);
            const double wm = 4.0*s*(1.0-s);
            const double w1 = s*(2.0*s-1.0);
            for (int c=0;c<ncontrols;c++) u[c] = w0*u0_[c] + wm*um_[c] + w1*u1_[c];
        }
        else if (linear_u) {
            for (int c=0;c<ncontrols;c++) u[c] = (1.0-s)*u0_[c] + s*u1_[c];
        }
        else {
            for (int c=0;c<ncontrols;c++) u[c] = u0_[c];
        }
    };

    const bool want_cost = ( Lint != NULL );
    if (want_cost) *Lint = 0.0;

    // Where the interior samples fall, as step indices. They are placed at step boundaries so
    // that the state at a sample is one the integrator actually produced, rather than an
    // interpolation of states either side of it.
    const int nsamp = ( xsamp != NULL || usamp != NULL || tsamp != NULL )
                      ? algorithm.ms_path_samples : 0;
    std::vector<int> sample_step(nsamp>0 ? nsamp : 1, 0);
    for (int q = 0; q < nsamp; q++) {
        int st = (int) ( ( (double)(q+1) * (double) nsteps )/((double)(nsamp+1)) + 0.5 );
        if ( st < 1 )        st = 1;
        if ( st > nsteps-1 ) st = nsteps-1;
        sample_step[q] = st;
    }

    // A held control is read once for the whole segment. The varying forms are read at every
    // stage, which is the only place the stage count reaches the control at all.
    if ( !varying_u ) eval_u(0.0);

    adouble t = tk;
    for (int s = 0; s < nsteps; s++) {

        for (int p = 0; p < nstg; p++) {

            // Where this stage sits: as a fraction of the SEGMENT, for the control, and as a
            // time, for the dynamics. The control's coordinate is a fixed double whatever the
            // segment's physical length turns out to be, which is what keeps a flexible
            // partition from changing the control's shape as well as its span.
            const double sp = ( (double) s + ctab[p] )/((double) nsteps);
            if ( varying_u ) eval_u(sp);

            adouble tp = t + ctab[p]*dt;

            // The stage state. Stage zero has no coefficients at all, so it reads the step's
            // own state rather than copying it; zero coefficients elsewhere are skipped, which
            // matters because the eighth-order table is mostly zeros and every skipped term is
            // a multiplication and an addition that never reach the tape.
            adouble* xin = xw;
            if ( p > 0 ) {
                for (int j=0;j<nstates;j++) xstg[j] = xw[j];
                for (int m = 0; m < p; m++) {
                    const double a = Atab[p*nstg + m];
                    if ( a == 0.0 ) continue;
                    adouble dta = a*dt;
                    const adouble* Km = K + m*nstates;
                    for (int j=0;j<nstates;j++) xstg[j] = xstg[j] + dta*Km[j];
                }
                xin = xstg;
            }

            problem.dae(K + p*nstates, pscr, xin, u, parameters, tp, xad, iphase, workspace);
            if (want_cost && problem.integrand_cost)
                Lstg[p] = problem.integrand_cost(xin, u, parameters, tp, xad, iphase, workspace);
        }

        if (workspace->enable_nlp_counters)
            workspace->solution->mesh_stats[ workspace->current_mesh_refinement_iteration-1 ].n_ode_rhs_evals += nstg;

        for (int p = 0; p < nstg; p++) {
            const double bp = btab[p];
            if ( bp == 0.0 ) continue;
            adouble dtb = bp*dt;
            const adouble* Kp = K + p*nstates;
            for (int j=0;j<nstates;j++) xw[j] = xw[j] + dtb*Kp[j];
            if (want_cost && problem.integrand_cost) *Lint = *Lint + dtb*Lstg[p];
        }

        t = t + dt;

        for (int q = 0; q < nsamp; q++) {
            if ( sample_step[q] != s+1 ) continue;
            if (xsamp) for (int j=0;j<nstates;j++) xsamp[q*nstates+j] = xw[j];
            if (tsamp) tsamp[q] = t;
            // The control at the end of this step. The last stage of both tableaux has c = 1,
            // so u already holds it; deriving it again here is what would let the sampled
            // control and the integrated control drift apart under a scheme added later, so it
            // is asked for rather than assumed.
            if (usamp && ncontrols > 0) {
                if ( varying_u ) eval_u( ((double) s + 1.0)/((double) nsteps) );
                for (int c=0;c<ncontrols;c++) usamp[q*ncontrols+c] = u[c];
            }
        }
    }

    for (int j=0;j<nstates;j++) xend[j] = xw[j];
}
