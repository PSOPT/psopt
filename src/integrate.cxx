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

// ===========================================================================================
// The two IMPLICIT tableaux, DERIVED rather than transcribed.
//
// The standing rule is that a mistyped Runge-Kutta coefficient does not fail -- it gives a
// scheme of lower order that still converges to the right answer -- so these are computed here
// from the conditions that define them, in the same order a derivation would take them, and
// were checked before use against properties none of the derivation used: row sums against c,
// the linear conditions B(k) holding to round-off up to the order and failing beyond it, the
// stage-order condition C(2), stiff accuracy (b equal to the last row of A), the stability
// function tending to zero at minus infinity (L-stability) and staying inside the unit disc on
// the imaginary axis (A-stability), the expansion of R(z) against the exponential series, and
// the observed order on a stiff nonlinear non-autonomous system in fifty-digit arithmetic.
//
// TR-BDF2 (Bank, Coughran, Fichtner, Grosse, Rose and Smith, 1985): a trapezoidal step to
// gamma followed by a BDF2 step to the end, which as an ESDIRK is three stages of order two.
// Its coefficients are short enough to check by hand, which is why it is here as well as the
// third-order table.
//
// ESDIRK3: four stages, order three. gamma is the root of x^3 - 3x^2 + 3x/2 - 1/6 in (1/3, 1),
// the value that makes a singly diagonally implicit scheme third order and A-stable; c2 = 2
// gamma is forced by C(2) at the second stage; c3 = 3/5 is the one free parameter; and a32,
// a31 and b then FOLLOW from C(2) and B(1..3), which is how they are obtained below.
// ===========================================================================================
static double ms_trbdf2_c[3], ms_trbdf2_b[3], ms_trbdf2_A[9];
static double ms_esdirk3_c[4], ms_esdirk3_b[4], ms_esdirk3_A[16];
static double ms_trbdf2_gamma = 0.0, ms_esdirk3_gamma = 0.0;
static bool   ms_implicit_built = false;

static void ms_build_implicit(void)
{
    if (ms_implicit_built) return;

    // ---- TR-BDF2 -------------------------------------------------------------------------
    {
        const double g = 2.0 - std::sqrt(2.0);          // the second abscissa
        const double d = g/2.0;                         // the repeated diagonal
        for (int q = 0; q < 9; q++) ms_trbdf2_A[q] = 0.0;
        ms_trbdf2_c[0] = 0.0;  ms_trbdf2_c[1] = g;  ms_trbdf2_c[2] = 1.0;
        ms_trbdf2_A[3*1 + 0] = d;  ms_trbdf2_A[3*1 + 1] = d;
        const double w = std::sqrt(2.0)/4.0;
        ms_trbdf2_A[3*2 + 0] = w;  ms_trbdf2_A[3*2 + 1] = w;  ms_trbdf2_A[3*2 + 2] = d;
        // Stiff accuracy IS the definition of b here: the last stage is the step.
        for (int q = 0; q < 3; q++) ms_trbdf2_b[q] = ms_trbdf2_A[3*2 + q];
        ms_trbdf2_gamma = d;
    }

    // ---- ESDIRK3 -------------------------------------------------------------------------
    {
        // gamma by Newton on its own defining cubic, from the standard bracket.
        double g = 0.435;
        for (int it = 0; it < 200; it++) {
            const double p  = g*g*g - 3.0*g*g + 1.5*g - 1.0/6.0;
            const double dp = 3.0*g*g - 6.0*g + 1.5;
            const double dg = p/dp;
            g -= dg;
            if ( fabs(dg) < 1.0e-16 ) break;
        }
        const double c3 = 0.6;
        for (int q = 0; q < 16; q++) ms_esdirk3_A[q] = 0.0;
        ms_esdirk3_c[0] = 0.0;  ms_esdirk3_c[1] = 2.0*g;
        ms_esdirk3_c[2] = c3;   ms_esdirk3_c[3] = 1.0;
        ms_esdirk3_A[4*1 + 0] = g;  ms_esdirk3_A[4*1 + 1] = g;
        // C(2) at the third stage: a32 c2 + gamma c3 = c3^2/2.
        const double a32 = ( 0.5*c3*c3 - g*c3 )/( 2.0*g );
        const double a31 = c3 - a32 - g;
        ms_esdirk3_A[4*2 + 0] = a31; ms_esdirk3_A[4*2 + 1] = a32; ms_esdirk3_A[4*2 + 2] = g;
        // B(1), B(2), B(3) for b1, b2, b3 with b4 = gamma: a 3x3 Vandermonde in (c1,c2,c3),
        // solved in closed form because c1 = 0 makes it small enough to write out.
        const double c2 = 2.0*g;
        const double r1 = 1.0 - g;
        const double r2 = 0.5 - g;
        const double r3 = 1.0/3.0 - g;
        // b2 c2 + b3 c3 = r2 ;  b2 c2^2 + b3 c3^2 = r3
        const double det = c2*c3*(c3 - c2);
        const double b2  = ( r2*c3*c3 - r3*c3 )/det;
        const double b3  = ( r3*c2 - r2*c2*c2 )/det;
        const double b1  = r1 - b2 - b3;
        ms_esdirk3_A[4*3 + 0] = b1; ms_esdirk3_A[4*3 + 1] = b2;
        ms_esdirk3_A[4*3 + 2] = b3; ms_esdirk3_A[4*3 + 3] = g;
        for (int q = 0; q < 4; q++) ms_esdirk3_b[q] = ms_esdirk3_A[4*3 + q];
        ms_esdirk3_gamma = g;
    }

    ms_implicit_built = true;
}


void ms_propagate_segment(adouble* xend, adouble* Lint, int k, adouble* xad, int iphase,
                          adouble& t0, adouble& tf, adouble* parameters, Workspace* workspace,
                          int nsteps_override,
                          adouble* xsamp, adouble* usamp, adouble* tsamp,
                          std::vector<adouble>* tau,
                          int nsamp_override)
{
    Prob& problem   = *workspace->problem;
    Alg&  algorithm = *workspace->algorithm;
    const int i         = iphase-1;
    const int nstates   = problem.phase[i].nstates;
    const int ncontrols = problem.phase[i].ncontrols;

    // How many steps this SEGMENT takes. Normally the user's ms_steps_per_segment, and under
    // ms_adaptive_steps whatever the step driver chose for this segment between solves --
    // frozen for the whole of this solve either way, which is the point of it. The override is
    // the error estimator asking for the same segment at half the step.
    int nsteps = ( nsteps_override > 0 ) ? nsteps_override
                                         : ms_segment_steps(i, k, workspace);
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
    double        gam_impl = 0.0;
    const bool    implicit = ms_implicit_integrator(algorithm);
    if      ( ms_rk8(algorithm) )    { ms_build_rk8(); Atab = ms_rk8_A; btab = ms_rk8_b; ctab = ms_rk8_c; }
    else if ( ms_trbdf2(algorithm) ) { ms_build_implicit(); Atab = ms_trbdf2_A; btab = ms_trbdf2_b;
                                       ctab = ms_trbdf2_c; gam_impl = ms_trbdf2_gamma; }
    else if ( ms_esdirk3(algorithm) ){ ms_build_implicit(); Atab = ms_esdirk3_A; btab = ms_esdirk3_b;
                                       ctab = ms_esdirk3_c; gam_impl = ms_esdirk3_gamma; }
    else                             { Atab = ms_rk4_A; btab = ms_rk4_b; ctab = ms_rk4_c; }

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

    // ===================================================================================
    // The algebraic block, when this phase declares one: the half-explicit scheme.
    //
    // For a semi-explicit index-1 system xdot = f(x,z,u,t), 0 = g(x,z,u,t), the differential
    // part of a stage is EXPLICIT -- the stage state is known from stages already taken --
    // so the only thing left to determine is z, and it is determined by n_z equations in n_z
    // unknowns and nothing else. That is the whole of the method.
    //
    // It keeps the tableau's own order, and that is a theorem rather than a hope. Index 1
    // means dg/dz is nonsingular, so the algebraic relation defines z = G(x,u,t) locally, and
    // an explicit Runge-Kutta applied to the REDUCED ordinary system xdot = f(x,G(x,u,t),u,t)
    // produces exactly the stage sequence written below. The two are the same computation,
    // not an approximation of one another, so RK8 carries a DAE at eighth order on the day it
    // carries an ODE at eighth order. Measured, in fifty-digit arithmetic outside PSOPT:
    // 8.40, 8.30, 8.18, 8.10 for successive halvings against the tableau's 8, and 4.16, 4.08,
    // 4.04, 4.02 against RK4's 4.
    //
    // And there is no drift. The algebraic relation is not an invariant that the integrator
    // is asked to preserve, it is an equation solved wherever z is defined at all, so |g| is
    // the inner iteration's residual and not a quantity that grows along the trajectory.
    // That is the structural difference from index reduction, which is the other way to reach
    // this class and leaves |g| to the integrator's own error.
    // ===================================================================================
    const int nalg    = ms_algebraic_vars(problem, i, algorithm);
    const int nfree_u = ncontrols - nalg;          // the controls the optimiser actually chooses
    int malg = algorithm.ms_algebraic_iterations;
    if ( malg < 1 ) malg = 1;

    std::vector<adouble> gres_( (nalg>0) ? nalg : 1 );
    std::vector<adouble> gnew_( (nalg>0) ? nalg : 1 );
    std::vector<adouble> Jb_  ( (nalg>0) ? nalg*nalg : 1 );
    std::vector<adouble> Jw_  ( (nalg>0) ? nalg*(nalg+1) : 1 );
    std::vector<adouble> dz_  ( (nalg>0) ? nalg : 1 );
    std::vector<adouble> Kscr_( (nalg>0 || implicit) ? nstates : 1 );
    std::vector<double>  gtar_( (nalg>0) ? nalg : 1, 0.0 );
    // The target is the path component's own bound rather than zero, so that a user who
    // writes the algebraic relation as g = c rather than g = 0 gets the equation solved
    // and not a different one. validate has already required lower == upper on these.
    for (int j = 0; j < nalg; j++) gtar_[j] = (problem.phase[i].bounds.lower.path)(j);
    bool alg_seeded = false;

    // The residual of the algebraic equations at the current (xin, u, tp). Every evaluation
    // costs a full call of the user's dae, because the user's dae computes the derivatives
    // and the path together; that is the price of the method and it is why the iteration
    // count matters. The derivative buffer is scratch here -- the stage's own derivatives are
    // taken from the LAST call, made after z has converged.
    auto alg_residual = [&](adouble* xin, adouble& tp, adouble* r) {
        problem.dae(Kscr_.data(), pscr, xin, u, parameters, tp, xad, iphase, workspace);
        for (int j = 0; j < nalg; j++) r[j] = pscr[j] - gtar_[j];
    };

    // Dense solve of Jb dz = r. Gaussian elimination with partial pivoting, the pivot chosen
    // on the taping-point values -- the same standing caveat any branch in a user's dae
    // carries, and harmless here because the pivot order only reorders exact arithmetic.
    auto alg_solve = [&](const adouble* J, const adouble* r, adouble* dz) {
        if ( nalg == 1 ) { dz[0] = r[0]/J[0]; return; }
        const int n = nalg;
        adouble* M = Jw_.data();
        for (int p = 0; p < n; p++) {
            for (int q = 0; q < n; q++) M[p*(n+1)+q] = J[p*n+q];
            M[p*(n+1)+n] = r[p];
        }
        for (int col = 0; col < n; col++) {
            int piv = col; double best = fabs( M[col*(n+1)+col].value() );
            for (int rr = col+1; rr < n; rr++) {
                const double v = fabs( M[rr*(n+1)+col].value() );
                if ( v > best ) { best = v; piv = rr; }
            }
            if ( piv != col )
                for (int cc = col; cc <= n; cc++) {
                    adouble tmp = M[col*(n+1)+cc];
                    M[col*(n+1)+cc] = M[piv*(n+1)+cc];
                    M[piv*(n+1)+cc] = tmp;
                }
            for (int rr = col+1; rr < n; rr++) {
                adouble fct = M[rr*(n+1)+col]/M[col*(n+1)+col];
                for (int cc = col; cc <= n; cc++)
                    M[rr*(n+1)+cc] = M[rr*(n+1)+cc] - fct*M[col*(n+1)+cc];
            }
        }
        for (int rr = n-1; rr >= 0; rr--) {
            adouble s = M[rr*(n+1)+n];
            for (int cc = rr+1; cc < n; cc++) s = s - M[rr*(n+1)+cc]*dz[cc];
            dz[rr] = s/M[rr*(n+1)+rr];
        }
    };

    // The stage solve: a FIXED, unrolled number of Broyden iterations.
    //
    // Fixed is the point. A loop whose length depends on the values cannot be taped, and --
    // the deeper objection -- would make the constraint function non-smooth in the decision
    // variables, because a convergence test flipping as the iterate moves changes the discrete
    // map the matching condition is written on. What is taped here is a fixed sequence of
    // arithmetic, so the derivative that comes back is the derivative of what was computed.
    // That is Bock's internal numerical differentiation applied one level further in.
    //
    // Broyden rather than Newton because Broyden needs no derivative of g: it builds its slope
    // from function values alone, so nothing inside the iteration has to be differentiated and
    // the nesting that a DAE capability is usually said to require never arises. The slope is
    // seeded once per segment by finite differences and updated thereafter; the seed's quality
    // does not matter -- measured, an identity seed gives the same answer to every digit --
    // because the updates correct it within an iteration.
    //
    // THE UPDATE IS DAMPED, and the damping is not a nicety. Once the iteration has converged
    // the step and the residual difference are both at round-off, and the update divides one
    // by the square of the other: undamped it replaces a good slope with noise. The symptom is
    // a count that appears to need to grow with the number of algebraic components --
    // measured, two coupled components stalled at four iterations and looked as though they
    // wanted five. Damped, four iterations give the full eighth order at one, two and three
    // components alike. It is written as arithmetic rather than a branch for a reason that is
    // specific to a taped computation; see the note at the update itself.
    auto solve_algebraic = [&](adouble* xin, adouble& tp) {
        if ( nalg <= 0 ) return;
        adouble* const r  = gres_.data();
        adouble* const rn = gnew_.data();
        adouble* const J  = Jb_.data();
        adouble* const dz = dz_.data();

        alg_residual(xin, tp, r);

        if ( !alg_seeded ) {
            for (int q = 0; q < nalg; q++) {
                const double zq = u[nfree_u+q].value();
                const double dl = 1.0e-7*(1.0 + fabs(zq));
                u[nfree_u+q] = u[nfree_u+q] + dl;
                alg_residual(xin, tp, rn);
                u[nfree_u+q] = u[nfree_u+q] - dl;
                for (int p = 0; p < nalg; p++) J[p*nalg+q] = (rn[p] - r[p])/dl;
            }
            alg_seeded = true;
            if (workspace->enable_nlp_counters)
                workspace->solution->mesh_stats[ workspace->current_mesh_refinement_iteration-1 ]
                    .n_ode_rhs_evals += nalg;
        }

        for (int it = 0; it < malg; it++) {
            alg_solve(J, r, dz);
            for (int j = 0; j < nalg; j++) u[nfree_u+j] = u[nfree_u+j] - dz[j];
            alg_residual(xin, tp, rn);

            // Broyden's first update with Delta z = -dz:
            //   J <- J + ((rn - r) + J dz)(-dz)^T/(||dz||^2 + eps)
            //
            // The eps is the whole of the guard, and it is arithmetic rather than a branch on
            // purpose. Once the iteration has converged, dz and (rn - r) are both at round-off
            // and the quotient is noise of order one, which replaces a good slope with a bad
            // one. A BRANCH cannot guard that: the condition is false at the taping point,
            // where the iterate is far from any solution, and true later -- so the tape would
            // record the unguarded division and then perform it at exactly the iterates where
            // it is unsafe. Measured, that is not a small effect: the cost stopped converging
            // in the step count and moved in the wrong direction.
            //
            // eps = (1e-8)^2 (1 + ||z||^2) leaves the update untouched wherever ||dz|| is
            // larger than about 1e-8 times the variable's own size, which covers every
            // iteration that is still making progress, and damps it smoothly to nothing below
            // that. The slope then stops improving at a relative accuracy of about 1e-8, which
            // is enough to carry the remaining iterations to machine precision.
            adouble s2 = 1.0e-16, z2 = 0.0;
            for (int j = 0; j < nalg; j++) {
                s2 = s2 + dz[j]*dz[j];
                z2 = z2 + u[nfree_u+j]*u[nfree_u+j];
            }
            s2 = s2 + 1.0e-16*z2;
            for (int p = 0; p < nalg; p++) {
                adouble Jd = 0.0;
                for (int q = 0; q < nalg; q++) Jd = Jd + J[p*nalg+q]*dz[q];
                adouble num = rn[p] - r[p] + Jd;
                for (int q = 0; q < nalg; q++)
                    J[p*nalg+q] = J[p*nalg+q] - num*dz[q]/s2;
            }
            for (int j = 0; j < nalg; j++) r[j] = rn[j];
        }

        if (workspace->enable_nlp_counters)
            workspace->solution->mesh_stats[ workspace->current_mesh_refinement_iteration-1 ]
                .n_ode_rhs_evals += malg + 1;
    };

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
    //
    // The algebraic components are NOT parameterised. They are the last nfree_u..ncontrols-1
    // slots, they are solved at every stage, and interpolating them between the node values
    // would be exactly the approximation this method exists to remove -- an algebraic variable
    // held or ramped across a segment while the state moves under it is the path-constraint
    // formulation, drifting at first order. So eval_u writes the free controls only, and the
    // algebraic slots carry the last solved value forward as the next stage's warm start.
    const int nparam_u = ( nalg > 0 ) ? nfree_u : ncontrols;
    auto eval_u = [&](double s) {
        if (nparam_u <= 0) return;
        if (quad_u) {
            const double w0 = (2.0*s-1.0)*(s-1.0);
            const double wm = 4.0*s*(1.0-s);
            const double w1 = s*(2.0*s-1.0);
            for (int c=0;c<nparam_u;c++) u[c] = w0*u0_[c] + wm*um_[c] + w1*u1_[c];
        }
        else if (linear_u) {
            for (int c=0;c<nparam_u;c++) u[c] = (1.0-s)*u0_[c] + s*u1_[c];
        }
        else {
            for (int c=0;c<nparam_u;c++) u[c] = u0_[c];
        }
    };

    const bool want_cost = ( Lint != NULL );
    if (want_cost) *Lint = 0.0;

    // Where the interior samples fall, as step indices. They are placed at step boundaries so
    // that the state at a sample is one the integrator actually produced, rather than an
    // interpolation of states either side of it.
    // How many interior samples to capture. Normally the user's ms_path_samples, because the
    // samples the constraints are imposed at and the samples captured have to be the same
    // points. The override exists for the refinement indicator, which needs to look inside a
    // segment at points where nothing is being enforced -- that being the whole of what it
    // measures -- without changing what the NLP sees.
    int nsamp = ( xsamp != NULL || usamp != NULL || tsamp != NULL )
                ? algorithm.ms_path_samples : 0;
    if ( nsamp_override > 0 && ( xsamp != NULL || usamp != NULL || tsamp != NULL ) )
        nsamp = nsamp_override;
    if ( nsamp > nsteps - 1 ) nsamp = nsteps - 1;
    if ( nsamp < 0 )          nsamp = 0;
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

    // The algebraic variables start from the segment's own start node, where they are decision
    // variables whose path row IS their algebraic equation -- so at a solution the warm start
    // is exact, and the finite-difference seed for the slope is taken at a point where the
    // residual is zero.
    for (int j = 0; j < nalg; j++) u[nfree_u+j] = u0_[nfree_u+j];

    // =======================================================================================
    // The implicit path: an ESDIRK stage solve, and the one factorisation that serves a step.
    //
    // For stage i the unknown is Y = (X, Z) with
    //
    //     X - S_i - h gamma f(X, Z, u, t_i) = 0,      g(X, Z, u, t_i) = 0,
    //
    // S_i being everything the earlier stages already fixed. The algebraic block, when there
    // is one, is folded into the SAME system rather than solved separately: for a stiff DAE
    // the differential and algebraic parts are coupled through the stiffness and solving them
    // in turn would be solving a different problem.
    //
    // W = I - h gamma J is formed by finite differences once per STEP, factorised once, and
    // reused by every stage. That is what a single repeated diagonal is for, and it is what a
    // stiff solver has always done: the stages differ in S and in the time, not in the matrix.
    // Nothing in it is differentiated to build it, so no automatic differentiation is nested;
    // the finite differences, the factorisation and the fixed iteration count are all ordinary
    // arithmetic and the tape records them as such.
    //
    // Stiff accuracy is why the step's state is the LAST STAGE rather than a b-weighted sum.
    // The two are equal in exact arithmetic, b being the last row of A, but the stage value is
    // the one that satisfies the algebraic constraint -- so on a DAE the constraint holds at
    // the step end, which is where the matching conditions live.
    const int nd_impl = implicit ? ( nstates + nalg ) : 1;
    std::vector<adouble> Wlu_( implicit ? nd_impl*nd_impl : 1 );
    std::vector<adouble> Rv_ ( nd_impl ), dY_( nd_impl ), Yv_( nd_impl ), Sv_( nstates );
    std::vector<int>     piv_( nd_impl, 0 );

    // The residual of the stage system at the trial Y. One call of the user's dae, which
    // returns the derivatives and the path together.
    auto impl_residual = [&](const adouble* Y, const adouble* S, adouble& hg, adouble& tp,
                             adouble* R) {
        for (int j = 0; j < nstates; j++) xstg[j] = Y[j];
        for (int q = 0; q < nalg;    q++) u[nfree_u+q] = Y[nstates+q];
        problem.dae(Kscr_.data(), pscr, xstg, u, parameters, tp, xad, iphase, workspace);
        for (int j = 0; j < nstates; j++) R[j] = Y[j] - S[j] - hg*Kscr_[j];
        for (int q = 0; q < nalg;    q++) R[nstates+q] = pscr[q] - gtar_[q];
    };

    // LU with partial pivoting, the pivot order decided on the taping-point values. The same
    // standing caveat as any branch inside a taped computation, and harmless here: a pivot
    // order reorders exact arithmetic, and W is nonsingular for any step small enough for the
    // scheme to be worth using.
    auto impl_factor = [&](adouble* W) {
        const int n = nd_impl;
        for (int q = 0; q < n; q++) piv_[q] = q;
        for (int col = 0; col < n; col++) {
            int    pv = col;
            double best = fabs( W[col*n+col].value() );
            for (int r = col+1; r < n; r++) {
                const double v = fabs( W[r*n+col].value() );
                if ( v > best ) { best = v; pv = r; }
            }
            if ( pv != col ) {
                for (int q = 0; q < n; q++) { adouble tmp = W[col*n+q];
                    W[col*n+q] = W[pv*n+q]; W[pv*n+q] = tmp; }
                const int ti = piv_[col]; piv_[col] = piv_[pv]; piv_[pv] = ti;
            }
            for (int r = col+1; r < n; r++) {
                adouble m = W[r*n+col]/W[col*n+col];
                W[r*n+col] = m;
                for (int q = col+1; q < n; q++) W[r*n+q] = W[r*n+q] - m*W[col*n+q];
            }
        }
    };
    auto impl_solve = [&](const adouble* W, const adouble* R, adouble* y) {
        const int n = nd_impl;
        for (int r = 0; r < n; r++) {
            adouble sum = R[ piv_[r] ];
            for (int q = 0; q < r; q++) sum = sum - W[r*n+q]*y[q];
            y[r] = sum;
        }
        for (int r = n-1; r >= 0; r--) {
            adouble sum = y[r];
            for (int q = r+1; q < n; q++) sum = sum - W[r*n+q]*y[q];
            y[r] = sum/W[r*n+r];
        }
    };

    int mimp = algorithm.ms_implicit_iterations;
    if ( mimp < 1 ) mimp = 1;

    adouble t = tk;
    for (int s = 0; s < nsteps; s++) {

      if ( implicit ) {

        // One factorisation for the whole step. W = I - h gamma J, with J formed by forward
        // differences on the same residual the stages will solve, at the step's own start --
        // where it is cheapest to be accurate, because every stage begins within O(h) of it.
        {
            adouble hg0 = gam_impl*dt;
            adouble t0s = t;
            for (int j = 0; j < nstates; j++) { Yv_[j] = xw[j]; Sv_[j] = xw[j]; }
            for (int q = 0; q < nalg;    q++) Yv_[nstates+q] = u[nfree_u+q];
            impl_residual(Yv_.data(), Sv_.data(), hg0, t0s, Rv_.data());
            for (int c2 = 0; c2 < nd_impl; c2++) {
                const double y0 = Yv_[c2].value();
                const double dl = 1.0e-7*(1.0 + fabs(y0));
                Yv_[c2] = Yv_[c2] + dl;
                impl_residual(Yv_.data(), Sv_.data(), hg0, t0s, dY_.data());
                Yv_[c2] = Yv_[c2] - dl;
                for (int r2 = 0; r2 < nd_impl; r2++)
                    Wlu_[r2*nd_impl + c2] = (dY_[r2] - Rv_[r2])/dl;
            }
            impl_factor(Wlu_.data());
            if (workspace->enable_nlp_counters)
                workspace->solution->mesh_stats[ workspace->current_mesh_refinement_iteration-1 ]
                    .n_ode_rhs_evals += nd_impl + 1;
        }

        // Stage zero is EXPLICIT -- the E of ESDIRK -- so it is the derivative at the step's
        // own state, with the algebraic variables it already carries. At the segment start
        // those come from the node, whose path row is their algebraic equation; afterwards
        // they come from the previous step's last stage, where stiff accuracy left them
        // satisfying the constraint.
        {
            const double sp0 = ( (double) s + ctab[0] )/((double) nsteps);
            if ( varying_u ) eval_u(sp0);
            adouble tp0 = t + ctab[0]*dt;
            problem.dae(K, pscr, xw, u, parameters, tp0, xad, iphase, workspace);
            if (want_cost && problem.integrand_cost)
                Lstg[0] = problem.integrand_cost(xw, u, parameters, tp0, xad, iphase, workspace);
        }

        for (int p = 1; p < nstg; p++) {
            const double sp = ( (double) s + ctab[p] )/((double) nsteps);
            if ( varying_u ) eval_u(sp);
            adouble tp = t + ctab[p]*dt;
            adouble hg = gam_impl*dt;

            // Everything the earlier stages already fixed.
            for (int j = 0; j < nstates; j++) Sv_[j] = xw[j];
            for (int m = 0; m < p; m++) {
                const double a = Atab[p*nstg + m];
                if ( a == 0.0 ) continue;
                adouble dta = a*dt;
                const adouble* Km = K + m*nstates;
                for (int j = 0; j < nstates; j++) Sv_[j] = Sv_[j] + dta*Km[j];
            }

            // Warm-started from the previous stage, which is O(h) away -- the reason a fixed
            // count works at all. Modified Newton with a Jacobian formed at the step's start
            // contracts by O(h) an iteration, so m iterations leave O(h^(m+1)).
            for (int j = 0; j < nstates; j++) Yv_[j] = xstg[j];
            for (int q = 0; q < nalg;    q++) Yv_[nstates+q] = u[nfree_u+q];
            if ( p == 1 ) for (int j = 0; j < nstates; j++) Yv_[j] = xw[j];

            for (int it = 0; it < mimp; it++) {
                impl_residual(Yv_.data(), Sv_.data(), hg, tp, Rv_.data());
                impl_solve(Wlu_.data(), Rv_.data(), dY_.data());
                for (int q = 0; q < nd_impl; q++) Yv_[q] = Yv_[q] - dY_[q];
            }

            // The stage's derivative, at the converged stage value.
            for (int j = 0; j < nstates; j++) xstg[j] = Yv_[j];
            for (int q = 0; q < nalg;    q++) u[nfree_u+q] = Yv_[nstates+q];
            problem.dae(K + p*nstates, pscr, xstg, u, parameters, tp, xad, iphase, workspace);
            if (want_cost && problem.integrand_cost)
                Lstg[p] = problem.integrand_cost(xstg, u, parameters, tp, xad, iphase, workspace);

            if (workspace->enable_nlp_counters)
                workspace->solution->mesh_stats[ workspace->current_mesh_refinement_iteration-1 ]
                    .n_ode_rhs_evals += mimp + 1;
        }

        // Stiff accuracy: the step IS the last stage. Not the b-weighted sum, which is the
        // same number in exact arithmetic and is not the vector that satisfies the algebraic
        // constraint. The running cost still takes the b-weighted sum, that being a quadrature
        // and not a state.
        for (int j = 0; j < nstates; j++) xw[j] = xstg[j];
        if (want_cost && problem.integrand_cost)
            for (int p = 0; p < nstg; p++) {
                const double bp = btab[p];
                if ( bp == 0.0 ) continue;
                *Lint = *Lint + (bp*dt)*Lstg[p];
            }

      }
      else {

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

            // The algebraic variables first, then the derivatives at the point they define.
            // The order matters: f is evaluated at the z that satisfies g, never the other
            // way round, which is what makes the stage a point of the reduced system.
            if ( nalg > 0 ) solve_algebraic(xin, tp);

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
                // The algebraic variables at a sample belong to the state at the sample, and
                // the state at the END of a step is not any stage's state -- it is the
                // b-weighted combination of them -- so the last stage's z does not answer for
                // it. Solved again here, at the point the sample actually is. A path
                // constraint imposed at a sample would otherwise be imposed at a z that
                // satisfies g nowhere.
                // Under an implicit scheme there is nothing to re-solve: stiff accuracy
                // means the step end IS the last stage, so u already carries the algebraic
                // variables that satisfy the constraint there.
                if ( nalg > 0 && !implicit ) solve_algebraic(xw, t);
                for (int c=0;c<ncontrols;c++) usamp[q*ncontrols+c] = u[c];
            }
        }
    }

    for (int j=0;j<nstates;j++) xend[j] = xw[j];
}
