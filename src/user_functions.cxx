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
e-mail:    v.m.becerra@ieee.org

**********************************************************************************************/


#include "psopt.h"

// Bring std names into this translation unit (formerly leaked via psopt.h).
using namespace std;




// ---------------------------------------------------------------------------
// Initial history of a phase with delayed states or controls.
//
// For t in [t0-d, t0) the delayed accessors below are asked for a trajectory value at a
// time before the phase begins. A delay differential equation is not well posed without
// an initial history function there, so PSOPT needs one of two things: either the user
// supplies problem.initial_history, or PSOPT falls back on the constant history
// phi(t) = x(t0), psi(t) = u(t0), which is what every earlier release did silently and
// which keeps the state continuous at t0.
//
// Returns true when the user's history supplied the value, in which case *value holds it.
// ---------------------------------------------------------------------------

static bool evaluate_initial_history(adouble* value, int index, bool want_state,
                                     int iphase, adouble& delayed_time, adouble* xad,
                                     Workspace* workspace)
{
   Prob& problem = *workspace->problem;

   if ( problem.initial_history == NULL ) return false;

   int i = iphase-1;

   adouble* hstates   = workspace->history_states[i].get();
   adouble* hcontrols = workspace->history_controls[i].get();
   adouble* params    = workspace->parameters[i].get();

   get_parameters( params, xad, iphase, workspace );

   for (int k=0; k<problem.phase[i].nstates;   k++) hstates[k]   = 0.0;
   for (int k=0; k<problem.phase[i].ncontrols; k++) hcontrols[k] = 0.0;

   problem.initial_history( hstates, hcontrols, params, delayed_time, xad, iphase, workspace );

   *value = want_state ? hstates[index] : hcontrols[index];

   return true;
}

// A one-off notice when the initial time of a phase using delayed variables is free. The
// comparison that detects t-d < t0 is taken on a taped quantity, so it is decided once,
// when the tape is recorded, and stays decided as t0 moves during the solve. With t0
// fixed -- as it is in every delayed example distributed with PSOPT -- this is harmless.

static void warn_if_initial_time_is_free(int iphase, Workspace* workspace)
{
   Prob& problem   = *workspace->problem;
   Alg&  algorithm = *workspace->algorithm;
   int   i         = iphase-1;

   if ( workspace->delay_free_time_warned ) return;
   if ( !useAutomaticDifferentiation(algorithm) ) return;

   double lo = problem.phase[i].bounds.lower.StartTime;
   double up = problem.phase[i].bounds.upper.StartTime;

   if ( lo == up ) return;

   workspace->delay_free_time_warned = true;

   snprintf(workspace->text,sizeof(workspace->text),
      "\n>>> Warning: phase %d uses delayed variables and has a free initial time."
      "\n    The test that decides whether a delayed time falls before t0 is recorded on"
      "\n    the derivative tape, so it is fixed at the value t0 had when the tape was"
      "\n    taken. If t0 moves across a delay boundary during the solve the derivatives"
      "\n    will be silently wrong. Fix t0, or use algorithm.derivatives = \"numerical\".",
      iphase);
   psopt_print(workspace,workspace->text);
}

void get_delayed_control(adouble* delayed_control, int control_index, int iphase, adouble& time, double delay, adouble* xad, Workspace* workspace)
{

 int k;
 int i = iphase-1;
 Prob& problem = *workspace->problem;
 int norder = problem.phase[i].current_number_of_intervals;
 adouble t0, tf;
 adouble delayed_time;
 adouble* time_array = workspace->time_array_tmp.get();;
 adouble* single_control_traj = workspace->single_trajectory_tmp.get();

 // control_index counts from zero, as every other index of the user interface does.
 // Releases before 2026 documented it as counting from one, so code written to that
 // description reads the wrong control, silently, whenever the index happens to be in
 // range. The check below catches the case where it is not.
 if ( control_index < 0 || control_index >= problem.phase[i].ncontrols ) {
     char msg[256];
     snprintf(msg,sizeof(msg),
        "get_delayed_control: control index %d is out of range for phase %d, which has %d "
        "controls. The index counts from zero.", control_index, iphase, problem.phase[i].ncontrols);
     error_message(msg);
 }

 get_individual_control_trajectory(single_control_traj, control_index, iphase, xad, workspace);
 get_times( &t0, &tf, xad, iphase, workspace);
 for (k=0; k<norder+1; k++) { // EIGEN_UPDATE
	time_array[k]  =  convert_to_original_time_ad( (workspace->snodes[i])(k), t0, tf );
 }
 // The test is taken on a taped quantity; see warn_if_initial_time_is_free.
 warn_if_initial_time_is_free(iphase, workspace);

 // The history is defined on [t0-d, t0); at a delayed time of exactly t0 the value
 // wanted is the trajectory's own, u(t0). The distinction is invisible with the default
 // constant history, whose value there is u(t0) anyway, but a user history that steps at
 // t0 would otherwise be applied at one node too many -- worth h/6 of the jump in a
 // Hermite-Simpson quadrature, which is far larger than the discretization error.
 if ( time-delay >= t0 ) {
        delayed_time = time-delay;
 }
 else {
        delayed_time = time-delay;
        if ( evaluate_initial_history( delayed_control, control_index, false,
                                       iphase, delayed_time, xad, workspace ) ) return;
        delayed_time = t0;   // default: the constant history psi(t) = u(t0)
 }

 spline_interpolation( delayed_control, delayed_time, time_array, single_control_traj, norder+1, workspace);

}

void get_delayed_state(adouble* delayed_state, int state_index, int iphase, adouble& time, double delay, adouble* xad, Workspace* workspace)
{

 int k;
 int i= iphase-1;
 Prob& problem = *workspace->problem;
 Alg&  algorithm=*workspace->algorithm;
 int norder = problem.phase[i].current_number_of_intervals;
 adouble t0, tf;
 double ts;
 adouble delayed_time;
 adouble* time_array = workspace->time_array_tmp.get();
 adouble* single_state_traj =  workspace->single_trajectory_tmp.get();

 // state_index counts from zero; see the note in get_delayed_control.
 if ( state_index < 0 || state_index >= problem.phase[i].nstates ) {
     char msg[256];
     snprintf(msg,sizeof(msg),
        "get_delayed_state: state index %d is out of range for phase %d, which has %d "
        "states. The index counts from zero.", state_index, iphase, problem.phase[i].nstates);
     error_message(msg);
 }

 get_individual_state_trajectory(single_state_traj, state_index, iphase, xad, workspace);
 get_times( &t0, &tf, xad, iphase, workspace);
 for (k=0; k<norder+1; k++) { // EIGEN_UPDATE
        ts = (workspace->snodes[i])(k);
	     time_array[k]  =  convert_to_original_time_ad( ts, t0, tf );
 }
 warn_if_initial_time_is_free(iphase, workspace);

 // See the note in get_delayed_control on the boundary case.
 if ( time-delay >= t0 ) {
        delayed_time = time-delay;
 }
 else {
        delayed_time = time-delay;
        if ( evaluate_initial_history( delayed_state, state_index, true,
                                       iphase, delayed_time, xad, workspace ) ) return;
        delayed_time = t0;   // default: the constant history phi(t) = x(t0)
 }

 if ( use_global_collocation(algorithm) &&  norder<100  ) {
 	lagrange_interpolation_ad( delayed_state, delayed_time, time_array, single_state_traj, norder+1, workspace);
 }
 else  {
	spline_interpolation( delayed_state, delayed_time, time_array, single_state_traj, norder+1, workspace);
 }

}

void get_interpolated_state(adouble* interp_state, int state_index, int iphase, adouble& time, adouble* xad, Workspace* workspace)
{

 int k;
 int i = iphase-1;
 Prob& problem = *workspace->problem;
 Alg&  algorithm=*workspace->algorithm;
 int norder = problem.phase[i].current_number_of_intervals;
 adouble t0, tf;
 double ts;

 adouble* time_array = workspace->time_array_tmp.get();
 adouble* single_state_traj =  workspace->single_trajectory_tmp.get();
 get_individual_state_trajectory(single_state_traj, state_index, iphase, xad, workspace);
 get_times( &t0, &tf, xad, iphase, workspace);
 for (k=0; k<norder+1; k++) { // EIGEN_UPDATE
        ts = (workspace->snodes[i])(k);
	     time_array[k]  =  convert_to_original_time_ad( ts, t0, tf );
 }

 if (  use_global_collocation(algorithm) && hp_mesh_active(problem.phase[i]) ) {
        // hp multi-interval mesh: a global Lagrange interpolant over the clustered composite
        // nodes is ill-conditioned and blows up the error estimator. Interpolate only on the
        // local interval containing `time`. Because this is per-interval (low local order) it
        // is well-conditioned at any total node count, so no norder cap is needed here.
        bool   gauss = ( algorithm.collocation_method == "Gauss" );
        int    K  = hp_num_intervals(problem.phase[i]);
        double tq = time.value();
        int    s  = 0;
        for (int jint=0; jint<K; jint++) {
            int ord = hp_interval_order(problem.phase[i], jint);
            // Routing boundary and interpolation node set differ for Gauss. The LG state on an
            // interval is the degree-(ord) polynomial through its support nodes [left breakpoint,
            // ord Gauss points] - the right breakpoint is a SEPARATE variable fixed by the
            // quadrature defining constraint, NOT a support node, so it must be excluded from the
            // interpolant (including it forces a degree ord+1 fit that oscillates to pass through
            // the quadrature-derived value, inflating the residual near interfaces). The closing
            // segment [last Gauss point, right breakpoint] is then a forward evaluation of this
            // same polynomial, so routing still advances at the right breakpoint.
            int en  = s + ord;                       // last interpolation node (degree-ord poly)
            int rb  = s + ord + (gauss ? 1 : 0);     // routing boundary: right breakpoint (Gauss)
            if ( rb > norder ) rb = norder;          // last Gauss interval: x_f is not stored here
            if ( en > norder ) en = norder;
            if ( tq <= time_array[rb].value() || jint == K-1 ) {
                lagrange_interpolation_ad( interp_state, time, time_array+s,
                                           single_state_traj+s, (en-s)+1, workspace );
                return;
            }
            s = rb;                                  // advance to the right breakpoint
        }
 	lagrange_interpolation_ad( interp_state, time, time_array, single_state_traj, norder+1, workspace);
 }
 else if (  use_global_collocation(algorithm) && norder<100 ) {
 	lagrange_interpolation_ad( interp_state, time, time_array, single_state_traj, norder+1, workspace);
 }
 else  {
	spline_interpolation( interp_state, time, time_array, single_state_traj, norder+1, workspace);
 }

}


void get_interpolated_control(adouble* interp_control, int control_index, int iphase, adouble& time, adouble* xad, Workspace* workspace)
{

 int k;
 int i = iphase-1;
 Prob& problem = *workspace->problem;
 Alg&  algorithm = *workspace->algorithm;
 int norder = problem.phase[i].current_number_of_intervals;
 adouble t0, tf;
 double ts;

 adouble* time_array = workspace->time_array_tmp.get();
 adouble* single_control_traj =  workspace->single_trajectory_tmp.get();
 get_individual_control_trajectory(single_control_traj, control_index, iphase, xad, workspace);
 get_times( &t0, &tf, xad, iphase, workspace);
 for (k=0; k<norder+1; k++) { // EIGEN_UPDATE
        ts = (workspace->snodes[i])(k);
	     time_array[k]  =  convert_to_original_time_ad( ts, t0, tf );
 }

 // For a Gauss hp mesh the control is collocated only at the interior Gauss points; the
 // (non-collocated) breakpoint storage slots carry control variables that enter no collocation
 // condition and are left arbitrary by the NLP. A global spline through them corrupts the
 // interpolated control near every interface and badly inflates the error estimate there.
 // Interpolate instead within the interval containing `time`, using only that interval's
 // Gauss-point controls. Radau (breakpoints collocated) keeps the global spline: bit-identical.
 if ( use_global_collocation(algorithm)
      && hp_mesh_active(problem.phase[i]) && algorithm.collocation_method=="Gauss" ) {
     int K = hp_num_intervals(problem.phase[i]); double tq = time.value(); int s = 0;
     for (int jint=0; jint<K; jint++) {
         int ord = hp_interval_order(problem.phase[i], jint);
         int rb  = s + ord + 1; if (rb > norder) rb = norder;   // right breakpoint (x_f if last)
         if ( tq <= time_array[rb].value() || jint == K-1 ) {
             // this interval's Gauss-point controls occupy storage indices [s+1 .. s+ord]
             lagrange_interpolation_ad( interp_control, time, time_array+(s+1),
                                        single_control_traj+(s+1), ord, workspace );
             return;
         }
         s = rb;
     }
 }

 // Hermite-Simpson carries a control variable of its own at the midpoint of every
 // interval, and the control representation the transcription assumes is the quadratic
 // through (u_k, ubar_k, u_{k+1}) on interval k -- that is the control the Simpson
 // defect integrates. A spline through the node values alone therefore ignores half of
 // the control variables. Wherever the node and midpoint branches separate -- on
 // bang-bang, singular and sliding arcs the separation is systematic and O(1), because
 // there the defects pin only a weighted combination of the two and leave the split
 // free -- such a spline describes a control the solver never used, and the residual
 // computed from it is an artefact of the interpolation rather than a property of the
 // solution. Since this routine feeds the discretization error estimate, and that
 // estimate drives mesh refinement, the artefact is not merely cosmetic: it refines
 // where there is nothing to refine and refuses to stop. Interpolate on the interval
 // containing `time` instead, using that interval's midpoint control.
 if ( need_midpoint_controls(algorithm, workspace) && norder >= 1 ) {
     int    nstates   = problem.phase[i].nstates;
     int    ncontrols = problem.phase[i].ncontrols;
     int    nparam    = problem.phase[i].nparameters;
     int    iphase_offset = get_iphase_offset(problem, iphase, workspace);
     int    bar_offset    = (nstates+ncontrols)*(norder+1)+nparam;
     MatrixXd& control_scaling = problem.phase[i].scale.controls;
     double tq = time.value();
     int    kk = 0;
     while ( kk < norder-1 && tq > time_array[kk+1].value() ) kk++;
     adouble ubar = xad[iphase_offset+bar_offset+(kk)*ncontrols+control_index]
                    /control_scaling(control_index);
     adouble hk   = time_array[kk+1] - time_array[kk];
     adouble s    = (time - time_array[kk])/hk;
     *interp_control =   (2.0*s-1.0)*(s-1.0)*single_control_traj[kk]
                       + 4.0*s*(1.0-s)*ubar
                       + s*(2.0*s-1.0)*single_control_traj[kk+1];
     return;
 }

 spline_interpolation( interp_control, time, time_array, single_control_traj, norder+1, workspace);

}

void get_control_derivative(adouble* control_derivative, int control_index, int iphase, adouble& time, adouble* xad, Workspace* workspace)
{
// This function computes an approximation to the time derivative of a specified control variable, based
// on an interpolated finite difference.
     adouble t0, tf;
     adouble pcontrol;
     adouble control;
     adouble ptime;

     get_times(&t0, &tf, xad, iphase, workspace);

     double h = sqrt(PSOPT_extras::GetEPS());

     if ( time == tf ) {
	     h = -h;
     }

     ptime = time + h;

     get_interpolated_control(&control, control_index, iphase, time, xad, workspace);

     get_interpolated_control(&pcontrol, control_index, iphase, ptime, xad, workspace);

     *control_derivative = (pcontrol - control)/h;

}

void get_state_derivative(adouble* state_derivative, int state_index, int iphase, adouble& time, adouble* xad, Workspace* workspace)
{
// This function computes an approximation to the time derivative of a specified state variable, based
// on an interpolated finite difference.
     adouble t0, tf;
     adouble pstate;
     adouble state;
     adouble ptime;

     get_times(&t0, &tf, xad, iphase, workspace);

     double h = sqrt(PSOPT_extras::GetEPS());

     if ( time == tf ) {
	     h = -h;
     }

     ptime = time + h;

     get_interpolated_state(&state, state_index, iphase, time, xad, workspace);

     get_interpolated_state(&pstate, state_index, iphase, ptime, xad, workspace);

     *state_derivative = (pstate - state)/h;

}

