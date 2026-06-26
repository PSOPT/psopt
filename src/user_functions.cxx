/*********************************************************************************************

This file is part of the PSOPT library, a software tool for computational optimal control

Copyright (C) 2009-2020 Victor M. Becerra

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
           School of Energy and Electronic Engineering
           Portsmouth PO1 3DJ
           United Kingdom
e-mail:    v.m.becerra@ieee.org

**********************************************************************************************/


#include "psopt.h"

// Bring std names into this translation unit (formerly leaked via psopt.h).
using namespace std;



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
 get_individual_control_trajectory(single_control_traj, control_index, iphase, xad, workspace);
 get_times( &t0, &tf, xad, iphase, workspace);
 for (k=0; k<norder+1; k++) { // EIGEN_UPDATE
	time_array[k]  =  convert_to_original_time_ad( (workspace->snodes[i])(k), t0, tf );
 }
 if ( time-delay>t0 ) // Careful because this if-then statement may not be differentiable
        delayed_time= time-delay;
 else
          delayed_time=t0; // what is best to do here?

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
 get_individual_state_trajectory(single_state_traj, state_index, iphase, xad, workspace);
 get_times( &t0, &tf, xad, iphase, workspace);
 for (k=0; k<norder+1; k++) { // EIGEN_UPDATE
        ts = (workspace->snodes[i])(k);
	     time_array[k]  =  convert_to_original_time_ad( ts, t0, tf );
 }
 if ( time-delay>t0 )
        delayed_time= time-delay;
 else
          delayed_time=t0; // nothing best to do here...

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

