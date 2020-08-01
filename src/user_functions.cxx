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



void get_delayed_control(adouble* delayed_control, int control_index, int iphase, adouble& time, double delay, adouble* xad, Workspace* workspace)
{

 int k;
 int i = iphase-1;
 Prob& problem = *workspace->problem;
 int norder = problem.phase[i].current_number_of_intervals;
 adouble t0, tf;
 adouble delayed_time;
 adouble* time_array = workspace->time_array_tmp;;
 adouble* single_control_traj = workspace->single_trajectory_tmp;
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
 adouble* time_array = workspace->time_array_tmp;
 adouble* single_state_traj =  workspace->single_trajectory_tmp;
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

 adouble* time_array = workspace->time_array_tmp;
 adouble* single_state_traj =  workspace->single_trajectory_tmp;
 get_individual_state_trajectory(single_state_traj, state_index, iphase, xad, workspace);
 get_times( &t0, &tf, xad, iphase, workspace);
 for (k=0; k<norder+1; k++) { // EIGEN_UPDATE
        ts = (workspace->snodes[i])(k);
	     time_array[k]  =  convert_to_original_time_ad( ts, t0, tf );
 }

 if (  use_global_collocation(algorithm) && norder<100 ) {
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
 int norder = problem.phase[i].current_number_of_intervals;
 adouble t0, tf;
 double ts;

 adouble* time_array = workspace->time_array_tmp;
 adouble* single_control_traj =  workspace->single_trajectory_tmp;
 get_individual_control_trajectory(single_control_traj, control_index, iphase, xad, workspace);
 get_times( &t0, &tf, xad, iphase, workspace);
 for (k=0; k<norder+1; k++) { // EIGEN_UPDATE
        ts = (workspace->snodes[i])(k);
	     time_array[k]  =  convert_to_original_time_ad( ts, t0, tf );
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

