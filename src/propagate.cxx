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

void euler_propagate( void (*dae)(adouble* derivatives, adouble* path, adouble* states,
         adouble* controls, adouble* parameters, adouble& time,
        adouble* xad, int iphase, Workspace* workspace),
        MatrixXd& control_trajectory,
        MatrixXd& time_vector,
        Prob & problem,
        MatrixXd& initial_state,
	MatrixXd& parameters,
        MatrixXd& state_trajectory,
        int iphase, Workspace* workspace)
{
        int nsteps = length(time_vector) - 1;

	int nstates = problem.phases(iphase).nstates;

	int ncontrols = problem.phases(iphase).ncontrols;

	int nparam    = problem.phases(iphase).nparameters;

	int npath = problem.phases(iphase).npath;

	adouble *states =  new adouble[nstates];

	adouble *controls = new adouble[ncontrols];

	adouble *derivatives = new adouble[nstates];

	adouble *path = new adouble[npath];

	adouble *param = new adouble[nparam];

	adouble *xad = NULL;

	adouble time;

	int i, k;

	for(i=1;i<=nparam;i++)  param[i-1] = parameters(i);


   state_trajectory.col(0) = initial_state;

	for (k=0; k< nsteps; k++) // EIGEN_UPDATE
	{

	     for(i=0;i<nstates;i++)  states[i] = state_trajectory(i,k);

	     for(i=0;i<ncontrols;i++)  controls[i] = control_trajectory(i,k);

	     time = time_vector(k);

	     dae( derivatives, path, states, controls, param, time, xad, iphase, workspace);

	     double h = time_vector(k+1)-time_vector(k);

	     for(i=0;i<nstates;i++) {
	        state_trajectory(i,k+1) = state_trajectory(i,k) + h*derivatives[i].value();
	     }
	}
}

void rk4_propagate( void (*dae)(adouble* derivatives, adouble* path, adouble* states,
        adouble* controls, adouble* parameters, adouble& time,
        adouble* xad, int iphase, Workspace* workspace),
        MatrixXd& control_trajectory,
        MatrixXd& time_vector,
        MatrixXd& initial_state,
	     MatrixXd& parameters,
        Prob & problem,
        int iphase,
        MatrixXd& state_trajectory, Workspace* workspace)
{
        // 4th order Runge Kutta algorithm with given time vector and control sequence
   
   int nsteps = length(time_vector) - 1;

	int nstates = problem.phases(iphase).nstates;

	int ncontrols = problem.phases(iphase).ncontrols;

	int nparam    = problem.phases(iphase).nparameters;

	int npath = problem.phases(iphase).npath;

	adouble *states =  new adouble[nstates];

	adouble *controls = new adouble[ncontrols];

	adouble *derivatives = new adouble[nstates];

	adouble *path = new adouble[npath];

	adouble *param = new adouble[nparam];

	adouble *xad = NULL;

	adouble time;

	state_trajectory.resize(nstates,nsteps+1);

	MatrixXd K1(nstates,1), K2(nstates,1), K3(nstates,1), K4(nstates,1);

	int i, k;

	for(i=0;i<nparam;i++)  param[i-1] = parameters(i); // EIGEN_UPDATE


   state_trajectory.col(0)= initial_state;

	for (k=0; k< nsteps; k++)
	{

	     for(i=0;i<nstates;i++)  states[i] = state_trajectory(i,k);

	     for(i=0;i<ncontrols;i++)  controls[i] = control_trajectory(i,k);

	     time = time_vector(k);

	     dae( derivatives, path, states, controls, param, time, xad, iphase, workspace);

	     double h = time_vector(k+1)-time_vector(k);

	     for(i=0;i<nstates;i++) { // EIGEN_UPDATE
	        K1(i,0) =  h*derivatives[i].value();
	     }

        time = time + h/2.0;

	     for(i=0;i<nstates;i++)  states[i] = state_trajectory(i,k) + K1(i)/2.0;

	     for(i=0;i<ncontrols;i++)  controls[i] = ( control_trajectory(i,k) + control_trajectory(i,k+1) )/2.0;

	     dae( derivatives, path, states, controls, param, time, xad, iphase, workspace);

	     for(i=0;i<nstates;i++) {
	        K2(i,0) = h*derivatives[i].value();
	     }

	     for(i=0;i<nstates;i++)  states[i] = state_trajectory(i,k) + K2(i)/2.0;

	     dae( derivatives, path, states, controls, param, time, xad, iphase, workspace);

	     for(i=0;i<nstates;i++) {
	        K3(i,0) = h*derivatives[i].value();
	     }

        time = time + h/2.0;

	     for(i=0;i<nstates;i++)  states[i] = state_trajectory(i,k) + K3(i);

	     for(i=0;i<ncontrols;i++)  controls[i] = control_trajectory(i,k+1);

	     dae( derivatives, path, states, controls, param, time, xad, iphase, workspace);

	     for(i=0;i<nstates;i++) {
	        K4(i,0) = h*derivatives[i].value();
	     }

        K1 = (K1 + 2*K2 + 2*K3 + K4)/6.0;

	     for(i=0;i<nstates;i++) {
	        state_trajectory(i,k+1) = state_trajectory(i,k)  + K1(i);
	     }


	}

	delete[] states;
	delete[] controls;
	delete[] path;
	delete[] param;
	delete[] derivatives;
}




void rkf_propagate( void (*dae)(adouble* derivatives, adouble* path, adouble* states,
        adouble* controls, adouble* parameters, adouble& time,
        adouble* xad, int iphase, Workspace* workspace),
        MatrixXd& control_trajectory,
        MatrixXd& time_vector,
        MatrixXd& initial_state,
	     MatrixXd& parameters,
        double tolerance,
        double hmin,
	     double hmax,
        Prob & problem,
        int iphase,
        MatrixXd& state_trajectory,
        MatrixXd& new_time_vector,
	     MatrixXd& new_control_trajectory, Workspace* workspace)
{
// Runge-Kutta-Fehlberg method with variable step size with local truncation error within a given tolerance.
// Reference: Burden (2005) "Numerical Analysis", page 287.

	int nstates = problem.phases(iphase).nstates;

	int ncontrols = problem.phases(iphase).ncontrols;

	int nparam    = problem.phases(iphase).nparameters;

	int npath = problem.phases(iphase).npath;

	adouble *states =  new adouble[nstates];

	adouble *controls = new adouble[ncontrols];

	adouble *derivatives = new adouble[nstates];

	adouble *path = new adouble[npath];

	adouble *param = new adouble[nparam];

	adouble *xad = NULL;

	adouble time, timep;

	int flag = 1;

	MatrixXd K1(nstates,1), K2(nstates,1), K3(nstates,1), K4(nstates,1), K5(nstates,1), K6(nstates,1);

	MatrixXd controlp(ncontrols,1);

	int i, k;

	for(i=0;i<nparam;i++)  param[i] = parameters(i); // EIGEN_UPDATE

   state_trajectory.col(0) = initial_state;

	time = time_vector(0);

	double h = hmax;

	double tf = time_vector(length(time_vector)-1);

	k = 1;

	double maxR, delta;

	MatrixXd R;

   new_time_vector(0) = time_vector(0); // Line added 16 Feb. 2011.


	while (flag==1) {

	     //////////// K1 CALCULATIONS ////////////////////////
	     for(i=0;i<nstates;i++)  states[i] = state_trajectory(i,k);

	     timep = time;

	     if ( ncontrols>0 ) linear_interpolation(controlp, timep.value(), time_vector, control_trajectory, length(time_vector));

	     for(i=0;i<ncontrols;i++)  controls[i] = controlp(i);

	     dae( derivatives, path, states, controls, param, timep, xad, iphase, workspace);

	     for(i=0;i<nstates;i++) {
	        K1(i) =  h*derivatives[i].value();
	     }

        /////////////// K2 CALCULATIONS ////////////////////
             timep = time + h/4.0;

	     for(i=0;i<nstates;i++)  states[i] = state_trajectory(i,k) + K1(i)/4.0;

	     if ( ncontrols>0 ) linear_interpolation(controlp, timep.value(), time_vector, control_trajectory, length(time_vector));


	     for(i=0;i<ncontrols;i++)  controls[i] = controlp(i);

	     dae( derivatives, path, states, controls, param, timep, xad, iphase, workspace);

	     for(i=0;i<nstates;i++) {
	        K2(i) = h*derivatives[i].value();
	     }

             /////////////// K3 CALCULATIONS //////////////////////////
             timep = time + 3.0*h/8.0;

	     for(i=0;i<nstates;i++)  states[i] = state_trajectory(i,k) + (3.0/32.0)*K1(i) + (9.0/32.0)*K2(i);

	     if ( ncontrols>0 ) linear_interpolation(controlp, timep.value(), time_vector, control_trajectory, length(time_vector));

	     for(i=0;i<ncontrols;i++)  controls[i] = controlp(i);

	     dae( derivatives, path, states, controls, param, timep, xad, iphase, workspace);

	     for(i=0;i<nstates;i++) {
	        K3(i) = h*derivatives[i].value();
	     }

	     /////////////// K4 CALCULATIONS //////////////////////////
        timep = time + 12.0*h/13.0;

	     for(i=0;i<nstates;i++)  states[i] = state_trajectory(i,k) + (1932.0/2197.0)*K1(i) - (7200.0/2197.0)*K2(i) + (7296.0/2197.0)*K3(i);

	     if ( ncontrols>0 ) linear_interpolation(controlp, timep.value(), time_vector, control_trajectory, length(time_vector));

	     for(i=0;i<ncontrols;i++)  controls[i] = controlp(i);

	     dae( derivatives, path, states, controls, param, timep, xad, iphase, workspace);

	     for(i=0;i<nstates;i++) {
	        K4(i) = h*derivatives[i].value();
	     }


	     /////////////// K5 CALCULATIONS //////////////////////////
        timep = time + h;

	     for(i=0;i<nstates;i++)  states[i] = state_trajectory(i,k) + (439.0/216.0)*K1(i) - (8.0)*K2(i) + (3680.0/513.0)*K3(i) - (845.0/4104.0)*K4(i);

	     if ( ncontrols>0 ) linear_interpolation(controlp, timep.value(), time_vector, control_trajectory, length(time_vector));

	     for(i=0;i<ncontrols;i++)  controls[i] = controlp(i);

	     dae( derivatives, path, states, controls, param, timep, xad, iphase, workspace);

	     for(i=0;i<nstates;i++) {
	        K5(i) = h*derivatives[i].value();
	     }

	     /////////////// K6 CALCULATIONS //////////////////////////
        timep = time + h/2.0;

	     for(i=0;i<nstates;i++)  states[i] = state_trajectory(i,k) - (8.0/27.0)*K1(i) + (2.0)*K2(i) - (3544.0/2565.0)*K3(i) + (1859.0/4104.0)*K4(i) - (11.0/40.0)*K5(i);

	     if ( ncontrols>0 ) linear_interpolation(controlp, timep.value(), time_vector, control_trajectory, length(time_vector));

	     for(i=0;i<ncontrols;i++)  controls[i] = controlp(i);

	     dae( derivatives, path, states, controls, param, timep, xad, iphase, workspace);

	     for(i=0;i<nstates;i++) {
	        K6(i) = h*derivatives[i].value();
	     }

	     // COMPUTE R

	     R = (1.0/h)*(  (1.0/360.0)*K1 - (128.0/4275.0)*K3 - (2197.0/75240.0)*K4 + (1.0/50.0)*K5 + (2.0/55.0)*K6  ).cwiseAbs();

	     for (i=0;i<nstates;i++) {
	        double scale =  fabs(states[i].value());
		if ( scale < 1.e-6  ) scale = 1.0;
		R(i)=R(i)/ scale;
	     }

	     maxR = Max(R);

	     if ( maxR <= tolerance ) { // Approximation accepted
	         time = time + h;

            state_trajectory.col(k+1) = state_trajectory.col(k) + (25.0/216.0)*K1 + (1408.0/2565.0)*K3 + (2197.0/4104.0)*K4 - (1.0/5.0)*K5;
		 new_time_vector(k+1) = time.value();
		 timep = new_time_vector(k+1);
		 if (ncontrols>0) {
		    linear_interpolation(controlp, timep.value(), time_vector, control_trajectory, length(time_vector));

		    new_control_trajectory.col(k+1) = controlp;
		 }
		 k= k+1;

	     }

	     delta = 0.84*pow(tolerance/maxR, 1.0/4.0 );

	     if ( delta <= 0.1 ) {
	          h = 0.1*h;
	     }
	     else {
	       if (  delta >= 4.0 ) {
	           h = 4*h;
	       }
	       else {
		   h = delta*h;
	       }
	     }

	     if ( h > hmax ) h = hmax;

	     if ( time.value() >= tf ) {
	       flag = 0;
	     }
	     else {
	         if ( time.value() + h > tf )
		    h = tf - time.value();
		 else if (h<hmin) {
		   flag = 0;
		   fprintf(stderr,"\nh=%e, hmin=%e", h, hmin);
		   error_message("\n Warning: minimum step size exceeded in rkf_propagate( )");
		 }

	     }

	} // End while loop

	state_trajectory.resize(nstates,k);
	new_control_trajectory.resize(ncontrols,k);
	new_time_vector.resize(1,k);

	delete[] states;
	delete[] controls;
	delete[] path;
	delete[] param;
	delete[] derivatives;
}


