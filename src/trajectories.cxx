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



void get_individual_control_trajectory(adouble *control_traj, int control_index, int iphase, adouble* xad, Workspace* workspace)
{
    int k;
    int i = iphase-1;
    Prob& problem = *workspace->problem;
    int ncontrols = problem.phase[i].ncontrols;
    int norder    = problem.phase[i].current_number_of_intervals;
    int iphase_offset = get_iphase_offset(problem,iphase, workspace);
    MatrixXd& control_scaling = problem.phase[i].scale.controls;

    for(k=0;k<norder+1;k++) { // EIGEN_UPDATE: k index shifted by -1

          control_traj[k] = xad[iphase_offset+(k)*ncontrols+control_index]/control_scaling(control_index);
    }

}

void get_individual_state_trajectory(adouble *state_traj, int state_index, int iphase, adouble* xad, Workspace* workspace)
{
    int k;
    int i = iphase-1;
    Prob& problem = *workspace->problem;
    int ncontrols = problem.phase[i].ncontrols;
    int nstates   = problem.phase[i].nstates;
    int norder    = problem.phase[i].current_number_of_intervals;
    int iphase_offset = get_iphase_offset(problem,iphase, workspace);
    int offset1   = ncontrols*(norder+1);
    MatrixXd& state_scaling = problem.phase[i].scale.states;

    for(k=0;k<norder+1;k++) { // EIGEN_UPDATE: k index shifted by -1

          state_traj[k] = xad[iphase_offset+offset1+(k)*nstates+state_index]/state_scaling(state_index); 
    }

}




void compute_derivatives_trajectory( MatrixXd& Xdot, Prob& problem, Sol& solution,  int iphase, Workspace* workspace )
{

	adouble *derivatives;
    	adouble *controls;
    	adouble *parameters;
    	adouble *states;
        adouble *path;
	adouble time;

        int iph;

        int  i = iphase-1;

	if ( problem.multi_segment_flag || workspace->auto_linked_flag ) {
	  iph = 1;
	}
	else {
	  iph = iphase;
	}

	states        = workspace->states[i];
	derivatives   = workspace->derivatives[i];
        controls      = workspace->controls[i];
        parameters    = workspace->parameters[iph-1];
        path          = workspace->path[i];

	int j, k, l;

	int norder    = problem.phase[i].current_number_of_intervals;

	int nstates   = problem.phase[i].nstates;

        int ncontrols = problem.phase[i].ncontrols;

        int nparam    = problem.phase[i].nparameters;


        for(l=0;l<nparam;  l++) parameters[l]  = (solution.parameters[i])(l); // EIGEN_UPDATE: Index l shifted by -1

	for(k=0; k<norder+1; k++)  // EIGEN_UPDATE: Index k shifted by -1
        {

            for(l=0;l<ncontrols;l++) controls[l] = (solution.controls[i])(l,k); // EIGEN_UPDATE: Index l shifted by -1

            for(l=0;l<nstates;  l++) states[l]   = (solution.states[i])(l,k);   // EIGEN_UPDATE: Index l shifted by -1


            time = (solution.nodes[i])(k);

             problem.dae(derivatives, path,  states, controls, parameters, time, solution.xad, iphase,workspace);


            for(j=0; j< nstates; j++) {  // EIGEN_UPDATE: j index shifted by -1
                 Xdot(j,k) = derivatives[j].value();
            }

        }

}



