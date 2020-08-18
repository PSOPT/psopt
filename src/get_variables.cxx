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


void get_controls(adouble* controls, adouble* xad, int iphase, int k, Workspace* workspace)
{
        int i = iphase-1;
        Prob& problem = *workspace->problem;
	     MatrixXd& control_scaling = problem.phase[i].scale.controls;

	     int j;

        int iphase_offset= get_iphase_offset(problem,iphase, workspace);

        // get controls

       int ncontrols = problem.phase[i].ncontrols;



        for(j=0;j<ncontrols;j++) {

           controls[j] =  xad[iphase_offset+(k)*ncontrols+j]/control_scaling(j);
        }

}

void get_controls_bar(adouble* controls_bar, adouble* xad, int iphase, int k, Workspace* workspace)
{
   int i = iphase-1;
   Prob& problem = *workspace->problem;
	MatrixXd& control_scaling = problem.phase[i].scale.controls;

	int j;

   int iphase_offset= get_iphase_offset(problem,iphase, workspace);

	int norder    = problem.phase[i].current_number_of_intervals;
	int ncontrols = problem.phase[i].ncontrols;
	int nstates   = problem.phase[i].nstates;
        int nparam    = problem.phase[i].nparameters;

        int offset = (nstates+ncontrols)*(norder+1)+nparam;

        for(j=0;j<ncontrols;j++) {
             controls_bar[j] =  xad[iphase_offset+offset+(k)*ncontrols+j]/control_scaling(j);
        }
}


void get_final_controls(adouble* controls, adouble* xad, int iphase, Workspace* workspace)
{
        int i = iphase-1;
        Prob& problem = *workspace->problem;
        int k = problem.phase[i].current_number_of_intervals;
        get_controls(controls, xad, iphase, k, workspace);
}

void get_initial_controls(adouble* controls, adouble* xad, int iphase, Workspace* workspace)
{
        get_controls(controls, xad, iphase, 0, workspace);
}

void get_states(adouble* states, adouble* xad, int iphase, int k, Workspace* workspace)
{
        int i = iphase-1;
        Prob& problem            = *workspace->problem;
	MatrixXd& state_scaling   = problem.phase[i].scale.states;


	int j;

        int iphase_offset= get_iphase_offset(problem, iphase, workspace);


        int nstates = problem.phase[i].nstates;
        int ncontrols=problem.phase[i].ncontrols;
        int norder   =problem.phase[i].current_number_of_intervals;
	int offset1   = ncontrols*(norder+1);
        // get states
        for(j=0;j<nstates;j++) {
           states[j] =  xad[iphase_offset+offset1+(k)*nstates+j]/state_scaling(j);
        }

}

void get_final_states(adouble* states, adouble* xad, int iphase, Workspace* workspace)
{
        Prob& problem = *workspace->problem;
        int k = problem.phase[iphase-1].current_number_of_intervals;
        get_states(states, xad, iphase, k, workspace);
}

void get_initial_states(adouble* states, adouble* xad, int iphase, Workspace* workspace)
{
        get_states(states, xad, iphase, 0, workspace);
}



void get_parameters(adouble* parameters, adouble* xad, int iphase, Workspace* workspace)
{
        Prob& problem = *workspace->problem;

	int iph;

	if ( problem.multi_segment_flag || workspace->auto_linked_flag ) {
	  iph = 1;
	}
	else {
	  iph = iphase;
	}


        int i = iph-1;
        MatrixXd& param_scaling   = problem.phase[i].scale.parameters;


	int j;

	int norder    = problem.phase[i].current_number_of_intervals;
	int ncontrols = problem.phase[i].ncontrols;
	int nstates   = problem.phase[i].nstates;
        int nparam    = problem.phase[i].nparameters;
        int offset2   = (ncontrols+nstates)*(norder+1);



        int iphase_offset = get_iphase_offset(problem,iph, workspace);

        // get parameters
        for(j=0;j<nparam;j++) {
             parameters[j] =  xad[iphase_offset+offset2+j]/param_scaling(j);
        }

}

void get_times(adouble *t0, adouble *tf, adouble* xad, int iphase, Workspace* workspace)
{
        int i = iphase-1;
        Prob& problem = *workspace->problem;
        double   time_scaling    =  problem.phase[i].scale.time;


	     int nvars_phase_i = get_nvars_phase_i(problem,i, workspace);

        int iphase_offset = get_iphase_offset(problem, iphase, workspace);

	     *t0  = xad[iphase_offset + nvars_phase_i-2]/time_scaling;
	     *tf  = xad[iphase_offset + nvars_phase_i-1]/time_scaling;

}

adouble get_initial_time(adouble* xad, int iphase, Workspace* workspace)
{
        int i = iphase-1;
        Prob& problem = *workspace->problem;
        double   time_scaling    =  problem.phase[i].scale.time;
        adouble t0;



	     int nvars_phase_i = get_nvars_phase_i(problem,i, workspace);

        int iphase_offset = get_iphase_offset(problem,iphase, workspace);

        t0  = xad[iphase_offset + nvars_phase_i-2]/time_scaling;

        return (t0);
}

adouble get_final_time(adouble* xad, int iphase, Workspace* workspace)
{
        int i = iphase-1;
        Prob& problem = *workspace->problem;
        double   time_scaling    =  problem.phase[i].scale.time;
        adouble tf;



	     int nvars_phase_i = get_nvars_phase_i(problem,i, workspace);

        int iphase_offset = get_iphase_offset(problem,iphase, workspace);

		  tf  = xad[iphase_offset + nvars_phase_i-1]/time_scaling;

        return (tf);
}

void get_scaled_decision_variables_and_bounds(MatrixXd& x, MatrixXd& xlb, MatrixXd& xub, Workspace* workspace)
{

    int i;

    Prob & problem = *(workspace->problem);

    adouble* xad = workspace->xad;


    int nvar = get_number_nlp_vars(problem, workspace);

    for(i=0;i<nvar;i++){  // EIGEN_UPDATE: index i shifted by -1
      x(i) = xad[i].value();
      xlb(i)= (*workspace->xlb)(i);
      xub(i)= (*workspace->xub)(i);
    }

}


