/*********************************************************************************************

TThis file is part of the PSOPT library, a software tool for computational optimal control

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



void  define_nlp_bounds(MatrixXd& xlb, MatrixXd& xub, Prob& problem, Alg& algorithm, Workspace* workspace)
{
  // This function defines the NLP bounds given the information provided by
  // the user.

   int i, k;

   int x_phase_offset = 0;

   xlb = zeros(workspace->nvars,1);
   xub = zeros(workspace->nvars,1);


   for (i=0; i< problem.nphases; i++)
   {
   	int norder    = problem.phase[i].current_number_of_intervals;
   	int ncontrols = problem.phase[i].ncontrols;
    int nparam    = problem.phase[i].nparameters;
   	int nstates   = problem.phase[i].nstates;
   	int offset1   = ncontrols*(norder+1);
    int offset2   = (ncontrols+nstates)*(norder+1);

	MatrixXd& control_scaling = problem.phase[i].scale.controls;
	MatrixXd& state_scaling   = problem.phase[i].scale.states;
    MatrixXd& param_scaling   = problem.phase[i].scale.parameters;
    double time_scaling      = problem.phase[i].scale.time;

	int nvars_phase_i = get_nvars_phase_i(problem,i, workspace);


	for (k=0; k<norder+1; k++) {   // EIGEN_UPDATE: k index shifted by -1.
                if (ncontrols>0) {

	              xlb.block(x_phase_offset+k*ncontrols, 0, ncontrols, 1 ) = elemProduct((problem.phase[i].bounds.lower.controls),control_scaling);
                }

		xlb.block(x_phase_offset+(k)*nstates+offset1,0, nstates, 1) = elemProduct((problem.phase[i].bounds.lower.states),state_scaling);

                if (ncontrols>0) {

		     xub.block(x_phase_offset+(k)*ncontrols,0, ncontrols, 1) = elemProduct((problem.phase[i].bounds.upper.controls),control_scaling);

                }

        xub.block(x_phase_offset+(k)*nstates+offset1,0, nstates, 1)=elemProduct((problem.phase[i].bounds.upper.states),state_scaling);

	}

        offset1 = (nstates+ncontrols)*(norder+1);

        if (nparam>=1) {


             xlb.block(x_phase_offset+offset2,0,nparam,1)= elemProduct((problem.phase[i].bounds.lower.parameters),param_scaling);



             xub.block(x_phase_offset+offset2,0,nparam,1)= elemProduct((problem.phase[i].bounds.upper.parameters),param_scaling);

        }

        offset1 = (nstates+ncontrols)*(norder+1) + nparam;

        if ( need_midpoint_controls(*workspace->algorithm, workspace) ) {
	   for (k=0; k<norder; k++) { // EIGEN_UPDATE: K index shifted by -1
                if (ncontrols>0) {

		          xlb.block(x_phase_offset+offset1+(k)*ncontrols,0,ncontrols,1) = elemProduct((problem.phase[i].bounds.lower.controls),control_scaling);


                  xub.block(x_phase_offset+offset1+(k)*ncontrols,0,ncontrols,1) = elemProduct((problem.phase[i].bounds.upper.controls),control_scaling);

                }
	   }
        }

	xlb(x_phase_offset+nvars_phase_i-2)   = problem.phase[i].bounds.lower.StartTime*time_scaling; //EIGEN_UPDATE
	xub(x_phase_offset+nvars_phase_i-2)   = problem.phase[i].bounds.upper.StartTime*time_scaling; //EIGEN_UPDATE

	xlb(x_phase_offset+nvars_phase_i-1) = problem.phase[i].bounds.lower.EndTime*time_scaling;  //EIGEN_UPDATE
	xub(x_phase_offset+nvars_phase_i-1) = problem.phase[i].bounds.upper.EndTime*time_scaling;  //EIGEN_UPDATE

        x_phase_offset += nvars_phase_i;

   }

}





void get_constraint_bounds(double* g_l, double* g_u, Workspace* workspace)
{

  Prob *problem = workspace->problem;
  Alg  *algorithm= workspace->algorithm;

  int lam_phase_offset=0;

  MatrixXd& constraint_scaling = *workspace->constraint_scaling;

  int offset, k, l, i, j;

  double path_sc, event_sc;

  for(i=0; i<problem->nphases; i++)
  {
  	int nstates = problem->phase[i].nstates;
  	int norder  = problem->phase[i].current_number_of_intervals;
  	int nevents = problem->phase[i].nevents;
  	int npath   = problem->phase[i].npath;
  	MatrixXd& path_scaling    = problem->phase[i].scale.path;
  	MatrixXd& event_scaling   = problem->phase[i].scale.events;

        int ncons_phase_i = get_ncons_phase_i(*problem,i, workspace);
	// lower and upper bounds on g(x)
	for (k=0;k<nstates*(norder+1);k++) {
		g_l[lam_phase_offset+k] = 0.0;
		g_u[lam_phase_offset+k] = 0.0;
	}

	offset = lam_phase_offset+nstates*(norder+1);

	for (k=0;k<nevents;k++) { // EIGEN_UPDATE: index k shifted by -1
		j = offset + k;
		if( algorithm->scaling=="user" )
		        event_sc = event_scaling(k);
 		else
			event_sc = constraint_scaling(j);

		g_l[j] = (problem->phase[i].bounds.lower.events)(k)*event_sc;
		g_u[j] = (problem->phase[i].bounds.upper.events)(k)*event_sc;
	}

	offset = offset + nevents;

	for (k=0; k<npath; k++) // EIGEN_UPDATE: index k shifted by -1
	{
		for (l=0;l<(norder + 1);l++) {   // EIGEN_UPDATE: index l shifted by -1.
		    j = offset + (l)*npath + k;
		    if( algorithm->scaling=="user" )
		       path_sc = path_scaling(k);
		    else
		       path_sc = constraint_scaling(j);

		    g_l[j] = (problem->phase[i].bounds.lower.path)(k)*path_sc;
		    g_u[j] = (problem->phase[i].bounds.upper.path)(k)*path_sc;
		}
	}


        if ( need_midpoint_controls(*workspace->algorithm, workspace) ) {
		offset = offset + npath*(norder+1);

		for (k=0; k<npath; k++)  // EIGEN_UPDATE: index k shifted by -1.
		{
			for (l=0;l<norder;l++) { // EIGEN_UPDATE: index l shifted by -1.
		    		j = offset + (l)*npath + k;
		    		if( algorithm->scaling=="user" )
		       			path_sc = path_scaling(k);
		    		else
		       			path_sc = constraint_scaling(j);

		    		g_l[j] = (problem->phase[i].bounds.lower.path)(k)*path_sc;
		    		g_u[j] = (problem->phase[i].bounds.upper.path)(k)*path_sc;
			}
		}
	}


        lam_phase_offset += ncons_phase_i;

        // Bounds for t0 <= tf constraint

        double diff_t0Min_tfMax= (problem->phase[i].bounds.lower.StartTime-problem->phase[i].bounds.upper.EndTime);
        diff_t0Min_tfMax *= problem->phase[i].scale.time;
        if (algorithm->scaling == "automatic")  {
           diff_t0Min_tfMax *= constraint_scaling(lam_phase_offset-1);
        }

        g_l[ lam_phase_offset - 1 ] = diff_t0Min_tfMax;
        g_u[ lam_phase_offset - 1 ] = 0.0;

  }

  for (k=0;k<problem->nlinkages; k++)
  {
        j = lam_phase_offset + k;
        g_l[j] = problem->bounds.lower.linkage(k); // EIGEN_UPDATE
        g_u[j] = problem->bounds.upper.linkage(k); // EIGEN_UPDATE

  }


  return;
}

