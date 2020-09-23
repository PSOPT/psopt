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

using namespace Eigen;


void  define_initial_nlp_guess(MatrixXd& x0, MatrixXd& lambda, Sol& solution, Prob& problem, Alg& algorithm, Workspace* workspace)
{
   // This function determines the initial guess for the NLP decision vector x0

   int i;
   int x_phase_offset =0;

   for(i=0; i<problem.nphases; i++)
   {

	MatrixXd& control_scaling = problem.phase[i].scale.controls;
	MatrixXd& state_scaling   = problem.phase[i].scale.states;
        MatrixXd& param_scaling   = problem.phase[i].scale.parameters;
        double&   time_scaling   = problem.phase[i].scale.time;

	int norder    = problem.phase[i].current_number_of_intervals;
	int ncontrols = problem.phase[i].ncontrols;
	int nstates   = problem.phase[i].nstates;
        int nparam    = problem.phase[i].nparameters;
	int offset1   = ncontrols*(norder+1);
        int offset2   = (ncontrols+nstates)*(norder+1);
	int k;
	double t00, tf0;

	MatrixXd umean, xmean;
	MatrixXd Xdot(nstates,norder+1);
	MatrixXd up, xp, un, xn;
	MatrixXd time_guess;


	int nvars_phase_i = get_nvars_phase_i(problem,i, workspace);

	if ( !isEmpty(problem.phase[i].guess.time) ) {
	  t00 = problem.phase[i].guess.time(0);    // EIGEN_UPDATE
	  tf0 = problem.phase[i].guess.time( length(problem.phase[i].guess.time) -1 ); //EIGEN_UPDATE
	}
	else {
		t00 = (problem.phase[i].bounds.lower.StartTime + problem.phase[i].bounds.upper.StartTime)/2.0;
		tf0 = (problem.phase[i].bounds.lower.EndTime + problem.phase[i].bounds.upper.EndTime)/2.0;
	}

        if (nparam > 0) {
               if( !isEmpty(problem.phase[i].guess.parameters) )  {
                  solution.parameters[i] = problem.phase[i].guess.parameters;
               }
               else {
                  solution.parameters[i] = 0.5*((problem.phase[i].bounds.upper.parameters) + (problem.phase[i].bounds.lower.parameters));
               }
        }


	if ( !isEmpty(problem.phase[i].guess.time)) {
		time_guess = problem.phase[i].guess.time;
	} else {
		time_guess = linspace(t00,tf0,norder+1);
	}

	for (k=0; k<norder+1; k++) { //EIGEN_UPDATE
		(solution.nodes[i])(0,k)  =  convert_to_original_time( (workspace->snodes[i])(k), t00, tf0 );
	}


	if ( !isEmpty(problem.phase[i].guess.controls ) && ncontrols>0 ) {

		for (k=0;k<ncontrols;k++) { // EIGEN_UPDATE

         up = problem.phase[i].guess.controls.row(k);

			linear_interpolation(un,solution.nodes[i],time_guess, up, length(up) );

         (solution.controls[i]).row(k)  = un;
		}


	}
	else {
		if ( ncontrols>0 )   {
		   umean = ((problem.phase[i].bounds.lower.controls)+(problem.phase[i].bounds.upper.controls))/2.0;
	      for(k=0;k<norder+1;k++) // EIGEN_UPDATE
	      {
               (solution.controls[i]).col(k) = umean;
	      } 
	   } 
	}


	if ( !isEmpty(problem.phase[i].guess.states) ) {
		for (k=0;k<nstates;k++) { // EIGEN_UPDATE

         xp = problem.phase[i].guess.states.row(k);

            linear_interpolation(xn,solution.nodes[i],time_guess, xp, length(xp));

         (solution.states[i]).row(k)  = xn;
		}

	}

	// Set initial states at mean of feasible region for the moment...

	if ( isEmpty(problem.phase[i].guess.states) )
	{
		xmean = ((problem.phase[i].bounds.lower.states)+(problem.phase[i].bounds.upper.states))/2.0;

      (solution.states[i]).col(0) = xmean;
	}


	// Take the state trajectory as constant and equal to the initial values found above.

	if ( isEmpty(problem.phase[i].guess.states) ) {
		for (k=0;k<norder;k++) // EIGEN_UPDATE
		{

         (solution.states[i]).col(k+1)=(solution.states[i]).col(0);
		}
	}

	// Determine the scaling factors to be used

	determine_scaling_factors_for_variables(solution,problem,algorithm);

	// Now assign the (scaled) variables to the initial sqp decision vector:

	x0(x_phase_offset+nvars_phase_i-2)   = t00*time_scaling;  //EIGEN_UPDATE
	x0(x_phase_offset+nvars_phase_i-1)   = tf0*time_scaling;  // EIGEN_UPDATE


	for (k=0; k<norder+1; k++) { // EIGEN_UPDATE
           if (ncontrols>0) {

             x0.block(x_phase_offset+(k)*ncontrols, 0, ncontrols,1) = elemProduct((solution.controls[i]).col(k), control_scaling);
           }

           x0.block(x_phase_offset+(k)*nstates+offset1,0, nstates,1)=elemProduct((solution.states[i]).col(k), state_scaling);
	}

        if (nparam >= 1) {

             x0.block(x_phase_offset+offset2, 0, nparam , 1)=elemProduct(solution.parameters[i],param_scaling);
        }

       if (need_midpoint_controls(*workspace->algorithm, workspace)) {
	  for (k=0; k<norder; k++) {  // EIGEN_UPDATE
            if (ncontrols>0) {

             x0.block( x_phase_offset+offset2+nparam+(k)*ncontrols, 0, ncontrols , 1) = elemProduct((solution.controls[i]).col(k), control_scaling);
            }

	  }
       }

        x_phase_offset += nvars_phase_i;

  }

  determine_objective_scaling(x0,solution,problem,algorithm, workspace);

  determine_constraint_scaling_factors(x0, solution, problem, algorithm, workspace);

  // Assign zeros to the vector of lagrange multipliers:

  lambda.resize(workspace->ncons, 1);
  lambda.setZero();


}



void hot_start_nlp_guess(MatrixXd& x0,MatrixXd& lambda, Sol& solution,Prob& problem,Alg& algorithm, MatrixXd* prev_states, MatrixXd* prev_controls, MatrixXd* prev_costates, MatrixXd* prev_path, MatrixXd* prev_nodes, MatrixXd* prev_param, MatrixXd& prev_t0, MatrixXd& prev_tf, Workspace* workspace )
{

     int i, k;

     int x_phase_offset   =0;
     int lam_phase_offset = 0;


     sprintf(workspace->text,"\nHot starting solution\n");
     psopt_print(workspace,workspace->text);

     x0.resize(workspace->nvars,1);

     x0.setZero();

     lambda.resize(workspace->ncons,1);

     lambda.setZero();


     for(i=0; i<problem.nphases;i++)
     {

	MatrixXd& control_scaling = problem.phase[i].scale.controls;
	MatrixXd& state_scaling   = problem.phase[i].scale.states;
   MatrixXd& param_scaling   = problem.phase[i].scale.parameters;
   double   time_scaling    =  problem.phase[i].scale.time;


	int norder    = problem.phase[i].current_number_of_intervals;
	int ncontrols = problem.phase[i].ncontrols;
	int npath     = problem.phase[i].npath;
	int nstates   = problem.phase[i].nstates;
        int nparam    = problem.phase[i].nparameters;
	int nevents   = problem.phase[i].nevents;
	int offset1   = ncontrols*(norder+1);
        int offset2   = (ncontrols+nstates)*(norder+1);
	int k;
	int offset;
	MatrixXd xn, xp;
	MatrixXd un, up;
	MatrixXd pn, pp;

	int nvars_phase_i = get_nvars_phase_i(problem,i, workspace);

   int ncons_phase_i = get_ncons_phase_i(problem,i, workspace);


	//     prev_nodes.Save("prev_nodes.dat");

	for (k=0; k<norder+1; k++) { // EIGEN_UPDATE
		(solution.nodes[i])(0,k)  =  convert_to_original_time( (workspace->snodes[i])(k), prev_t0(i), prev_tf(i) );
	}


	// Interpolate states into new nodes
	for (k=0;k<nstates;k++) {  // EIGEN_UPDATE


      xp = (prev_states[i]).row(k);
		if (!use_local_collocation(algorithm) ) {
		    lagrange_interpolation(xn,solution.nodes[i],prev_nodes[i], xp);
		}
		else {
		    linear_interpolation(xn,solution.nodes[i],prev_nodes[i], xp, length(xp));
		}

      (solution.states[i]).row(k)  = xn;
	}

	// Interpolate costates into new nodes
	(workspace->dual_costates[i]).resize(nstates,norder+1);
	for (k=0;k<nstates;k++) { // EIGEN_UPDATE

      xp = (prev_costates[i]).row(k);
		if (!use_local_collocation(algorithm)) {
		    lagrange_interpolation(xn,solution.nodes[i],prev_nodes[i], xp);
		}
		else {
		    linear_interpolation(xn,solution.nodes[i],prev_nodes[i], xp, length(xp));
		}

      (workspace->dual_costates[i]).row(k)  = xn;
	}

	// Interpolate controls into new nodes
	for (k=0;k<ncontrols;k++) { // EIGEN_UPDATE

      up = (prev_controls[i]).row(k);
		linear_interpolation(un,solution.nodes[i],prev_nodes[i], up, length(up) );

      (solution.controls[i]).row(k)  = un;
	}

	// Interpolate path constraint multipliers into new nodes
	if (npath) {
		(workspace->dual_path[i]).resize(npath,norder+1);
		for (k=0;k<npath;k++) { // EIGEN_UPDATE

         pp = (prev_path[i]).row(k);
			if (!use_local_collocation(algorithm)) {
			    lagrange_interpolation(pn,solution.nodes[i],prev_nodes[i], pp);
			}
			else {
			    linear_interpolation(pn,solution.nodes[i],prev_nodes[i], pp, length(pp));
			}

         (workspace->dual_path[i]).row(k)  = pn;
		}
	}

	// Now copy relevant variables into the decision vector

	for (k=0; k<norder+1; k++) {  // EIGEN_UPDATE
          if(ncontrols>0) {

        x0.block(x_phase_offset+(k)*ncontrols,0,ncontrols,1) = elemProduct( (solution.controls[i]).col(k), control_scaling);
          }

        x0.block(x_phase_offset+(k)*nstates+offset1,0,nstates,1)=elemProduct((solution.states[i]).col(k), state_scaling);
	}

        if (nparam>0) {

           x0.block(x_phase_offset+offset2, 0, nparam, 1 ) = elemProduct(prev_param[i], param_scaling);
        }

        if ( need_midpoint_controls(*workspace->algorithm, workspace) ) {

	  for (k=0; k<norder; k++) { // EIGEN_UPDATE
             if(ncontrols>0) {

	            x0.block(x_phase_offset+offset2+nparam+(k)*ncontrols,0,ncontrols,1) =   elemProduct( (solution.controls[i]).col(k), control_scaling);
             }
	  }
   }

	x0(x_phase_offset+ nvars_phase_i-2) = prev_t0(i)*time_scaling;
	x0(x_phase_offset+ nvars_phase_i-1) = prev_tf(i)*time_scaling;


	// And finally copy the lagrange multiplier variables into vector lambda


   for (k=0;k<(norder+1);k++) {
      for (int j=0; j<nstates;j++) {
          lambda(lam_phase_offset + k*nstates+j ) = (workspace->dual_costates[i])(j,k);
      }   
   }

	offset = lam_phase_offset+nstates*(norder+1);

	if (nevents>0)
        {

                 lambda.block(offset,0,nevents,1) = (workspace->dual_events[i]);
	       offset = offset + nevents;
        }
	if (npath>0) {


      for (k=0;k<(norder+1);k++) {
         for (int j=0; j<npath;j++) {
         	lambda(offset+k*npath+j) = (workspace->dual_path[i])(j,k);
         }
      }
   }
        x_phase_offset += nvars_phase_i;
        lam_phase_offset += ncons_phase_i;
  }

  // Now deal with the Lagrange multipliers of the linkage constraints

  if (problem.nlinkages)
  {
     for(k=0; k<problem.nlinkages; k++) // EIGEN_UPDATE
     {
          lambda(lam_phase_offset + k) = -(*solution.dual.linkages)(k);
     }

  }

  // Recalculate the scaling factors for objective and constraints

  determine_objective_scaling(x0,solution,problem,algorithm, workspace);

  determine_constraint_scaling_factors(x0, solution, problem, algorithm, workspace);

}


