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


void evaluate_differential_error_in_phase(MatrixXd& state_error, int iphase, adouble time, adouble* xad, Workspace* workspace)
{
     //   Computes  the differential error epsilon(t) = (xdot(t)-f(x,u,p,t)) within a phase

     Prob* problem = workspace->problem;
     adouble dot_state_j;
     adouble state_j;
     adouble control_j;
     adouble* states;
     adouble* controls;
     adouble* path;
     adouble* parameters;
     adouble* derivatives;
     int i = iphase-1;
     int nstates   = problem->phase[i].nstates;
     int ncontrols = problem->phase[i].ncontrols;
     int j, iph;

     if ( problem->multi_segment_flag || workspace->auto_linked_flag ) {
	  iph = 1;
     }
     else {
	  iph = iphase;
     }

     states        = workspace->states[i];
     controls      = workspace->controls[i];
     parameters    = workspace->parameters[iph-1];
     path          = workspace->path[i];
     derivatives   = workspace->derivatives[i];

     for (j=0;j<nstates;j++) {
          get_interpolated_state(&state_j, j, iphase, time, xad, workspace);    // EIGEN_UPDATE
          states[j]       = state_j;
     }

     for (j=0;j<ncontrols;j++) {
          get_interpolated_control(&control_j, j, iphase, time, xad, workspace); // EIGEN_UPDATE
          controls[j]       = control_j;
     }

     get_parameters(parameters, xad, iphase, workspace );
     problem->dae(derivatives, path, states, controls, parameters, time, xad, iphase, workspace);

     for (j=0;j<nstates;j++) {
          get_state_derivative(&dot_state_j,j,iphase,time,xad, workspace);  // EIGEN_UPDATE
          state_error(j) = dot_state_j.value() - derivatives[j].value();
     }

}




void evaluate_integral_of_differential_error(MatrixXd& eta, int iphase, adouble t1, adouble t2, adouble* xad, int n, Workspace* workspace)
{
// This function evaluates integral[t1,t2]{ |xdot-f(x,u,p,t)| } dt
// by using composite Simpson intergration with n steps.

   Prob* problem = workspace->problem;
   int nstates   = problem->phase[iphase-1].nstates;
	double h = (t2.value()-t1.value())/n;
	int j;

	MatrixXd R1(nstates,1);
	MatrixXd state_error1(nstates,1);
	MatrixXd state_error2(nstates,1);

   evaluate_differential_error_in_phase( state_error1, iphase, t1, xad, workspace );
   evaluate_differential_error_in_phase( state_error2, iphase, t2, xad, workspace );


   R1 = state_error1.cwiseAbs() + state_error2.cwiseAbs() ;

   int nover2 = (int) n/2;

	for (j=1; j<=nover2-1; j++) {

			evaluate_differential_error_in_phase( state_error1, iphase, t1 +2*j*h, xad, workspace );

         R1+= 2.0*state_error1.cwiseAbs();

	}

	for (j=1; j<=nover2; j++) {

			evaluate_differential_error_in_phase( state_error1, iphase, t1 +(2*j-1)*h, xad, workspace );

         R1 += 4.0*state_error1.cwiseAbs();
	}


	eta = (h/3.0)*R1;
}


void evaluate_integral_of_differential_error_L2(MatrixXd& eta, int iphase, adouble t1, adouble t2, adouble* xad, int n, Workspace* workspace)
{
// This function evaluates the L2 norm SQRT[ integral[t1,t2]{ |xdot-f(x,u,p,t)|^2 } dt ]
// by using composite Simpson intergration with n steps.

     	Prob* problem = workspace->problem;
     	int nstates   = problem->phase[iphase-1].nstates;

	double h = (t2.value()-t1.value())/n;
	int j;

	MatrixXd R1(nstates,1);
	MatrixXd state_error1(nstates,1);
	MatrixXd state_error2(nstates,1);

        evaluate_differential_error_in_phase( state_error1, iphase, t1, xad, workspace );
        evaluate_differential_error_in_phase( state_error2, iphase, t2, xad, workspace );

   R1 =  state_error1.cwiseAbs2() + state_error2.cwiseAbs2() ;

        int nover2 = (int) n/2;

	for (j=1; j<=nover2-1; j++) {

			evaluate_differential_error_in_phase( state_error1, iphase, t1 +2*j*h, xad, workspace );
         R1 += 2.0*( state_error1.cwiseAbs2() );

	}

	for (j=1; j<=nover2; j++) {

			evaluate_differential_error_in_phase( state_error1, iphase, t1 +(2*j-1)*h, xad, workspace );
         R1 += 4.0*( state_error1.cwiseAbs2() );
	}


   eta = ((h/3.0)*R1).cwiseSqrt();
}


void evaluate_matrix_of_integrated_errors_in_phase(MatrixXd& eta, int iphase, adouble* xad, int n, Workspace* workspace)
{
//	This function computes a matrix of integrated absolute differential errors, where element (i,j)
//	corresponds to state i and interval j within the phase.
	int k;
	adouble t1, t2;
	adouble t0, tf;
     	Prob* problem = workspace->problem;
        int norder    = problem->phase[iphase-1].current_number_of_intervals;
        int nnodes    = norder + 1;
     	int nstates   = problem->phase[iphase-1].nstates;
	MatrixXd eta_k(nstates,1);
        get_times(&t0, &tf, xad, iphase, workspace );

	for (k=0;k< nnodes-1;k++){  // EIGEN_UPDATE: k index shifted by -1
		t1 = convert_to_original_time_ad( (workspace->snodes[iphase-1])(k), t0, tf );
		t2 = convert_to_original_time_ad( (workspace->snodes[iphase-1])(k+1), t0, tf );;

      evaluate_integral_of_differential_error(eta_k,iphase,t1,t2,xad,n, workspace);
      eta.col(k) = eta_k;
	}
}




void evaluate_solution(Prob& problem,Alg& algorithm,Sol& solution, Workspace* workspace)
{
// This function evaluates the solution found by calculating the relative local error
// over the intervals of each phase. See Betts (2001), Section 4.7.2, and the PSOPT
// handbook for details on the method used.
//
	int nphases = problem.nphases;
	int iphase;
	int n = algorithm.nsteps_error_integration;
	adouble* xad = workspace->xad;
	MatrixXd eta;
	char msg[100];

	MatrixXd states;
	MatrixXd Xdot;
	MatrixXd w;
	MatrixXd states_i;
	MatrixXd Xdot_i;
	MatrixXd eta_i;
	int i;
	psopt_print(workspace,"\n>>> Evaluating the discretization (ODE) error...\n\n");
        for(iphase=1;iphase<=nphases;iphase++) {
	        MatrixXd& epsilon = solution.relative_errors[iphase-1];
		// store previous relative errors before calculating the new ones
		workspace->old_relative_errors[iphase-1] = epsilon;
		MatrixXd& w = workspace->error_scaling_weights[iphase-1];
		int norder = problem.phase[iphase-1].current_number_of_intervals;
		int nstates = problem.phase[iphase-1].nstates;
		if (workspace->current_mesh_refinement_iteration==1) {
			states = solution.get_states_in_phase(iphase);
			Xdot   = workspace->Xdot[iphase-1];
			w.resize(nstates,1);
			for (i=0; i<nstates;i++ ) { // EIGEN_UPDATE: i index shifted by -1

            states_i = states.row(i);

            Xdot_i   = Xdot.row(i);

            w(i) = max(  states_i.lpNorm<Infinity>() , Xdot_i.lpNorm<Infinity>() ) + 1.0;
			}
		}
		eta.resize(nstates,norder);
    		evaluate_matrix_of_integrated_errors_in_phase(eta,iphase,xad,n, workspace);
		for(i=0;i<norder;i++) { // EIGEN_UPDATE: index i shifted by -1

         eta.col(i) = eta.col(i).cwiseQuotient(w);

         eta_i = eta.col(i);

         epsilon(0,i)   = eta_i.maxCoeff();
		}

	}

	// Now print statistics to file

	double mv;

	psopt_print(workspace,"\n*******************************************************************************");
        sprintf(msg,"\n                 Evaluation of mesh refinement iteration %i                     \n", workspace->current_mesh_refinement_iteration );
	psopt_print(workspace,msg);
	sprintf(msg,"\n_____________________________Statistics per phase______________________________");
	psopt_print(workspace,msg);
	sprintf(msg,"\nPhase\t\tNodes\t\tMax ODE Error\tMin ODE error\tMean ODE Error", workspace->current_mesh_refinement_iteration );
	psopt_print(workspace,msg);

   solution.mesh_stats[ workspace->current_mesh_refinement_iteration-1 ].epsilon_max = 0;
	solution.mesh_stats[ workspace->current_mesh_refinement_iteration-1 ].nnodes = 0;

	if (use_local_collocation(algorithm) && workspace->differential_defects=="trapezoidal") {
             solution.mesh_stats[ workspace->current_mesh_refinement_iteration-1 ].method = "TRP";
        }

        else if (use_local_collocation(algorithm)&& workspace->differential_defects == "Hermite-Simpson") {
             solution.mesh_stats[ workspace->current_mesh_refinement_iteration-1 ].method = "H-S";
        }

	else if (use_global_collocation(algorithm)&&algorithm.collocation_method == "Legendre") {
             solution.mesh_stats[ workspace->current_mesh_refinement_iteration-1 ].method = "LGL";
        }

	else if (use_global_collocation(algorithm)&&algorithm.collocation_method == "Chebyshev") {
             solution.mesh_stats[ workspace->current_mesh_refinement_iteration-1 ].method = "CGL";
        }

	if ( solution.mesh_stats[ workspace->current_mesh_refinement_iteration-1 ].method == "LGL" || solution.mesh_stats[ workspace->current_mesh_refinement_iteration-1 ].method == "CGL" ) {
	    if (algorithm.diff_matrix == "standard")
	        solution.mesh_stats[ workspace->current_mesh_refinement_iteration-1 ].method += "-ST";
	    else if ( algorithm.diff_matrix == "reduced-roundoff" )
	        solution.mesh_stats[ workspace->current_mesh_refinement_iteration-1 ].method += "-RR";
	    else if ( algorithm.diff_matrix == "central-differences")
	        solution.mesh_stats[ workspace->current_mesh_refinement_iteration-1 ].method += "-CD";

	}


	for(iphase=1;iphase<=nphases;iphase++) {

	  	MatrixXd& emax_history = workspace->emax_history[iphase-1];

	   mv = solution.relative_errors[iphase-1].transpose().mean();

		emax_history( workspace->current_mesh_refinement_iteration-1, 0) = (double) solution.nodes[iphase-1].size();

		emax_history( workspace->current_mesh_refinement_iteration-1, 1) = solution.relative_errors[iphase-1].maxCoeff();
		sprintf(msg,"\n%i\t\t%li\t\t%e\t%e\t%e", iphase, length(solution.nodes[iphase-1]), Max(solution.relative_errors[iphase-1]),  Min(solution.relative_errors[iphase-1]), mv );
		psopt_print(workspace,msg);


      if ( emax_history( workspace->current_mesh_refinement_iteration-1, 1)>solution.mesh_stats[ workspace->current_mesh_refinement_iteration-1 ].epsilon_max )
		{

        solution.mesh_stats[ workspace->current_mesh_refinement_iteration-1 ].epsilon_max = emax_history( workspace->current_mesh_refinement_iteration-1, 1);
		}

      solution.mesh_stats[ workspace->current_mesh_refinement_iteration-1 ].nnodes += (int) solution.nodes[iphase-1].size();
	}

	solution.mesh_stats[ workspace->current_mesh_refinement_iteration-1 ].nvars = get_number_nlp_vars(problem, workspace);
	solution.mesh_stats[ workspace->current_mesh_refinement_iteration-1 ].ncons = get_number_nlp_constraints(problem, workspace );

	int jj = workspace->current_mesh_refinement_iteration-1;
	sprintf(msg,"\n\n_____________________________Overall Statistics________________________________");
	psopt_print(workspace,msg);
	sprintf(msg,"\n\nTotal Nodes\tNVARS\t\tNCONS\t\tN Obj Eval\tN Cons Eval");
	psopt_print(workspace,msg);

	sprintf(msg,"\n%i\t\t%i\t\t%i\t\t%i\t\t%i", solution.mesh_stats[jj].nnodes, solution.mesh_stats[jj].nvars,
		solution.mesh_stats[jj].ncons, solution.mesh_stats[jj].n_obj_evals, solution.mesh_stats[jj].n_con_evals
		);
	psopt_print(workspace,msg);

	sprintf(msg,"\n\nN Jac Eval\tN Hes Eval\tN ODE RHS\tMax ODE Error\tNLP CPU (sec)");
	psopt_print(workspace,msg);

	sprintf(msg,"\n%i\t\t%i\t\t%i\t\t%e\t%e",
		solution.mesh_stats[jj].n_jacobian_evals, solution.mesh_stats[jj].n_hessian_evals,
		solution.mesh_stats[jj].n_ode_rhs_evals, solution.mesh_stats[jj].epsilon_max,
		solution.mesh_stats[jj].CPU_time);
	psopt_print(workspace,msg);

        psopt_print(workspace,"\n*******************************************************************************\n\n");

}

