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

     states        = workspace->states[i].get();
     controls      = workspace->controls[i].get();
     parameters    = workspace->parameters[iph-1].get();
     path          = workspace->path[i].get();
     derivatives   = workspace->derivatives[i].get();

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




// ---------------------------------------------------------------------------------------
// The state representation of a local collocation method, used by the error estimator
//
// On the interval [t_k, t_{k+1}] the polynomial that trapezoidal and Hermite-Simpson both
// define is the cubic Hermite through (x_k, f_k) and (x_{k+1}, f_{k+1}). Both methods satisfy
// xdot(t_k) = f_k exactly at every node, so the residual of that polynomial measures what
// happens strictly BETWEEN the nodes, which is what the discretization error is.
//
// get_interpolated_state fits a global cubic spline through the node states instead. A spline
// does not interpolate the derivatives, so it carries residual at the nodes themselves, where
// the collocation is satisfied to solver tolerance, and it rings across any kink in the
// solution. Measured against the local error obtained by integrating each interval accurately
// with the transcription's own control, on the terrain example at 160, 320 and 640 uniform
// nodes the spline over-reports by median factors of 7.2, 8.5 and 15.1 -- a ratio that grows
// as the mesh is refined, so the estimate converges more slowly than the quantity it is
// estimating -- while the Hermite form holds 8.4, 5.0 and 4.6 and correlates far better with
// the interval-by-interval distribution (0.99, 0.99, 0.86 against 0.86, 0.95, 0.78).
//
// The polynomial is formed here, in the estimator, rather than in get_interpolated_state: that
// routine belongs to the user interface and is called from inside dae(), so it cannot evaluate
// f without recursion. Forming it here also disposes of get_state_derivative's forward
// difference with a step of sqrt(eps), which is first order and carries a cancellation error of
// order eps*|x|/sqrt(eps); the Hermite derivative is available in closed form.
// ---------------------------------------------------------------------------------------

static void local_node_state_and_derivative(double* xnode, double* fnode, int iphase, int k,
                                            adouble* xad, Workspace* workspace)
{
     Prob* problem = workspace->problem;
     int i         = iphase-1;
     int nstates   = problem->phase[i].nstates;
     int iph = ( problem->multi_segment_flag || workspace->auto_linked_flag ) ? 1 : iphase;

     adouble* states      = workspace->states[i].get();
     adouble* controls    = workspace->controls[i].get();
     adouble* parameters  = workspace->parameters[iph-1].get();
     adouble* path        = workspace->path[i].get();
     adouble* derivatives = workspace->derivatives[i].get();

     adouble t0, tf;
     get_times(&t0, &tf, xad, iphase, workspace);
     adouble time = convert_to_original_time_ad( (workspace->snodes[i])(k), t0, tf );

     get_states(   states,   xad, iphase, k, workspace );
     get_controls( controls, xad, iphase, k, workspace );
     get_parameters( parameters, xad, iphase, workspace );

     problem->dae(derivatives, path, states, controls, parameters, time, xad, iphase, workspace);

     for (int j=0;j<nstates;j++) {
          xnode[j] = states[j].value();
          fnode[j] = derivatives[j].value();
     }
}

static void evaluate_differential_error_hermite(MatrixXd& state_error, int iphase, adouble time,
                                                const double* x0, const double* f0,
                                                const double* x1, const double* f1,
                                                double tk, double hk,
                                                adouble* xad, Workspace* workspace)
{
     Prob* problem = workspace->problem;
     int i         = iphase-1;
     int nstates   = problem->phase[i].nstates;
     int ncontrols = problem->phase[i].ncontrols;
     int iph = ( problem->multi_segment_flag || workspace->auto_linked_flag ) ? 1 : iphase;

     adouble* states      = workspace->states[i].get();
     adouble* controls    = workspace->controls[i].get();
     adouble* parameters  = workspace->parameters[iph-1].get();
     adouble* path        = workspace->path[i].get();
     adouble* derivatives = workspace->derivatives[i].get();

     double s  = (time.value() - tk)/hk;
     double s2 = s*s, s3 = s2*s;

     double b00 =  2.0*s3 - 3.0*s2 + 1.0,  b10 =      s3 - 2.0*s2 + s;
     double b01 = -2.0*s3 + 3.0*s2,        b11 =      s3 -     s2;
     double d00 =  6.0*s2 - 6.0*s,         d10 =  3.0*s2 - 4.0*s + 1.0;
     double d01 = -6.0*s2 + 6.0*s,         d11 =  3.0*s2 - 2.0*s;

     for (int j=0;j<nstates;j++)
          states[j] = b00*x0[j] + hk*b10*f0[j] + b01*x1[j] + hk*b11*f1[j];

     adouble control_j;
     for (int j=0;j<ncontrols;j++) {
          get_interpolated_control(&control_j, j, iphase, time, xad, workspace);
          controls[j] = control_j;
     }

     get_parameters( parameters, xad, iphase, workspace );

     problem->dae(derivatives, path, states, controls, parameters, time, xad, iphase, workspace);

     for (int j=0;j<nstates;j++) {
          double xdot = ( d00*x0[j] + hk*d10*f0[j] + d01*x1[j] + hk*d11*f1[j] )/hk;
          state_error(j) = xdot - derivatives[j].value();
     }
}


static void evaluate_integral_of_differential_error_hermite(MatrixXd& eta, int iphase,
                 const double* x0, const double* f0, const double* x1, const double* f1,
                 double tk, double hk, adouble* xad, int n, Workspace* workspace)
{
// The same composite Simpson quadrature of |xdot-f| as below, taken over one interval of a
// local collocation mesh, with the state read off the interval's own cubic Hermite.

     Prob* problem = workspace->problem;
     int nstates   = problem->phase[iphase-1].nstates;
     double h      = hk/n;
     int j;

     MatrixXd R1(nstates,1);
     MatrixXd e1(nstates,1);
     MatrixXd e2(nstates,1);

     evaluate_differential_error_hermite( e1, iphase, (adouble) tk,      x0,f0,x1,f1, tk,hk, xad, workspace );
     evaluate_differential_error_hermite( e2, iphase, (adouble)(tk+hk),  x0,f0,x1,f1, tk,hk, xad, workspace );

     R1 = e1.cwiseAbs() + e2.cwiseAbs();

     int nover2 = (int) n/2;

     for (j=1; j<=nover2-1; j++) {
          evaluate_differential_error_hermite( e1, iphase, (adouble)(tk+2*j*h), x0,f0,x1,f1, tk,hk, xad, workspace );
          R1 += 2.0*e1.cwiseAbs();
     }

     for (j=1; j<=nover2; j++) {
          evaluate_differential_error_hermite( e1, iphase, (adouble)(tk+(2*j-1)*h), x0,f0,x1,f1, tk,hk, xad, workspace );
          R1 += 4.0*e1.cwiseAbs();
     }

     eta = (h/3.0)*R1;
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

   bool local = use_local_collocation(*workspace->algorithm)
                && ( workspace->differential_defects == "trapezoidal"
                     || workspace->differential_defects == "Hermite-Simpson" );

   vector<double> x0(nstates), f0(nstates), x1(nstates), f1(nstates);

   if ( local && nnodes > 1 )
        local_node_state_and_derivative( x0.data(), f0.data(), iphase, 0, xad, workspace );

	for (k=0;k< nnodes-1;k++){  // EIGEN_UPDATE: k index shifted by -1
		t1 = convert_to_original_time_ad( (workspace->snodes[iphase-1])(k), t0, tf );
		t2 = convert_to_original_time_ad( (workspace->snodes[iphase-1])(k+1), t0, tf );;

      if ( local ) {
           // Sweep the nodes once: the right-hand end of one interval is the left-hand end of
           // the next, so each node's derivative is evaluated exactly once over the phase.
           local_node_state_and_derivative( x1.data(), f1.data(), iphase, k+1, xad, workspace );

           double tk = t1.value();
           double hk = t2.value() - tk;

           if ( hk > 0.0 ) {
                evaluate_integral_of_differential_error_hermite( eta_k, iphase,
                      x0.data(), f0.data(), x1.data(), f1.data(), tk, hk, xad, n, workspace );
           }
           else {
                eta_k.setZero();     // defensive: a degenerate interval contributes nothing
           }

           x0 = x1;  f0 = f1;
      }
      else {
           evaluate_integral_of_differential_error(eta_k,iphase,t1,t2,xad,n, workspace);
      }
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
	adouble* xad = workspace->xad.get();
	MatrixXd eta;

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
        snprintf(workspace->text,sizeof(workspace->text),"\n                 Evaluation of mesh refinement iteration %i                     \n", workspace->current_mesh_refinement_iteration );
	psopt_print(workspace,workspace->text);
	snprintf(workspace->text,sizeof(workspace->text),"\n_____________________________Statistics per phase______________________________");
	psopt_print(workspace,workspace->text);
//	snprintf(workspace->text,sizeof(workspace->text),,"\nPhase\t\tNodes\t\tMax ODE Error\tMin ODE error\tMean ODE Error", workspace->current_mesh_refinement_iteration );
	snprintf(workspace->text,sizeof(workspace->text),"\nPhase\t\tNodes\t\tMax ODE Error\tMin ODE error\tMean ODE Error");
	psopt_print(workspace,workspace->text);

   solution.mesh_stats[ workspace->current_mesh_refinement_iteration-1 ].epsilon_max = 0;
	solution.mesh_stats[ workspace->current_mesh_refinement_iteration-1 ].nnodes = 0;

	// Discretisation-method abbreviation for the mesh-statistics reports (the text files and
	// the LaTeX table), recorded once per mesh-refinement iteration for every refinement mode
	// (hp-adaptive pseudospectral and local Betts alike). The integrated-residual transcription
	// is checked first, as it can run on top of any collocation method. The final branch is
	// defensive: validate_user_input rejects unknown collocation methods, so it should not be
	// reached, but it labels the row with the raw method name rather than leaving it blank,
	// guarding against this site drifting behind a future method.
	{
	    auto& ms = solution.mesh_stats[ workspace->current_mesh_refinement_iteration-1 ];
	    if (algorithm.transcription_method == "integrated-residual")
	        ms.method = "IR";
	    else if (use_local_collocation(algorithm) && workspace->differential_defects == "trapezoidal")
	        ms.method = "TRP";
	    else if (use_local_collocation(algorithm) && workspace->differential_defects == "Hermite-Simpson")
	        ms.method = "H-S";
	    else if (use_global_collocation(algorithm) && algorithm.collocation_method == "Legendre")
	        ms.method = "LGL";
	    else if (use_global_collocation(algorithm) && algorithm.collocation_method == "Chebyshev")
	        ms.method = "CGL";
	    else if (use_global_collocation(algorithm) && algorithm.collocation_method == "Radau")
	        ms.method = "LGR";
	    else if (use_global_collocation(algorithm) && algorithm.collocation_method == "Gauss")
	        ms.method = "LG";
	    else
	        ms.method = algorithm.collocation_method;   // defensive: never blank

	    if (ms.method == "LGL" || ms.method == "CGL") {
	        if (algorithm.diff_matrix == "standard")
	            ms.method += "-ST";
	        else if (algorithm.diff_matrix == "reduced-roundoff")
	            ms.method += "-RR";
	    }
	}


	for(iphase=1;iphase<=nphases;iphase++) {

	  	MatrixXd& emax_history = workspace->emax_history[iphase-1];

	   mv = solution.relative_errors[iphase-1].transpose().mean();

		emax_history( workspace->current_mesh_refinement_iteration-1, 0) = (double) solution.nodes[iphase-1].size();

		emax_history( workspace->current_mesh_refinement_iteration-1, 1) = solution.relative_errors[iphase-1].maxCoeff();
		snprintf(workspace->text,sizeof(workspace->text),"\n%i\t\t%li\t\t%e\t%e\t%e", iphase, length(solution.nodes[iphase-1]), Max(solution.relative_errors[iphase-1]),  Min(solution.relative_errors[iphase-1]), mv );
		psopt_print(workspace,workspace->text);


      if ( emax_history( workspace->current_mesh_refinement_iteration-1, 1)>solution.mesh_stats[ workspace->current_mesh_refinement_iteration-1 ].epsilon_max )
		{

        solution.mesh_stats[ workspace->current_mesh_refinement_iteration-1 ].epsilon_max = emax_history( workspace->current_mesh_refinement_iteration-1, 1);
		}

      solution.mesh_stats[ workspace->current_mesh_refinement_iteration-1 ].nnodes += (int) solution.nodes[iphase-1].size();
	}

	solution.mesh_stats[ workspace->current_mesh_refinement_iteration-1 ].nvars = get_number_nlp_vars(problem, workspace);
	solution.mesh_stats[ workspace->current_mesh_refinement_iteration-1 ].ncons = get_number_nlp_constraints(problem, workspace );

	int jj = workspace->current_mesh_refinement_iteration-1;
	snprintf(workspace->text,sizeof(workspace->text),"\n\n_____________________________Overall Statistics________________________________");
	psopt_print(workspace,workspace->text);
	snprintf(workspace->text,sizeof(workspace->text),"\n\nTotal Nodes\tNVARS\t\tNCONS\t\tN Obj Eval\tN Cons Eval");
	psopt_print(workspace,workspace->text);

	snprintf(workspace->text,sizeof(workspace->text),"\n%i\t\t%i\t\t%i\t\t%i\t\t%i", solution.mesh_stats[jj].nnodes, solution.mesh_stats[jj].nvars,
		solution.mesh_stats[jj].ncons, solution.mesh_stats[jj].n_obj_evals, solution.mesh_stats[jj].n_con_evals
		);
	psopt_print(workspace,workspace->text);

	snprintf(workspace->text,sizeof(workspace->text),"\n\nN Jac Eval\tN Hes Eval\tN ODE RHS\tMax ODE Error\tNLP CPU (sec)");
	psopt_print(workspace,workspace->text);

	snprintf(workspace->text,sizeof(workspace->text),"\n%i\t\t%i\t\t%i\t\t%e\t%e",
		solution.mesh_stats[jj].n_jacobian_evals, solution.mesh_stats[jj].n_hessian_evals,
		solution.mesh_stats[jj].n_ode_rhs_evals, solution.mesh_stats[jj].epsilon_max,
		solution.mesh_stats[jj].CPU_time);
	psopt_print(workspace,workspace->text);

        psopt_print(workspace,"\n*******************************************************************************\n\n");

}

