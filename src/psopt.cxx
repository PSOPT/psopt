//
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




int psopt(Sol& solution, Prob& problem, Alg& algorithm)
{
	
	 unique_ptr<Workspace> workspace_up{ new Workspace{problem, algorithm,solution} }; 


    initialize_solution(solution,problem,algorithm, workspace_up.get() );


    try {
           psopt_main(solution, problem, algorithm, workspace_up );
    }
    catch (ErrorHandler handler)
    {
           solution.error_msg = handler.error_message;
           solution.error_flag = 1;
    }
    
    return solution.error_flag;
}


void psopt_main(Sol& solution, Prob& problem, Alg& algorithm,  unique_ptr<Workspace>& workspace_up)
{
// PSOPT:  main algorithm

int MAX_STANDARD_PS_NODES = 200;


Workspace* workspace = workspace_up.get();

string startup_message= "\n *******************************************************************************\n * This is PSOPT, an optimal control solver based on pseudospectral and local  *\n * collocation methods, together with large scale nonlinear programming        *";

sprintf(workspace->text, "%s %s %s", "\n *******************************************************************************\n * PSOPT release number: ", PSOPT_RELEASE_STRING, "                                                 *");

string release_message= workspace->text;

string license_notice=  "\n * Copyright (C) 2010-2020  Victor M. Becerra.                                 *\n *                                                                             *\n * This library is free software; you can redistribute it and/or               *\n * modify it under the terms of the GNU Lesser General Public License          *\n * as published by the Free Software Foundation;  version 2.1.                 *\n * This library is distributed in the hope that it will be useful,             *\n * but WITHOUT ANY WARRANTY; without even the implied warranty of              *\n * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU           *\n * Lesser General Public License for more details.                             *\n * You should have received a copy of the GNU Lesser General Public            *\n * License along with this library;                                            *\n * If not please visit http://www.gnu.org/licenses                             *\n *                                                                             *";

string contact_notice=  "\n * The author can be contacted at his email address:    v.m.becerra@ieee.org   *\n *                                                                             *\n *******************************************************************************\n\n";


  PSOPT_extras::tic();

  get_local_time( solution.start_date_and_time );



  workspace->problem   = &problem;

  workspace->algorithm = &algorithm;

  workspace->solution  = &solution;

  int nlp_ncons;
  int nlp_neq;
  int offset;
  int nphases = problem.nphases;

  int number_of_mesh_refinement_iterations = get_number_of_mesh_refinement_iterations(problem,algorithm);


  int i, k;
  int iter_nodes;
  int hotflag = 0;
  int x_phase_offset = 0;
  int lam_phase_offset = 0;



  sprintf(workspace->text,"%s",startup_message.c_str());
  psopt_print(workspace,workspace->text);
  sprintf(workspace->text,"%s",release_message.c_str());
  psopt_print(workspace,workspace->text);
  sprintf(workspace->text,"%s",license_notice.c_str());
  psopt_print(workspace,workspace->text);
  sprintf(workspace->text,"%s",contact_notice.c_str());
  psopt_print(workspace,workspace->text);

  validate_user_input(problem,algorithm, workspace);

  PSOPT_extras::SetPrintLevel( algorithm.print_level );


  if (problem.integrand_cost == NULL )  {
      for(i=1;i<=problem.nphases;i++) {
	  problem.phase[i-1].zero_cost_integrand = true;
      }
  }
  
  for (iter_nodes=1; iter_nodes<= number_of_mesh_refinement_iterations; iter_nodes++) {


    workspace->current_mesh_refinement_iteration = iter_nodes;

    MatrixXd& x0     = *workspace->x0;
    MatrixXd& lambda = *workspace->lambda;
    MatrixXd& xlb    = *workspace->xlb;
    MatrixXd& xub    = *workspace->xub;



    if (algorithm.collocation_method=="trapezoidal") {
             workspace->differential_defects = "trapezoidal";
    }

    else if (algorithm.collocation_method == "Hermite-Simpson") {
             workspace->differential_defects = "Hermite-Simpson";
    }

    else if (use_global_collocation(algorithm)) {
          if (algorithm.diff_matrix=="standard") {
             workspace->differential_defects = "standard";
	  }

	  else if (algorithm.diff_matrix=="reduced-roundoff") {
             workspace->differential_defects = "reduced-roundoff";
	  }

	  else if (algorithm.diff_matrix=="central-differences") {
             workspace->differential_defects = "central-differences";
	  }
    }



    if (algorithm.mesh_refinement == "manual" )
    {
			for (i=0; i<nphases; i++)
			{
	  			problem.phase[i].current_number_of_intervals    = ( (int) problem.phase[i].nodes(iter_nodes-1)) -1;
			}
    }

    else  if (algorithm.mesh_refinement=="automatic" && use_global_collocation(algorithm)  ) {

         	if ( iter_nodes == 1 ) {
	    			for (i=0; i<nphases; i++)
					{
		    			problem.phase[i].current_number_of_intervals    = ( (int) problem.phase[i].nodes(iter_nodes-1)) -1;
					}
	  			}

	  			else if (iter_nodes <= algorithm.mr_min_extrapolation_points) {
	 				for (i=0; i<nphases; i++)
					{
		    			problem.phase[i].current_number_of_intervals    += algorithm.mr_initial_increment;
					}

          	}
    }

    else  if (algorithm.mesh_refinement=="automatic" && use_local_collocation(algorithm) ) {
          // Local mesh refinement algorithm by Betts (2001)

	  // Step 3: Estimate primary order for the new mesh (either trapezoidal or Hermite-Simpson)
         if ( iter_nodes == 1 ) {
	    			for (i=0; i<nphases; i++)
					{
		    			problem.phase[i].current_number_of_intervals    = ( (int) problem.phase[i].nodes(iter_nodes-1)) -1;
					}
					if ( workspace->differential_defects != "trapezoidal" && algorithm.switch_order > 0) {
		   	  		workspace->differential_defects = "trapezoidal";
					}
	  		}	

	 		else if (iter_nodes == 2) {

	      	bool equi_error = check_for_equidistributed_error(problem,algorithm,solution);

	      	if (equi_error && algorithm.switch_order>0) {
				 	  workspace->differential_defects = "Hermite-Simpson";
	      	}
	 		}


	 		else if (workspace->differential_defects == "trapezoidal" && iter_nodes> algorithm.switch_order && algorithm.switch_order>0) {

	           workspace->differential_defects = "Hermite-Simpson";

    		}

	  		// Step 4: Estimate order reduction for each interval.

	 		if (iter_nodes>2) {
	      		estimate_order_reduction(problem,algorithm,solution, workspace);
	 		}

	 		else {
	      		zero_order_reduction(problem,algorithm,solution, workspace);
	 		}		

	 		// Step 5: construct new mesh

	 		if (iter_nodes>=2) {

	      		construct_new_mesh(problem,algorithm,solution, workspace);

	 		}



    }

    workspace->trace_f_done = false;

    workspace->nvars     = get_number_nlp_vars(problem, workspace);

    nlp_ncons           = get_number_nlp_constraints(problem, workspace);

    workspace->ncons = nlp_ncons;

    resize_workspace_vars(problem,algorithm,solution, workspace);

    resize_solution(solution,problem,algorithm);

    // Compute the nodes for each phase

    if (algorithm.collocation_method == "Legendre") {
    	for(i=0; i<nphases; i++)
    	{
		 if (problem.phase[i].current_number_of_intervals > MAX_STANDARD_PS_NODES && algorithm.mesh_refinement=="automatic")
		 {
		    // When using automatic mesh refinement, switch to central differences when the order is too large to avoid numerical problems.
		    workspace->differential_defects="central-differences";
		 }
        	 lglnodes( problem.phase[i].current_number_of_intervals, workspace->snodes[i], workspace->w[i], workspace->P[i], workspace->D[i], workspace);

                sort_vector(workspace->snodes[i],workspace->sindex[i]);

                rearrange_vector(workspace->w[i], workspace->sindex[i] );

    	}
    }

    else if ( algorithm.collocation_method == "Chebyshev" ) {

	    for(i=0; i<nphases; i++)
    	    {
	         if (problem.phase[i].current_number_of_intervals > MAX_STANDARD_PS_NODES && algorithm.mesh_refinement=="automatic")
		 {
		    // When using automatic mesh refinement, switch to central differences when the order is too large to avoid numerical problems.
		    workspace->differential_defects="central-differences";
		 }
         	cglnodes( problem.phase[i].current_number_of_intervals, workspace->snodes[i], workspace->w[i], workspace->D[i], workspace );

                sort_vector(workspace->snodes[i],workspace->sindex[i]);

                rearrange_vector(workspace->w[i], workspace->sindex[i] );
            }

    }

    else if ( ( use_local_collocation(algorithm) && (iter_nodes==1)) || (use_local_collocation(algorithm) && (iter_nodes>1) && (algorithm.mesh_refinement=="manual") )  ) {

	    for(i=0; i<nphases; i++)
    	    {
	        workspace->snodes[i] = linspace(-1.0, 1.0, problem.phase[i].current_number_of_intervals+1);
            }

    }


    // Define initial NLP guess
    if (iter_nodes==1) {
       hotflag = 0;
       define_initial_nlp_guess(x0, lambda, solution, problem, algorithm, workspace);
    }
    else {
        	hotflag = 1;
			hot_start_nlp_guess(x0, lambda, solution,problem,algorithm, workspace->prev_states, workspace->prev_controls, workspace->prev_costates, workspace->prev_path, workspace->prev_nodes, workspace->prev_param, *workspace->prev_t0, *workspace->prev_tf, workspace);
    }

  	// Define NLP bounds on the decision vector

    define_nlp_bounds(*workspace->xlb, *workspace->xub, problem, algorithm, workspace);

    nlp_neq = 0; //   equality constraints

    // Note: the decision variables are:
    // x = [vec(controls)' vec(states)'  parameters' t0 tf ... (this is repeated for each phase)]';
    // where "controls" is a real matrix with dimensions [ncontrols x norder+1]
    //       "states" is a real matrix with dimensions [nstates x norder+1]
    //       "parameters" is a real vector with dimensions [nparameters x 1]
    //       "t0" and "tf" are real numbers
    //       and the vec(.) operator stacks the columns of a matrix one below the next


    sprintf(workspace->text,"\nProblem:\t\t\t\t\t\t%s", problem.name.c_str());
    psopt_print(workspace,workspace->text);


    sprintf(workspace->text, "\nThis is mesh refinement iteration:\t\t\t%i", iter_nodes);
    psopt_print(workspace,workspace->text);
    if ( use_global_collocation(algorithm) ) {
      sprintf(workspace->text, "\nCollocation method:\t\t\t\t\t%s", algorithm.collocation_method.c_str());
      psopt_print(workspace,workspace->text);
    }
    else {
      sprintf(workspace->text, "\nCollocation method:\t\t\t\t\t%s", workspace->differential_defects.c_str());
      psopt_print(workspace,workspace->text);
    }
    if ( use_global_collocation(algorithm) ) {
	  sprintf(workspace->text, "\nDifferentiation matrix:\t\t\t\t\t%s", algorithm.diff_matrix.c_str());
	  psopt_print(workspace,workspace->text);
    }

    sprintf(workspace->text, "\nNumber of NLP variables\t\t\t\t\t%i", workspace->nvars );
    psopt_print(workspace,workspace->text);
    sprintf(workspace->text, "\nNumber of NLP nonlinear constraints:\t\t\t%i", nlp_ncons);
    psopt_print(workspace,workspace->text);
    for(i=1;i<=problem.nphases;i++) {
      sprintf(workspace->text, "\nNumber of nodes phase %i:\t\t\t\t%i", i,problem.phases(i).current_number_of_intervals+1);
      psopt_print(workspace,workspace->text);
    }
    sprintf(workspace->text,"\n");
    psopt_print(workspace,workspace->text);

    workspace->enable_nlp_counters = true;

    chronometer_tic(workspace);

    NLP_interface( algorithm, &x0,  ff_num, gg_num, nlp_ncons,  nlp_neq , &xlb, &xub, &lambda, hotflag, 1, workspace, problem.user_data   );

    solution.mesh_stats[workspace->current_mesh_refinement_iteration-1].CPU_time = chronometer_toc(workspace);

    workspace->enable_nlp_counters = false;

    // Copy the resultant decision vector into the relevant solution variables.

    copy_decision_variables(solution, x0, problem, algorithm, workspace);

    solution.cost = ff_num(x0, workspace)/problem.scale.objective;

    sprintf(workspace->text,"\nReturned (unscaled) cost function value: %e", solution.cost);
    psopt_print(workspace,workspace->text);
    x_phase_offset   = 0;
    lam_phase_offset = 0;

    for(i=0; i<nphases; i++)
    {
        int nstates   = problem.phase[i].nstates;
        int ncontrols = problem.phase[i].ncontrols;
        int nparam    = problem.phase[i].nparameters;
        int norder    = problem.phase[i].current_number_of_intervals;
        int nevents   = problem.phase[i].nevents;
        int npath     = problem.phase[i].npath;
        int iphase = i+1;
        MatrixXd& deriv_scaling = problem.phase[i].scale.defects;
        MatrixXd pz(nstates,1);
        MatrixXd pint;
        MatrixXd tint(1,norder);
        MatrixXd ts;
        MatrixXd pl;
        MatrixXd pextra(1,norder+1);
        double   time_scaling    =  problem.phase[i].scale.time;
        double hk, t0, tf;

        int nvars_phase_i = get_nvars_phase_i(problem, i, workspace);

        int ncons_phase_i =  get_ncons_phase_i(problem,i, workspace);

        solution.dual.costates[i] = lambda.block(lam_phase_offset,0,nstates*(norder+1),1);

		  solution.dual.costates[i] = reshape(solution.dual.costates[i], nstates, norder+1).eval();

		  workspace->prev_costates[i]      = solution.dual.costates[i];

    	  t0 = (solution.nodes[i])(1);

	     tf = (solution.nodes[i])(0, solution.nodes[i].cols()-1); 

	     if ( algorithm.collocation_method=="Legendre" ) {
	    		for(k=0;k<norder+1;k++) { // EIGEN_UPDATE: Index k shifted to start at 0

                   (solution.dual.costates[i]).block(0,k,nstates,1) = (solution.dual.costates[i]).block(0,k,nstates,1)/(workspace->w[i])(k);  // See PhD thesis by Huntington (2006).
	    		}
	     }	

		  if (use_local_collocation(algorithm)) {
        		for(k=0;k<norder;k++) { // EIGEN_UPDATE: Index k shifted to start at 0
	         // For local collocation, this is an estimate of the costate.
	         	hk = (solution.nodes[i])(k+1)-(solution.nodes[i])(k);

	         	(solution.dual.costates[i]).block(0,k,nstates,1) = (solution.dual.costates[i]).block(0,k,nstates,1)/(2.0*hk); // This gives the costates at the midpoints between the collocation nodes.
               (solution.dual.costates[i]).block(0,k,nstates,1)= (solution.dual.costates[i]).block(0,k,nstates,1)*(tf-t0);
 	   		}
		  }		

    	  if ( algorithm.collocation_method == "Chebyshev" ) {
        		for(k=1;k<norder;k++) {  // EIGEN_UPDATE: Index k shifted by -1
                  double tk = (workspace->snodes[i])(k);

                    (solution.dual.costates[i]).block(0,k,nstates,1) = (solution.dual.costates[i]).block(0,k,nstates,1)/(workspace->w[i])(k);

                    (solution.dual.costates[i]).block(0,k,nstates,1) =(1.0/sqrt(1.0 - tk*tk))*(solution.dual.costates[i]).block(0,k,nstates,1); // See PhD thesis by Pietz (2003)


        		}
    	  }

    	  if ( algorithm.scaling=="user" ) {
	     		for (k=0;k<norder+1;k++) {  // EIGEN_UPDATE: Index k shifted by -1
 //                  (solution.dual.costates[i])(colon(),k)=(solution.dual.costates[i])(colon(),k) & deriv_scaling;
                     (solution.dual.costates[i]).block(0,k,nstates,1) =(solution.dual.costates[i]).block(0,k,nstates,1).cwiseProduct(deriv_scaling);
            }
		  }

   	  solution.dual.costates[i] = solution.dual.costates[i]/problem.scale.objective;

   	  if ( algorithm.scaling=="automatic" ) {
	  			solution.dual.costates[i] = reshape(solution.dual.costates[i], nstates*(norder+1), 1 );
          	solution.dual.costates[i] = solution.dual.costates[i].cwiseProduct((*workspace->constraint_scaling).block(lam_phase_offset,0, nstates*(norder+1),1) );
	  			solution.dual.costates[i] = reshape(solution.dual.costates[i], nstates, norder+1);
   	  }

   	  if ( algorithm.collocation_method == "Chebyshev") {
                // use linear extrapolation to approximate costate values at both ends (not perfect but better than nothing)

                MatrixXd Pl, Pl1, Pl2;
                double sl, sl1;

                  Pl1 = solution.dual.costates[i].block(0,norder-3, nstates, 1);
                  Pl2 = solution.dual.costates[i].block(0,norder-4, nstates, 1);
                  sl  = (workspace->snodes[i])(norder-2)-(workspace->snodes[i])(norder -3  );
                  sl1 = (workspace->snodes[i])(norder-3)-(workspace->snodes[i])(  norder-4 );

                Pl = Pl1 + (Pl1-Pl2)/sl1*sl;
               solution.dual.costates[i].block(0, norder-2, nstates, 1) = Pl;


                Pl1 = solution.dual.costates[i].block(0,4, nstates, 1);

                Pl2 = solution.dual.costates[i].block(0,3, nstates, 1);


                sl1  = (workspace->snodes[i])(4)-(workspace->snodes[i])(3);

                sl   = (workspace->snodes[i])(3)-(workspace->snodes[i])(2);

                Pl = Pl2 - (Pl1-Pl2)/sl1*sl;


                solution.dual.costates[i].block(0,2, nstates,1) = Pl;

                /*  */



                  Pl1 = solution.dual.costates[i].block(0,norder-2,nstates,1);

                  Pl2 = solution.dual.costates[i].block(0,norder-3,nstates,1);


                sl  = (workspace->snodes[i])(norder-1)-(workspace->snodes[i])(norder -2  );

                sl1 = (workspace->snodes[i])(norder-2)-(workspace->snodes[i])(  norder-3 );

                Pl = Pl1 + (Pl1-Pl2)/sl1*sl;


                solution.dual.costates[i].block(0,norder-1,nstates,1) = Pl;


                Pl1 = solution.dual.costates[i].block(0,3,nstates,1);

                Pl2 = solution.dual.costates[i].block(0,2,nstates,1);


                sl1  = (workspace->snodes[i])(3)-(workspace->snodes[i])(2);

                sl   = (workspace->snodes[i])(2)-(workspace->snodes[i])(1);

                Pl = Pl2 - (Pl1-Pl2)/sl1*sl;


                solution.dual.costates[i].block(0,1,nstates,1) = Pl;


                /*  */


                Pl1 = solution.dual.costates[i].block(0,norder-1,nstates,1);

                Pl2 = solution.dual.costates[i].block(0,norder-2,nstates,1);


                sl  = (workspace->snodes[i])(norder)-(workspace->snodes[i])(norder -1  );

                sl1 = (workspace->snodes[i])(norder-1)-(workspace->snodes[i])(  norder-2 );

                Pl = Pl1 + (Pl1-Pl2)/sl1*sl;


                solution.dual.costates[i].block(0,norder,nstates,1) = Pl;


                Pl1 = solution.dual.costates[i].block(0,2,nstates,1);


                Pl2 = solution.dual.costates[i].block(0,1,nstates,1);


                sl1  = (workspace->snodes[i])(2)-(workspace->snodes[i])(1);

                sl   = (workspace->snodes[i])(1)-(workspace->snodes[i])(0);

                Pl = Pl2 - (Pl1-Pl2)/sl1*sl;

                solution.dual.costates[i].block(0,0,nstates,1) = Pl;

        }

    	  if ( (algorithm.collocation_method == "Legendre" || algorithm.collocation_method == "Chebyshev") && algorithm.diff_matrix != "central-differences") {
                // Smooth the costates, as they can be noisy in the case of standard Legendre or Chebyshev collocation.
                // (see Farhroo and Ross "Costate estimation by a Legendre Pseudospectral Method", Journal of Guidance
                //    Control and Dynamics, 2001).

                pint.resize(nstates,norder+1);

                pint.block(0,0,nstates,1) = (solution.dual.costates[i].block(0,0,nstates,1) + solution.dual.costates[i].block(0,1,nstates,1))/2.0;

                for (k=1;k< norder;k++) {  // EIGEN_UPDATE: k index shifted by -1

                pint.block(0,k,nstates,1) = (0.25*solution.dual.costates[i].block(0,k-1,nstates,1) + 0.5*solution.dual.costates[i].block(0,k,nstates,1) + 0.25*solution.dual.costates[i].block(0,k+1,nstates,1) )/1.0;
                }


                pint.block(0,norder,nstates,1) = (solution.dual.costates[i].block(0,norder-1,nstates,1) + solution.dual.costates[i].block(0,norder,nstates,1))/2.0;

                solution.dual.costates[i] = pint;

         }


			if ( use_local_collocation(algorithm) ) {
        // use linear extrapolation to approximate the costate values at the collocation nodes...

         	pint = solution.dual.costates[i].block(0,0,nstates,norder);
				for(int k=0;k<norder;k++) {  // EIGEN_UPDATE: k index shifted by -1.
		   		tint(k) = ( (workspace->snodes[i])(k)+(workspace->snodes[i])(k+1) )/2.0;
				}
        		for (int l=0;l<nstates;l++) { // EIGEN_UPDATE: l index shifted by -1.

                   ts = workspace->snodes[i].transpose();

                   tint = tint.transpose().eval();

                   pl = pint.row(l);
                   linear_interpolation(pextra, ts, tint, pl,norder); // EIGEN_UPDATE - Check function

                   solution.dual.costates[i].row(l) = pextra;
       		}		

        // use linear extrapolation to approximate costate values at end point (not perfect but better than nothing)

                pint = solution.dual.costates[i].block(0  ,norder-2  ,nstates  , 2 );

                tint = workspace->snodes[i].block(0  ,norder-2  ,1  , 2 );
        		for (int l=0;l<nstates;l++) {  // EIGEN_UPDATE: l index shifted by -1.

                   double tss = (workspace->snodes[i])(norder);

                   tint = tint.transpose().eval();

                   long ncols = pint.cols();
                   pl = pint.block(l, 0, 1, ncols);
                   linear_interpolation(pextra, tss, tint, pl,2);

                   solution.dual.costates[i](l,norder) = pextra(0);
        		}	

         }


		offset = lam_phase_offset+nstates*(norder+1);
	
		if (   algorithm.nlp_method == "IPOPT"   ) {
	                solution.dual.costates[i] = -solution.dual.costates[i];
	    }
	

		solution.dual.events[i]  = -lambda.block(offset,0,nevents,1);
		workspace->dual_events[i] = -(solution.dual.events[i]);
	
		if (algorithm.scaling=="user") {

			   solution.dual.events[i] = solution.dual.events[i].cwiseProduct( problem.phase[i].scale.events );
		}
	
	    if (algorithm.scaling=="automatic" && nevents>0) {

		   solution.dual.events[i] = solution.dual.events[i].cwiseProduct(  (*workspace->constraint_scaling).block(offset,0,nevents,1) );
	    }
		solution.dual.events[i] /= problem.scale.objective;
	
		if (   algorithm.nlp_method == "IPOPT"   ) {
	                solution.dual.events[i] = -solution.dual.events[i];
	    }

		offset = offset+nevents;
		if (npath>0) {

	             solution.dual.path[i]      = -lambda.block(offset,0, npath*(norder+1), 1);
		     solution.dual.path[i]      = reshape(solution.dual.path[i], npath, norder+1);
		     workspace->prev_path[i] = -(solution.dual.path[i]);
		     if (use_local_collocation(algorithm)) {
	        	for (k=0;k<norder;k++) {  // EIGEN_UPDATE: k index shifted by -1
	        	    // Below are the path constraint adjoint estimates for trapezoidal discretization
	        	    // See Betts (2010), p. 176.
	        	    if ( algorithm.collocation_method == "trapezoidal") {
	                    double hk = (solution.nodes[i])(k+1)-(solution.nodes[i])(k);
	                    if (k==0) {

	                          (solution.dual.path[i]).block(0,k,npath,1) = -(solution.dual.path[i]).block(0,k,npath,1)*(2.0/hk);
	                    }
	                    else if (k==norder-1) {
	                          double hk_1 = (solution.nodes[i])(k)-(solution.nodes[i])(k-1);

	                          (solution.dual.path[i]).block(0,k,npath,1) = -(solution.dual.path[i]).block(0,k,npath,1)*(2.0/hk_1);
	                    }
	                    else {
	                        double hk_1 = (solution.nodes[i])(k)-(solution.nodes[i])(k-1);

	                          (solution.dual.path[i]).block(0,k,npath,1) = -(solution.dual.path[i]).block(0,k,npath,1)*(2.0/(hk+hk_1));
	        	    }
	                    }
	           	    if ( algorithm.collocation_method == "Hermite-Simpson") {
	           	        // These are the path constraint adjoint estimates for Hermite-Simpson discretization
	           	        // See Betts (2010), p. 177.
	           	        double hk = (solution.nodes[i])(k+1)-(solution.nodes[i])(k);
	                    if (k==0) {

	                          (solution.dual.path[i]).block(0,k,npath,1) = -(solution.dual.path[i]).block(0,k,npath,1)*(6.0/hk);
	                    }
	                    else if (k==norder-1) {
	                        double hk_1 = (solution.nodes[i])(k)-(solution.nodes[i])(k-1);

	                        (solution.dual.path[i]).block(0,k,npath,1) = -(solution.dual.path[i]).block(0,k,npath,1)*(6.0/hk_1);
	                    }
	                    else {
	                        double hk_1 = (solution.nodes[i])(k)-(solution.nodes[i])(k-1);

	                        (solution.dual.path[i]).block(0,k,npath,1) = -(solution.dual.path[i]).block(0,k,npath,1)*(6.0/(hk+hk_1));
	                    }
	           	    }
	
	
	        	}
             // use linear extrapolation to approximate multiplier values at end point (not perfect but better than nothing)

          	pint = solution.dual.path[i].block(0,norder-2,npath,2);

          	tint = workspace->snodes[i].block(0,norder-2,1,2);
             for (int l=0;l<npath;l++) {  //EIGEN_UPDATE: l index shifted by -1.

                double tss = (workspace->snodes[i])(norder);

                tint = tint.transpose().eval();

                pl   = pint.block(l,0,1,2);
                linear_interpolation(pextra, tss, tint, pl,2);

                  solution.dual.path[i](l,norder) = pextra(0);

             }

	     }

	     if ( algorithm.collocation_method == "Legendre") {
             for (k=0;k<norder+1;k++) {  // EIGEN_UPDATE: k index shifted by -1

		     		(solution.dual.path[i]).block(0,k,npath,1) = -(solution.dual.path[i]).block(0,k,npath,1)/(workspace->w[i])(k);  // See PhD thesis by Huntington (2006).

		     		(solution.dual.path[i]).block(0,k,npath,1) = (solution.dual.path[i]).block(0,k,npath,1)*(2.0/(tf-t0));
             }
	     }

	     if ( algorithm.collocation_method == "Chebyshev" ) {
             for (k=1; k< norder;k++) { // EIGEN_UPDATE: k index shifted by -1
                  double tk = (workspace->snodes[i])(k);

                    (solution.dual.path[i]).block(0,k,npath,1) = (solution.dual.path[i]).block(0,k,npath,1)/(workspace->w[i])(k);

                    (solution.dual.path[i]).block(0,k,npath,1) = -(1.0/sqrt(1.0 - tk*tk))*(solution.dual.path[i]).block(0,k,npath,1); // See PhD thesis by Pietz (2003)

                    (solution.dual.path[i]).block(0,k,npath,1) = (solution.dual.path[i]).block(0,k,npath,1)*(2.0/(tf-t0));

             }
            // use linear extrapolation to approximate multiplier values at both ends (not perfect but better than nothing)

                pint = solution.dual.path[i].block(0,2,npath,norder-3);


             tint.resize(1,norder-4);
             for (int ii=0; ii<(norder-4); ii++) {
                 tint(ii) = workspace->snodes[i](2+ii);
             }
            for (int l=0;l<npath;l++) {  //EIGEN_UPDATE: l index shifted by -1

                   ts = workspace->snodes[i].transpose();

                   tint = tint.transpose().eval();

                     long ncols = pint.cols();
                     pl = pint.block(l,0,1,ncols);
                   linear_interpolation(pextra, ts, tint, pl,norder-4);

                   for(k=0; k<length(pextra);k++) {
                	   solution.dual.path[i](l,k) = pextra(k);
                   }
            }
	     }

	     if ( algorithm.scaling == "user") {
	         for(k=0;k<norder+1;k++) { //EIGEN_UPDATE: k index shifted by -1

                      (solution.dual.path[i]).block(0,k,npath,1) = (solution.dual.path[i]).block(0,k,npath,1).cwiseProduct(problem.phase[i].scale.path);
	         }
         }

         if (algorithm.scaling=="automatic") {
             for (k=0;k<norder+1;k++) { //EIGEN_UPDATE: k index shifted by -1

                    solution.dual.path[i].block(0,k,npath,1)  =solution.dual.path[i].block(0,k,npath,1).cwiseProduct( (*workspace->constraint_scaling).block(offset+(k-1)*npath,0,npath,1) );
  	         }
         }
         solution.dual.path[i] /= problem.scale.objective;

	}

    offset = offset + npath*(norder+1);


    compute_derivatives_trajectory( workspace->Xdot[i], problem, solution, iphase, workspace );


       solution.dual.Hamiltonian[i] = (solution.integrand_cost[i]);
       MatrixXd Temp1 =  solution.dual.costates[i].cwiseProduct(workspace->Xdot[i]);
       MatrixXd Temp2 = sum_columns(Temp1);
       solution.dual.Hamiltonian[i]+= Temp2;



		workspace->prev_states[i]   =   solution.states[i];
		workspace->prev_controls[i] =   solution.controls[i];
		workspace->prev_nodes[i]    =   solution.nodes[i];
        if (problem.phase[i].nparameters) workspace->prev_param[i]    =   solution.parameters[i];

        (*workspace->prev_t0)(i) = x0(x_phase_offset+nvars_phase_i-2)/problem.phase[i].scale.time;

        (*workspace->prev_tf)(i) = x0(x_phase_offset+nvars_phase_i-1)/problem.phase[i].scale.time;

        x_phase_offset   += nvars_phase_i;
        lam_phase_offset += ncons_phase_i;

	// store scaled nodes
	workspace->old_snodes[i] = workspace->snodes[i];


    }

    if (problem.nlinkages) {

         *solution.dual.linkages = -lambda.block(lam_phase_offset, 0, problem.nlinkages, 1);
         if (algorithm.scaling=="user") {

             *solution.dual.linkages = (*solution.dual.linkages).cwiseProduct(problem.scale.linkages);

         }
         if (algorithm.scaling=="automatic") {

           *solution.dual.linkages =  (*solution.dual.linkages).cwiseProduct(   (*workspace->constraint_scaling).block(offset,0,problem.nlinkages,1)   );
         }
         *solution.dual.linkages /= problem.scale.objective;

         if (algorithm.nlp_method == "IPOPT") {
             *solution.dual.linkages = -(*solution.dual.linkages);
         }
    }

    if (!useAutomaticDifferentiation(algorithm) && algorithm.nlp_method=="IPOPT")  {
//          deleteIndexGroups( workspace->igroup, workspace->nvars );
    }

    evaluate_solution(problem, algorithm, solution, workspace);

    if ( algorithm.mesh_refinement == "automatic" ) {
       // Check satisfaction of mesh refinement tolerance
       int mr_phase_convergence_count = 0;
       for ( i=0; i< problem.nphases; i++ ) {
	    MatrixXd& emax_history = workspace->emax_history[i];

            if ( emax_history( iter_nodes-1, 1 ) <= algorithm.ode_tolerance )
	        mr_phase_convergence_count++;
       }



       if (mr_phase_convergence_count == problem.nphases ) {
	    psopt_print(workspace,"\n>>> PSOPT: automatic mesh refinement iterations converged as the maximum");
	    psopt_print(workspace,"\n>>> relative error in all phases is lower than algorithm.ode_tolerance\n");
	    break; // break the iterations.
       }

    }

    if (algorithm.mesh_refinement == "automatic" && iter_nodes>= algorithm.mr_min_extrapolation_points && use_global_collocation(algorithm) && iter_nodes<number_of_mesh_refinement_iterations  )
    {

	   // Calculate the next number of nodes for each phase
	      compute_next_mesh_size( problem, algorithm, solution, workspace );

    }


  } // End of mesh refinement iterations loop

  solution.cpu_time = PSOPT_extras::toc();

  get_local_time( solution.end_date_and_time );



  if (algorithm.print_level>0) {

    print_algorithm_summary(problem, algorithm, solution, workspace);

    print_solution_summary(problem, algorithm, solution, workspace);

    print_constraint_summary(problem, solution, workspace);

    print_iterations_summary(problem,algorithm,solution, workspace);

    print_iterations_summary_tex(problem,algorithm,solution, workspace);

    print_psopt_summary(problem, algorithm, solution, workspace);

  }


//  if (workspace) delete  workspace;
  
  return;

}


