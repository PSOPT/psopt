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


void validate_user_input(Prob& problem, Alg& algorithm, Workspace* workspace)
{
    int i;

    if (algorithm.nlp_method != "IPOPT" && algorithm.nlp_method != "SNOPT" )
       error_message("Incorrect NLP method specified. The only valid values are \"IPOPT\" or \"SNOPT\" ");
    if (algorithm.collocation_method != "Legendre" && algorithm.collocation_method!="Chebyshev" && algorithm.collocation_method!="trapezoidal" && algorithm.collocation_method!="Hermite-Simpson" && algorithm.collocation_method!="Radau" && algorithm.collocation_method!="Gauss")
       error_message("Incorrect collocation method specified. Valid options are \"Legendre\" , \"Chebyshev\", \"trapezoidal\", \"Hermite-Simpson\", \"Radau\", and \"Gauss\" ");
    if (algorithm.scaling != "automatic" && algorithm.scaling!="user")
       error_message("Incorrect scaling option specified. Valid options are \"automatic\" and \"user\" ");
    if (algorithm.transcription_method != "collocation" && algorithm.transcription_method != "integrated-residual")
       error_message("Incorrect transcription_method specified. Valid options are \"collocation\" and \"integrated-residual\" ");
    if (algorithm.transcription_method == "integrated-residual") {
       if (algorithm.collocation_method != "Hermite-Simpson")
          error_message("integrated-residual transcription currently requires collocation_method = \"Hermite-Simpson\" ");
       if (algorithm.ir_residual_nodes < 2)
          error_message("algorithm.ir_residual_nodes must be >= 2 for integrated-residual transcription ");
    }
    if (algorithm.ir_objective != "residual" && algorithm.ir_objective != "cost")
       error_message("Incorrect ir_objective specified. Valid options are \"residual\" and \"cost\" ");
    if (algorithm.ir_objective == "cost") {
       if (algorithm.transcription_method != "integrated-residual")
          error_message("ir_objective=\"cost\" requires transcription_method=\"integrated-residual\" (DAIR optimality step) ");
       if (algorithm.ir_regularization <= 0.0 && algorithm.ir_residual_bound < 0.0)
          error_message("ir_objective=\"cost\" needs the dynamics enforced: set ir_regularization>0 (penalty form) or ir_residual_bound>=0 (robust constraint form) ");
    }
    if (algorithm.ir_residual_bound >= 0.0) {
       if (algorithm.transcription_method != "integrated-residual" || algorithm.ir_objective != "cost")
          error_message("ir_residual_bound>=0 (robust-DAIR constraint form) requires transcription_method=\"integrated-residual\" and ir_objective=\"cost\" ");
    }
    if (algorithm.ir_dair) {
       if (algorithm.transcription_method != "integrated-residual")
          error_message("ir_dair=true requires transcription_method=\"integrated-residual\" ");
       if (algorithm.collocation_method != "Hermite-Simpson")
          error_message("ir_dair=true requires collocation_method=\"Hermite-Simpson\" ");
       if (algorithm.ir_dair_delta_factor <= 0.0)
          error_message("algorithm.ir_dair_delta_factor must be > 0 ");
    }
    if (algorithm.ir_local_order != 0) {
       if (algorithm.ir_local_order < 2)
          error_message("algorithm.ir_local_order must be 0 (legacy cubic-Hermite IR) or >= 2 (Nie-Kerrigan) ");
       if (algorithm.transcription_method != "integrated-residual")
          error_message("ir_local_order requires transcription_method=\"integrated-residual\" ");
       if (algorithm.ir_residual_nodes < algorithm.ir_local_order)
          error_message("ir_residual_nodes must be >= ir_local_order so the residual is well sampled ");
    }
    if (algorithm.ir_regularization < 0.0)
       error_message("algorithm.ir_regularization must be >= 0 ");
    if (algorithm.ir_regularization > 0.0) {
       if (algorithm.collocation_method != "Hermite-Simpson")
          error_message("integrated-residual regularization (ir_regularization>0) currently requires collocation_method = \"Hermite-Simpson\" ");
       if (algorithm.ir_residual_nodes < 2)
          error_message("algorithm.ir_residual_nodes must be >= 2 when ir_regularization>0 ");
    }
    if (algorithm.defect_scaling != "state-based" && algorithm.defect_scaling!="jacobian-based")
       error_message("Incorrect differential defect scaling option specified. Valid options are \"state-based\" and \"jacobian-based\" ");
    if (algorithm.derivatives != "automatic" && algorithm.derivatives!="numerical")
       error_message("Incorrect derivatives option specified. Valid options are \"automatic\" and \"numerical\" ");
    if (algorithm.hessian != "exact" && algorithm.hessian!="limited-memory" && algorithm.hessian!="numerical")
       error_message("Incorrect algorithm.hessian option specified. Valid options are \"limited-memory\", \"exact\" and \"numerical\" ");
    if (algorithm.on_error != "fail-fast" && algorithm.on_error != "fail-soft")
       error_message("Incorrect algorithm.on_error option specified. Valid options are \"fail-fast\" and \"fail-soft\" ");
    if ((algorithm.hessian == "exact" || algorithm.hessian == "numerical") && algorithm.nlp_method !="IPOPT") {
       snprintf(workspace->text,sizeof(workspace->text),"\n*** Warning: the '%s' algorithm.hessian option is only available with the IPOPT solver", algorithm.hessian.c_str());
       psopt_print(workspace,workspace->text);
    }
    if (algorithm.diff_matrix != "standard" && algorithm.diff_matrix!="diff_matrix" && algorithm.diff_matrix!="central-differences" &&  algorithm.diff_matrix!="reduced-roundoff" )
       error_message("Incorrect algorithm.diff_matrix option specified. Valid options are \"standard\", \"reduced-roundoff\" , \"central-differences\" ");


    if (algorithm.hessian == "exact" && algorithm.derivatives !="automatic") {
       snprintf(workspace->text,sizeof(workspace->text),"\n*** Warning: the 'exact' algorithm.hessian option is only available with automatic derivatives");
       psopt_print(workspace,workspace->text);
    }
    if (algorithm.nlp_tolerance <= 0)
       error_message("algorithm.nlp_tolerance must be positive");
    if (algorithm.nlp_iter_max <= 0)
       error_message("algorithm.iter_max must be positive");

    if (algorithm.nsteps_error_integration <= 0)
       error_message("algorithm.nsteps_error_integration must be positive");

    if (algorithm.ode_tolerance <= 0)
       error_message("algorithm.ode_tolerance must be positive");

    if (algorithm.mr_max_growth_factor <= 0 || algorithm.mr_max_growth_factor > 1 )
       error_message("algorithm.mr_max_growth_factor must be in the interval (0,1]");
    if (algorithm.mr_max_iterations <= 0)
       error_message("algorithm.mr_max_iterations must be positive");
    if (algorithm.mr_min_order < 2)
       error_message("algorithm.mr_min_order must be >= 2");
    if (algorithm.mr_max_order < algorithm.mr_min_order)
       error_message("algorithm.mr_max_order must be >= algorithm.mr_min_order");

    if (algorithm.mr_kappa <= 0 || algorithm.mr_kappa>1.0 )
       error_message("algorithm.mr_kappa must be in the interval (0,1]");

    if (algorithm.mr_M1 <= 0  )
       error_message("algorithm.mr_M1 must be positive");

    if (algorithm.switch_order < 0  )
       error_message("algorithm.switch_order must be >= 0");

    if (algorithm.mesh_refinement != "automatic" &&  algorithm.mesh_refinement != "manual" )
       error_message("algorithm.mesh_refinement must either \"manual\" or \"automatic\" ");



    for (i=0;i<problem.nphases;i++)
    {

         if (problem.phase[i].ncontrols>1) {
		      (problem.phase[i].bounds.lower.controls).resize(problem.phase[i].ncontrols,1);
         	(problem.phase[i].bounds.upper.controls).resize(problem.phase[i].ncontrols,1);
         }

         if (problem.phase[i].nstates>1) {
		      (problem.phase[i].bounds.lower.states).resize(problem.phase[i].nstates,1);
         	(problem.phase[i].bounds.upper.states).resize(problem.phase[i].nstates,1);
         }

         if (problem.phase[i].nevents>1) {
	       	(problem.phase[i].bounds.lower.events).resize(problem.phase[i].nevents,1); 
         	(problem.phase[i].bounds.upper.events).resize(problem.phase[i].nevents,1);
         }

         if (problem.phase[i].nparameters>1) {
		      (problem.phase[i].bounds.lower.parameters).resize(problem.phase[i].nparameters,1);
         	(problem.phase[i].bounds.upper.parameters).resize(problem.phase[i].nparameters,1);
         }



         if (problem.phase[i].ncontrols>0 && (problem.phase[i].bounds.lower.controls.array() >  problem.phase[i].bounds.upper.controls.array() ).any() )
         {
                snprintf(workspace->text,sizeof(workspace->text),"Infeasible control variable bounds supplied by the user in phase %i",i);
 		error_message(workspace->text);
         }

         if (problem.phase[i].nstates>0 && ( problem.phase[i].bounds.lower.states.array() >  problem.phase[i].bounds.upper.states.array() ).any() )
         {
                snprintf(workspace->text,sizeof(workspace->text),"Infeasible state variable bounds supplied by the user in phase %i",i);
 		error_message(workspace->text);
         }

         if (problem.phase[i].nevents>0 && ( problem.phase[i].bounds.lower.events.array() >  problem.phase[i].bounds.upper.events.array() ).any() )
         {
                snprintf(workspace->text,sizeof(workspace->text),"Infeasible event bounds supplied by the user in phase %i",i);
 		error_message(workspace->text);
         }

         if (problem.phase[i].nparameters>0 && ( problem.phase[i].bounds.lower.parameters.array() >  problem.phase[i].bounds.upper.parameters.array() ).any() )
         {
                snprintf(workspace->text,sizeof(workspace->text),"Infeasible static parameter bounds supplied by the user in phase %i",i);
 		error_message(workspace->text);
         }

         if ( problem.phase[i].bounds.lower.StartTime >  problem.phase[i].bounds.upper.StartTime  )
         {
                snprintf(workspace->text,sizeof(workspace->text),"Infeasible start time bounds supplied by the user in phase %i",i);
 		error_message(workspace->text);
         }

         if ( problem.phase[i].bounds.lower.EndTime >  problem.phase[i].bounds.upper.EndTime  )
         {
                snprintf(workspace->text,sizeof(workspace->text),"Infeasible end time bounds supplied by the user in phase %i",i);
 		error_message(workspace->text);
         }

         // hp-adaptive fixed mesh (Route B, increment 1). When a phase carries an
         // explicit multi-interval mesh via hp_orders / hp_breakpoints, validate it:
         // Radau-only for now; K orders and K-1 breaks; breaks strictly increasing in
         // (0,1); each interval order >= 2 (a single-point interval is degenerate).
         if ( hp_mesh_active(problem.phase[i]) )
         {
            if ( algorithm.collocation_method != "Radau" && algorithm.collocation_method != "Gauss" && algorithm.collocation_method != "Legendre" && algorithm.collocation_method != "Chebyshev" )
            {
               snprintf(workspace->text,sizeof(workspace->text),"hp-adaptive mesh (hp_orders) in phase %i currently requires collocation_method = \"Radau\", \"Gauss\", \"Legendre\" or \"Chebyshev\"",i+1);
               error_message(workspace->text);
            }
            int Khp = (int) problem.phase[i].hp_orders.size();
            if ( (int) problem.phase[i].hp_breakpoints.size() != Khp-1 )
            {
               snprintf(workspace->text,sizeof(workspace->text),"In phase %i, length(hp_breakpoints) must equal length(hp_orders)-1",i+1);
               error_message(workspace->text);
            }
            for (int kk=0; kk<Khp; kk++)
            {
               if ( problem.phase[i].hp_orders(kk) < 2 )
               {
                  snprintf(workspace->text,sizeof(workspace->text),"In phase %i, every hp_orders entry must be >= 2",i+1);
                  error_message(workspace->text);
               }
            }
            for (int kk=0; kk<Khp-1; kk++)
            {
               double bk = problem.phase[i].hp_breakpoints(kk);
               if ( bk <= 0.0 || bk >= 1.0 )
               {
                  snprintf(workspace->text,sizeof(workspace->text),"In phase %i, hp_breakpoints must lie strictly in the open interval (0,1)",i+1);
                  error_message(workspace->text);
               }
               if ( kk>0 && bk <= problem.phase[i].hp_breakpoints(kk-1) )
               {
                  snprintf(workspace->text,sizeof(workspace->text),"In phase %i, hp_breakpoints must be strictly increasing",i+1);
                  error_message(workspace->text);
               }
            }
         }

	 if ( problem.phase[i].nobserved >  0  )
         {

//	        if ( problem.phase[i].observation_nodes(1) != problem.phase[i].bounds.lower.StartTime && problem.phase[i].observation_nodes(1) != problem.phase[i].bounds.upper.StartTime )
//		{
//                 snprintf(workspace->text,sizeof(workspace->text),"Initial observation time must be equal to start time in phase %i",i);
// 		  error_message(workspace->text);
//		}

//	        if ( fabs( problem.phase[i].observation_nodes("end") - problem.phase[i].bounds.lower.EndTime ) > 0.0001 && fabs(problem.phase[i].observation_nodes("end") - problem.phase[i].bounds.upper.EndTime )>0.0001 )
//		{
//		  snprintf(workspace->text,sizeof(workspace->text),"Final observation time must be equal to end time in phase %i",i+1);
// 		  error_message(workspace->text);
//		}

	        if ( problem.phase[i].nsamples != length( problem.phase[i].observation_nodes)  )
		{
		  snprintf(workspace->text,sizeof(workspace->text),"Length of observation nodes vector in phase %i must be equal to problem.phases(%i).nsamples",i+1, i+1);
 		  error_message(workspace->text);
		}

		if (  isEmpty( problem.phase[i].residual_weights ) ) {
                    problem.phase[i].residual_weights = ones( problem.phase[i].nobserved, problem.phase[i].nsamples );
		}

		if ( problem.phase[i].nsamples !=  problem.phase[i].residual_weights.cols() )
		{
		  snprintf(workspace->text,sizeof(workspace->text),"The number of columns of the residual weight vector in phase %i must be equal to problem.phases(%i).nsamples",i+1, i+1);
 		  error_message(workspace->text);
		}

		if ( problem.phase[i].nobserved !=  problem.phase[i].residual_weights.rows() )
		{
		  snprintf(workspace->text,sizeof(workspace->text),"The number of rows of the residual weight vector in phase %i must be equal to the number of observed variables", i+1 );
 		  error_message(workspace->text);
		}

		if (  isEmpty(problem.phase[i].covariance)  ) {
                    problem.phase[i].covariance = eye( problem.phase[i].nobserved );
		}

		if ( problem.phase[i].nobserved !=  problem.phase[i].covariance.rows() && problem.phase[i].nobserved !=  problem.phase[i].covariance.cols()  )
		{
		  snprintf(workspace->text,sizeof(workspace->text),"The number of rows and columns of matrix problem.phases(%i).covariance must be equal to problem.phases(%i).nobserved",i+1, i+1);
 		  error_message(workspace->text);
		}

		if ( !isSymmetric(problem.phase[i].covariance)  )
		{
		  snprintf(workspace->text,sizeof(workspace->text),"Matrix problem.phases(%i).covariance must be symmetric",i+1);
 		  error_message(workspace->text);
		}

		if ( problem.phase[i].regularization_factor< 0  )
		{
		  snprintf(workspace->text,sizeof(workspace->text),"problem.phases(%i).regularization_factor must be positive",i+1);
 		  error_message(workspace->text);
		}

		problem.integrand_cost 	= NULL;
                problem.endpoint_cost 	= &endpoint_cost_for_parameter_estimation;
                problem.phase[i].zero_cost_integrand = true;



         }

    }

   if (problem.nlinkages>0 && ( problem.bounds.lower.linkage.array() >  problem.bounds.upper.linkage.array() ).any() )
   {
         snprintf(workspace->text,sizeof(workspace->text),"Infeasible phase linkage bounds supplied by the user");
         error_message(workspace->text);
   }

   if ( length(problem.bounds.lower.times) !=  length(problem.bounds.upper.times) || (!isEmpty(problem.bounds.lower.times) && length(problem.bounds.lower.times)!=problem.nphases+1) )
   {
         snprintf(workspace->text,sizeof(workspace->text),"Incorrect length of problem.bounds.lower.times or problem.bounds.upper.times");
         error_message(workspace->text);
   }

   if ( !isEmpty(problem.bounds.lower.times) ) {
     for (i=0;i<problem.nphases;i++) { //EIGEN_UPDATE
	    problem.phase[i].bounds.lower.StartTime = problem.bounds.lower.times(i);
	    problem.phase[i].bounds.upper.StartTime = problem.bounds.upper.times(i);
	    problem.phase[i].bounds.lower.EndTime   = problem.bounds.lower.times(i+1);
	    problem.phase[i].bounds.upper.EndTime   = problem.bounds.upper.times(i+1);
     }
   }

}
