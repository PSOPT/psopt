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


void validate_user_input(Prob& problem, Alg& algorithm, Workspace* workspace)
{
    char m[500];
    int i;

    if (algorithm.nlp_method != "IPOPT" && algorithm.nlp_method != "SNOPT" )
       error_message("Incorrect NLP method specified. The only valid values are \"IPOPT\" or \"SNOPT\" ");
    if (algorithm.collocation_method != "Legendre" && algorithm.collocation_method!="Chebyshev" && algorithm.collocation_method!="trapezoidal" && algorithm.collocation_method!="Hermite-Simpson")
       error_message("Incorrect collocation method specified. Valid options are \"Legendre\" , \"Chebyshev\", \"trapezoidal\", and \"Hermite-Simpson\" ");
    if (algorithm.scaling != "automatic" && algorithm.scaling!="user")
       error_message("Incorrect scaling option specified. Valid options are \"automatic\" and \"user\" ");
    if (algorithm.defect_scaling != "state-based" && algorithm.defect_scaling!="jacobian-based")
       error_message("Incorrect differential defect scaling option specified. Valid options are \"state-based\" and \"jacobian-based\" ");
    if (algorithm.derivatives != "automatic" && algorithm.derivatives!="numerical")
       error_message("Incorrect derivatives option specified. Valid options are \"automatic\" and \"numerical\" ");
    if (algorithm.hessian != "exact" && algorithm.hessian!="limited-memory")
       error_message("Incorrect algorithm.hessian option specified. Valid options are \"limited-memory\" and \"exact\" ");
    if (algorithm.hessian == "exact" && algorithm.nlp_method !="IPOPT") {
       sprintf(workspace->text,"\n*** Warning: the 'exact' algorithm.hessian option is only available with the IPOPT solver");
       psopt_print(workspace,workspace->text);
    }
    if (algorithm.diff_matrix != "standard" && algorithm.diff_matrix!="diff_matrix" && algorithm.diff_matrix!="central-differences" &&  algorithm.diff_matrix!="reduced-roundoff" )
       error_message("Incorrect algorithm.diff_matrix option specified. Valid options are \"standard\", \"reduced-roundoff\" , \"central-differences\" ");


    if (algorithm.hessian == "exact" && algorithm.derivatives !="automatic") {
       sprintf(workspace->text,"\n*** Warning: the 'exact' algorithm.hessian option is only available with automatic derivatives");
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

    if (algorithm.mr_max_increment_factor <= 0 || algorithm.mr_max_increment_factor > 1 )
       error_message("algorithm.mr_max_increment_factor must be in the interval (0,1]");
    if (algorithm.mr_max_iterations <= 0)
       error_message("algorithm.mr_max_iterations must be positive");
    if (algorithm.mr_min_extrapolation_points < 2)
       error_message("algorithm.mr_min_extrapolation_points must be >= 2");

    if (algorithm.mr_initial_increment <= 0)
       error_message("algorithm.mr_initial_increment must be positive");

    if (algorithm.mr_kappa <= 0 || algorithm.mr_kappa>1.0 )
       error_message("algorithm.mr_kappa must be in the interval (0,1]");

    if (algorithm.mr_M1 <= 0  )
       error_message("algorithm.mr_M1 must be positive");

    if (algorithm.switch_order < 0  )
       error_message("algorithm.mr_M1 must be >= 0");

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
                sprintf(m,"Infeasible control variable bounds supplied by the user in phase %i",i);
 		error_message(m);
         }

         if (problem.phase[i].nstates>0 && ( problem.phase[i].bounds.lower.states.array() >  problem.phase[i].bounds.upper.states.array() ).any() )
         {
                sprintf(m,"Infeasible state variable bounds supplied by the user in phase %i",i);
 		error_message(m);
         }

         if (problem.phase[i].nevents>0 && ( problem.phase[i].bounds.lower.events.array() >  problem.phase[i].bounds.upper.events.array() ).any() )
         {
                sprintf(m,"Infeasible event bounds supplied by the user in phase %i",i);
 		error_message(m);
         }

         if (problem.phase[i].nparameters>0 && ( problem.phase[i].bounds.lower.parameters.array() >  problem.phase[i].bounds.upper.parameters.array() ).any() )
         {
                sprintf(m,"Infeasible static parameter bounds supplied by the user in phase %i",i);
 		error_message(m);
         }

         if ( problem.phase[i].bounds.lower.StartTime >  problem.phase[i].bounds.upper.StartTime  )
         {
                sprintf(m,"Infeasible start time bounds supplied by the user in phase %i",i);
 		error_message(m);
         }

         if ( problem.phase[i].bounds.lower.EndTime >  problem.phase[i].bounds.upper.EndTime  )
         {
                sprintf(m,"Infeasible end time bounds supplied by the user in phase %i",i);
 		error_message(m);
         }

	 if ( problem.phase[i].nobserved >  0  )
         {

//	        if ( problem.phase[i].observation_nodes(1) != problem.phase[i].bounds.lower.StartTime && problem.phase[i].observation_nodes(1) != problem.phase[i].bounds.upper.StartTime )
//		{
//                 sprintf(m,"Initial observation time must be equal to start time in phase %i",i);
// 		  error_message(m);
//		}

//	        if ( fabs( problem.phase[i].observation_nodes("end") - problem.phase[i].bounds.lower.EndTime ) > 0.0001 && fabs(problem.phase[i].observation_nodes("end") - problem.phase[i].bounds.upper.EndTime )>0.0001 )
//		{
//		  sprintf(m,"Final observation time must be equal to end time in phase %i",i+1);
// 		  error_message(m);
//		}

	        if ( problem.phase[i].nsamples != length( problem.phase[i].observation_nodes)  )
		{
		  sprintf(m,"Length of observation nodes vector in phase %i must be equal to problem.phases(%i).nsamples",i+1, i+1);
 		  error_message(m);
		}

		if (  isEmpty( problem.phase[i].residual_weights ) ) {
                    problem.phase[i].residual_weights = ones( problem.phase[i].nobserved, problem.phase[i].nsamples );
		}

		if ( problem.phase[i].nsamples !=  problem.phase[i].residual_weights.cols() )
		{
		  sprintf(m,"The number of columns of the residual weight vector in phase %i must be equal to problem.phases(%i).nsamples",i+1, i+1);
 		  error_message(m);
		}

		if ( problem.phase[i].nobserved !=  problem.phase[i].residual_weights.rows() )
		{
		  sprintf(m,"The number of rows of the residual weight vector in phase %i must be equal to the number of observed variables", i+1 );
 		  error_message(m);
		}

		if (  isEmpty(problem.phase[i].covariance)  ) {
                    problem.phase[i].covariance = eye( problem.phase[i].nobserved );
		}

		if ( problem.phase[i].nobserved !=  problem.phase[i].covariance.rows() && problem.phase[i].nobserved !=  problem.phase[i].covariance.cols()  )
		{
		  sprintf(m,"The number of rows and columns of matrix problem.phases(%i).covariance must be equal to problem.phases(%i).nobserved",i+1, i+1);
 		  error_message(m);
		}

		if ( !isSymmetric(problem.phase[i].covariance)  )
		{
		  sprintf(m,"Matrix problem.phases(%i).covariance must be symmetric",i+1);
 		  error_message(m);
		}

		if ( problem.phase[i].regularization_factor< 0  )
		{
		  sprintf(m,"problem.phases(%i).regularization_factor must be positive",i+1);
 		  error_message(m);
		}

		problem.integrand_cost 	= NULL;
                problem.endpoint_cost 	= &endpoint_cost_for_parameter_estimation;
                problem.phase[i].zero_cost_integrand = true;



         }

    }

   if (problem.nlinkages>0 && ( problem.bounds.lower.linkage.array() >  problem.bounds.upper.linkage.array() ).any() )
   {
         sprintf(m,"Infeasible phase linkage bounds supplied by the user");
         error_message(m);
   }

   if ( length(problem.bounds.lower.times) !=  length(problem.bounds.upper.times) || (!isEmpty(problem.bounds.lower.times) && length(problem.bounds.lower.times)!=problem.nphases+1) )
   {
         sprintf(m,"Incorrect length of problem.bounds.lower.times or problem.bounds.upper.times");
         error_message(m);
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
