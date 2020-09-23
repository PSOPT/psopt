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


adouble ff_ad(adouble* xad, Workspace* workspace)
{
    // This function implements the NLP cost function for automatic differentiation

    adouble retval=0;
    adouble *states;
    adouble *states_next;
    adouble *controls;
    adouble *parameters;
    adouble *initial_states;
    adouble time;
    adouble t0;
    adouble tf;
    adouble sum_cost;
    adouble tmp1;
    adouble integrand_cost;
    adouble endpoint_cost;
    adouble phase_sum_cost;

    Sol& solution = *workspace->solution;


    int i,k, iph;

    Prob& problem = *workspace->problem;

    Alg& algorithm = *workspace->algorithm;

    sum_cost = 0.0;

    for(i=0;i<problem.nphases;i++)
    {
        int iphase = i+1;
	MatrixXd& w = workspace->w[i];

        int norder    = problem.phase[i].current_number_of_intervals;


        phase_sum_cost = 0.0;

	if ( problem.multi_segment_flag || workspace->auto_linked_flag ) {
	  iph = 1;
	}
	else {
	  iph = iphase;
	}

	states        = workspace->states[i];
	states_next   = workspace->states_next[i];
        controls      = workspace->controls[i];
        parameters    = workspace->parameters[iph-1];
        initial_states= workspace->initial_states[i];

        get_parameters(parameters, xad, iphase, workspace);

        get_times(&t0, &tf, xad, iphase, workspace);

	if (problem.phase[i].zero_cost_integrand == true) {
	     phase_sum_cost = 0.0;
	}
	else {

	      if ( !use_local_collocation(algorithm) ) {

		for(k=0; k<norder+1; k++)  // EIGEN_UPDATE: k index shifted by -1
		{

		    get_controls(controls, xad, iphase, k, workspace);

		    get_states(states, xad, iphase, k, workspace);

		    time = convert_to_original_time_ad( (workspace->snodes[i])(k), t0, tf );

		    adouble stime = (workspace->snodes[i])(k);

		    integrand_cost = problem.integrand_cost(states,controls,parameters,time,xad,iphase,workspace);

		    if (workspace->algorithm->collocation_method=="Chebyshev") {
			// Multiply by the reciprocal of the Chebyshev weighting function to evaluate the
			// correct integral.
			integrand_cost *= sqrt(1.0-stime*stime);
		    }

		    (solution.integrand_cost[i])(k) = integrand_cost.value();

		    phase_sum_cost += ((tf-t0)/2.0)*integrand_cost*w(k);

		}

	    }

	    else {

		  for (k=0; k<norder;k++) { // EIGEN_UPDATE: k index shifted by -1
		      // Uses trapezoidal integration to integrate the cost
		      int l;


		      adouble interval_cost = 0.0;

		      adouble integrand;

		      get_controls(controls, xad, iphase, k, workspace);
		      get_states(states, xad, iphase, k, workspace);

		      adouble tk = convert_to_original_time_ad( (workspace->snodes[i])(k),   t0, tf );
		      adouble tk1= convert_to_original_time_ad( (workspace->snodes[i])(k+1), t0, tf );

		      adouble h = tk1-tk;

		      interval_cost = problem.integrand_cost(states,controls,parameters,tk,xad,iphase,workspace);

		      (solution.integrand_cost[i])(k) = interval_cost.value();


		      get_controls(controls, xad, iphase,k+1, workspace );
		      get_states(states_next, xad, iphase, k+1, workspace);

		      integrand = problem.integrand_cost(states_next,controls,parameters,tk1,xad,iphase,workspace);

		      interval_cost += integrand;

              if(k==norder-1) {
                   (solution.integrand_cost[i])(k+1) = integrand.value();
              }

		      if ( need_midpoint_controls(algorithm, workspace) ) {

			  adouble tmiddle = (tk+tk1)/2.0;

			  get_controls_bar(controls,xad,iphase,k, workspace);

			  for( l =0; l< problem.phase[i].nstates; l++ ) {

			          states[l] = 0.5*(states[l]+states_next[l]);



			  }

			  interval_cost += 4.0*problem.integrand_cost(states,controls,parameters,tmiddle,xad,iphase,workspace);


			  interval_cost *= h/6.0;

		      }

		      else {
		         interval_cost *= h/2.0;
		      }

		      phase_sum_cost += interval_cost;

		  }

	    }

	} // End if-else (zero_cost_integrand)

        sum_cost += phase_sum_cost;

        solution.integrated_cost[i] = phase_sum_cost.value();

        get_states(initial_states, xad, iphase, 0, workspace); // EIGEN_UPDATE 1 changed to 0

        get_states(states, xad, iphase, norder, workspace);  // EIGEN_UPDATE norder+1 changed to norder.

        endpoint_cost = problem.endpoint_cost(initial_states,states,parameters,t0,tf,xad,iphase, workspace);

        solution.endpoint_cost[i] = endpoint_cost.value();

	sum_cost += endpoint_cost;


    }

    if (problem.scale.objective != -1)
    {
	retval = sum_cost*problem.scale.objective;
    }
    else {
	retval = sum_cost;
    }

    if (workspace->enable_nlp_counters) {
         solution.mesh_stats[  workspace->current_mesh_refinement_iteration-1 ].n_obj_evals++;
    }

    return (retval);
}



double ff_num(MatrixXd& x, Workspace* workspace)
{
   // This function implements the NLP cost function for numerical differentiation

   int j;

   adouble retval;

   adouble* xad = workspace->xad;


   for(j=0; j<workspace->nvars; j++)
   {
        xad[j] = x(j); //EIGEN_UPDATE
   }

   retval = ff_ad( xad, workspace );

   return (retval.value());

}

