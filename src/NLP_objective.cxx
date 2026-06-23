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


// Integral over one phase of ||xdot - f||^2, using the per-interval cubic Hermite state
// representation (knots x_k,x_{k+1} with endpoint slopes f_k,f_{k+1}) and the quadratic
// control through (u_k, ubar_k, u_{k+1}), evaluated on the Gauss-Legendre residual grid
// (workspace->ir_nodes/ir_weights on [0,1], distinct from the knots/midpoint where the
// cubic makes the residual vanish). Shared by:
//   * the integrated-residual transcription  (increment 1: this IS the phase objective);
//   * integrated-residual regularisation      (increment 2: added to the user objective
//     with weight algorithm.ir_regularization, defects retained so the costates survive).
static adouble integrated_residual_phase(int i, int iphase, adouble* xad,
        adouble t0, adouble tf, adouble* parameters, Workspace* workspace)
{
    Prob& problem = *workspace->problem;
    int norder    = problem.phase[i].current_number_of_intervals;
    int nstates   = problem.phase[i].nstates;
    int ncontrols = problem.phase[i].ncontrols;
    adouble* xk    = workspace->states[i].get();
    adouble* xk1   = workspace->states_next[i].get();
    adouble* fk    = workspace->derivatives[i].get();
    adouble* fk1   = workspace->derivatives_next[i].get();
    adouble* uk    = workspace->controls[i].get();
    adouble* uk1   = workspace->controls_next[i].get();
    adouble* ubar  = workspace->controls_bar[i].get();
    adouble* xq    = workspace->states_bar[i].get();        // x(tau) scratch
    adouble* fq    = workspace->derivatives_bar[i].get();   // f(x(tau),u(tau)) scratch
    adouble* pscr  = workspace->path_bar[i].get();          // dae path output (discarded)
    std::vector<adouble> uq( (ncontrols>0)? ncontrols : 1 );
    adouble* uptr  = (ncontrols>0) ? uq.data() : uk;
    MatrixXd& tau  = workspace->ir_nodes;
    MatrixXd& wq   = workspace->ir_weights;
    int m = workspace->ir_m;
    int j, k;
    adouble R = 0.0;

    for (k=0; k<norder; k++) {
        adouble tk  = convert_to_original_time_ad( (workspace->snodes[i])(k),   t0, tf );
        adouble tk1 = convert_to_original_time_ad( (workspace->snodes[i])(k+1), t0, tf );
        adouble hk  = tk1 - tk;

        get_states(xk,  xad, iphase, k,   workspace);
        get_states(xk1, xad, iphase, k+1, workspace);
        if (ncontrols>0) {
            get_controls(uk,   xad, iphase, k,   workspace);
            get_controls(uk1,  xad, iphase, k+1, workspace);
            get_controls_bar(ubar, xad, iphase, k, workspace);
        }
        problem.dae(fk,  pscr, xk,  uk,  parameters, tk,  xad, iphase, workspace);
        problem.dae(fk1, pscr, xk1, uk1, parameters, tk1, xad, iphase, workspace);

        for (int q=0; q<m; q++) {
            double s   = tau(q);                 // local coordinate in [0,1]
            double h00 =  2*s*s*s - 3*s*s + 1;
            double h10 =      s*s*s - 2*s*s + s;
            double h01 = -2*s*s*s + 3*s*s;
            double h11 =      s*s*s -   s*s;
            double g00 =  6*s*s - 6*s;           // d/ds of the basis (for xdot)
            double g10 =  3*s*s - 4*s + 1;
            double g01 = -6*s*s + 6*s;
            double g11 =  3*s*s - 2*s;
            for (j=0; j<nstates; j++)
                xq[j] = h00*xk[j] + h10*hk*fk[j] + h01*xk1[j] + h11*hk*fk1[j];
            if (ncontrols>0) {
                double L0 =  2.0*(s-0.5)*(s-1.0);  // quadratic through u_k, ubar_k, u_{k+1}
                double Lm = -4.0*s*(s-1.0);
                double L1 =  2.0*s*(s-0.5);
                for (int c=0;c<ncontrols;c++) uq[c] = L0*uk[c] + Lm*ubar[c] + L1*uk1[c];
            }
            adouble tq = tk + s*hk;
            problem.dae(fq, pscr, xq, uptr, parameters, tq, xad, iphase, workspace);

            adouble rsq = 0.0;
            for (j=0; j<nstates; j++) {
                adouble xdot = (g00*xk[j] + g01*xk1[j])/hk + g10*fk[j] + g11*fk1[j];
                adouble rj   = xdot - fq[j];
                rsq += rj*rj;
            }
            R += hk * wq(q) * rsq;
        }
    }
    return R;
}


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

	states        = workspace->states[i].get();
	states_next   = workspace->states_next[i].get();
        controls      = workspace->controls[i].get();
        parameters    = workspace->parameters[iph-1].get();
        initial_states= workspace->initial_states[i].get();

        get_parameters(parameters, xad, iphase, workspace);

        get_times(&t0, &tf, xad, iphase, workspace);

	if ( workspace->transcription_method == "integrated-residual" && algorithm.ir_objective != "cost" ) {
	    // Integrated-residual transcription, feasibility step (increment 1 / DAIR
	    // feasibility): the whole phase objective is the integral of ||xdot - f||^2;
	    // no user cost, no endpoint cost. (ir_objective=="cost" instead minimises the
	    // user cost J with the residual as a penalty -- see below -- for the DAIR
	    // optimality step, increment 3.)
	    phase_sum_cost = integrated_residual_phase(i, iphase, xad, t0, tf, parameters, workspace);
	}
	else if (problem.phase[i].zero_cost_integrand == true) {
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

        // Integrated-residual penalty term:
        //  * increment 2 (collocation): J + rho*R with the defects retained as hard
        //    constraints (costates survive);
        //  * increment 3 DAIR optimality step (transcription=="integrated-residual",
        //    ir_objective=="cost"): J + rho*R with the defects dropped, so the dynamics
        //    are enforced only through the residual penalty and the alternating loop.
        // Skipped only for the pure-residual feasibility step (increment 1), whose
        // objective is already R.
        bool pure_residual = ( workspace->transcription_method == "integrated-residual"
                               && algorithm.ir_objective != "cost" );
        if ( !pure_residual && algorithm.ir_regularization > 0.0 ) {
            phase_sum_cost += algorithm.ir_regularization
                            * integrated_residual_phase(i, iphase, xad, t0, tf, parameters, workspace);
        }

        sum_cost += phase_sum_cost;

        solution.integrated_cost[i] = phase_sum_cost.value();

        get_states(initial_states, xad, iphase, 0, workspace); // EIGEN_UPDATE 1 changed to 0

        get_states(states, xad, iphase, norder, workspace);  // EIGEN_UPDATE norder+1 changed to norder.

        if ( workspace->transcription_method == "integrated-residual" && algorithm.ir_objective != "cost" ) {
            // Pure least-squares feasibility objective (increment 1 / DAIR feasibility):
            // no endpoint cost. The optimality step (ir_objective=="cost") keeps it.
            endpoint_cost = 0.0;
        }
        else {
            endpoint_cost = problem.endpoint_cost(initial_states,states,parameters,t0,tf,xad,iphase, workspace);
        }

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

   adouble* xad = workspace->xad.get();


   for(j=0; j<workspace->nvars; j++)
   {
        xad[j] = x(j); //EIGEN_UPDATE
   }

   retval = ff_ad( xad, workspace );

   return (retval.value());

}

