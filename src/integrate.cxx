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



adouble integrate( adouble (*integrand)(adouble*,adouble*,adouble*,adouble&,adouble*,int, Workspace* workspace), adouble* xad, int iphase, Workspace* workspace )
{
// Evaluates the integral of a user supplied function over the time
// span of phase iphase
    	   Prob& problem = *workspace->problem;
	     Alg&  algorithm = *workspace->algorithm;
        int i = iphase-1;

        int k;

        adouble* controls;
        adouble* states;
	     adouble* states_next;
        adouble* parameters;
        adouble time;
        adouble t0;
        adouble tf;
        adouble ieval;
        adouble retval = 0.0;

	     MatrixXd& w = workspace->w[i];

        int norder    = problem.phase[i].current_number_of_intervals;

	     states        = workspace->states[i].get();
	     states_next   = workspace->states_next[i].get();
        controls      = workspace->controls[i].get();
        parameters    = workspace->parameters[i].get();

        get_parameters(parameters, xad, iphase, workspace);

        get_times(&t0, &tf, xad, iphase, workspace);

	     if ( !use_local_collocation(algorithm) ) {

	        for(k=0; k<norder+1; k++)  // EIGEN_UPDATE: k index shifted by -1
	        {

		         get_controls(controls, xad, iphase, k, workspace);

		         get_states(states, xad, iphase, k, workspace);

		         time = convert_to_original_time_ad( (workspace->snodes[i])(k), t0, tf );

		         ieval = (*integrand)(states,controls,parameters,time,xad,iphase, workspace);

		         retval += ((tf-t0)/2.0)*ieval*w(k);
	       }

	}

	else if ( workspace->algorithm->ir_local_order >= 2
	          && norder >= workspace->algorithm->ir_local_order
	          && (norder % workspace->algorithm->ir_local_order) == 0 ) {

		  // Nie-Kerrigan local representation: the same per-element Lobatto quadrature
		  // that phase_running_cost uses for the objective, and for the same reason. The
		  // state and control on an element are degree-d polynomials through that
		  // element's own d+1 nodes, and the Hermite-Simpson midpoint of the branch below
		  // -- the cubic Hermite state and the midpoint control VARIABLE -- belongs to a
		  // different representation. The midpoint control is not read by anything in this
		  // transcription, so integrating against it returned the integral of a trajectory
		  // the solver never chose: on the test problem min int u^2 with xdot = u over
		  // [0,1] from 0 to 1, whose answer is 1, this routine returned about a half,
		  // because the midpoint control variables sit wherever the initial guess left
		  // them. The snodes within element e are that element's d+1 LGL nodes; the rule
		  // is exact to degree 2d-1 and so matched to the representation.
		  int d = workspace->algorithm->ir_local_order;
		  int M = norder / d;
		  MatrixXd& wl = workspace->ir_lgl_w;         // d+1 LGL weights on [-1,1], sum 2
		  for (int e=0; e<M; e++) {
		      int base = e*d;
		      adouble te0 = convert_to_original_time_ad( (workspace->snodes[i])(base),   t0, tf );
		      adouble te1 = convert_to_original_time_ad( (workspace->snodes[i])(base+d), t0, tf );
		      adouble he  = te1 - te0;
		      for (int r=0; r<=d; r++) {
		          int gk = base + r;
		          get_controls(controls, xad, iphase, gk, workspace);
		          get_states(states,     xad, iphase, gk, workspace);
		          adouble tnode = convert_to_original_time_ad( (workspace->snodes[i])(gk), t0, tf );
		          retval += (he/2.0) * wl(r)
		                    * (*integrand)(states,controls,parameters,tnode,xad,iphase, workspace);
		      }
		  }

	}

	else {


		  // Under Hermite-Simpson the midpoint state of Simpson's rule is the cubic
		  // Hermite value the transcription defines, not the arithmetic mean; see the
		  // note in phase_running_cost (NLP_objective.cxx), which had the same defect.
		  const bool     hs_int    = need_midpoint_controls(algorithm, workspace);
		  const int      ns_int    = problem.phase[i].nstates;
		  adouble* const derivs    = workspace->derivatives[i].get();
		  adouble* const derivs_nx = workspace->derivatives_next[i].get();
		  adouble* const path_scr  = workspace->path[i].get();
		  adouble* const path_scr2 = workspace->path_next[i].get();
		  adouble* const states_bar= workspace->states_bar[i].get();

		  for (k=0; k<norder;k++) {  // EIGEN_UPDATE: k index shifted by -1
		      int l;

		      adouble interval_value = 0.0;

		      get_controls(controls, xad, iphase, k, workspace);
		      get_states(states, xad, iphase, k, workspace);

		      adouble tk = convert_to_original_time_ad( (workspace->snodes[i])(k),   t0, tf );
		      adouble tk1= convert_to_original_time_ad( (workspace->snodes[i])(k+1), t0, tf );

		      adouble h = tk1-tk;

		      interval_value = (*integrand)(states,controls,parameters,tk,xad,iphase, workspace);

		      if (hs_int) {
		          if (k == 0)
		              problem.dae(derivs, path_scr, states, controls, parameters, tk, xad, iphase, workspace);
		          else
		              for (l=0; l<ns_int; l++) derivs[l] = derivs_nx[l];
		      }

		      get_controls(controls, xad, iphase,k+1 , workspace);
		      get_states(states_next, xad, iphase, k+1, workspace);


		      interval_value += (*integrand)(states_next,controls,parameters,tk1,xad,iphase, workspace);

		      if (hs_int) {

			         problem.dae(derivs_nx, path_scr2, states_next, controls, parameters, tk1, xad, iphase, workspace);

			         adouble tmiddle = (tk+tk1)/2.0;

			         for( l =0; l< ns_int; l++ ) {

			             states_bar[l] = 0.5*(states[l]+states_next[l])
			                             + h*(derivs[l]-derivs_nx[l])/8.0;

			        }

			         get_controls_bar(controls,xad,iphase,k, workspace);

			        interval_value += 4.0*(*integrand)(states_bar,controls,parameters,tmiddle,xad,iphase, workspace);


			        interval_value *= h/6.0;

		     }

		     else {
		           interval_value *= h/2.0;
		      }

		      retval+= interval_value;

		  }

       }

       return retval;

}

