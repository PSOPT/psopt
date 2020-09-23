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

	     states        = workspace->states[i];
	     states_next   = workspace->states_next[i];
        controls      = workspace->controls[i];
        parameters    = workspace->parameters[i];

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

	else {


		  for (k=0; k<norder;k++) {  // EIGEN_UPDATE: k index shifted by -1
		      int l;

		      adouble interval_value = 0.0;

		      get_controls(controls, xad, iphase, k, workspace);
		      get_states(states, xad, iphase, k, workspace);

		      adouble tk = convert_to_original_time_ad( (workspace->snodes[i])(k),   t0, tf );
		      adouble tk1= convert_to_original_time_ad( (workspace->snodes[i])(k+1), t0, tf );

		      adouble h = tk1-tk;

		      interval_value = (*integrand)(states,controls,parameters,tk,xad,iphase, workspace);


		      get_controls(controls, xad, iphase,k+1 , workspace);
		      get_states(states_next, xad, iphase, k+1, workspace);


		      interval_value += (*integrand)(states_next,controls,parameters,tk1,xad,iphase, workspace);

		      if (need_midpoint_controls(algorithm, workspace)) {

			         adouble tmiddle = (tk+tk1)/2.0;

			         get_controls_bar(controls,xad,iphase,k, workspace);

			         for( l =0; l< problem.phase[i].nstates; l++ ) {

			             states[l] = 0.5*(states[l]+states_next[l]);

			        }

			        interval_value += 4.0*(*integrand)(states,controls,parameters,tmiddle,xad,iphase, workspace);


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

