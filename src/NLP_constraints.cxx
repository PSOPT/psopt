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




void gg_num( MatrixXd& x, MatrixXd* g, Workspace*  workspace )
{
   // This function implements the NLP inequality  constraints for numerical differentiation

   int j;

   adouble* xad = workspace->xad;

   adouble* gad = workspace->gad;



   for(j=0; j<workspace->nvars; j++)
   {
        xad[j] = x(j); // EIGEN_UPDATE
   }

   gg_ad( xad, gad, workspace );


   for(j=0; j<workspace->ncons; j++)
   {
        (*g)(j) = gad[j].value();  //EIGEN_UPDATE
   }

}




void gg_ad( adouble* xad, adouble* gad, Workspace* workspace )
{
    // This function implements the NLP inequality  constraints for automatic differentiation

    Prob* problem = workspace->problem;

    Alg* algorithm = workspace->algorithm;

    MatrixXd& constraint_scaling = *workspace->constraint_scaling;

    MatrixXd& linkage_scaling = problem->scale.linkages;

    adouble retval=0;
    adouble *states;
    adouble *resid;
    adouble *derivatives;
    adouble *controls;
    adouble *parameters;
    adouble *initial_states;
    adouble *final_states;
    adouble *events;
    adouble *path;
    adouble *linkages;
    adouble time;
    adouble t0;
    adouble tf;
    adouble sum;
    adouble *states_traj;
    adouble *derivs_traj;


    int i, j, iph;

    int phase_offset  = 0;

    linkages = workspace->linkages;

    for(i=0;i< problem->nphases; i++)
    {
        int iphase = i+1;
	     MatrixXd& D               = workspace->D[i];
	     MatrixXd& deriv_scaling   = problem->phase[i].scale.defects;
	     MatrixXd& path_scaling    = problem->phase[i].scale.path;
	     MatrixXd& event_scaling   = problem->phase[i].scale.events;
        double   time_scaling     = problem->phase[i].scale.time;

	if ( problem->multi_segment_flag || workspace->auto_linked_flag ) {
	  iph = 1;
	}
	else {
	  iph = iphase;
	}

	states        = workspace->states[i];
	resid         = workspace->resid[i];
	derivatives   = workspace->derivatives[i];
   controls      = workspace->controls[i];
   parameters    = workspace->parameters[iph-1];
	initial_states= workspace->initial_states[i];
	final_states  = workspace->final_states[i];
	events        = workspace->events[i];
	path          = workspace->path[i];
   states_traj   = workspace->states_traj[i];
   derivs_traj   = workspace->derivs_traj[i];

	int j, k,  l;

   int ncons_phase_i;

	int norder    = problem->phase[i].current_number_of_intervals;

	int nstates   = problem->phase[i].nstates;

	int nevents   = problem->phase[i].nevents;

	int npath     = problem->phase[i].npath;

	int offset;

   ncons_phase_i = get_ncons_phase_i(*problem,i, workspace);

   int path_offset = phase_offset+nstates*(norder+1)+nevents;

   get_parameters(parameters, xad, iphase, workspace );

   get_times(&t0, &tf, xad, iphase, workspace);

	for(k=0; k<norder+1; k++) // EIGEN_UPDATE: k index shifted by -1
   {
             get_states(states, xad, iphase, k, workspace);
             for(j=0;j<nstates;j++) {
                 states_traj[(k)*nstates+j] = states[j];
             }
   }

   if ( workspace->differential_defects != "Hermite-Simpson" && workspace->differential_defects != "trapezoidal") {
	            mtrx_mul_trans(states_traj,&D(0), derivs_traj,nstates, norder+1,norder+1,norder+1); //EIGEN_UPDATE
	}

	for(k=0; k<norder+1; k++)  // EIGEN_UPDATE: index k shifted by -1
   {

            get_controls(controls, xad, iphase, k, workspace);

            get_states(states, xad, iphase, k, workspace);

            if (k==0) {  // EIGEN_UPDATE
               for(j=0;j<nstates;j++)
                    initial_states[j] = states[j];
            }

            if (k==(norder)) { // EIGEN_UPDATE
               for(j=0;j<nstates;j++)
                    final_states[j] = states[j];
            }

            time = convert_to_original_time_ad( (workspace->snodes[i])(k), t0, tf );
            problem->dae(derivatives, path, states, controls, parameters, time, xad, iphase,workspace);
	         if (workspace->enable_nlp_counters) {
		           workspace->solution->mesh_stats[  workspace->current_mesh_refinement_iteration-1 ].n_ode_rhs_evals++;
	         }

            if (workspace->differential_defects != "Hermite-Simpson" && workspace->differential_defects != "trapezoidal"  ) {
                // Differentiation matrix based defects

               for (j=0; j<nstates; j++) {
                     resid[j] = derivs_traj[(k)*nstates+j] - (tf-t0)/2.0*derivatives[j];  // EIGEN_UPDATE
	   	            l = phase_offset+(k)*nstates+j;  // EIGEN_UPDATE
		               gad[l] = resid[j];

		            if ( algorithm->scaling=="user" )
				         gad[l] *=deriv_scaling(j);  // EIGEN_UPDATE
               }

            }
            else if (workspace->differential_defects == "trapezoidal") {
            // Trapezoidal method
                if (k!=(norder)) { // EIGEN_UPDATE
                    adouble* states_next      = workspace->states_next[i];
                    adouble* controls_next    = workspace->controls_next[i];
                    adouble* derivatives_next = workspace->derivatives_next[i];
                    adouble* path_next        = workspace->path_next[i];
                    adouble  time_next        = convert_to_original_time_ad( (workspace->snodes[i])(k+1), t0, tf );
                    adouble  hk               = time_next-time;
                    get_states(states_next, xad, iphase, k+1, workspace);
                    get_controls(controls_next, xad, iphase, k+1, workspace);
                    problem->dae(derivatives_next,path_next,states_next,controls_next,parameters,time_next,xad, iphase,workspace);
		              if (workspace->enable_nlp_counters) {
			               workspace->solution->mesh_stats[  workspace->current_mesh_refinement_iteration-1 ].n_ode_rhs_evals++;
		              }
                    for (j=0; j<nstates; j++) {
                          resid[j] = states_next[j]-states[j]-hk*(derivatives[j]+derivatives_next[j])/2.0;
	   	                 l = phase_offset+(k)*nstates+j; // EIGEN_UPDATE
		                    gad[l] = resid[j]*(tf-t0)/(2.0*hk);
		                   if ( algorithm->scaling=="user" )
		   		                gad[l] *=deriv_scaling(j);  // EIGEN_UPDATE
                         }
                }
                else {
                    for (j=0; j<nstates; j++) {
	   	                 l = phase_offset+(k)*nstates+j; // EIGEN_UPDATE
		                    gad[l] = 0.0;
                    }
                }

            }
            else if (workspace->differential_defects == "Hermite-Simpson") {
              // Hermite Simpson defects
              if (k!=(norder)) { // EIGEN_UPDATE
                    adouble* states_next      = workspace->states_next[i];
                    adouble* controls_next    = workspace->controls_next[i];
                    adouble* derivatives_next = workspace->derivatives_next[i];
                    adouble* path_next        = workspace->path_next[i];
                    adouble* path_bar         = workspace->path_bar[i];
                    adouble* states_bar       = workspace->states_bar[i];
                    adouble* controls_bar     = workspace->controls_bar[i];
                    adouble* derivatives_bar  = workspace->derivatives_bar[i];
                    adouble  time_next        = convert_to_original_time_ad( (workspace->snodes[i])(k+1), t0, tf );
                    adouble  hk               = time_next-time;
                    adouble  time_bar         = time + 0.5*hk;
                    int path_bar_offset = phase_offset+nstates*(norder+1)+nevents+npath*(norder+1);
                    get_controls_bar(controls_bar,xad,iphase,k, workspace);
                    get_states(states_next, xad, iphase, k+1, workspace);
                    get_controls(controls_next, xad, iphase, k+1, workspace);
                    problem->dae(derivatives_next,path_next,states_next,controls_next,parameters,time_next,xad, iphase,workspace);
		              if (workspace->enable_nlp_counters) {
			                workspace->solution->mesh_stats[  workspace->current_mesh_refinement_iteration-1 ].n_ode_rhs_evals++;
		              }
                    for (j=0;j<nstates;j++) {
                        states_bar[j] = 0.5*(states[j]+states_next[j])+hk*(derivatives[j]-derivatives_next[j])/8.0;
                    }

                    problem->dae(derivatives_bar,path_bar,states_bar,controls_bar,parameters,time_bar,xad,iphase,workspace);
		              if (workspace->enable_nlp_counters) {
			                workspace->solution->mesh_stats[  workspace->current_mesh_refinement_iteration-1 ].n_ode_rhs_evals++;
		              }

                    for (j=0; j<nstates; j++) {
                        resid[j] = states_next[j]-states[j]-hk*(derivatives[j]+4.0*derivatives_bar[j]+derivatives_next[j] )/6.0;

	   	               l = phase_offset+(k)*nstates+j;  // EIGEN_UPDATE
		                  gad[l] = resid[j]*(tf-t0)/(2.0*hk);

		                  if ( algorithm->scaling=="user" ) {
  		   		             gad[l] *=deriv_scaling(j); // EIGEN_UPDATE
				                constraint_scaling(l)= deriv_scaling(j); // EIGEN_UPDATE
				            }
                    }
	   	           for (j=0; j<npath; j++)
                    {
				            l = path_bar_offset + (k)*npath + j;  // EIGEN_UPDATE
				            gad[l] = path_bar[j];
				            if ( algorithm->scaling=="user" ) {
				               gad[l] *= path_scaling(j);  // EIGEN_UPDATE
    		                  constraint_scaling(l)= path_scaling(j);  // EIGEN_UPDATE
    		               }
                    }   



                }
                else {
                      for (j=0; j<nstates; j++) {
	   	                   l = phase_offset+(k)*nstates+j;  // EIGEN_UPDATE
		                      gad[l] = 0.0;
                      }
                }

            }


	      for (j=0; j<npath; j++)
	      {
		          l      = path_offset + (k)*npath + j;  // EIGEN_UPDATE
		          gad[l] = path[j];
		          if ( algorithm->scaling=="user" ) {
  			               gad[l] *= path_scaling(j);  // EIGEN_UPDATE
			               constraint_scaling(l)= path_scaling(j);  // EIGEN_UPDATE
			       }

         }



     } // end for( k...)


	offset = phase_offset+nstates*(norder+1);

   problem->events(events, initial_states, final_states, parameters, t0, tf, xad, iphase,workspace);


	// Define nevents constraints functions related to the event inequalities
	for (k=0; k<nevents;k++) {
		j = offset + k;
		gad[j] =  events[k];
		if ( algorithm->scaling=="user" ) {
			  gad[j] *= event_scaling(k);  // EIGEN_UPDATE
		          constraint_scaling(j)= event_scaling(k); // EIGEN_UPDATE
		    }

   	}

        // Add tf >= t0 constraint [ t0MIN-tfMAX <= t0-tf <= 0 ]

      gad[ phase_offset + ncons_phase_i - 1] =  (t0 - tf)*time_scaling;

	   if ( algorithm->scaling=="user" ) {
//	         gad[ phase_offset + ncons_phase_i-1] *= time_scaling;
            constraint_scaling( phase_offset + ncons_phase_i )= time_scaling;
        }

        phase_offset += ncons_phase_i;

  }

  // Now include the phase linkage constraints into the constraint vector

  if (problem->nlinkages) {


     if ( problem->multi_segment_flag) {
  	       auto_link_multiple(linkages, xad, problem->nphases, workspace);
     }
     else {
          problem->linkages( linkages, xad, workspace );
     }

     for(j=0;j<problem->nlinkages;j++)
     {
     	int l = phase_offset+j;
      	gad[l] = linkages[j];
        if ( algorithm->scaling=="user" ) {
  	          gad[l] *= linkage_scaling(j);  // EIGEN_UPDATE
	          constraint_scaling(l)= linkage_scaling(j);  // EIGEN_UPDATE
	    }
     }

  }

  if ( algorithm->scaling=="automatic" || algorithm->scaling!="user" )
  {
	// Scale the constraints using automatic scaling
	if ( workspace->use_constraint_scaling )
	{
		for (j=0;j<workspace->ncons;j++) {
			gad[j] *= constraint_scaling(j);  // EIGEN_UPDATE
		}
	}
  }

  if (workspace->enable_nlp_counters) {
         workspace->solution->mesh_stats[  workspace->current_mesh_refinement_iteration-1 ].n_con_evals++;
  }
	
}

