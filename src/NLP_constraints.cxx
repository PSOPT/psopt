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




void gg_num( MatrixXd& x, MatrixXd* g, Workspace*  workspace )
{
   // This function implements the NLP inequality  constraints for numerical differentiation

   int j;

   adouble* xad = workspace->xad.get();

   adouble* gad = workspace->gad.get();



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

    linkages = workspace->linkages.get();

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

	states        = workspace->states[i].get();
	resid         = workspace->resid[i].get();
	derivatives   = workspace->derivatives[i].get();
   controls      = workspace->controls[i].get();
   parameters    = workspace->parameters[iph-1].get();
	initial_states= workspace->initial_states[i].get();
	final_states  = workspace->final_states[i].get();
	events        = workspace->events[i].get();
	path          = workspace->path[i].get();
   states_traj   = workspace->states_traj[i].get();
   derivs_traj   = workspace->derivs_traj[i].get();

	int j, k,  l;

   int ncons_phase_i;

	int norder    = problem->phase[i].current_number_of_intervals;

	int nstates   = problem->phase[i].nstates;

	int nevents   = problem->phase[i].nevents;

	int npath     = problem->phase[i].npath;

	int offset;

   // Multi-interval Gauss bookkeeping. Gauss collocates strictly interior points, so the
   // K stored breakpoints (tau_0..tau_{K-1}) interleaved with each interval's Gauss points
   // are all non-collocated. Mark the breakpoint storage indices (zero-padded defect rows)
   // and record per-interval Gauss orders/left-breakpoint indices for the K defining
   // constraints. Single-block Gauss => K=1, one breakpoint at node 0 (the initial node).
   std::vector<char> gauss_is_bp;
   std::vector<int>  gauss_bp_idx;     // storage index of each interval's left breakpoint
   std::vector<int>  gauss_n;          // per-interval Gauss order
   int gauss_K = 1;
   if ( workspace->differential_defects == "Gauss" ) {
       if ( hp_mesh_active(problem->phase[i]) ) {
           gauss_K = hp_num_intervals(problem->phase[i]);
           gauss_n.resize(gauss_K);
           for (int jj=0; jj<gauss_K; jj++) gauss_n[jj] = hp_interval_order(problem->phase[i], jj);
       } else {
           gauss_n.assign(1, norder);                 // single block: norder Gauss points
       }
       gauss_bp_idx.resize(gauss_K);
       gauss_is_bp.assign(norder+1, 0);
       int gidx = 0;
       for (int jj=0; jj<gauss_K; jj++) {
           gauss_bp_idx[jj] = gidx;                   // left breakpoint (non-collocated)
           gauss_is_bp[gidx] = 1;
           gidx += 1 + gauss_n[jj];                   // breakpoint + this interval's Gauss pts
       }
       // gidx == norder+1; the terminal tau_K is the appended x_f, not stored here.
   }

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

            if (workspace->differential_defects != "Hermite-Simpson" && workspace->differential_defects != "trapezoidal" && workspace->differential_defects != "Radau" && workspace->differential_defects != "Gauss"  ) {
                // Differentiation matrix based defects

               for (j=0; j<nstates; j++) {
                     resid[j] = derivs_traj[(k)*nstates+j] - (tf-t0)/2.0*derivatives[j];  // EIGEN_UPDATE
	   	            l = phase_offset+(k)*nstates+j;  // EIGEN_UPDATE
		               gad[l] = resid[j];

		            if ( algorithm->scaling=="user" )
				         gad[l] *=deriv_scaling(j);  // EIGEN_UPDATE
               }

            }
            else if (workspace->differential_defects == "Radau") {
                // Radau pseudospectral: rectangular differentiation matrix occupies the
                // top norder rows of D (collocation defects); the terminal node (k==norder)
                // is the non-collocated boundary -> zero-padded (free terminal state),
                // exactly like the trapezoidal/Hermite-Simpson final row.
               for (j=0; j<nstates; j++) {
                     l = phase_offset+(k)*nstates+j;
                     if (k != norder) {
                         resid[j] = derivs_traj[(k)*nstates+j] - (tf-t0)/2.0*derivatives[j];
                         gad[l] = resid[j];
                         if ( algorithm->scaling=="user" )
                                 gad[l] *= deriv_scaling(j);
                     }
                     else {
                         gad[l] = 0.0;
                     }
               }

            }
            else if (workspace->differential_defects == "Gauss") {
                // Legendre-Gauss pseudospectral: collocation defects at the interior Gauss
                // points; every stored breakpoint (single-block: just the initial node k=0;
                // multi-interval: tau_0..tau_{K-1}) is non-collocated -> zero-padded. The
                // breakpoint/terminal states are fixed by the K Gauss-quadrature defining
                // constraints assembled after this loop.
               for (j=0; j<nstates; j++) {
                     l = phase_offset+(k)*nstates+j;
                     if ( !gauss_is_bp[k] ) {
                         resid[j] = derivs_traj[(k)*nstates+j] - (tf-t0)/2.0*derivatives[j];
                         gad[l] = resid[j];
                         if ( algorithm->scaling=="user" )
                                 gad[l] *= deriv_scaling(j);
                     }
                     else {
                         gad[l] = 0.0;
                     }
               }

            }
            else if (workspace->differential_defects == "trapezoidal") {
            // Trapezoidal method
                if (k!=(norder)) { // EIGEN_UPDATE
                    adouble* states_next      = workspace->states_next[i].get();
                    adouble* controls_next    = workspace->controls_next[i].get();
                    adouble* derivatives_next = workspace->derivatives_next[i].get();
                    adouble* path_next        = workspace->path_next[i].get();
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
              if ( workspace->transcription_method == "integrated-residual" ) {
                  // Integrated-residual transcription: the dynamics are enforced via the
                  // integrated-residual objective (see ff_ad), not by collocation defects.
                  // Keep the defect rows but zero them (same zero-padding convention as the
                  // final HS/Radau/Gauss row) so the constraint layout, bounds and the
                  // event/path/costate offsets are completely undisturbed.
                  for (j=0; j<nstates; j++) {
                      l = phase_offset+(k)*nstates+j;
                      gad[l] = 0.0;
                  }
              }
              else
              // Hermite Simpson defects
              if (k!=(norder)) { // EIGEN_UPDATE
                    adouble* states_next      = workspace->states_next[i].get();
                    adouble* controls_next    = workspace->controls_next[i].get();
                    adouble* derivatives_next = workspace->derivatives_next[i].get();
                    adouble* path_next        = workspace->path_next[i].get();
                    adouble* path_bar         = workspace->path_bar[i].get();
                    adouble* states_bar       = workspace->states_bar[i].get();
                    adouble* controls_bar     = workspace->controls_bar[i].get();
                    adouble* derivatives_bar  = workspace->derivatives_bar[i].get();
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

   // Gauss: final_states is the appended terminal-state variable, not the last stored node.
   if ( workspace->differential_defects == "Gauss" ) {
       get_gauss_terminal_states(final_states, xad, iphase, workspace);
   }

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

        // Radau: in-solve terminal-control interpolation pin. Pin the non-collocated
        // terminal control U(norder) to the Lagrange interpolant of the collocation
        // controls at tau = +1, as ncontrols equality constraints placed just before
        // the t0<=tf constraint (which stays the last slot of the phase block).
        //
        // Single block: the interpolant uses all norder collocation controls
        // U(0..norder-1). hp multi-interval mesh: tau=+1 belongs to the LAST interval,
        // so the interpolant must be that interval's LOCAL polynomial, built only on
        // its collocation nodes U(m0..norder-1) with m0 = N_eff - n_last. Using a
        // global interpolant over the clustered multi-interval node set would be both
        // wrong (not the hp local extrapolant) and numerically unstable. For a single
        // interval m0 = 0, so this reduces exactly to the legacy pin (bit-identical).
        if ( workspace->differential_defects == "Radau" ) {
            int ncontrols = problem->phase[i].ncontrols;
            int pin_base  = phase_offset + nstates*(norder+1) + nevents + npath*(norder+1);
            MatrixXd& sn  = workspace->snodes[i];
            int m0 = 0;                                              // first node of the interpolating interval
            if ( hp_mesh_active(problem->phase[i]) )
                m0 = norder - hp_interval_order(problem->phase[i], hp_num_intervals(problem->phase[i]) - 1);
            double xe = sn(norder);                                 // terminal node, tau = +1
            get_controls(controls, xad, iphase, norder, workspace); // terminal control
            for (int l2=0; l2<ncontrols; l2++) gad[pin_base+l2] = controls[l2];
            for (int m=m0; m<norder; m++) {                         // subtract local interpolant
                double Lwm = 1.0;
                for (int jj=m0; jj<norder; jj++) if (jj!=m) Lwm *= (xe - sn(jj))/(sn(m)-sn(jj));
                get_controls(controls, xad, iphase, m, workspace);
                for (int l2=0; l2<ncontrols; l2++) gad[pin_base+l2] -= Lwm*controls[l2];
            }
        }

        // Gauss: terminal/interface defining (quadrature) constraints, one per interval:
        //   x(tau_{j+1}) - x(tau_j) - (tf-t0)/2 * sum_{k in Gauss_j} w_k f_k = 0   (nstates each)
        // x(tau_K) is the appended terminal x_f; interior breakpoint states tau_1..tau_{K-1}
        // are ordinary stored nodes (shared once => C0 continuity). Composite weights already
        // carry the sub-interval scaling Dtau_j/2. Placed just before the t0<=tf constraint.
        // Single block (K=1) reduces to the original x_f - x_0 - (tf-t0)/2 sum w_k f_k.
        if ( workspace->differential_defects == "Gauss" ) {
            int quad_base = phase_offset + nstates*(norder+1) + nevents + npath*(norder+1);
            for (int jint=0; jint<gauss_K; jint++) {
                int cbase = quad_base + jint*nstates;
                // + right boundary state (interior breakpoint node, or x_f for the last interval)
                if (jint < gauss_K-1) get_states(states, xad, iphase, gauss_bp_idx[jint+1], workspace);
                else                  get_gauss_terminal_states(states, xad, iphase, workspace);
                for (j=0; j<nstates; j++) gad[cbase+j] = states[j];
                // - left boundary state
                get_states(states, xad, iphase, gauss_bp_idx[jint], workspace);
                for (j=0; j<nstates; j++) gad[cbase+j] -= states[j];
                // - (tf-t0)/2 sum_{k in interval jint} w_k f_k
                int kf = gauss_bp_idx[jint] + 1;
                int kl = gauss_bp_idx[jint] + gauss_n[jint];
                for (k=kf; k<=kl; k++) {
                    adouble tk = convert_to_original_time_ad( (workspace->snodes[i])(k), t0, tf );
                    get_states(states, xad, iphase, k, workspace);
                    get_controls(controls, xad, iphase, k, workspace);
                    problem->dae(derivatives, path, states, controls, parameters, tk, xad, iphase, workspace);
                    double wk = (workspace->w[i])(k);
                    for (j=0; j<nstates; j++) gad[cbase+j] -= (tf-t0)/2.0 * wk * derivatives[j];
                }
            }
        }

        // LGL hp: K-1 interface defects at the interior breakpoints. Each interior
        // breakpoint is collocated from both sides; the LEFT interval's terminal-row defect
        // is the primary one (in the square top-left block of D, applied in the per-node
        // loop above), and the RIGHT interval's initial-node derivative is enforced here:
        //   D.col(M+e) . X  -  (tf-t0)/2 f(x(bp_e)) = 0,   bp_e = c_{e+1} = sum_{m<=e} n_m.
        // D.col(M+e) is the interface coefficient vector stored by lgl_nodes_multi (nonzero
        // only over interval (e+1)'s nodes). Placed in the same constraint slot that the
        // Radau pin / Gauss quadrature defining constraints occupy. K=1 => no interface rows.
        if ( ( algorithm->collocation_method == "Legendre" || algorithm->collocation_method == "Chebyshev" ) && hp_mesh_active(problem->phase[i]) ) {
            int M  = norder + 1;
            int Kl = hp_num_intervals(problem->phase[i]);
            int iface_base = phase_offset + nstates*(norder+1) + nevents + npath*(norder+1);
            int cbp = 0;
            for (int e = 0; e < Kl-1; e++) {
                cbp += hp_interval_order(problem->phase[i], e);     // bp_e = c_{e+1}
                int cbase = iface_base + e*nstates;
                adouble tk = convert_to_original_time_ad( (workspace->snodes[i])(cbp), t0, tf );
                get_states(states, xad, iphase, cbp, workspace);
                get_controls(controls, xad, iphase, cbp, workspace);
                problem->dae(derivatives, path, states, controls, parameters, tk, xad, iphase, workspace);
                for (j=0; j<nstates; j++) {
                    adouble acc = 0.0;
                    for (int nd=0; nd<M; nd++) {
                        double dcoef = D(nd, M + e);
                        if (dcoef != 0.0) acc += dcoef * states_traj[nd*nstates + j];
                    }
                    gad[cbase+j] = acc - (tf-t0)/2.0 * derivatives[j];
                    if ( algorithm->scaling=="user" ) gad[cbase+j] *= deriv_scaling(j);
                }
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

  // Robust-DAIR optimality step (Option B): append the residual-box constraint block
  //   |(xdot - f)_{k,q,j}| <= ir_residual_bound
  // for every raw residual component, immediately after the linkages. Bounding the residual
  // components (each with an O(1) gradient that does not vanish near feasibility) keeps the
  // trajectory determined and well-conditioned, unlike the single scalar integral. Defects
  // remain dropped (non-binding); the objective stays the pure user cost J. Rows left
  // unscaled so the bound applies to the raw residual.
  if ( algorithm->transcription_method == "integrated-residual"
       && ( algorithm->ir_dair
            || ( algorithm->ir_objective == "cost"
                 && algorithm->ir_residual_bound >= 0.0 ) ) ) {
      adouble t0r, tfr;
      int rb = phase_offset + problem->nlinkages;   // start of the residual-box block
      int m  = workspace->ir_m;
      for (int ip=0; ip<problem->nphases; ip++) {
          int iphr   = ip+1;
          int iphpar = ( problem->multi_segment_flag || workspace->auto_linked_flag ) ? 1 : iphr;
          adouble* params = workspace->parameters[iphpar-1].get();
          get_parameters(params, xad, iphr, workspace);
          get_times(&t0r, &tfr, xad, iphr, workspace);
          int cnt = ir_box_rows( problem->phase[ip].current_number_of_intervals,
                                 problem->phase[ip].nstates, m,
                                 workspace->algorithm->ir_local_order );
          integrated_residual_phase(ip, iphr, xad, t0r, tfr, params, workspace, &gad[rb]);
          for (int t=0; t<cnt; t++) constraint_scaling(rb+t) = 1.0;
          rb += cnt;
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

