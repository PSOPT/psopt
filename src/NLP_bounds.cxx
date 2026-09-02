/*********************************************************************************************

TThis file is part of the PSOPT library, a software tool for computational optimal control

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
#include <vector>

// Bring std names into this translation unit (formerly leaked via psopt.h).
using namespace std;



void  define_nlp_bounds(MatrixXd& xlb, MatrixXd& xub, Prob& problem, Alg& algorithm, Workspace* workspace)
{
  // This function defines the NLP bounds given the information provided by
  // the user.

   int i, k;

   int x_phase_offset = 0;

   xlb = zeros(workspace->nvars,1);
   xub = zeros(workspace->nvars,1);


   for (i=0; i< problem.nphases; i++)
   {
   	int norder    = problem.phase[i].current_number_of_intervals;
   	int ncontrols = problem.phase[i].ncontrols;
    int nparam    = problem.phase[i].nparameters;
   	int nstates   = problem.phase[i].nstates;
   	int offset1   = ncontrols*(norder+1);
    int offset2   = (ncontrols+nstates)*(norder+1);

	MatrixXd& control_scaling = problem.phase[i].scale.controls;
	MatrixXd& state_scaling   = problem.phase[i].scale.states;
    MatrixXd& param_scaling   = problem.phase[i].scale.parameters;
    double time_scaling      = problem.phase[i].scale.time;

	int nvars_phase_i = get_nvars_phase_i(problem,i, workspace);


	for (k=0; k<norder+1; k++) {   // EIGEN_UPDATE: k index shifted by -1.
                if (ncontrols>0) {

	              xlb.block(x_phase_offset+k*ncontrols, 0, ncontrols, 1 ) = elemProduct((problem.phase[i].bounds.lower.controls),control_scaling);
                }

		xlb.block(x_phase_offset+(k)*nstates+offset1,0, nstates, 1) = elemProduct((problem.phase[i].bounds.lower.states),state_scaling);

                if (ncontrols>0) {

		     xub.block(x_phase_offset+(k)*ncontrols,0, ncontrols, 1) = elemProduct((problem.phase[i].bounds.upper.controls),control_scaling);

                }

        xub.block(x_phase_offset+(k)*nstates+offset1,0, nstates, 1)=elemProduct((problem.phase[i].bounds.upper.states),state_scaling);

	}

        offset1 = (nstates+ncontrols)*(norder+1);

        if (nparam>=1) {


             xlb.block(x_phase_offset+offset2,0,nparam,1)= elemProduct((problem.phase[i].bounds.lower.parameters),param_scaling);



             xub.block(x_phase_offset+offset2,0,nparam,1)= elemProduct((problem.phase[i].bounds.upper.parameters),param_scaling);

        }

        offset1 = (nstates+ncontrols)*(norder+1) + nparam;

        if ( midpoint_control_vars(*workspace->algorithm, workspace) ) {
	   for (k=0; k<norder; k++) { // EIGEN_UPDATE: K index shifted by -1
                if (ncontrols>0) {

		          xlb.block(x_phase_offset+offset1+(k)*ncontrols,0,ncontrols,1) = elemProduct((problem.phase[i].bounds.lower.controls),control_scaling);


                  xub.block(x_phase_offset+offset1+(k)*ncontrols,0,ncontrols,1) = elemProduct((problem.phase[i].bounds.upper.controls),control_scaling);

                }
	   }
        }

        // Gauss: appended terminal-state variable takes the state bounds.
        if ( algorithm.collocation_method == "Gauss" ) {
            int xf_off = (nstates+ncontrols)*(norder+1) + nparam;
            xlb.block(x_phase_offset+xf_off,0,nstates,1) = elemProduct((problem.phase[i].bounds.lower.states),state_scaling);
            xub.block(x_phase_offset+xf_off,0,nstates,1) = elemProduct((problem.phase[i].bounds.upper.states),state_scaling);
        }

	xlb(x_phase_offset+nvars_phase_i-2)   = problem.phase[i].bounds.lower.StartTime*time_scaling; //EIGEN_UPDATE
	xub(x_phase_offset+nvars_phase_i-2)   = problem.phase[i].bounds.upper.StartTime*time_scaling; //EIGEN_UPDATE

	xlb(x_phase_offset+nvars_phase_i-1) = problem.phase[i].bounds.lower.EndTime*time_scaling;  //EIGEN_UPDATE
	xub(x_phase_offset+nvars_phase_i-1) = problem.phase[i].bounds.upper.EndTime*time_scaling;  //EIGEN_UPDATE

        x_phase_offset += nvars_phase_i;

   }

   // The cost variables carry no bound of their own: what they may be is decided by the
   // equality that ties each to its phase's quadrature.
   {
       const int nm = mayer_extra_vars(problem, algorithm, workspace);
       for (int q = 0; q < nm; q++) {
           xlb(x_phase_offset+q) = -PSOPT::inf;
           xub(x_phase_offset+q) =  PSOPT::inf;
       }
   }

}





void get_constraint_bounds(double* g_l, double* g_u, Workspace* workspace)
{

  Prob *problem = workspace->problem;
  Alg  *algorithm= workspace->algorithm;

  int lam_phase_offset=0;

  MatrixXd& constraint_scaling = *workspace->constraint_scaling;

  int offset, k, l, i, j;

  double path_sc, event_sc;

  for(i=0; i<problem->nphases; i++)
  {
  	int nstates = problem->phase[i].nstates;
  	int norder  = problem->phase[i].current_number_of_intervals;
  	int nevents = problem->phase[i].nevents;
  	int npath   = problem->phase[i].npath;
  	MatrixXd& path_scaling    = problem->phase[i].scale.path;
  	MatrixXd& event_scaling   = problem->phase[i].scale.events;

        int ncons_phase_i = get_ncons_phase_i(*problem,i, workspace);
	// lower and upper bounds on g(x)
	if ( workspace->transcription_method == "integrated-residual" ) {
		// The (zeroed) defect rows carry no dynamics under integrated-residual
		// transcription; make them non-binding range constraints so the active-set
		// Jacobian is not rank-deficient (degenerate 0=0 equalities break IPOPT's
		// multiplier computation). Dynamics are enforced by the IR objective instead.
		for (k=0;k<nstates*(norder+1);k++) {
			g_l[lam_phase_offset+k] = -1.0e20;
			g_u[lam_phase_offset+k] =  1.0e20;
		}
	}
	else {
	for (k=0;k<nstates*(norder+1);k++) {
		g_l[lam_phase_offset+k] = 0.0;
		g_u[lam_phase_offset+k] = 0.0;
	}
	}

	offset = lam_phase_offset+nstates*(norder+1);

	for (k=0;k<nevents;k++) { // EIGEN_UPDATE: index k shifted by -1
		j = offset + k;
		if( algorithm->scaling=="user" )
		        event_sc = event_scaling(k);
 		else
			event_sc = constraint_scaling(j);

		g_l[j] = (problem->phase[i].bounds.lower.events)(k)*event_sc;
		g_u[j] = (problem->phase[i].bounds.upper.events)(k)*event_sc;
	}

	offset = offset + nevents;

	// Path constraints that have been folded into the integrated residual are removed from
	// the pointwise set, following Neuenhofen and Kerrigan, in whose formulation the
	// algebraic equations are enforced through the residual and not additionally as
	// pointwise constraints. The rows are left in place with free bounds rather than
	// physically deleted, so that the constraint layout, the Jacobian structure and the
	// path-multiplier arrays returned to the user keep their declared shapes; the
	// multipliers of the freed rows are zero, which is the correct report for a constraint
	// that is no longer imposed here.
	std::vector<int>    ir_alg_idx;
	std::vector<double> ir_alg_tgt;
	if ( ir_algebraic_rows(*problem, *algorithm, i) > 0 )
	    ir_algebraic_index(*problem, i, ir_alg_idx, ir_alg_tgt);
	std::vector<bool> folded( npath > 0 ? npath : 1, false );
	for (size_t a = 0; a < ir_alg_idx.size(); a++) folded[ ir_alg_idx[a] ] = true;

	for (k=0; k<npath; k++) // EIGEN_UPDATE: index k shifted by -1
	{
		for (l=0;l<(norder + 1);l++) {   // EIGEN_UPDATE: index l shifted by -1.
		    j = offset + (l)*npath + k;
		    if ( folded[k] ) { g_l[j] = -PSOPT::inf; g_u[j] = PSOPT::inf; continue; }
		    if( algorithm->scaling=="user" )
		       path_sc = path_scaling(k);
		    else
		       path_sc = constraint_scaling(j);

		    g_l[j] = (problem->phase[i].bounds.lower.path)(k)*path_sc;
		    g_u[j] = (problem->phase[i].bounds.upper.path)(k)*path_sc;
		}
	}


        if ( need_midpoint_controls(*workspace->algorithm, workspace) ) {
		offset = offset + npath*(norder+1);

		for (k=0; k<npath; k++)  // EIGEN_UPDATE: index k shifted by -1.
		{
			for (l=0;l<norder;l++) { // EIGEN_UPDATE: index l shifted by -1.
		    		j = offset + (l)*npath + k;
		    		if ( folded[k] ) { g_l[j] = -PSOPT::inf; g_u[j] = PSOPT::inf; continue; }
		    		if( algorithm->scaling=="user" )
		       			path_sc = path_scaling(k);
		    		else
		       			path_sc = constraint_scaling(j);

		    		g_l[j] = (problem->phase[i].bounds.lower.path)(k)*path_sc;
		    		g_u[j] = (problem->phase[i].bounds.upper.path)(k)*path_sc;
			}
		}
	}


        // Radau: terminal-control interpolation pin constraints are equalities (=0).
        if ( algorithm->collocation_method == "Radau" ) {
            int ncontrols = problem->phase[i].ncontrols;
            int pin_base  = lam_phase_offset + nstates*(norder+1) + nevents + npath*(norder+1);
            for (int l2=0; l2<ncontrols; l2++) { g_l[pin_base+l2] = 0.0; g_u[pin_base+l2] = 0.0; }
        }

        // Gauss: the K Gauss-quadrature defining constraints (one per interval) are equalities (=0).
        if ( algorithm->collocation_method == "Gauss" ) {
            int Kg = hp_mesh_active(problem->phase[i]) ? hp_num_intervals(problem->phase[i]) : 1;
            int quad_base = lam_phase_offset + nstates*(norder+1) + nevents + npath*(norder+1);
            for (int l2=0; l2<Kg*nstates; l2++) { g_l[quad_base+l2] = 0.0; g_u[quad_base+l2] = 0.0; }
        }

        // Legendre hp: the K-1 LGL interface defects (interior breakpoints) are equalities (=0).
        if ( ( algorithm->collocation_method == "Legendre" || algorithm->collocation_method == "Chebyshev" ) && hp_mesh_active(problem->phase[i]) ) {
            int Kl = hp_num_intervals(problem->phase[i]);
            int iface_base = lam_phase_offset + nstates*(norder+1) + nevents + npath*(norder+1);
            for (int l2=0; l2<(Kl-1)*nstates; l2++) { g_l[iface_base+l2] = 0.0; g_u[iface_base+l2] = 0.0; }
        }

        lam_phase_offset += ncons_phase_i;

        // Bounds for t0 <= tf constraint

        double diff_t0Min_tfMax= (problem->phase[i].bounds.lower.StartTime-problem->phase[i].bounds.upper.EndTime);
        diff_t0Min_tfMax *= problem->phase[i].scale.time;
        if (algorithm->scaling == "automatic")  {
           diff_t0Min_tfMax *= constraint_scaling(lam_phase_offset-1);
        }

        g_l[ lam_phase_offset - 1 ] = diff_t0Min_tfMax;
        g_u[ lam_phase_offset - 1 ] = 0.0;

  }

  for (k=0;k<problem->nlinkages; k++)
  {
        j = lam_phase_offset + k;
        g_l[j] = problem->bounds.lower.linkage(k); // EIGEN_UPDATE
        g_u[j] = problem->bounds.upper.linkage(k); // EIGEN_UPDATE

  }

  // J_i - (quadrature of phase i) = 0, one per converted phase, immediately after the
  // linkages. This block and the robust-DAIR block below are mutually exclusive:
  // mayer_extra_vars declines the integrated-residual transcriptions outright.
  {
      const int nm = mayer_extra_vars(*problem, *algorithm, workspace);
      for (k = 0; k < nm; k++) {
          j = lam_phase_offset + problem->nlinkages + k;
          g_l[j] = 0.0;
          g_u[j] = 0.0;
      }
  }

  // Robust-DAIR (Option B): box bounds on every residual component  |r_{k,q,j}| <= delta.
  // Under ir_dair the per-phase tolerances ir_delta_phase[p] = K*h_p^2 (set by the DAIR
  // optimality sub-solve) override the scalar ir_residual_bound; when ir_delta_phase is not
  // populated (the feasibility sub-solve, and direct-box mode) the scalar is used for all phases.
  if ( algorithm->transcription_method == "integrated-residual"
       && ( algorithm->ir_dair
            || ( algorithm->ir_objective == "cost"
                 && algorithm->ir_residual_bound >= 0.0 ) ) ) {
      bool use_phase = ( (int) workspace->ir_delta_phase.size() == problem->nphases );
      int m  = algorithm->ir_residual_nodes;
      int rb = lam_phase_offset + problem->nlinkages;
      for (i=0; i<problem->nphases; i++) {
          double delta = use_phase ? workspace->ir_delta_phase[i] : algorithm->ir_residual_bound;
          int cnt = ir_box_rows( problem->phase[i].current_number_of_intervals,
                                 problem->phase[i].nstates, m, algorithm->ir_local_order,
                                 ir_algebraic_rows(*problem, *algorithm, i) );
          for (int t=0; t<cnt; t++) { g_l[rb+t] = -delta; g_u[rb+t] = delta; }
          rb += cnt;
      }
  }


  return;
}

