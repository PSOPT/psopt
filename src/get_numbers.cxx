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

int get_number_of_mesh_refinement_iterations(Prob& problem, Alg& algorithm)
{
  int number_of_mesh_refinement_iterations;

  if (algorithm.mesh_refinement == "manual") {	  	
        int n1 = problem.phases(1).nodes.cols();
        for (int iphase=2;iphase<=problem.nphases;iphase++) {
             if (problem.phases(iphase).nodes.cols() != n1) {
                error_message("Number of manual mesh refinement iterations must be the same in all phases");             
             }       
        }
  	
        number_of_mesh_refinement_iterations = n1;
  }
  else
  {
        number_of_mesh_refinement_iterations = algorithm.mr_max_iterations;
  }

  return number_of_mesh_refinement_iterations;
}


int get_number_nlp_vars(Prob& problem, Workspace* workspace)
{
   int i;
   int nlp_vars = 0;


   for (i=0;i< problem.nphases;i++)
   {
   	int    nstates   = problem.phase[i].nstates;
   	int    ncontrols = problem.phase[i].ncontrols;
        int    nparam    = problem.phase[i].nparameters;
   	int    nodes     = problem.phase[i].current_number_of_intervals;
      	nlp_vars += (ncontrols+nstates)*(nodes+1)+nparam+2;
        if ( need_midpoint_controls(*workspace->algorithm, workspace) ) {
            nlp_vars += (ncontrols)*(nodes);
        }
        if ( workspace->algorithm->collocation_method == "Gauss" ) {
            nlp_vars += nstates;   // appended terminal-state variable (before t0,tf)
        }
   }

   return nlp_vars;
}

int get_max_number_nlp_vars(Prob& problem, Alg& algorithm)
{
   int i;
   int nlp_vars = 0;


   for (i=0;i< problem.nphases;i++)
   {
   	int    nstates   = problem.phase[i].nstates;
   	int    ncontrols = problem.phase[i].ncontrols;
        int    nparam    = problem.phase[i].nparameters;
   	int    max_nodes = get_max_nodes( problem, i+1, &algorithm );
      	nlp_vars += (ncontrols+nstates)*(max_nodes+1)+nparam+2;
        nlp_vars += (ncontrols)*(max_nodes);
        if ( algorithm.collocation_method == "Gauss" ) nlp_vars += nstates;

   }

   return nlp_vars;
}


int get_number_nlp_constraints(Prob& problem,Workspace* workspace)
{

    int i;
    int nlp_ncons = 0;

    for(i=0; i<problem.nphases; i++)
    {

       nlp_ncons  += get_ncons_phase_i(problem,i, workspace);

    }

   nlp_ncons +=  problem.nlinkages;

   // Robust-DAIR (Option B): a box constraint |(xdot-f)_{g,q,j}| <= ir_residual_bound on every
   // raw residual component (group g, GL point q, state j) of every phase. The number of groups
   // per phase is norder (cubic-Hermite) or M=norder/ir_local_order (Nie-Kerrigan); see ir_box_rows.
   if ( workspace->algorithm->transcription_method == "integrated-residual"
        && ( workspace->algorithm->ir_dair
             || ( workspace->algorithm->ir_objective == "cost"
                  && workspace->algorithm->ir_residual_bound >= 0.0 ) ) ) {
       int m = workspace->algorithm->ir_residual_nodes;
       for (i=0; i<problem.nphases; i++)
          nlp_ncons += ir_box_rows( problem.phase[i].current_number_of_intervals,
                                    problem.phase[i].nstates, m,
                                    workspace->algorithm->ir_local_order );
   }

   return nlp_ncons;
}

int get_max_number_nlp_constraints(Prob& problem, Alg& algorithm)
{

    int i;
    int nlp_ncons = 0;

    for(i=0; i<problem.nphases; i++)
    {

       int nstates   = problem.phase[i].nstates;
       int nevents   = problem.phase[i].nevents;
       int npath     = problem.phase[i].npath;
       int max_nodes = get_max_nodes(problem, i+1, &algorithm);


       nlp_ncons  += nstates*(max_nodes+1)+ (nevents) + npath*(max_nodes+1) + 1;


       nlp_ncons += npath*(max_nodes);

       if ( algorithm.collocation_method == "Radau" ) nlp_ncons += problem.phase[i].ncontrols;
       if ( algorithm.collocation_method == "Gauss" ) {
           int Kg;
           if ( hp_auto_active(algorithm) ) {
               // The ph driver grows the interval count across refinement, so the workspace
               // must be sized for the worst-case K rather than the seed mesh. With storage
               // (Nc + K - 1) bounded by the node ceiling R = max_nodes and every interval
               // carrying at least one collocation node, K <= (R+1)/2; size for that bound.
               Kg = (max_nodes + 1) / 2 + 1;
           } else {
               Kg = hp_mesh_active(problem.phase[i]) ? (int) problem.phase[i].hp_orders.size() : 1;
           }
           nlp_ncons += Kg * problem.phase[i].nstates;   // K Gauss-quadrature defining constraints (one per interval; K=1 single-block)
       }
       if ( ( algorithm.collocation_method == "Legendre" || algorithm.collocation_method == "Chebyshev" ) && ( hp_auto_active(algorithm) || hp_mesh_active(problem.phase[i]) ) ) {
           int Kl;
           if ( hp_auto_active(algorithm) ) {
               // ph driver grows K across refinement: size for worst-case. LGL shares
               // collocated breakpoints (storage M = N_eff+1 <= max_nodes) and each interval
               // carries order >= 2, so K <= N_eff/2; (max_nodes+1)/2+1 bounds it with slack.
               Kl = (max_nodes + 1) / 2 + 1;
           } else {
               Kl = (int) problem.phase[i].hp_orders.size();
           }
           nlp_ncons += (Kl - 1) * problem.phase[i].nstates;   // K-1 LGL interface defects
       }



   }

   nlp_ncons +=  problem.nlinkages;

   // Robust-DAIR (Option B): reserve the residual-box block (upper bound via max nodes).
   if ( algorithm.transcription_method == "integrated-residual"
        && ( algorithm.ir_dair
             || ( algorithm.ir_objective == "cost"
                  && algorithm.ir_residual_bound >= 0.0 ) ) ) {
       int m = algorithm.ir_residual_nodes;
       for (i=0; i<problem.nphases; i++)
          nlp_ncons += ir_box_rows( get_max_nodes(problem, i+1, &algorithm),
                                    problem.phase[i].nstates, m, algorithm.ir_local_order );
   }

   return nlp_ncons;
}


int get_max_nodes(Prob& problem,int iphase, Alg* algorithm)
{

    int i;

    // hp-adaptive automatic driver (Route B, Liu-Hager-Rao ph): the mesh grows across
    // refinement iterations, so the workspace must be sized for an a-priori ceiling on
    // N_eff, not the (as-yet-unseeded) current mesh. Checked before hp_mesh_active so the
    // ceiling wins even if the user supplied an initial hp mesh that the driver will grow.
    if ( hp_auto_active(*algorithm) )
        return hp_node_ceiling(problem, *algorithm, iphase);

    // hp-adaptive fixed mesh (Route B): when a phase carries an explicit
    // multi-interval Radau mesh (hp_orders populated), its effective single-block
    // order is N_eff = sum(hp_orders). The Workspace allocation and every
    // max-count function key off this, so report it directly. Phases without an
    // hp mesh fall through to the legacy schedule below (bit-identical).
    if ( hp_mesh_active(problem.phase[iphase-1]) ) {
        int Nc = problem.phase[iphase-1].hp_orders.sum();
        // Gauss collocates strictly interior points, so every interface breakpoint is a
        // non-collocated stored node: storage = Nc + K and norder = Nc + K - 1. Radau (and
        // the other endpoint-collocating schemes) share the breakpoint, so N_eff = Nc.
        if ( algorithm->collocation_method == "Gauss" )
            return Nc + (int) problem.phase[iphase-1].hp_orders.size() - 1;
        return Nc;
    }

    // Default to the manual-mode value so that retval is always initialized,
    // even if mesh_refinement holds an unrecognized value. Returning an
    // uninitialized retval here is dangerous: it propagates into
    // get_max_nodes_in_all_phases() and sizes the entire Workspace allocation.
    long length_nodes = problem.phase[iphase-1].nodes.size(); // EIGEN_UPDATE
    int retval = problem.phase[iphase-1].nodes(length_nodes-1);

    if (algorithm->mesh_refinement == "manual") {
         // retval already holds the manual-mode value computed above.
    }

    else if (algorithm->mesh_refinement == "automatic" && !use_local_collocation(*algorithm) ) {

         retval = problem.phase[iphase-1].nodes(0) + (algorithm->mr_min_extrapolation_points-1)*(algorithm->mr_initial_increment);
         int count = retval;
         for (i=1;i<=(algorithm->mr_max_iterations-2);i++) {
//                int increment = algorithm->mr_max_increment_factor*count;
// 		          count += increment;
// The above two lines have been deleted has they led to inconsistent mesh refinement iterations with different maximum number of mesh refinement iterations.
// Thanks to Emmanuel Schneider for pointing out this issue.             
            count += (int) algorithm->mr_max_increment_factor * count; 
	      }
	      retval += count;
    }

    else if (algorithm->mesh_refinement == "automatic" && use_local_collocation(*algorithm) ) {
         int M = problem.phase[iphase-1].nodes(0);
	      int mcount = M;
	      for (i=1; i<= algorithm->mr_max_iterations;i++) {
	         mcount += (int) mcount*algorithm->mr_max_increment_factor;
	      }
	      retval = mcount;
    }

    return retval;
}



int get_max_nodes_in_all_phases(Prob& problem, Alg& algorithm)
{
    int iphase;

    int max_nodes = 0;

    for (iphase=1;iphase<=problem.nphases;iphase++) {
        if ( get_max_nodes(problem,iphase,&algorithm) > max_nodes )
	     max_nodes = get_max_nodes(problem,iphase, &algorithm);
    }

    return (max_nodes);
}




int get_nvars_phase_i(Prob& problem, int i, Workspace* workspace)
{
	int norder    = problem.phase[i].current_number_of_intervals;
	int ncontrols = problem.phase[i].ncontrols;
	int nstates   = problem.phase[i].nstates;
        int nparam    = problem.phase[i].nparameters;


	int nvars_phase_i = (nstates+ncontrols)*(norder+1)+nparam;

        if ( need_midpoint_controls(*workspace->algorithm, workspace) ) {
                    nvars_phase_i += ncontrols*norder;
        }

        if ( workspace->algorithm->collocation_method == "Gauss" ) {
                    nvars_phase_i += nstates;   // appended terminal-state variable (before t0,tf)
        }

        nvars_phase_i += 2;

        return nvars_phase_i;

}

int get_ncons_phase_i(Prob& problem, int i, Workspace* workspace)
{
	int norder    = problem.phase[i].current_number_of_intervals;
	int nstates   = problem.phase[i].nstates;
        int nevents   = problem.phase[i].nevents;
        int npath     = problem.phase[i].npath;

        int ncons_phase_i = nstates*(norder+1) + nevents + npath*(norder+1)+1;

        if ( need_midpoint_controls(*workspace->algorithm, workspace) ) {
                    ncons_phase_i += npath*norder;
        }

        if ( workspace->algorithm->collocation_method == "Radau" ) {
                    ncons_phase_i += problem.phase[i].ncontrols;  // terminal-control interpolation pin
        }

        if ( workspace->algorithm->collocation_method == "Gauss" ) {
                    int Kg = hp_mesh_active(problem.phase[i]) ? (int) problem.phase[i].hp_orders.size() : 1;
                    ncons_phase_i += Kg * problem.phase[i].nstates;    // K Gauss-quadrature defining constraints (one per interval; K=1 single-block)
        }

        if ( ( workspace->algorithm->collocation_method == "Legendre" || workspace->algorithm->collocation_method == "Chebyshev" ) && hp_mesh_active(problem.phase[i]) ) {
                    int Kl = (int) problem.phase[i].hp_orders.size();
                    ncons_phase_i += (Kl - 1) * problem.phase[i].nstates;  // K-1 LGL interface defects (interior breakpoints collocated from both sides)
        }

        return ncons_phase_i;

}


int get_iphase_offset(Prob& problem, int iphase, Workspace* workspace)
{
        int ii;

        int iphase_offset=0;

        for(ii=0;ii< iphase-1;ii++) {

		int nvars_phase_i = get_nvars_phase_i(problem,ii, workspace);

        	iphase_offset +=nvars_phase_i;

       }

       return (iphase_offset);

}



int get_number_of_controls(Prob& problem, int iphase)
{

	return problem.phase[iphase-1].ncontrols;
}

int get_number_of_states(Prob& problem, int iphase)
{
	return problem.phase[iphase-1].nstates;
}

int get_number_of_events(Prob& problem, int iphase)
{
	return problem.phase[iphase-1].nevents;
}

int get_number_of_path_constraints(Prob& problem, int iphase)
{
	return  problem.phase[iphase-1].npath;
}

int get_number_of_nodes(Prob& problem, int iphase)
{
	return problem.phase[iphase-1].current_number_of_intervals;
}



int get_number_of_parameters(Prob& problem, int iphase)
{
	return problem.phase[iphase-1].nparameters;
}

int get_number_of_linkages(Prob& problem)
{
	return problem.nlinkages;
}

int get_number_of_phases(Prob& problem)
{
	return problem.nphases;
}


int get_total_number_of_parameters(Prob& problem)
{
	 int pcount=0;
	 for (int iphase=1;iphase<=problem.nphases;iphase++) {
        pcount+=get_number_of_parameters(problem,iphase);
    }
	 return pcount;
}
