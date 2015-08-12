/*********************************************************************************************

This file is part of the PSOPT library, a software tool for computational optimal control

Copyright (C) 2009-2015 Victor M. Becerra

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
           University of Reading
           School of Systems Engineering
           P.O. Box 225, Reading RG6 6AY
           United Kingdom
           e-mail: vmbecerra99@gmail.com

**********************************************************************************************/


#include "psopt.h"


DMatrix& Sol::get_states_in_phase(int iphase)
{
  //   if (iphase <1 || iphase > workspace->problem->nphases)
  //        error_message("incorrect phase index in Prob::phases()");
     return states[iphase-1];
}

DMatrix& Sol::get_parameters_in_phase(int iphase)
{
  //   if (iphase <1 || iphase > workspace->problem->nphases)
  //        error_message("incorrect phase index in Prob::phases()");
     return parameters[iphase-1];
}



DMatrix& Sol::get_controls_in_phase(int iphase)
{
     if (iphase <1 || iphase > problem->nphases) {
          error_message("incorrect phase index in Sol::get_controls_in_phase()");
     }
     return controls[iphase-1];
}

DMatrix& Sol::get_time_in_phase(int iphase)
{
     if (iphase <1 || iphase > problem->nphases)
          error_message("incorrect phase index in Sol::get_time_in_phase()");
     return nodes[iphase-1];
}

DMatrix& Sol::get_dual_costates_in_phase(int iphase)
{
     if (iphase <1 || iphase > problem->nphases)
          error_message("incorrect phase index in Sol::get_dual_costates_in_phase()");
     return dual.costates[iphase-1];
}

DMatrix& Sol::get_dual_hamiltonian_in_phase(int iphase)
{
     if (iphase <1 || iphase > problem->nphases)
          error_message("incorrect phase index in Sol::get_dual_hamiltonian_in_phase()");
     return dual.Hamiltonian[iphase-1];
}

DMatrix& Sol::get_dual_path_in_phase(int iphase)
{
     if (iphase <1 || iphase > problem->nphases)
          error_message("incorrect phase index in Sol::get_dual_path_in_phase()");
     return dual.path[iphase-1];
}

DMatrix& Sol::get_dual_events_in_phase(int iphase)
{
     if (iphase <1 || iphase > problem->nphases)
          error_message("incorrect phase index in Sol::get_dual_events_in_phase()");
     return dual.events[iphase-1];
}

DMatrix& Sol::get_dual_linkages()
{
     return *dual.linkages;
}

DMatrix& Sol::get_relative_local_error_in_phase(int iphase)
{
     if (iphase <1 || iphase > problem->nphases)
          error_message("incorrect phase index in Prob::phases()");
     return relative_errors[iphase-1];
}

void initialize_solution(Sol& solution, Prob& problem, Alg& algorithm, Workspace* workspace)
{
   int nphases = problem.nphases;
   int nparam;
   int i;

   solution.states      = new DMatrix[nphases];
   solution.controls    = new DMatrix[nphases];
   solution.nodes       = new DMatrix[nphases];
   solution.integrand_cost= new DMatrix[nphases];
   solution.parameters  = new DMatrix[nphases];
   solution.relative_errors     = new DMatrix[nphases];

   solution.dual.costates = new DMatrix[nphases];
   solution.dual.path     = new DMatrix[nphases];
   solution.dual.events   = new DMatrix[nphases];
   solution.dual.Hamiltonian = new DMatrix[nphases];
   solution.dual.linkages    = new DMatrix;
   solution.endpoint_cost    = new double[nphases];
   solution.integrated_cost  = new double[nphases];
   solution.problem     = &problem;


   for (i=0;i<nphases; i++)
   {
      nparam = problem.phase[i].nparameters;
      solution.parameters[i].Resize(nparam,1);
   }

   solution.error_flag = false;
   solution.error_msg = "";

   solution.mesh_stats = new MeshStats[ get_number_of_mesh_refinement_iterations(problem,algorithm)];
   for (i=0;i<get_number_of_mesh_refinement_iterations(problem,algorithm) ; i++)
   {
      solution.mesh_stats[i].n_con_evals = 0;
      solution.mesh_stats[i].n_obj_evals = 0;
      solution.mesh_stats[i].n_ode_rhs_evals = 0;
      solution.mesh_stats[i].n_jacobian_evals = 0;
      solution.mesh_stats[i].n_hessian_evals = 0;
   }

   return;
}

void resize_solution(Sol& solution, Prob& problem, Alg& algorithm)
{

  int i;

  for(i=0; i<problem.nphases; i++)
  {
        int nstates       = problem.phase[i].nstates;
        int ncontrols     = problem.phase[i].ncontrols;
        int current_number_of_intervals = problem.phase[i].current_number_of_intervals;
        int npath         = problem.phase[i].npath;

  	(solution.states[i]).Resize( nstates, current_number_of_intervals+1);
   	(solution.controls[i]).Resize(ncontrols, current_number_of_intervals+1);
   	(solution.nodes[i]).Resize(1, current_number_of_intervals+1);
   	(solution.integrand_cost[i]).Resize(1, current_number_of_intervals+1);
   	(solution.dual.costates[i]).Resize(nstates, current_number_of_intervals+1);
   	(solution.dual.Hamiltonian[i]).Resize(1, current_number_of_intervals+1);
	(solution.relative_errors[i]).Resize(1, current_number_of_intervals);
   	if (npath) {
     		(solution.dual.path[i]).Resize(npath, current_number_of_intervals+1);
   	}
  }
  return;
}


