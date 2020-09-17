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


MatrixXd& Sol::get_states_in_phase(int iphase)
{
  //   if (iphase <1 || iphase > workspace->problem->nphases)
  //        error_message("incorrect phase index in Prob::phases()");
     return states[iphase-1];
}

MatrixXd& Sol::get_parameters_in_phase(int iphase)
{
  //   if (iphase <1 || iphase > workspace->problem->nphases)
  //        error_message("incorrect phase index in Prob::phases()");
     return parameters[iphase-1];
}



MatrixXd& Sol::get_controls_in_phase(int iphase)
{
     if (iphase <1 || iphase > problem->nphases) {
          error_message("incorrect phase index in Sol::get_controls_in_phase()");
     }
     return controls[iphase-1];
}

MatrixXd& Sol::get_time_in_phase(int iphase)
{
     if (iphase <1 || iphase > problem->nphases)
          error_message("incorrect phase index in Sol::get_time_in_phase()");
     return nodes[iphase-1];
}

MatrixXd& Sol::get_dual_costates_in_phase(int iphase)
{
     if (iphase <1 || iphase > problem->nphases)
          error_message("incorrect phase index in Sol::get_dual_costates_in_phase()");
     return dual.costates[iphase-1];
}

MatrixXd& Sol::get_dual_hamiltonian_in_phase(int iphase)
{
     if (iphase <1 || iphase > problem->nphases)
          error_message("incorrect phase index in Sol::get_dual_hamiltonian_in_phase()");
     return dual.Hamiltonian[iphase-1];
}

MatrixXd& Sol::get_dual_path_in_phase(int iphase)
{
     if (iphase <1 || iphase > problem->nphases)
          error_message("incorrect phase index in Sol::get_dual_path_in_phase()");
     return dual.path[iphase-1];
}

MatrixXd& Sol::get_dual_events_in_phase(int iphase)
{
     if (iphase <1 || iphase > problem->nphases)
          error_message("incorrect phase index in Sol::get_dual_events_in_phase()");
     return dual.events[iphase-1];
}

MatrixXd& Sol::get_dual_linkages()
{
     return *dual.linkages;
}

MatrixXd& Sol::get_relative_local_error_in_phase(int iphase)
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

   solution.states      = new MatrixXd[nphases];           
   solution.controls    = new MatrixXd[nphases];           
   solution.nodes       = new MatrixXd[nphases];           
   solution.integrand_cost= new MatrixXd[nphases];         
   solution.parameters  = new MatrixXd[nphases];           
   solution.relative_errors     = new MatrixXd[nphases];   

   solution.dual.costates = new MatrixXd[nphases];         
   solution.dual.path     = new MatrixXd[nphases];         
   solution.dual.events   = new MatrixXd[nphases];         
   solution.dual.Hamiltonian = new MatrixXd[nphases];      
   solution.dual.linkages    = new MatrixXd;               
   solution.endpoint_cost    = new double[nphases];
   solution.integrated_cost  = new double[nphases];
   solution.problem     = &problem;


   for (i=0;i<nphases; i++)
   {
      nparam = problem.phase[i].nparameters;
      solution.parameters[i].resize(nparam,1);
   }

   solution.error_flag = 0;
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

  	(solution.states[i]).resize( nstates, current_number_of_intervals+1);
   	(solution.controls[i]).resize(ncontrols, current_number_of_intervals+1);
   	(solution.nodes[i]).resize(1, current_number_of_intervals+1);
   	(solution.integrand_cost[i]).resize(1, current_number_of_intervals+1);
   	(solution.dual.costates[i]).resize(nstates, current_number_of_intervals+1);
   	(solution.dual.Hamiltonian[i]).resize(1, current_number_of_intervals+1);
	(solution.relative_errors[i]).resize(1, current_number_of_intervals);
   	if (npath) {
     		(solution.dual.path[i]).resize(npath, current_number_of_intervals+1);
   	}
  }
  return;
}


