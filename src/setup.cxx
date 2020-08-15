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

void psopt_level1_setup(Prob& problem)
{
   int i;

   int nphases = problem.nphases;
   problem.phase   = new Phases[nphases];
   for(i=0;i<nphases;i++) {
       problem.phase[i].nparameters = 0;
       problem.phase[i].nobserved   = 0;
       problem.phase[i].nsamples    = 0;
       problem.phase[i].zero_cost_integrand = false;
       problem.phase[i].regularization_factor = 0.0;
       problem.phase[i].nodes.resize(1);
   }



}


void psopt_level2_setup(Prob& problem, Alg& algorithm)
{

  int i;
  int nlinkages = problem.nlinkages;

  for(i=0; i< problem.nphases; i++)
  {
        int nstates   = problem.phase[i].nstates;
        int ncontrols = problem.phase[i].ncontrols;
        int nevents   = problem.phase[i].nevents;
        int npath     = problem.phase[i].npath;
        int nparam    = problem.phase[i].nparameters;
        int nnodes    = problem.phase[i].nodes(0);
        int nobserved = problem.phase[i].nobserved;
        int nsamples  = problem.phase[i].nsamples;

   problem.phase[i].observation_nodes.resize(nobserved,nnodes);
   problem.phase[i].observations.resize(nobserved,nsamples);

	problem.phase[i].scale.controls.resize(ncontrols,1);
	problem.phase[i].scale.states.resize(nstates,1);
	problem.phase[i].scale.defects.resize(nstates,1);
	problem.phase[i].scale.events.resize(nevents,1);
	problem.phase[i].scale.path.resize(npath,1);
   problem.phase[i].scale.parameters.resize(nparam,1);

	problem.phase[i].bounds.upper.states.resize(nstates,1);
	problem.phase[i].bounds.lower.states.resize(nstates,1);
	problem.phase[i].bounds.lower.controls.resize(ncontrols,1);
	problem.phase[i].bounds.upper.controls.resize(ncontrols,1);
	problem.phase[i].bounds.upper.events.resize(nevents,1);
	problem.phase[i].bounds.lower.events.resize(nevents,1);
	problem.phase[i].bounds.upper.parameters.resize(nparam,1);
	problem.phase[i].bounds.lower.parameters.resize(nparam,1);
	problem.phase[i].bounds.lower.path.resize(npath,1);
	problem.phase[i].bounds.upper.path.resize(npath,1);



   problem.scale.linkages.resize(nlinkages,1);

   problem.phase[i].guess.controls.resize(0,0);
   problem.phase[i].guess.states.resize(0,0);
   problem.phase[i].guess.time.resize(0,0);
   problem.phase[i].guess.parameters.resize(nparam,1);

   problem.phase[i].name.states_ = new string[nstates];
   problem.phase[i].name.controls_ = new string[ncontrols];
   problem.phase[i].name.parameters_ = new string[nparam];
   problem.phase[i].units.states_  = new string[nstates];
   problem.phase[i].units.controls_ = new string[ncontrols];
   problem.phase[i].units.parameters_ = new string[nparam];



  }

  //problem.bounds.lower.times.resize(problem.nphases+1,1);
  //problem.bounds.upper.times.resize(problem.nphases+1,1);

  problem.bounds.lower.linkage.resize(nlinkages,1);
  problem.bounds.upper.linkage.resize(nlinkages,1);

  problem.bounds.lower.linkage.setZero();
  problem.bounds.lower.linkage.setZero();


  problem.endpoint_cost               = NULL;
  problem.integrand_cost              = NULL;
  problem.dae                         = NULL;
  problem.events                      = NULL;
  problem.linkages                    = NULL;
  problem.observation_function        = NULL;


  // Set default values for some parameters

  algorithm.nlp_method                  = "IPOPT";
  algorithm.scaling                     = "automatic";
  algorithm.defect_scaling              = "state-based";
  algorithm.derivatives                 = "automatic";
  algorithm.constraint_scaling          = "automatic";
  algorithm.nlp_iter_max                = 1000;
  algorithm.nlp_tolerance               = 1.e-6;
  algorithm.jac_sparsity_ratio  	       = 0.5;
  algorithm.hess_sparsity_ratio 	       = 0.2;
  algorithm.hessian                     = "limited-memory";
  algorithm.collocation_method          = "Legendre";
  algorithm.diff_matrix                 = "standard";
  algorithm.ipopt_linear_solver         = "mumps";
  algorithm.print_level                 = 1;
  algorithm.save_sparsity_pattern       = 0;
  algorithm.nsteps_error_integration    = 10;
  algorithm.ode_tolerance               = 1.e-3;
  algorithm.mr_max_increment_factor     = 0.4;
  algorithm.mr_max_iterations		       = 7;
  algorithm.mr_min_extrapolation_points = 2;
  algorithm.mr_initial_increment        = 5;
  algorithm.mr_kappa                    = 0.1;
  algorithm.mr_M1                       = 5;
  algorithm.mesh_refinement 		       = "manual";
  algorithm.switch_order                = 2;
  algorithm.parameter_statistics        = "yes";
  algorithm.parameter_estimation_norm   = 2;
  algorithm.ipopt_max_cpu_time          = 3600.0;
  problem.multi_segment_flag            = false;


  return;
}

