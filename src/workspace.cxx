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



void initialize_workspace_vars(Prob& problem, Alg& algorithm, Sol& solution, Workspace* workspace)
{

  int nphases = problem.nphases;
  int i;

  int max_nvars = get_max_number_nlp_vars(problem, algorithm);
  int max_ncons = get_max_number_nlp_constraints(problem, algorithm);

  int max_nodes = get_max_nodes_in_all_phases(problem, algorithm);

  workspace->P         = make_unique<MatrixXd[]>(nphases);
  workspace->sindex    = make_unique<RowVectorXi[]>(nphases);
  workspace->w         = make_unique<MatrixXd[]>(nphases);
  workspace->D         = make_unique<MatrixXd[]>(nphases);
  workspace->snodes    = make_unique<MatrixXd[]>(nphases);
  workspace->old_snodes= make_unique<MatrixXd[]>(nphases);
  workspace->xlb       = make_unique<MatrixXd>();
  workspace->xub       = make_unique<MatrixXd>();
  workspace->x0        = make_unique<MatrixXd>();
  workspace->lambda    = make_unique<MatrixXd>();
  workspace->dual_costates = make_unique<MatrixXd[]>(nphases);
  workspace->dual_events   = make_unique<MatrixXd[]>(nphases);
  workspace->dual_path     = make_unique<MatrixXd[]>(nphases);
  workspace->linkage    = make_unique<MatrixXd>();
  workspace->constraint_scaling = make_unique<MatrixXd>();
  workspace->prev_states  = make_unique<MatrixXd[]>(nphases);
  workspace->prev_costates= make_unique<MatrixXd[]>(nphases);
  workspace->prev_controls= make_unique<MatrixXd[]>(nphases);
  workspace->prev_path    = make_unique<MatrixXd[]>(nphases);
  workspace->prev_param   = make_unique<MatrixXd[]>(nphases);
  workspace->prev_nodes   = make_unique<MatrixXd[]>(nphases);
  workspace->Ax           = make_unique<TripletSparseMatrix>();
  workspace->Gsp          = make_unique<TripletSparseMatrix>();
  workspace->Xsnopt       = make_unique<MatrixXd>();
  workspace->gsnopt       = make_unique<MatrixXd>();
  workspace->Xip          = make_unique<MatrixXd>();
  workspace->JacRow       = make_unique<MatrixXd>();
  workspace->Gip          = make_unique<MatrixXd>();
  workspace->GFip         = make_unique<MatrixXd>();
  workspace->Xdot         = make_unique<MatrixXd[]>(nphases);
  workspace->DerivResid   = make_unique<MatrixXd[]>(nphases);
  workspace->Xdotgg       = make_unique<MatrixXd[]>(nphases);
  workspace->e            = make_unique<MatrixXd[]>(nphases);
  workspace->hgg          = make_unique<MatrixXd[]>(nphases);
  workspace->prev_t0      = make_unique<MatrixXd>(nphases,1);
  workspace->prev_tf      = make_unique<MatrixXd>(nphases,1);
  workspace->JacCol1      = make_unique<MatrixXd>();
  workspace->JacCol2      = make_unique<MatrixXd>();
  workspace->JacCol3      = make_unique<MatrixXd>();
  workspace->xp           = make_unique<MatrixXd>();
  workspace->emax_history = make_unique<MatrixXd[]>(nphases);
  workspace->order_reduction=make_unique<MatrixXd[]>(nphases);
  workspace->old_relative_errors = make_unique<MatrixXd[]>(nphases);
  workspace->error_scaling_weights = make_unique<MatrixXd[]>(nphases);


  (*workspace->linkage).resize(problem.nlinkages,1);

  workspace->grw = make_unique<GRWORK>();

  workspace->grw->dfdx_j= make_unique<MatrixXd>();
  workspace->grw->F1    = make_unique<MatrixXd>();
  workspace->grw->F2    = make_unique<MatrixXd>();
  workspace->grw->F3    = make_unique<MatrixXd>();
  workspace->grw->F4    = make_unique<MatrixXd>();



  if (algorithm.nlp_method=="IPOPT") {
	workspace->iArow     = make_unique<int[]>((int) (algorithm.jac_sparsity_ratio*max_nvars*max_ncons));
	workspace->jAcol     = make_unique<int[]>((int) (algorithm.jac_sparsity_ratio*max_nvars*max_ncons));
	workspace->iGrow     = make_unique<int[]>((int) (algorithm.jac_sparsity_ratio*max_nvars*max_ncons));
	workspace->jGcol     = make_unique<int[]>((int) (algorithm.jac_sparsity_ratio*max_nvars*max_ncons));
	workspace->jac_Aij   = make_unique<double[]>((int) (algorithm.jac_sparsity_ratio*max_nvars*max_ncons));
	workspace->jac_Gij   = make_unique<double[]>((int) (algorithm.jac_sparsity_ratio*max_nvars*max_ncons));
	if (algorithm.hessian == "exact" || algorithm.hessian == "numerical" ) {
		workspace->hess_nnz_capacity = (int) (algorithm.hess_sparsity_ratio*max_nvars*max_nvars);
		workspace->hess_ir   = new unsigned int[workspace->hess_nnz_capacity];
		workspace->hess_jc   = new unsigned int[workspace->hess_nnz_capacity];
		workspace->lambda_d  = make_unique<double[]>(max_ncons);
		if (algorithm.hessian == "numerical")
			workspace->hess_col_group = make_unique<int[]>(max_nvars);
	}
	else{
      workspace->hess_ir   = NULL;
		workspace->hess_jc   = NULL;
		workspace->hess_nnz_capacity = 0;
		workspace->lambda_d  = NULL;
	}
  }
  else {
    workspace->iArow     = NULL;
	 workspace->jAcol     = NULL;
	 workspace->iGrow     = NULL;
	 workspace->jGcol     = NULL;
	 workspace->jac_Aij   = NULL;
	 workspace->jac_Gij   = NULL;
    workspace->hess_ir   = NULL;
    workspace->hess_jc   = NULL;
    workspace->hess_nnz_capacity = 0;
    workspace->lambda_d  = NULL;
  }
  workspace->hess_verify_done = false;
  workspace->hess_maps_built  = false;
  
  if ( algorithm.nlp_method == "SNOPT") {
  	workspace->iGfun     = new unsigned int[(int) (algorithm.jac_sparsity_ratio*max_nvars*(max_ncons+1))];
  	workspace->jGvar     = new unsigned int[(int) (algorithm.jac_sparsity_ratio*max_nvars*(max_ncons+1))];
  	workspace->iGfun1    = make_unique<int[]>((int) (algorithm.jac_sparsity_ratio*max_nvars*(max_ncons+1)));
  	workspace->jGvar1    = make_unique<int[]>((int) (algorithm.jac_sparsity_ratio*max_nvars*(max_ncons+1)));
  }
  else {
  	workspace->iGfun     = NULL;
  	workspace->jGvar     = NULL;
  	workspace->iGfun1    = NULL;
  	workspace->jGvar1    = NULL;
  }
  workspace->xad       = make_unique<adouble[]>(max_nvars);
  workspace->gad       = make_unique<adouble[]>(max_ncons);
  workspace->fgad      = make_unique<adouble[]>(max_ncons+1);
  workspace->fg        = make_unique<double[]>(max_ncons+1);
  workspace->nrm_row   = make_unique<double[]>(max_ncons+1);

  workspace->states    = make_unique<unique_ptr<adouble[]>[]>(nphases);
  workspace->controls  = make_unique<unique_ptr<adouble[]>[]>(nphases);
  workspace->parameters= make_unique<unique_ptr<adouble[]>[]>(nphases);
  workspace->resid     = make_unique<unique_ptr<adouble[]>[]>(nphases);
  workspace->derivatives     = make_unique<unique_ptr<adouble[]>[]>(nphases);
  workspace->initial_states  = make_unique<unique_ptr<adouble[]>[]>(nphases);
  workspace->final_states    = make_unique<unique_ptr<adouble[]>[]>(nphases);
  workspace->initial_controls= make_unique<unique_ptr<adouble[]>[]>(nphases);
  workspace->final_controls  = make_unique<unique_ptr<adouble[]>[]>(nphases);
  workspace->events          = make_unique<unique_ptr<adouble[]>[]>(nphases);
  workspace->path            = make_unique<unique_ptr<adouble[]>[]>(nphases);
  workspace->states_traj     = make_unique<unique_ptr<adouble[]>[]>(nphases);
  workspace->derivs_traj     = make_unique<unique_ptr<adouble[]>[]>(nphases);
  workspace->linkages        = make_unique<adouble[]>(problem.nlinkages);
  workspace->states_next     = make_unique<unique_ptr<adouble[]>[]>(nphases);
  workspace->controls_next   = make_unique<unique_ptr<adouble[]>[]>(nphases);
  workspace->derivatives_next   = make_unique<unique_ptr<adouble[]>[]>(nphases);
  workspace->path_next          = make_unique<unique_ptr<adouble[]>[]>(nphases);
  workspace->states_bar         = make_unique<unique_ptr<adouble[]>[]>(nphases);
  workspace->controls_bar       = make_unique<unique_ptr<adouble[]>[]>(nphases);
  workspace->derivatives_bar    = make_unique<unique_ptr<adouble[]>[]>(nphases);
  workspace->path_bar           = make_unique<unique_ptr<adouble[]>[]>(nphases);
  workspace->observed_variable  = make_unique<unique_ptr<adouble[]>[]>(nphases);
  workspace->observed_residual  = make_unique<unique_ptr<adouble[]>[]>(nphases);
  workspace->interp_states_pe   = make_unique<unique_ptr<adouble[]>[]>(nphases);
  workspace->interp_controls_pe = make_unique<unique_ptr<adouble[]>[]>(nphases);
  workspace->lam_resid  = make_unique<unique_ptr<adouble[]>[]>(nphases);



  workspace->trace_f_done    = false;

  workspace->time_array_tmp = make_unique<adouble[]>(max_nodes +1);
  workspace->single_trajectory_tmp = make_unique<adouble[]>(max_nodes +1);
  workspace->L_ad_tmp = make_unique<adouble[]>(max_nodes +1);
  workspace->u_spline   = make_unique<adouble[]>(max_nodes +1);
  workspace->z_spline   = make_unique<adouble[]>(max_nodes +1);
  workspace->y2a_spline = make_unique<adouble[]>(max_nodes +1);


 for(i=0; i< problem.nphases; i++)
  {

        int nevents   = problem.phase[i].nevents;
        int npath     = problem.phase[i].npath;
        int nparam    = problem.phase[i].nparameters;
        int nstates   = problem.phase[i].nstates;
        int ncontrols = problem.phase[i].ncontrols;
	     int nobserved = problem.phase[i].nobserved;

        int max_nodes = get_max_nodes(problem,i+1, &algorithm);


        workspace->dual_events[i].resize(nevents,1);

        workspace->e[i].resize(nevents,1);

        if (nparam>=1) {
          workspace->prev_param[i].resize(nparam,1);
        }

        workspace->states[i]= make_unique<adouble[]>(nstates);
        workspace->controls[i] = make_unique<adouble[]>(ncontrols);
        workspace->parameters[i] = make_unique<adouble[]>(nparam);
        workspace->resid[i]= make_unique<adouble[]>(nstates);
        workspace->derivatives[i]= make_unique<adouble[]>(nstates);
        workspace->initial_states[i]= make_unique<adouble[]>(nstates);
        workspace->final_states[i]= make_unique<adouble[]>(nstates);
        workspace->initial_controls[i]= make_unique<adouble[]>(ncontrols);
        workspace->final_controls[i]= make_unique<adouble[]>(ncontrols);
        workspace->events[i]= make_unique<adouble[]>(nevents);
        workspace->path[i]= make_unique<adouble[]>(npath);

        workspace->states_next[i]     = make_unique<adouble[]>(nstates);
        workspace->controls_next[i]   = make_unique<adouble[]>(ncontrols);
        workspace->derivatives_next[i]= make_unique<adouble[]>(nstates);
        workspace->path_next[i]       = make_unique<adouble[]>(npath);
        workspace->states_bar[i]      = make_unique<adouble[]>(nstates);
        workspace->controls_bar[i]    = make_unique<adouble[]>(ncontrols);
        workspace->derivatives_bar[i] = make_unique<adouble[]>(nstates);

        workspace->path_bar[i]        = make_unique<adouble[]>(npath);

   	  workspace->observed_variable[i] = make_unique<adouble[]>(nobserved);
	     workspace->observed_residual[i] = make_unique<adouble[]>(nobserved);
  	     workspace->lam_resid[i]              = make_unique<adouble[]>(nobserved);

   	  workspace->interp_states_pe[i]   = make_unique<adouble[]>(nstates);
	     workspace->interp_controls_pe[i] = make_unique<adouble[]>(ncontrols);

        workspace->states_traj[i]= make_unique<adouble[]>(problem.phase[i].nstates*(max_nodes +1));
        workspace->derivs_traj[i]= make_unique<adouble[]>(problem.phase[i].nstates*(max_nodes +1));
        workspace->order_reduction[i].resize(1,max_nodes+1);


  }

  int dotindex = problem.outfilename.find_first_of(".");

  workspace->igroup = make_unique<IGroup>();

  string fname = "psopt_solution_" + problem.outfilename.substr(0,dotindex) + ".txt";

  workspace->psopt_solution_summary_file = fopen(fname.c_str(),"w");

  workspace->mesh_statistics = fopen("mesh_statistics.txt","w");

  fname = "mesh_statistics_" + problem.outfilename.substr(0,dotindex) + ".tex";

  workspace->mesh_statistics_tex = fopen( fname.c_str(),"w");

  if (workspace->psopt_solution_summary_file == NULL) error_message("Error opening \"psopt_solution_summary.txt\" file");

  if (workspace->mesh_statistics == NULL) error_message("Error opening \"mesh_statistics.txt\" file");

  if (workspace->mesh_statistics_tex == NULL) error_message("Error opening \"mesh_statistics.tex\" file");

  for(i=0;i<problem.nphases;i++) {
     workspace->emax_history[i].resize(50,2);
  }

  workspace->enable_nlp_counters = false;

  workspace->auto_linked_flag = false;

// Initialise tape tags to be used by ADOL_C

  workspace->ad_f.tag    = 1;
  workspace->ad_g.tag    = 2;
  workspace->ad_hess.tag = 3;
  workspace->ad_fg.tag   = 4;
  workspace->ad_gc.tag   = 5;

  workspace->user_data = problem.user_data;
  
  workspace->nS           = 0;

}

void resize_workspace_vars(Prob& problem, Alg& algorithm, Sol& solution, Workspace* workspace)
{
  int i;

  int nlp_ncons = get_number_nlp_constraints(problem, workspace );
  int nvars     = get_number_nlp_vars(problem, workspace);

  workspace->Xsnopt->resize(nvars, 1);
  workspace->gsnopt->resize(nlp_ncons, 1);

  workspace->Xip->resize(nvars,1);
  workspace->JacRow->resize(1,nvars);
  workspace->Gip->resize(nlp_ncons,1);
  workspace->GFip->resize(nvars,1);

  (*workspace->grw->dfdx_j).resize( nlp_ncons,1 );
  (*workspace->grw->F1).resize( nlp_ncons, 1 );
  (*workspace->grw->F2).resize( nlp_ncons, 1 );
  (*workspace->grw->F3).resize( nlp_ncons, 1 );
  (*workspace->grw->F4).resize( nlp_ncons, 1 );
  if ( algorithm.nlp_method == "SNOPT") {
    workspace->JacCol1->resize(nlp_ncons+1,1);
    workspace->JacCol2->resize(nlp_ncons+1,1);
    workspace->JacCol3->resize(nlp_ncons+1,1);
  } 
  else {
    workspace->JacCol1->resize(nlp_ncons,1);
    workspace->JacCol2->resize(nlp_ncons,1);
    workspace->JacCol3->resize(nlp_ncons,1);  
  }
  workspace->xp->resize(nvars,1);
  workspace->constraint_scaling->resize(nlp_ncons,1);


  solution.xad = workspace->xad.get();

  workspace->x0->resize(nvars,1);
  workspace->lambda->resize(nlp_ncons,1);

  // workspace->xlb->resize(nvars,1);
  // workspace->xub->resize(nvars,1);
 *workspace->xlb = MatrixXd::Zero(nvars,1);
 *workspace->xub = MatrixXd::Zero(nvars,1);

  workspace->nphases = problem.nphases;

  for(i=0; i< problem.nphases; i++)
  {
        int nstates   = problem.phase[i].nstates;
        int npath     = problem.phase[i].npath;
        int norder    = problem.phase[i].current_number_of_intervals;

        workspace->w[i].resize(norder+1,1);

        if ( algorithm.collocation_method == "Legendre") {
	        workspace->P[i].resize(norder+1,norder+1);
	}
	if ( use_global_collocation(algorithm) ) {
                workspace->D[i].resize(norder+1,norder+1);
	}

        workspace->snodes[i].resize(1,norder+1);
        workspace->sindex[i].resize(norder+1);



        workspace->dual_costates[i].resize(nstates,norder+1);
        workspace->dual_path[i].resize(npath,norder+1);
        workspace->DerivResid[i].resize(nstates,norder+1);
        workspace->Xdot[i].resize(nstates,norder+1);
        workspace->DerivResid[i].resize(nstates,norder+1);
        workspace->Xdotgg[i].resize(nstates,norder+1);
        workspace->hgg[i].resize(npath,norder+1);
	     workspace->order_reduction[i].resize(1,norder);

  }

}


work_str::work_str(Prob& problem, Alg& algorithm, Sol& solution)
{
    initialize_workspace_vars(problem, algorithm, solution, this);
}


work_str::~work_str()
{
  // Per-row arrays (states, controls, ... ) are unique_ptr<unique_ptr<adouble[]>[]>
  // and free themselves; no manual per-row delete loop is needed.
  if (this->hess_ir) delete [] this->hess_ir;
  if (this->hess_jc) delete [] this->hess_jc;
  if (this->iGfun) delete [] this->iGfun;
  if (this->jGvar) delete [] this->jGvar;

  // states, controls, ... (the 2-D arrays) are unique_ptr-owned and free themselves.



  // xlb, xub are unique_ptr<MatrixXd> and free themselves.



}

