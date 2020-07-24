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


void initialize_workspace_vars(Prob& problem, Alg& algorithm, Sol& solution, Workspace* workspace)
{

  int nphases = problem.nphases;
  int i;

  int max_nvars = get_max_number_nlp_vars(problem, algorithm);
  int max_ncons = get_max_number_nlp_constraints(problem, algorithm);

  int max_nodes = get_max_nodes_in_all_phases(problem, algorithm);

  workspace->P         = new DMatrix[nphases];
  workspace->sindex    = new DMatrix[nphases];
  workspace->w         = new DMatrix[nphases];
  workspace->D         = new DMatrix[nphases];
  workspace->snodes    = new DMatrix[nphases];
  workspace->old_snodes= new DMatrix[nphases];
  workspace->xlb       = new DMatrix;
  workspace->xub       = new DMatrix;
  workspace->x0        = new DMatrix;
  workspace->lambda    = new DMatrix;
  workspace->dual_costates = new DMatrix[nphases];
  workspace->dual_events   = new DMatrix[nphases];
  workspace->dual_path     = new DMatrix[nphases];
  workspace->linkage    = new DMatrix;
  workspace->constraint_scaling = new DMatrix;
  workspace->prev_states  = new DMatrix[nphases];
  workspace->prev_costates= new DMatrix[nphases];
  workspace->prev_controls= new DMatrix[nphases];
  workspace->prev_path    = new DMatrix[nphases];
  workspace->prev_param   = new DMatrix[nphases];
  workspace->prev_nodes   = new DMatrix[nphases];
  workspace->Ax           = new SparseMatrix;
  workspace->Gsp          = new SparseMatrix;
  workspace->Xsnopt       = new DMatrix;
  workspace->gsnopt       = new DMatrix;
  workspace->Xip          = new DMatrix;
  workspace->JacRow       = new DMatrix;
  workspace->Gip          = new DMatrix;
  workspace->GFip         = new DMatrix;
  workspace->Xdot         = new DMatrix[nphases];
  workspace->DerivResid   = new DMatrix[nphases];
  workspace->Xdotgg       = new DMatrix[nphases];
  workspace->e            = new DMatrix[nphases];
  workspace->hgg          = new DMatrix[nphases];
  workspace->prev_t0      = new DMatrix(nphases,1);
  workspace->prev_tf      = new DMatrix(nphases,1);
  workspace->JacCol1      = new DMatrix;
  workspace->JacCol2      = new DMatrix;
  workspace->JacCol3      = new DMatrix;
  workspace->xp           = new DMatrix;
  workspace->emax_history = new DMatrix[nphases];
  workspace->order_reduction=new DMatrix[nphases];
  workspace->old_relative_errors = new DMatrix[nphases];
  workspace->error_scaling_weights = new DMatrix[nphases];


  (*workspace->linkage).Resize(problem.nlinkages,1);

  workspace->grw = new GRWORK;

  workspace->grw->dfdx_j= new DMatrix;
  workspace->grw->F1    = new DMatrix;
  workspace->grw->F2    = new DMatrix;
  workspace->grw->F3    = new DMatrix;
  workspace->grw->F4    = new DMatrix;



  if (algorithm.nlp_method=="IPOPT") {
	workspace->iArow     = new int[(int) (algorithm.jac_sparsity_ratio*max_nvars*max_ncons)];
	workspace->jAcol     = new int[(int) (algorithm.jac_sparsity_ratio*max_nvars*max_ncons)];
	workspace->iGrow     = new int[(int) (algorithm.jac_sparsity_ratio*max_nvars*max_ncons)];
	workspace->jGcol     = new int[(int) (algorithm.jac_sparsity_ratio*max_nvars*max_ncons)];
	workspace->jac_Aij   = new double[(int) (algorithm.jac_sparsity_ratio*max_nvars*max_ncons)];
	workspace->jac_Gij   = new double[(int) (algorithm.jac_sparsity_ratio*max_nvars*max_ncons)];
	if (algorithm.hessian == "exact" ) {
		workspace->hess_ir   = new unsigned int[(int) (algorithm.hess_sparsity_ratio*max_nvars*max_nvars)];
		workspace->hess_jc   = new unsigned int[(int) (algorithm.hess_sparsity_ratio*max_nvars*max_nvars)];
		workspace->lambda_d  = new double [max_ncons];
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
    workspace->lambda_d  = NULL;
  }
  
  if ( algorithm.nlp_method == "SNOPT") {
  	workspace->iGfun     = new unsigned int[(int) (algorithm.jac_sparsity_ratio*max_nvars*(max_ncons+1))];
  	workspace->jGvar     = new unsigned int[(int) (algorithm.jac_sparsity_ratio*max_nvars*(max_ncons+1))];
  	workspace->iGfun1    = new unsigned int[(int) (algorithm.jac_sparsity_ratio*max_nvars*(max_ncons+1))];
  	workspace->jGvar1    = new unsigned int[(int) (algorithm.jac_sparsity_ratio*max_nvars*(max_ncons+1))];
  	workspace->iGfun2    = new unsigned int[(int) (algorithm.jac_sparsity_ratio*max_nvars*(max_ncons+1))];
  	workspace->jGvar2    = new unsigned int[(int) (algorithm.jac_sparsity_ratio*max_nvars*(max_ncons+1))];
  	workspace->G2        = new double[(int) (algorithm.jac_sparsity_ratio*max_nvars*(max_ncons+1))];
  }
  else {
  	workspace->iGfun     = NULL;
  	workspace->jGvar     = NULL;
  	workspace->iGfun1    = NULL;
  	workspace->jGvar1    = NULL;
  	workspace->iGfun2    = NULL;
  	workspace->jGvar2    = NULL;
  	workspace->G2        = NULL;
  }
  workspace->xad       = new adouble[max_nvars];
  workspace->gad       = new adouble[max_ncons];
  workspace->fgad      = new adouble[max_ncons+1];
  workspace->fg        = new double[max_ncons+1];
  workspace->nrm_row   = new double[max_ncons+1];

  workspace->states    = new adouble*[nphases];
  workspace->controls  = new adouble*[nphases];
  workspace->parameters= new adouble*[nphases];
  workspace->resid     = new adouble*[nphases];
  workspace->derivatives     = new adouble*[nphases];
  workspace->initial_states  = new adouble*[nphases];
  workspace->final_states    = new adouble*[nphases];
  workspace->initial_controls= new adouble*[nphases];
  workspace->final_controls  = new adouble*[nphases];
  workspace->events          = new adouble*[nphases];
  workspace->path            = new adouble*[nphases];
  workspace->states_traj     = new adouble*[nphases];
  workspace->derivs_traj     = new adouble*[nphases];
  workspace->linkages        = new adouble[problem.nlinkages];
  workspace->states_next     = new adouble*[nphases];
  workspace->controls_next   = new adouble*[nphases];
  workspace->derivatives_next   = new adouble*[nphases];
  workspace->path_next          = new adouble*[nphases];
  workspace->states_bar         = new adouble*[nphases];
  workspace->controls_bar       = new adouble*[nphases];
  workspace->derivatives_bar    = new adouble*[nphases];
  workspace->path_bar           = new adouble*[nphases];
  workspace->observed_variable  = new adouble*[nphases];
  workspace->observed_residual  = new adouble*[nphases];
  workspace->interp_states_pe   = new adouble*[nphases];
  workspace->interp_controls_pe = new adouble*[nphases];
  workspace->lam_resid  = new adouble*[nphases];



  workspace->trace_f_done    = false;

  workspace->time_array_tmp = new adouble[max_nodes +1];
  workspace->single_trajectory_tmp = new adouble[max_nodes +1];
  workspace->L_ad_tmp = new adouble[max_nodes +1];
  workspace->u_spline   = new adouble[max_nodes +1];
  workspace->z_spline   = new adouble[max_nodes +1];
  workspace->y2a_spline = new adouble[max_nodes +1];


 for(i=0; i< problem.nphases; i++)
  {

        int nevents   = problem.phase[i].nevents;
        int npath     = problem.phase[i].npath;
        int nparam    = problem.phase[i].nparameters;
        int nstates   = problem.phase[i].nstates;
        int ncontrols = problem.phase[i].ncontrols;
	    int nobserved = problem.phase[i].nobserved;

        int max_nodes = get_max_nodes(problem,i+1, &algorithm);


        workspace->dual_events[i].Resize(nevents,1);

        workspace->e[i].Resize(nevents,1);

        if (nparam>=1) {
          workspace->prev_param[i].Resize(nparam,1);
        }

        workspace->states[i]= new adouble[nstates];
        workspace->controls[i] = new adouble[ncontrols];
        workspace->parameters[i] = new adouble[nparam];
        workspace->resid[i]= new adouble[nstates];
        workspace->derivatives[i]= new adouble[nstates];
        workspace->initial_states[i]= new adouble[nstates];
        workspace->final_states[i]= new adouble[nstates];
        workspace->initial_controls[i]= new adouble[ncontrols];
        workspace->final_controls[i]= new adouble[ncontrols];
        workspace->events[i]= new adouble[nevents];
        workspace->path[i]= new adouble[npath];

        workspace->states_next[i]     = new adouble[nstates];
        workspace->controls_next[i]   = new adouble[ncontrols];
        workspace->derivatives_next[i]= new adouble[nstates];
        workspace->path_next[i]       = new adouble[npath];
        workspace->states_bar[i]      = new adouble[nstates];
        workspace->controls_bar[i]    = new adouble[ncontrols];
        workspace->derivatives_bar[i] = new adouble[nstates];

        workspace->path_bar[i]        = new adouble[npath];

   	workspace->observed_variable[i] = new adouble[nobserved];
	workspace->observed_residual[i] = new adouble[nobserved];
  	workspace->lam_resid[i]              = new adouble[nobserved];

   	workspace->interp_states_pe[i]   = new adouble[nstates];
	workspace->interp_controls_pe[i] = new adouble[ncontrols];

        workspace->states_traj[i]= new adouble[problem.phase[i].nstates*(max_nodes +1)];
        workspace->derivs_traj[i]= new adouble[problem.phase[i].nstates*(max_nodes +1)];


  }

  int dotindex = problem.outfilename.find_first_of(".");

  workspace->igroup = new IGroup;

  string fname = "psopt_solution_" + problem.outfilename.substr(0,dotindex) + ".txt";

  workspace->psopt_solution_summary_file = fopen(fname.c_str(),"w");

  workspace->mesh_statistics = fopen("mesh_statistics.txt","w");

  fname = "mesh_statistics_" + problem.outfilename.substr(0,dotindex) + ".tex";

  workspace->mesh_statistics_tex = fopen( fname.c_str(),"w");

  if (workspace->psopt_solution_summary_file == NULL) error_message("Error opening \"psopt_solution_summary.txt\" file");

  if (workspace->mesh_statistics == NULL) error_message("Error opening \"mesh_statistics.txt\" file");

  if (workspace->mesh_statistics_tex == NULL) error_message("Error opening \"mesh_statistics.tex\" file");

  for(i=0;i<problem.nphases;i++) {
     workspace->emax_history[i].Resize(50,2);
  }

  workspace->enable_nlp_counters = false;

  workspace->auto_linked_flag = false;

// Initialise tape tags to be used by ADOL_C

  workspace->tag_f        = 1;
  workspace->tag_g 	     = 2;
  workspace->tag_hess     = 3;
  workspace->tag_fg 	     = 4;
  workspace->tag_gc       = 5;

  workspace->user_data = problem.user_data;

}

void resize_workspace_vars(Prob& problem, Alg& algorithm, Sol& solution, Workspace* workspace)
{
  int i;

  int nlp_ncons = get_number_nlp_constraints(problem, workspace );
  int nvars     = get_number_nlp_vars(problem, workspace);

  workspace->Xsnopt->Resize(nvars, 1);
  workspace->gsnopt->Resize(nlp_ncons, 1);

  workspace->Xip->Resize(nvars,1);
  workspace->JacRow->Resize(1,nvars);
  workspace->Gip->Resize(nlp_ncons,1);
  workspace->GFip->Resize(nvars,1);

  (*workspace->grw->dfdx_j).Resize( nlp_ncons,1 );
  (*workspace->grw->F1).Resize( nlp_ncons, 1 );
  (*workspace->grw->F2).Resize( nlp_ncons, 1 );
  (*workspace->grw->F3).Resize( nlp_ncons, 1 );
  (*workspace->grw->F4).Resize( nlp_ncons, 1 );

  workspace->JacCol1->Resize(nlp_ncons,1);
  workspace->JacCol2->Resize(nlp_ncons,1);
  workspace->JacCol3->Resize(nlp_ncons,1);
  workspace->xp->Resize(nvars,1);
  workspace->constraint_scaling->Resize(nlp_ncons,1);


  solution.xad = workspace->xad;

  workspace->x0->Resize(nvars,1);
  workspace->lambda->Resize(nlp_ncons,1);

  workspace->xlb->Resize(nvars,1);
  workspace->xub->Resize(nvars,1);
  
  workspace->nphases = problem.nphases;

  for(i=0; i< problem.nphases; i++)
  {
        int nstates   = problem.phase[i].nstates;
        int npath     = problem.phase[i].npath;
        int norder    = problem.phase[i].current_number_of_intervals;

        workspace->w[i].Resize(norder+1,1);

        if ( algorithm.collocation_method == "Legendre") {
	        workspace->P[i].Resize(norder+1,norder+1);
	}
	if ( use_global_collocation(algorithm) ) {
                workspace->D[i].Resize(norder+1,norder+1);
	}

        workspace->snodes[i].Resize(1,norder+1);



        workspace->dual_costates[i].Resize(nstates,norder+1);
        workspace->dual_path[i].Resize(npath,norder+1);
        workspace->DerivResid[i].Resize(nstates,norder+1);
        workspace->Xdot[i].Resize(nstates,norder+1);
        workspace->DerivResid[i].Resize(nstates,norder+1);
        workspace->Xdotgg[i].Resize(nstates,norder+1);
        workspace->hgg[i].Resize(npath,norder+1);
	workspace->order_reduction[i].Resize(1,norder);

  }

}



work_str::~work_str()
{
  for(long unsigned int i=0; i< this->nphases; i++)
  {
    delete [] this->states[i];
    delete [] this->controls[i];
    delete [] this->parameters[i];
    delete [] this->resid[i];
    delete [] this->derivatives[i];
    delete [] this->initial_states[i];
    delete [] this->final_states[i];
    delete [] this->initial_controls[i];
    delete [] this->final_controls[i];
    delete [] this->events[i];
    delete [] this->path[i];

    delete [] this->states_next[i];
    delete [] this->controls_next[i];
    delete [] this->derivatives_next[i];
    delete [] this->path_next[i];
    delete [] this->states_bar[i];
    delete [] this->controls_bar[i];
    delete [] this->derivatives_bar[i];

    delete [] this->path_bar[i];

    delete [] this->observed_variable[i];
    delete [] this->observed_residual[i];
    delete [] this->lam_resid[i];

    delete [] this->interp_states_pe[i];
    delete [] this->interp_controls_pe[i];

    delete [] this->states_traj[i];
    delete [] this->derivs_traj[i];
  }

  if (this->G2) delete [] this->G2;
  if (this->hess_ir) delete [] this->hess_ir;
  if (this->hess_jc) delete [] this->hess_jc;
  if (this->iArow) delete [] this->iArow;
  if (this->iGfun1) delete [] this->iGfun1;
  if (this->iGfun2) delete [] this->iGfun2;
  if (this->iGfun) delete [] this->iGfun;
  if (this->iGrow) delete [] this->iGrow;
  if (this->jac_Aij) delete [] this->jac_Aij;
  if (this->jac_Gij) delete [] this->jac_Gij;
  if (this->jAcol) delete [] this->jAcol;
  if (this->jGcol) delete [] this->jGcol;
  if (this->jGvar1) delete [] this->jGvar1;
  if (this->jGvar2) delete [] this->jGvar2;
  if (this->jGvar) delete [] this->jGvar;
  if (this->lambda_d) delete [] this->lambda_d;

  delete [] this->xad;
  delete [] this->gad;
  delete [] this->fgad;
  delete [] this->fg;
  delete [] this->nrm_row;

  delete [] this->states;
  delete [] this->controls;
  delete [] this->parameters;
  delete [] this->resid;
  delete [] this->derivatives;
  delete [] this->initial_states;
  delete [] this->final_states;
  delete [] this->initial_controls;
  delete [] this->final_controls;
  delete [] this->events;
  delete [] this->path;
  delete [] this->states_traj;
  delete [] this->derivs_traj;
  delete [] this->linkages;
  delete [] this->states_next;
  delete [] this->controls_next;
  delete [] this->derivatives_next;
  delete [] this->path_next;
  delete [] this->states_bar;
  delete [] this->controls_bar;
  delete [] this->derivatives_bar;
  delete [] this->path_bar;
  delete [] this->observed_variable;
  delete [] this->observed_residual;
  delete [] this->interp_states_pe;
  delete [] this->interp_controls_pe;
  delete [] this->lam_resid;

  delete [] this->time_array_tmp;
  delete [] this->single_trajectory_tmp;
  delete [] this->L_ad_tmp;
  delete [] this->u_spline;
  delete [] this->z_spline;
  delete [] this->y2a_spline;

  delete this->igroup;

  delete [] this->P;
  delete [] this->sindex;
  delete [] this->w;
  delete [] this->D;
  delete [] this->snodes;
  delete [] this->old_snodes;
  delete    this->xlb;
  delete    this->xub;
  delete    this->x0;
  delete    this->lambda;
  delete [] this->dual_costates;
  delete [] this->dual_events;
  delete [] this->dual_path;
  delete    this->linkage;
  delete    this->constraint_scaling;
  delete [] this->prev_states;
  delete [] this->prev_costates;
  delete [] this->prev_controls;
  delete [] this->prev_path;
  delete [] this->prev_param;
  delete [] this->prev_nodes;
  delete    this->Ax;
  delete    this->Gsp;
  delete    this->Xsnopt;
  delete    this->gsnopt;
  delete    this->Xip;
  delete    this->JacRow;
  delete    this->Gip;
  delete    this->GFip;
  delete [] this->Xdot;
  delete [] this->DerivResid;
  delete [] this->Xdotgg;
  delete [] this->e;
  delete [] this->hgg;
  delete    this->prev_t0;
  delete    this->prev_tf;
  delete    this->JacCol1;
  delete    this->JacCol2;
  delete    this->JacCol3;
  delete    this->xp;
  delete [] this->emax_history;
  delete [] this->order_reduction;
  delete [] this->old_relative_errors;
  delete [] this->error_scaling_weights;

  delete this->grw->dfdx_j;
  delete this->grw->F1;
  delete this->grw->F2;
  delete this->grw->F3;
  delete this->grw->F4;

  delete this->grw;

}

