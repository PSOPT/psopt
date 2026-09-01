/*********************************************************************************************

This file is part of the PSOPT library, a software tool for computational optimal control

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
#include "psopt_sqp.hpp"

// Bring std/Ipopt names into this translation unit (formerly leaked via psopt.h).
using namespace std;
using namespace Ipopt;

#ifdef PSOPT_ALLOW_ENV_OVERRIDES
// ---------------------------------------------------------------------------------
// Overriding algorithm settings from the environment.
//
// Built only when PSOPT_ALLOW_ENV_OVERRIDES is defined, which CMake leaves off. It
// exists for benchmarking: comparing two QP backends, or two of Betts's strategies,
// across the whole example set means running one binary many ways, and the alternative
// is editing sixty-six example sources or -- as was done here for a while -- patching
// this file with a script before each sweep. That patch is a tracked source file
// modified out of band, which has now twice caused trouble: once a sweep whose examples
// had not all been rebuilt and quietly reported IPOPT's results as the SQP's, and once
// a commit that nearly shipped the patch itself.
//
// Every override is announced. A program that silently disregards the algorithm its own
// source asked for is a thoroughly confusing thing to debug, and the confusion is worse
// than the convenience is worth.
// ---------------------------------------------------------------------------------
static void psopt_env_override(const char* var, string& field, Workspace* workspace)
{
    const char* v = getenv(var);
    if (v == NULL || field == v) return;
    snprintf(workspace->text, sizeof(workspace->text),
             ">>> %s overrides the algorithm setting in the source: \"%s\" -> \"%s\"\n",
             var, field.c_str(), v);
    psopt_print(workspace, workspace->text);
    field = v;
}

static void psopt_env_override_int(const char* var, int& field, Workspace* workspace)
{
    const char* v = getenv(var);
    if (v == NULL) return;
    const int nv = atoi(v);
    if (nv == field) return;
    snprintf(workspace->text, sizeof(workspace->text),
             ">>> %s overrides the algorithm setting in the source: %d -> %d\n",
             var, field, nv);
    psopt_print(workspace, workspace->text);
    field = nv;
}

static void psopt_env_override_double(const char* var, double& field, Workspace* workspace)
{
    const char* v = getenv(var);
    if (v == NULL) return;
    const double nv = atof(v);
    if (nv == field) return;
    snprintf(workspace->text, sizeof(workspace->text),
             ">>> %s overrides the algorithm setting in the source: %g -> %g\n",
             var, field, nv);
    psopt_print(workspace, workspace->text);
    field = nv;
}

static void psopt_apply_environment_overrides(Alg& algorithm, Workspace* workspace)
{
    psopt_env_override("PSOPT_NLP_METHOD",       algorithm.nlp_method,      workspace);
    psopt_env_override("PSOPT_DERIVATIVES",      algorithm.derivatives,     workspace);
    psopt_env_override("PSOPT_HESSIAN",          algorithm.hessian,         workspace);
    psopt_env_override("PSOPT_QP_SOLVER",        algorithm.qp_solver,       workspace);
    psopt_env_override("PSOPT_QP_RESTORATION",   algorithm.qp_restoration,  workspace);
    psopt_env_override("PSOPT_ELASTIC_PENALTY",  algorithm.elastic_penalty, workspace);
    psopt_env_override("PSOPT_SQP_STRATEGY",     algorithm.sqp_strategy,    workspace);
    psopt_env_override_int("PSOPT_QP_ITER_MAX",  algorithm.qp_iter_max,     workspace);
    // Asking for a tolerance the solver cannot reach is how the acceptable-level
    // termination is exercised deliberately rather than waited for.
    psopt_env_override_double("PSOPT_NLP_TOLERANCE", algorithm.nlp_tolerance, workspace);
    psopt_env_override_int("PSOPT_NLP_ITER_MAX",  algorithm.nlp_iter_max,    workspace);

    // The Hessian setting is read again through the workspace further down, so it has to
    // be carried across as well; the others are taken from this Alg.
    if (getenv("PSOPT_HESSIAN") && workspace->algorithm)
        workspace->algorithm->hessian = algorithm.hessian;
}
// Applied before the mesh loop rather than with the rest, and it has to be. The others
// are read once per NLP solve, which is where they are applied; mesh_refinement is read
// by psopt_main *before* the first solve, to fix how many mesh iterations the loop will
// run. Applied at the same point as the others it changes the option after that bound
// has been taken, and a manual mesh with one node entry is then indexed a second time --
// out of range, an unchecked Eigen read, a node count of -1 and a bad_alloc from the
// workspace resize. That is how this was found.
//
// Turning mesh refinement off is what makes a comparison of two NLP solvers a controlled
// one: with it on, each solver's answer steers its own next mesh, so after the first
// refinement the two are not solving the same problem and their objectives are not
// comparable. examples/wheat and examples/mpec turned out to be demonstrating exactly
// that.
// Applied before the Workspace exists, so it cannot report through it. The Workspace is
// sized from algorithm.mesh_refinement -- get_max_nodes asks a different question of the
// node schedule under "manual" than under "automatic" -- so an override applied after the
// allocation leaves every array in the Workspace sized for the mesh schedule that was NOT
// used. Writing past them is silent on glibc and traps on macOS; see the note in psopt.cxx.
void psopt_apply_mesh_environment_override(Alg& algorithm)
{
    const char* v = getenv("PSOPT_MESH_REFINEMENT");
    if (v != NULL && algorithm.mesh_refinement != v) {
        if (algorithm.print_level)
            fprintf(stderr, ">>> PSOPT_MESH_REFINEMENT overrides the algorithm setting in "
                            "the source: \"%s\" -> \"%s\"\n",
                    algorithm.mesh_refinement.c_str(), v);
        algorithm.mesh_refinement = v;
    }
    // Scaling belongs here rather than with the per-solve overrides for the same reason:
    // the factors are computed once, before the mesh loop, so an override applied at the
    // NLP interface arrives after they have been fixed and does nothing at all. "user"
    // with no factors set is unit scaling, since psopt_level1_setup and psopt_level2_setup
    // initialise every factor to one, which makes a scaling study over the whole example
    // set a sweep rather than sixty-six edited sources.
    const char* sc = getenv("PSOPT_SCALING");
    if (sc != NULL && algorithm.scaling != sc) {
        if (algorithm.print_level)
            fprintf(stderr, ">>> PSOPT_SCALING overrides the algorithm setting in the "
                            "source: \"%s\" -> \"%s\"\n", algorithm.scaling.c_str(), sc);
        algorithm.scaling = sc;
    }

    const char* w = getenv("PSOPT_MR_SWITCH_DETECTION");
    if (w != NULL && atoi(w) != algorithm.mr_switch_detection) {
        if (algorithm.print_level)
            fprintf(stderr, ">>> PSOPT_MR_SWITCH_DETECTION overrides the algorithm setting "
                            "in the source: %d -> %d\n",
                    algorithm.mr_switch_detection, atoi(w));
        algorithm.mr_switch_detection = atoi(w);
    }
}

#endif  // PSOPT_ALLOW_ENV_OVERRIDES





int NLP_interface(
         Alg& algorithm,
         MatrixXd* x0,
         double (*f)(MatrixXd&, Workspace*),
	 void (*g)(MatrixXd&,MatrixXd*, Workspace*),
	 int nlp_ncons,
	 int nlp_neq  ,
	 MatrixXd* xlb,
	 MatrixXd* xub,
         MatrixXd* lambda,
         int hotflag,
         int iprint,
         Workspace* workspace,
         void* user_data
         )
{

    Sol*  solution= workspace->solution;

#ifdef PSOPT_ALLOW_ENV_OVERRIDES
    psopt_apply_environment_overrides(algorithm, workspace);
#endif
    if ( algorithm.nlp_method=="SQP" )
    {
        // PSOPT's own sequential quadratic programming solver. Dense at this stage:
        // see include/psopt_sqp.hpp for what that costs and what the sparse stage
        // will change.
        return SQP_interface(algorithm, x0, f, g, nlp_ncons, nlp_neq,
                             xlb, xub, lambda, hotflag, iprint, workspace, user_data);
    }

    else if ( algorithm.nlp_method=="IPOPT" )
    {


  // Create a new instance of nlp
  SmartPtr<TNLP> mynlp = new IPOPT_PSOPT(workspace, user_data);

  // Create a new instance of IpoptApplication
  SmartPtr<IpoptApplication> app = new IpoptApplication();

  // Change some options
  app->Options()->SetNumericValue("tol", workspace->algorithm->nlp_tolerance );
  app->Options()->SetStringValue("mu_strategy", "adaptive");
  app->Options()->SetStringValue("output_file", "ipopt.out");
  app->Options()->SetStringValue("nlp_scaling_method","gradient-based");
  app->Options()->SetStringValue("linear_solver",workspace->algorithm->ipopt_linear_solver);
  app->Options()->SetNumericValue("max_cpu_time", workspace->algorithm->ipopt_max_cpu_time );


  if ( ( useAutomaticDifferentiation(algorithm) && algorithm.hessian=="exact" )
       || algorithm.hessian=="numerical" ) {
     app->Options()->SetStringValue("hessian_approximation", "exact");
  }
  else {
     app->Options()->SetStringValue("hessian_approximation", "limited-memory");
  }

  if (algorithm.print_level==0) {
      app->Options()->SetIntegerValue("print_level", 0);
  }
  else if (algorithm.print_level>0) {
      app->Options()->SetIntegerValue("print_level", 5);
  }

  app->Options()->SetIntegerValue("max_iter", workspace->algorithm->nlp_iter_max);
  if (hotflag) {
     app->Options()->SetStringValue("warm_start_init_point", "yes");
  }
  else {
	app->Options()->SetStringValue("warm_start_init_point", "no");
  }

  // Intialize the IpoptApplication and process the options
  ApplicationReturnStatus status;
  status = app->Initialize();
  if (status != Solve_Succeeded) {
    printf("\n\n*** Error during initialization!\n");
    return (int) status;
  }

  // Ask Ipopt to solve the problem
  status = app->OptimizeTNLP(mynlp);

  if (psopt_ipopt_solved((int) status)) {
    psopt_print(workspace,"\n\n*** The problem has been solved!\n");
  }
  else {
    psopt_print(workspace,"\n\n*** The problem FAILED!\n");
  }

  // As the SmartPtrs go out of scope, the reference count
  // will be decremented and the objects will automatically
  // be deleted.

  solution->nlp_return_code = (int) status;

  return (int) status;


    }

    else {
        snprintf(workspace->text,sizeof(workspace->text),"\n Incorrect NLP method has been specified");
        error_message(workspace->text);
    }

    return 0;

}


