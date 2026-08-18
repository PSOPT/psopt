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

static void psopt_apply_environment_overrides(Alg& algorithm, Workspace* workspace)
{
    psopt_env_override("PSOPT_NLP_METHOD",       algorithm.nlp_method,      workspace);
    psopt_env_override("PSOPT_HESSIAN",          algorithm.hessian,         workspace);
    psopt_env_override("PSOPT_QP_SOLVER",        algorithm.qp_solver,       workspace);
    psopt_env_override("PSOPT_QP_RESTORATION",   algorithm.qp_restoration,  workspace);
    psopt_env_override("PSOPT_ELASTIC_PENALTY",  algorithm.elastic_penalty, workspace);
    psopt_env_override("PSOPT_SQP_STRATEGY",     algorithm.sqp_strategy,    workspace);
    psopt_env_override_int("PSOPT_QP_ITER_MAX",  algorithm.qp_iter_max,     workspace);

    // The Hessian setting is read again through the workspace further down, so it has to
    // be carried across as well; the others are taken from this Alg.
    if (getenv("PSOPT_HESSIAN") && workspace->algorithm)
        workspace->algorithm->hessian = algorithm.hessian;
}
#endif  // PSOPT_ALLOW_ENV_OVERRIDES



#ifdef USE_SNOPT
#include "snoptProblem.hpp"
#include "snopt_psopt.h"
// Workspace* tempsnoptworkspace;

Workspace* snoptProbLocal::workspace   = NULL;
void*      snoptProbLocal::_user_data  = NULL;

#endif




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
    [[maybe_unused]] int use_sparse_jac_function = 1; // referenced only in the USE_SNOPT path

#ifdef PSOPT_ALLOW_ENV_OVERRIDES
    psopt_apply_environment_overrides(algorithm, workspace);
#endif
    if ( algorithm.nlp_method=="SNOPT" )
    {


#ifdef USE_SNOPT
  // C++ interface to SNOPT.
  // Create a new instance of snoptProbLocal
  
  snoptProbLocal snprob(workspace, user_data);



  // Allocate and initialize. 
  int n     =  length(*x0);
  int neF   =  nlp_ncons+1;

  int lenA  =  (nlp_ncons+1)*n;

  int *iAfun = new int[lenA];
  int *jAvar = new int[lenA];
  double *A  = new double[lenA];

  int lenG   = lenA;
  int *iGfun = new int[lenG];
  int *jGvar = new int[lenG];

  double *x      = new double[n];
  double *xlow   = new double[n];
  double *xupp   = new double[n];
  double *xmul   = new double[n];
  int    *xstate = new int[n];

  double *F      = new double[neF];
  double *Flow   = new double[neF];
  double *Fupp   = new double[neF];
  double *Fmul   = new double[neF];
  int    *Fstate = new int[neF];



  int    ObjRow = 0;
  double ObjAdd = 0;

  int Cold = 0;
  int Basis = 1;
  int Warm = 2;
  int StartOption;

  double InfValue = 1.0e20;
  
  int nInf;
  double sInf;
  
        

  memcpy(xlow, &(*xlb)(0), n*sizeof(double) );

  memcpy(xupp, &(*xub)(0), n*sizeof(double) );
  
  for (int ix = 0; ix < n; ++ix) {
      x[ix]=0.0;
      xstate[ix]=0;
  }

  Flow[0] = -InfValue;
  Fupp[0] =  InfValue;
  Fmul[0] =  0.0;

  get_constraint_bounds( &Flow[1], &Fupp[1], workspace );


  for (int iF = 1; iF < neF; ++iF) {
     Fmul[iF] = (*lambda)(iF-1);
  }

  for (int iF = 0; iF < neF; ++ iF) {
     Fstate[iF] = 0; //  why is this 0 and xstate is 3?
  }

  memcpy(x, &(*x0)(0), n*sizeof(double) );
  
  
  snprob.initialize    ("", 1);      
  snprob.setPrintFile  ("snopt.out"); 
  snprob.setProbName   ("psopt");  




  int neA = -1;
  int neG = -1;


  workspace->jac_done = 0;
  
  if ( useAutomaticDifferentiation(algorithm) )
  {
  	

     psopt_ad::ad_record(workspace->ad_fg, n, neF, x,
         [&](const adouble* xin, adouble* yout){ fg_ad(const_cast<adouble*>(xin), yout, workspace); });
     psopt_ad::SparseTriplet Jfg = psopt_ad::ad_sparse_jacobian(workspace->ad_fg, x, /*reuse=*/false);
     workspace->F_nnz = Jfg.nnz();
     workspace->iGfun2.assign(Jfg.row.begin(), Jfg.row.end());
     workspace->jGvar2.assign(Jfg.col.begin(), Jfg.col.end());
     workspace->G2.assign(Jfg.val.begin(), Jfg.val.end());


     snprintf(workspace->text,sizeof(workspace->text),"\nJacobian sparsity detected using ADOLC:");
     psopt_print(workspace,workspace->text);

     double jsratio = (double) ((double)  workspace->F_nnz/((double) (n*neF)));

     if (jsratio > workspace->algorithm->jac_sparsity_ratio) {
           snprintf(workspace->text,sizeof(workspace->text), "increase algorithm.jac_sparsity_ratio to just above %f", jsratio);
           error_message(workspace->text);
     }

     snprintf(workspace->text,sizeof(workspace->text),"\n%i nonzero elements out of %li [ratio=%f]\n", workspace->F_nnz, n*neF, jsratio);
     psopt_print(workspace,workspace->text);

     if (!use_sparse_jac_function) {
        DetectJacobianSparsityAD(fg_num, *x0, neF,  &neA,  iAfun, jAvar, A,
                                            &neG,  iGfun, jGvar, 
                                            workspace->grw, workspace );
     }
     else {
       snprob.computeJac(neF, n, snoptProbLocal::snPSOPTusrf_,
		  x, xlow, xupp,
		  iAfun, jAvar, A, neA,
		  iGfun, jGvar, neG);     	                                          
     }                                         

  }
  
  else {
   

   if (use_sparse_jac_function) {
      snprob.computeJac(neF, n, snoptProbLocal::snPSOPTusrf_,
		  x, xlow, xupp,
		  iAfun, jAvar, A, neA,
		  iGfun, jGvar, neG);
		  
	}
	else {

     DetectJacobianSparsity(fg_num, *x0, neF,  &neA,  iAfun, jAvar, A,
                                            &neG,  iGfun, jGvar, 
                                            workspace->grw, workspace ); 
                                            

   }                          

  }

  int * iGfuni = new int[neG];
  int * jGvari = new int[neG];  
  
  for (int iG = 0; iG < neG; iG++) {
  	  if (use_sparse_jac_function) {
        workspace->iGfun1[iG] = iGfun[iG]-1;
        workspace->jGvar1[iG] = jGvar[iG]-1;
        iGfuni[iG] = workspace->iGfun1[iG];
        jGvari[iG] = workspace->jGvar1[iG];
     } else {
        iGfuni[iG] = iGfun[iG];
        jGvari[iG] = jGvar[iG];
     }
  }
  
  // Create and store matrix A as a triplet based sparse matrix.
  int * iAfuni = new int[neA];
  int * jAvari = new int[neA];
  

  for (int iA = 0; iA < neA; iA++) {
  	 if (use_sparse_jac_function) {
       iAfuni[iA] = iAfun[iA]-1;
       jAvari[iA] = jAvar[iA]-1;
    } else {
       iAfuni[iA] = iAfun[iA];
       jAvari[iA] = jAvar[iA];
    }

  }

  TripletSparseMatrix As(A, neF, n, neA, iAfuni, jAvari);
  


  workspace->As = &As;
  
//  (*workspace->Gsp).resize(neF, n, neG);


  workspace->jac_done=1;
  

  // Set optimizer options.
  int derivative_option;
  if ( useAutomaticDifferentiation(algorithm) ) {
       derivative_option = 1;
  }
  else {
       derivative_option = 0;
  }
    
  snprob.setIntParameter("Derivative option", derivative_option);

  snprob.setIntParameter("Verify level ", 3);
  snprob.setIntParameter("Major iterations limit",
		  workspace->algorithm->nlp_iter_max);
  snprob.setIntParameter("Minor iterations limit",
		  workspace->algorithm->nlp_iter_max);   
  snprob.setIntParameter("Iterations limit",
		  50 * workspace->algorithm->nlp_iter_max);	  		  
  snprob.setRealParameter("Major optimality tolerance",
		  workspace->algorithm->nlp_tolerance);
  if (!algorithm.print_level) {
    snprob.setIntParameter("Major print level", 0);
    snprob.setIntParameter("Minor print level", 0);
  }
  snprob.setParameter("LU complete pivoting");

  if (hotflag) StartOption=Warm;
  else StartOption=Cold;

  // ---------
  // Go for it
  // ---------

  // int inform = snprob.solve(StartOption);
  int inform;
  
  if (useAutomaticDifferentiation(algorithm)) {
  
		
		inform = snprob.solve(StartOption, neF, n, ObjAdd, ObjRow, snoptProbLocal::snPSOPTusrf_,
		iAfuni, jAvari, A, neA,
		iGfuni, jGvari, neG,
		xlow, xupp, Flow, Fupp,
		x, xstate, xmul,
		F, Fstate, Fmul,
		workspace->nS, nInf, sInf);		
		
		
  }
  else {
        
     inform = snprob.solve(StartOption, neF, n, ObjAdd, ObjRow, snoptProbLocal::snPSOPTusrf_,
              xlow, xupp, Flow, Fupp,
              x, xstate, xmul, F, Fstate, Fmul,
              workspace->nS, nInf, sInf);
              
  }

  solution->nlp_return_code = int (inform/10);

  // Copy results.

  memcpy( &(*x0)(0), x, n*sizeof(double));

  for (int ix = 0; ix < n; ++ix) solution->xad[ix] = x[ix];

  for (int iF = 1; iF < neF; ++iF) {
    (*lambda)(iF-1) = Fmul[iF];
  }

  // Delete allocated variables.
  delete [] iAfun;
  delete [] jAvar;
  delete [] A;

  delete [] iGfun;
  delete [] jGvar;
  
  delete [] x;
  delete [] xlow;
  delete [] xupp;
  delete [] xmul;
  delete [] xstate;

  delete [] F;
  delete [] Flow;
  delete [] Fupp;
  delete [] Fmul;
  delete [] Fstate;

  return 0;

#else
        snprintf(workspace->text,sizeof(workspace->text),"\nSNOPT method has been specified but not linked");
        error_message(workspace->text);
#endif

    }

    else if ( algorithm.nlp_method=="SQP" )
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

  if (status == Solve_Succeeded) {
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


