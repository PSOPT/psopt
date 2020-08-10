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

#ifdef USE_SNOPT
#include "snoptProblem.hpp"
Workspace* tempsnoptworkspace;
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


    if ( algorithm.nlp_method=="SNOPT" )
    {

#ifdef USE_SNOPT
  // C++ interface to SNOPT.
  snoptProblemA snprob;

  tempsnoptworkspace = workspace;

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

  // TODO: the below is for the old f2c interface and
  // I don't know how to port it over to the C++ interface.
  // integer nxnames = 1;
  // integer nFnames = 1;
  // char *xnames = (char*) my_calloc(nxnames*8, sizeof(char));
  // char *Fnames = (char*) my_calloc(nFnames*8, sizeof(char));

  int    ObjRow = 0;
  double ObjAdd = 0;

  int Cold = 0;
  int Basis = 1;
  int Warm = 2;
  int StartOption;

  double InfValue = 1.0e20;

  // TODO could be improved.
//  memcpy(xlow, xlb->GetPr(), n*sizeof(double) );
  memcpy(xlow, &(*xlb)(0), n*sizeof(double) );
//  memcpy(xupp, xub->GetPr(), n*sizeof(double) );
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
     Fmul[iF] = (*lambda)(iF);
  }

  for (int iF = 0; iF < neF; ++ iF) {
     Fstate[iF] = 0; // TODO why is this 0 and xstate is 3?
  }

  memcpy(x, x0->GetPr(), n*sizeof(double) );

  snprob.setProbName(problem->name.c_str());
  snprob.setPrintFile("snopt.out");

  snprob.setProblemSize(n, neF);
  snprob.setObjective(ObjRow, ObjAdd);
  snprob.setX(x, xlow, xupp, xmul, xstate);
  snprob.setF(F, Flow, Fupp, Fmul, Fstate);
  snprob.setUserFun(snPSOPTusrf_);

  // snopta will compute the Jacobian by finite-differences.
  int neA = -1;
  int neG = -1;
  snprob.setA(lenA, neA, iAfun, jAvar, A);
  snprob.setG(lenG, neG, iGfun, jGvar);

  workspace->jac_done = 0;

  // computeJac will determine sparsity and will set neA and neG.  
  snprob.computeJac(neA, neG);

  for (int iG = 0; iG < neG; ++iG) {
    workspace->iGfun[iG] = (unsigned int) iGfun[iG];
    workspace->jGvar[iG] = (unsigned int) jGvar[iG];
  }

  (*workspace->Gsp).Resize(neF, n, neG);

  if ( useAutomaticDifferentiation(algorithm) )
  {
     adouble* xad  = workspace->xad;
     adouble* fgad = workspace->fgad;
     double*    fg = workspace->fg;

     /* Tracing of function fg() */
     trace_on(workspace->tag_fg);
     for(i=0;i<n;i++)
		xad[i] <<= x[i];

     fg_ad(xad, fgad, workspace);

     for(i=0;i<neF;i++)
	     fgad[i] >>= fg[i];
     trace_off();



    int options[4];
    options[0]=0; options[1]=0; options[2]=0;options[3]=0;
    sparse_jac(workspace->tag_fg, neF, n, 0, x, &workspace->F_nnz, &workspace->iGfun2, &workspace->jGvar2, &workspace->G2, options);


     sprintf(workspace->text,"\nJacobian sparsity detected using ADOLC:");
     psopt_print(workspace,workspace->text);

     double jsratio = (double) ((double)  workspace->F_nnz/((double) (n*neF)));

     if (jsratio > workspace->algorithm->jac_sparsity_ratio) {
           sprintf(workspace->text, "increase algorithm.jac_sparsity_ratio to just above %f", jsratio);
           error_message(workspace->text);
     }

     sprintf(workspace->text,"\n%i nonzero elements out of %li [ratio=%f]\n", workspace->F_nnz, n*neF, jsratio);
     psopt_print(workspace,workspace->text);


      for (i=0;i<workspace->F_nnz;i++) {
           workspace->iGfun1[i] = workspace->iGfun2[i];  
        	  workspace->jGvar1[i] = workspace->jGvar2[i];
      }

  }


  workspace->jac_done=1;

  // Create and store matrix A (sparse).
  int * iAfuni = new int[neA];
  int * jAvari = new int[neA];

  for (int iA = 0; iA < neA; ++iA) {
       iAfuni[iA] = iAfun[iA];
       jAvari[iA] = jAvar[iA];
  }

  TripletSparseMatrix As(A, neF, n, neA, iAfuni, jAvari);

//  As.SaveSparsityPattern("SNOPT_Linear_pattern.txt");

  workspace->As = &As;

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

  int inform = snprob.solve(StartOption);

  solution->nlp_return_code = int (inform/10);

  // Copy results.
//  memcpy(x0->GetPr(), x, n*sizeof(double));
  memcpy( &(*x0)(0), x, n*sizeof(double));

  for (int ix = 0; ix < n; ++ix) solution->xad[ix] = x[ix];

  for (int iF = 1; iF < neF; ++iF) {
    (*lambda)(iF) = Fmul[iF];
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

  // TODO
  tempsnoptworkspace = NULL;

  return 0;

#else
        sprintf(workspace->text,"\nSNOPT method has been specified but not linked");
        error_message(workspace->text);
#endif

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


  if ( useAutomaticDifferentiation(algorithm) && algorithm.hessian=="exact" ) {
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
        sprintf(workspace->text,"\n Incorrect NLP method has been specified");
        error_message(workspace->text);
    }

    return 0;

}


