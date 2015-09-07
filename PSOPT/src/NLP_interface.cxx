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

#ifdef USE_SNOPT
#include "snoptProblem.hpp"
Workspace* tempsnoptworkspace;
#endif


int NLP_interface(
         Alg& algorithm,
         DMatrix* x0,
         double (*f)(DMatrix&, Workspace*),
	 void (*g)(DMatrix&,DMatrix*, Workspace*),
	 int nlp_ncons,
	 int nlp_neq  ,
	 DMatrix* xlb,
	 DMatrix* xub,
         DMatrix* lambda,
         int hotflag,
         int iprint,
         Workspace* workspace,
         void* user_data
         )
{
    Prob *problem = workspace->problem;
    Sol*  solution= workspace->solution;

    int  i;

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
  memcpy(xlow, xlb->GetPr(), n*sizeof(double) );
  memcpy(xupp, xub->GetPr(), n*sizeof(double) );
  for (int ix = 0; ix < n; ++ix) {
      x[ix]=0.0;
      xstate[ix]=3;
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

#ifdef ADOLC_VERSION_1
     sparse_jac(workspace->tag_fg, neF, n, 0, x, &workspace->F_nnz, &workspace->iGfun2, &workspace->jGvar2, &workspace->G2);
#endif


#ifdef ADOLC_VERSION_2
    int options[4];
    options[0]=0; options[1]=0; options[2]=0;options[3]=0;
    sparse_jac(workspace->tag_fg, neF, n, 0, x, &workspace->F_nnz, &workspace->iGfun2, &workspace->jGvar2, &workspace->G2, options);
#endif

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
          	  workspace->iGfun1[i] = workspace->iGfun2[i]+1;
        	  workspace->jGvar1[i] = workspace->jGvar2[i]+1;
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

  SparseMatrix As(A, neF, n, neA, iAfuni, jAvari);

//  As.SaveSparsityPattern("SNOPT_Linear_pattern.txt");

  workspace->As = &As;

  // Set optimizer options.
  snprob.setIntParameter("Derivative option", 0);
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
  memcpy(x0->GetPr(), x, n*sizeof(double));

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

/*
// ************* C INTERFACE ************************
  integer Cold = 0, Basis = 1, Warm = 2;
  integer StartOption;
  integer    iSpecs = 4,  spec_len;
  integer    iSumm  = 6;
  integer    iPrint = 1,  prnt_len;

  integer    INFO;
  integer    neA, neG;

  integer    nxname = 1, nFname = 1, npname;
  char       Prob[200];

  integer    minrw, miniw, mincw;
  integer    lenrw = 4000*(neF+n), leniw = 2000*(neF+n), lencw = 1000;
  doublereal *rw  = (doublereal*) my_calloc(lenrw, sizeof(doublereal));
  integer    *iw  = (integer*)    my_calloc(leniw, sizeof(integer));

  char       cw[8*500];

  char       printname[200];
  char       specname[200];

  integer    nS, nInf;
  doublereal sInf;
  integer    DerOpt, Major, iSum, iPrt, strOpt_len;
  char       strOpt[200];
  integer    IterOpt;

  if ( !algorithm.print_level )
  {
  	iSumm  = 0;
	iPrint = 0;
  }

  // open output files using snfilewrappers.[ch]
  sprintf(specname ,   "%s", "psopt.spc");   spec_len = strlen(specname);
  sprintf(printname,   "%s", "psopt.out");   prnt_len = strlen(printname);

  // Open the print file, fortran style
  snopenappend_
    ( &iPrint, printname,   &INFO, prnt_len );

  sninit_
    ( &iPrint, &iSumm, cw, &lencw, iw, &leniw, rw, &lenrw, 8*500 );

  sprintf(Prob,"%s",problem->name.c_str());

  workspace->jac_done = 0;

  snjac_
    ( &INFO, &neF, &n,  snPSOPTusrf_,
      iAfun, jAvar, &lenA, &neA, A,
      iGfun, jGvar, &lenG, &neG,
      x, xlow, xupp, &mincw, &miniw, &minrw,
      cw, &lencw, iw, &leniw, rw, &lenrw,
      cw, &lencw, iw, &leniw, rw, &lenrw,
      8*500, 8*500 );

  // Store co-ordinates of nonlinear Jacobian

  for(k=0; k<neG; k++) {
        workspace->iGfun[k] = (unsigned int) iGfun[k];
        workspace->jGvar[k] = (unsigned int) jGvar[k];
  }

  (*workspace->Gsp).Resize(neF,n,neG);

  if ( useAutomaticDifferentiation(algorithm) )
  {
     adouble* xad  = workspace->xad;
     adouble* fgad = workspace->fgad;
     double*    fg = workspace->fg;

     // Tracing of function fg()
     trace_on(workspace->tag_fg);
     for(i=0;i<n;i++)
		xad[i] <<= x[i];

     fg_ad(xad, fgad, workspace);

     for(i=0;i<neF;i++)
	fgad[i] >>= fg[i];
     trace_off();

#ifdef ADOLC_VERSION_1
     sparse_jac(workspace->tag_fg, neF, n, 0, x, &workspace->F_nnz, &workspace->iGfun2, &workspace->jGvar2, &workspace->G2);
#endif


#ifdef ADOLC_VERSION_2
    int options[4];
    options[0]=0; options[1]=0; options[2]=0;options[3]=0;
    sparse_jac(workspace->tag_fg, neF, n, 0, x, &workspace->F_nnz, &workspace->iGfun2, &workspace->jGvar2, &workspace->G2, options);
#endif

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
          	  workspace->iGfun1[i] = workspace->iGfun2[i]+1;
        	  workspace->jGvar1[i] = workspace->jGvar2[i]+1;
      }

  }


  workspace->jac_done=1;

  int * iAfuni = (int*) my_calloc(neA, sizeof(int) );
  int * jAvari = (int*) my_calloc(neA, sizeof(int) );

  for (i=0; i< neA; i++) {
       iAfuni[i] = iAfun[i];
       jAvari[i] = jAvar[i];
  }

  SparseMatrix As(A, neF, n, neA, iAfuni, jAvari);

//  As.SaveSparsityPattern("SNOPT_Linear_pattern.txt");

  workspace->As = &As;


  DerOpt = 0;
  iPrt   = 0;
  iSum   = 0;
  sprintf(strOpt,"%s","Derivative option");
  strOpt_len = strlen(strOpt);
  snseti_
    ( strOpt, &DerOpt, &iPrt, &iSum, &INFO,
      cw, &lencw, iw, &leniw, rw, &lenrw, strOpt_len, 8*500 );



  IterOpt = workspace->algorithm->nlp_iter_max;
  iPrt   = 0;
  iSum   = 0;
  sprintf(strOpt,"%s","Major iterations limit");
  strOpt_len = strlen(strOpt);
  snseti_
    ( strOpt, &IterOpt, &iPrt, &iSum, &INFO,
      cw, &lencw, iw, &leniw, rw, &lenrw, strOpt_len, 8*500 );
  sprintf(strOpt,"%s","Minor iterations limit");
  strOpt_len = strlen(strOpt);
  snseti_
    ( strOpt, &IterOpt, &iPrt, &iSum, &INFO,
      cw, &lencw, iw, &leniw, rw, &lenrw, strOpt_len, 8*500 );
  sprintf(strOpt,"%s","Iterations limit");
  IterOpt = 50*IterOpt;
  strOpt_len = strlen(strOpt);
  snseti_
    ( strOpt, &IterOpt, &iPrt, &iSum, &INFO,
      cw, &lencw, iw, &leniw, rw, &lenrw, strOpt_len, 8*500 );
  if (!algorithm.print_level) {
  	sprintf(strOpt,"%s","Major print level 0");
  	IterOpt = 0;
  	strOpt_len = strlen(strOpt);
        snset_
        ( strOpt, &iPrt, &iSum, &INFO,
            cw, &lencw, iw, &leniw, rw, &lenrw, strOpt_len, 8*500 );
  	sprintf(strOpt,"%s","Minor print level 0");
  	IterOpt = 0;
  	strOpt_len = strlen(strOpt);
        snset_
        ( strOpt, &iPrt, &iSum, &INFO,
            cw, &lencw, iw, &leniw, rw, &lenrw, strOpt_len, 8*500 );
  }

  // Set some extra options

  sprintf(strOpt,"%s","LU complete pivoting");
  strOpt_len = strlen(strOpt);
  snset_
    ( strOpt, &iPrt, &iSum, &INFO,
      cw, &lencw, iw, &leniw, rw, &lenrw, strOpt_len, 8*500 );



  //     ------------------------------------------------------------------ 
  //     Go for it                                                          
  //     ------------------------------------------------------------------ 

  if (hotflag) StartOption=Warm;
  else StartOption=Cold;

  snopta_
    ( &StartOption, &neF, &n, &nxname, &nFname,
      &ObjAdd, &ObjRow, Prob,  snPSOPTusrf_,
      iAfun, jAvar, &lenA, &neA, A,
      iGfun, jGvar, &lenG, &neG,
      xlow, xupp, xnames, Flow, Fupp, Fnames,
      x, xstate, xmul, F, Fstate, Fmul,
      &INFO, &mincw, &miniw, &minrw,
      &nS, &nInf, &sInf,
      cw, &lencw, iw, &leniw, rw, &lenrw,

      cw, &lencw, iw, &leniw, rw, &lenrw,
      npname, 8*nxname, 8*nFname,
      8*500, 8*500);

   solution->nlp_return_code = int (INFO/10);

// *************************************************


// Copy results...

  memcpy( x0->GetPr(), x, n*sizeof(double) );

  for(ii=0;ii<n;ii++) solution->xad[ii] = x[ii];


  for (ii=1;ii<neF;ii++) {
       (*lambda)(ii) = Fmul[ii];
  }


*/

  return 0;

#else
        sprintf(workspace->text,"\nSNOPT method has been specified but not linked");
        error_message(workspace->text);
#endif

    }

    else if ( algorithm.nlp_method=="IPOPT" )
    {

#ifdef USE_IPOPT


  // Create a new instance of nlp
  SmartPtr<TNLP> mynlp = new IPOPT_PSOPT(workspace, user_data);

  // Create a new instance of IpoptApplication
  SmartPtr<IpoptApplication> app = new IpoptApplication();

  // Change some options
  app->Options()->SetNumericValue("tol", workspace->algorithm->nlp_tolerance );
  app->Options()->SetStringValue("mu_strategy", "adaptive");
  app->Options()->SetStringValue("output_file", "ipopt.out");
  app->Options()->SetStringValue("nlp_scaling_method","gradient-based");
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
    psopt_print(workspace,"\n\n*** The problem solved!\n");
  }
  else {
    psopt_print(workspace,"\n\n*** The problem FAILED!\n");
  }

  // As the SmartPtrs go out of scope, the reference count
  // will be decremented and the objects will automatically
  // be deleted.

  solution->nlp_return_code = (int) status;

  return (int) status;


#else
        sprintf(workspace->text,"\nIPOPT method has been specified but not linked");
        error_message(workspace->text);
#endif


    }

    else {
        sprintf(workspace->text,"\n Incorrect NLP method has been specified");
        error_message(workspace->text);
    }

    return 0;

}


