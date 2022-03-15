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



// constructor
IPOPT_PSOPT::IPOPT_PSOPT(Workspace *pr, void *user_data)
{
    workspace       = pr;
    _user_data      = user_data;
}

//destructor
IPOPT_PSOPT::~IPOPT_PSOPT()
{}

adouble Lagrangian_ad(adouble* xad, double* lambda, double& obj_factor, Index& m, Workspace* workspace)
{
	adouble L;
	adouble f;
	adouble *g = workspace->gad;
	Index i;

	L = obj_factor*ff_ad(xad, workspace);

	gg_ad(xad, g, workspace);

	for(i=0; i<m ; i++) {
		L += lambda[ i ]*g[ i ];
	}

	return L;
}


bool check_no_cancel(void *user_data)
{
#ifdef WIN32
    if (user_data)
    {
        HANDLE *handle = (HANDLE *)user_data;
        int no_cancel = (WAIT_OBJECT_0!=WaitForSingleObject(*handle, 0));
        if (!no_cancel)
            fprintf(stderr, "\n --- User cancel event received! ---");
        return no_cancel;
    }
#endif
    return true;
}


bool IPOPT_PSOPT::intermediate_callback(AlgorithmMode mode,
                                       Index iter, Number obj_value,
                                       Number inf_pr, Number inf_du,
                                       Number mu, Number d_norm,
                                       Number regularization_size,
                                       Number alpha_du, Number alpha_pr,
                                       Index ls_trials,
                                       const IpoptData* ip_data,
                                       IpoptCalculatedQuantities* ip_cq)
{
    return check_no_cancel(_user_data);
}

// returns the size of the problem
bool IPOPT_PSOPT::get_nlp_info(Index& n, Index& m, Index& nnz_jac_g,
                             Index& nnz_h_lag, IndexStyleEnum& index_style)
{
  int nnz;
  int nnzA;
  int nnzG;
  int i;
  double jsratio;


  // Number of variables
  n = workspace->nvars;

  // Number of constraints in g(x)
  m = workspace->ncons;

  MatrixXd *X0 = workspace->x0;

  double  *x  = &(*X0)(0);

  if( !useAutomaticDifferentiation(*workspace->algorithm) ) {


     DetectJacobianSparsity(gg_num, *X0, m,  &nnzA,  workspace->iArow, workspace->jAcol, workspace->jac_Aij,
                                             &nnzG,  workspace->iGrow, workspace->jGcol,
                                             workspace->grw, workspace );

     nnz = nnzA+nnzG;

     jsratio = (double) ((double)  nnz/((double) (n*m)));

     if (jsratio > workspace->algorithm->jac_sparsity_ratio)
     {
           sprintf(workspace->text, "increase algorithm.jac_sparsity_ratio to just above %f", jsratio);
           error_message(workspace->text);
     }

     sprintf(workspace->text,"\nJacobian sparsity detected numerically:");
     psopt_print(workspace,workspace->text);
     sprintf(workspace->text,"\n*** %i nonzero elements out of %i [ratio=%f]", nnz, n*m, jsratio );
     psopt_print(workspace,workspace->text);
     sprintf(workspace->text,"\n*** %i nonzero elements are constant", nnzA );
     psopt_print(workspace,workspace->text);
     sprintf(workspace->text,"\n*** %i nonzero elements are not constant", nnzG );
     psopt_print(workspace,workspace->text);


  }


  if( useAutomaticDifferentiation(*workspace->algorithm) ) {

	unsigned int *jac_rind  = NULL;
	unsigned int *jac_cind  = NULL;
	double       *jac_values = NULL;

	adouble *xad = workspace->xad;
	adouble *gad = workspace->gad;
	double  *g   = workspace->fg;

	/* Tracing of function gg() */
	trace_on(workspace->tag_g);
	for(i=0;i<n;i++)
		xad[i] <<= x[i];

	gg_ad(xad, gad, workspace);

	for(i=0;i<m;i++)
		gad[i] >>= g[i];
	trace_off();


	/* Entries in row-compressed format using sparse_jac: */


// ASSUMING ADOL-C - VERSION 2.X.X
    int options[4];
    options[0]=0; options[1]=0; options[2]=0;options[3]=0;
    sparse_jac(workspace->tag_g, m, n, 0, x, &nnz, &jac_rind, &jac_cind, &jac_values, options);



	for(i=0;i<nnz;i++)
	{
		workspace->jGcol[i] = jac_cind[i];
		workspace->iGrow[i] = jac_rind[i];
	}

  ///Bug fix
  free(jac_rind); free(jac_cind); free(jac_values);

        sprintf(workspace->text,"\nJacobian sparsity detected using ADOLC:");
        psopt_print(workspace,workspace->text);

        jsratio = (double) ((double)  nnz/((double) (n*m)));

        if (jsratio > workspace->algorithm->jac_sparsity_ratio) {
           sprintf(workspace->text, "increase algorithm.jac_sparsity_ratio to just above %f", jsratio);
           error_message(workspace->text);
        }

        sprintf(workspace->text,"\n%i nonzero elements out of %i [ratio=%f]\n", nnz, n*m, jsratio);
        psopt_print(workspace,workspace->text);

  } // end if (autoderiv)

  int activate_hess;

  if (workspace->algorithm->hessian=="exact")
      activate_hess = 1;
  else
      activate_hess = 0;

  if( activate_hess*useAutomaticDifferentiation(*workspace->algorithm)  ) {

	double       *hess_values = NULL;
	adouble *xad = workspace->xad;
	adouble Lad;
	double  obj_factor = 1.0;

   double *lambda = &(*workspace->lambda)(0);
	double  L;
   int nnz_hess;
	/* Tracing of function Lagrangian_ad() */
	trace_on(workspace->tag_hess);
	for(i=0;i<n;i++)
		xad[i] <<= x[i];
	Lad = Lagrangian_ad(xad, lambda, obj_factor, m, workspace);
        Lad >>=L;
	trace_off();
	/* Entries in row-compressed format using sparse_hess: */

	unsigned int *hess_ir = NULL;
 	unsigned int *hess_jc = NULL;



// ASSUMING ADOL-C - VERSION 2.X.X
    int options[2];
    options[0]=1; options[1]=0;
    sparse_hess(workspace->tag_hess, n,0,x,&nnz_hess,&hess_ir, &hess_jc,&hess_values, options);


    for (i=0; i< nnz_hess; i++) {
		workspace->hess_ir[i] = hess_ir[i];
		workspace->hess_jc[i] = hess_jc[i];
    }

       sprintf(workspace->text,"\nHessian sparsity detected using ADOLC:");
       psopt_print(workspace,workspace->text);
       double hsratio = (double) ((double)  nnz_hess/((double) (n*n)));
       if (hsratio > workspace->algorithm->hess_sparsity_ratio) {
            sprintf(workspace->text, "increase algorithm.hess_sparsity_ratio to just above %f", hsratio);
            error_message(workspace->text);
       }

       sprintf(workspace->text,"\n%i nonzero elements out of %i [ratio = %f] \n", nnz_hess, n*n, hsratio );
       psopt_print(workspace,workspace->text);

       nnz_h_lag = nnz_hess;

  } // end if (autoderiv)

    nnz_jac_g = nnz;


  // the hessian is in this case assumed to be a square dense matrix but we
  // only need the lower left corner (since it is symmetric)
  if( !useAutomaticDifferentiation(*workspace->algorithm) || workspace->algorithm->hessian!="exact" )
        nnz_h_lag = (int) ((n*n)+n)/2;

/*   *
     * *
     * * *
     * * * *
     * * * * *
*/

  // use the C style indexing (0-based)
  index_style = TNLP::C_STYLE;

  return true;
}

// returns the variable bounds
bool IPOPT_PSOPT::get_bounds_info(Index n, Number* x_l, Number* x_u,
                                Index m, Number* g_l, Number* g_u)
{

  // here, the n and m we gave IPOPT in get_nlp_info are passed back to us.
  // If desired, we could assert to make sure they are what we think they are.
  assert(n == workspace->nvars);
  assert(m == workspace->ncons);


    double *xlb = &(*workspace->xlb)(0);

    double *xub = &(*workspace->xub)(0);


  int  j;


  // lower bounds on x
  for (j=0; j<workspace->nvars; j++) {
    x_l[j] = xlb[j];
  }

  // upper bounds on x
  for (j=0; j<workspace->nvars; j++) {
    x_u[j] = xub[j];
  }

  get_constraint_bounds(g_l, g_u, workspace);

  return true;
}

// returns the initial point for the problem
bool IPOPT_PSOPT::get_starting_point(Index n, bool init_x, Number* x,
                                   bool init_z, Number* z_L, Number* z_U,
                                   Index m, bool init_lambda,
                                   Number* lambda)
{
  // Here, we assume we only have starting values for x, if you code
  // your own NLP, you can provide starting values for the dual variables
  // if you wish
  assert(init_x == true);
//  assert(init_z == false);
//  assert(init_lambda == false);

  Index i;

//double *x0 = (workspace->x0)->GetPr();
  double *x0 = &(*workspace->x0)(0);

  // initialize to the given starting point

  for (i=0; i<workspace->nvars;i++)
  {
	  x[i] = x0[i];
  }


  return true;
}

// returns the value of the objective function
bool IPOPT_PSOPT::eval_f(Index n, const Number* x, bool new_x, Number& obj_value)
{
  assert(n == workspace->nvars);

  MatrixXd& X = *workspace->Xip;

  memcpy( &X(0), x, workspace->nvars*sizeof(double) );

  obj_value = ff_num(X, workspace);

  return true;
}

// return the gradient of the objective function grad_{x} f(x)
bool IPOPT_PSOPT::eval_grad_f(Index n, const Number* x, bool new_x, Number* grad_f)
{
  assert(n == workspace->nvars);

  MatrixXd& X = *workspace->Xip;

  MatrixXd& GF = *workspace->GFip;


  memcpy( &X(0), x, workspace->nvars*sizeof(double) );

  if(!useAutomaticDifferentiation(*workspace->algorithm))
     ScalarGradient( ff_num, X, &GF , workspace->grw, workspace );
  else
     ScalarGradientAD( ff_ad, X, &GF, &workspace->trace_f_done, workspace->tag_f, workspace );

  memcpy( grad_f, &GF(0), workspace->nvars*sizeof(double));

  return true;
}

// return the value of the constraints: g(x)
bool IPOPT_PSOPT::eval_g(Index n, const Number* x, bool new_x, Index m, Number* g)
{
  assert(n == workspace->nvars);
  assert(m == workspace->ncons);

  MatrixXd& X = *workspace->Xip;

  MatrixXd& G  = *workspace->Gip;

  memcpy( &X(0), x, workspace->nvars*sizeof(double) );

  gg_num(X, &G, workspace);

  memcpy( g, &G(0), workspace->ncons*sizeof(double) );

  return true;
}




void save_jacobian_sparsity_pattern(Index* rindex, Index* cindex, long nvars, long ncols, long nnz, Workspace* workspace)
{
    // Saves the jacobian structure in triplet form with 1 at the non-zeros.

    FILE* jac_file;

    jac_file = fopen("jacobian_pattern.dat", "w");

    if (jac_file == NULL) {
         error_message("save_jacobian_sparsity_pattern(): error opening jacobian pattern file");
    }

    double one  = 1.0;
    double zero = 0.0;

    fprintf(jac_file,"% li %li  %f", ncols, nvars, zero );

    for (int i = 1; i< nnz; i++ ) {
        fprintf(jac_file,"\n%li %li  %f", rindex[i]+1, cindex[i]+1, one ); // SPARSITY PATTERN SAVED USING 1-BASED INDICES (MATLAB-STYLE)
    }

    fclose(jac_file);

    workspace->algorithm->save_sparsity_pattern = 0;

}


// THE UPDATE OF THIS FUNCTION NEEDS TO BE FINALISED - EIGEN_UPDATE.

// return the structure or values of the jacobian
bool IPOPT_PSOPT::eval_jac_g(Index n, const Number* x, bool new_x,
                           Index m, Index nele_jac, Index* iRow, Index *jCol,
                           Number* values)
{
  MatrixXd& X = *workspace->Xip;

  double *xpr = &(*workspace->xp)(0); 

  int nnzA, nnzG, i;


  if (values == NULL) {
  // return the structure of the jacobian
    X = *workspace->x0;
    if (!useAutomaticDifferentiation(*workspace->algorithm))
    {


        	nnzA = workspace->jac_nnzA;
        	nnzG = workspace->jac_nnzG;


        	getIndexGroups( workspace->igroup, m, n, nnzG, workspace->iGrow, workspace->jGcol, workspace);

	     	for (i=0;i<nnzG;i++)
      	{
				iRow[i] = workspace->iGrow[i]; // EIGEN_UPDATE
				jCol[i] = workspace->jGcol[i]; // EIGEN_UPDATE
			}

			for (i=0;i<nnzA;i++)
			{
				iRow[i+nnzG] =workspace->iArow[i]; // EIGEN_UPDATE
				jCol[i+nnzG] =workspace->jAcol[i]; // EIGEN_UPDATE
			}

    }


    if (useAutomaticDifferentiation(*workspace->algorithm)) {


			for(i=0;i<nele_jac;i++)
			{
				iRow[i] = workspace->iGrow[i];
				jCol[i] = workspace->jGcol[i];
			}

    } // End if (autoderiv)
    
  }
  else {
    // return the values of the jacobian of the constraints
    if (!useAutomaticDifferentiation(*workspace->algorithm)) {

          memcpy( &X(0), x, workspace->nvars*sizeof(double) );
          
          // Compute by sparse finite differences only the non-constant Jacobian elements...
          EfficientlyComputeJacobianNonZeros(gg_num, X, m, values, workspace->jac_nnzG, workspace->iGrow,workspace->jGcol, workspace->igroup, workspace->grw, workspace );

          // Now include in array values[] the constant Jacobian elements calculated previously
          for(i=0;i<workspace->jac_nnzA;i++) {
                  values[workspace->jac_nnzG+i] = workspace->jac_Aij[i];
          }

    }

    /* find the jacobian values */

    if (useAutomaticDifferentiation(*workspace->algorithm)) {

    	int nnz = nele_jac;
    	unsigned int *jac_rind = NULL;
    	unsigned int *jac_cind = NULL;



	   for (i=0;i<n;i++) {
		   xpr[i] = x[i];
	   }

      double* jac_values = NULL;



// ASSUMING ADOL-C - VERSION 2.X.X
       int options[4];
       options[0]=0; options[1]=0; options[2]=0; options[3]=0;
	    sparse_jac(workspace->tag_g, m, n, 0, xpr, &nnz, &jac_rind, &jac_cind, &jac_values, options);


   	 for(i=0;i<nnz;i++) {

                values[i] = jac_values[i];
	    }

      ///Bug fix
      free(jac_rind); free(jac_cind); free(jac_values);



    }

    if (workspace->enable_nlp_counters) {
        workspace->solution->mesh_stats[ workspace->current_mesh_refinement_iteration-1 ].n_jacobian_evals++;
    }

  }

  if (workspace->algorithm->save_sparsity_pattern == 1) {
     save_jacobian_sparsity_pattern(iRow , jCol, n, m, nele_jac, workspace );
  }

  return true;
}


void IPOPT_PSOPT::finalize_solution(SolverReturn status,
                                  Index n, const Number* x, const Number* z_L, const Number* z_U,
                                  Index m, const Number* g, const Number* lambda,
                                  Number obj_value,
				  const IpoptData* ip_data,
				  IpoptCalculatedQuantities* ip_cq)
{
  // here is where we would store the solution to variables, or write to a file, etc
  // so we could use the solution.

  Sol* solution = workspace->solution;


    memcpy( &(*workspace->x0)(0), x, n*sizeof(double) );


    memcpy( &(*workspace->lambda)(0), lambda, m*sizeof(double) );

  for(int ii=0;ii<n;ii++) solution->xad[ii]=x[ii];

}






//return the structure or values of the hessian
bool IPOPT_PSOPT::eval_h(Index n, const Number* x, bool new_x,
                       Number obj_factor, Index m, const Number* lambda,
                       bool new_lambda, Index nele_hess, Index* iRow,
                       Index* jCol, Number* values)
{

  int i;

  if (workspace->algorithm->hessian!="exact")
    return false;

 if (!useAutomaticDifferentiation(*workspace->algorithm) ) return false;

  if (values == NULL) {
    if (useAutomaticDifferentiation(*workspace->algorithm)) {

	for(i=0;i<nele_hess;i++)
	{
		iRow[i] = workspace->hess_ir[i];
		jCol[i] = workspace->hess_jc[i];
	}

    } // End if
  }
  else {
    // return the values of the Hessian

    if (useAutomaticDifferentiation(*workspace->algorithm) && nele_hess>0) {
//    	double *xpr = workspace->Xsnopt->GetPr();
       double *xpr = &(*workspace->Xsnopt)(0);


// *******************************************************************
	adouble *xad = workspace->xad;
	adouble Lad;
	double  obj_factor_d = obj_factor;
	double*  lambda_d     = workspace->lambda_d;
	double  L;
	/* Tracing of Lagrangian function. It needs to be repeated because obj_factor and lambda change  */
	trace_on(workspace->tag_hess);
	for(i=0;i<n;i++)
		xad[i] <<= x[i];

	for(i=0;i<m;i++)
		lambda_d[i] = lambda[i];

	Lad = Lagrangian_ad(xad, lambda_d, obj_factor_d, m, workspace);
        Lad >>=L;
	trace_off();
// *******************************************************************

	for (i=0;i<n;i++) {
		xpr[i] = x[i];
	}

        unsigned int* hess_ir = NULL;
        unsigned int* hess_jc = NULL;

        double* hess_values = NULL;



// ASSUMING ADOL-C - VERSION 2.X.X
    int options[2];
    options[0]=1; options[1]=0;
    sparse_hess(workspace->tag_hess, n, 0, xpr, &nele_hess, &hess_ir, &hess_jc, &hess_values, options);


        for(i=0;i<nele_hess;i++) {
             values[i] = hess_values[i];
        }

	if (workspace->enable_nlp_counters) {
	    workspace->solution->mesh_stats[ workspace->current_mesh_refinement_iteration-1 ].n_hessian_evals++;
	}

    }

  }

  return true;

}



