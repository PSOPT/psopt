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

// Numerical Gradient Functions


void JacobianColumn( void fun(MatrixXd& x, MatrixXd* f, Workspace* ), MatrixXd& x, MatrixXd& xlb, MatrixXd& xub, int jCol,
                MatrixXd* JacColumn, GRWORK* grw, Workspace* workspace )
{
  // Computes only one column of the Jacobian matrix

  int  j, k;
  double delj;
  double sqreps;
  double xs;
  long nf  = JacColumn->rows();

  MatrixXd *dfdx_j = grw->dfdx_j;
  MatrixXd *F1   = grw->F1;
  MatrixXd *F2   = grw->F2;
  MatrixXd *F3   = grw->F3;
  int nvar = x.rows();


  dfdx_j->resize( nf, 1 );
  F1->resize(   nf, 1 );
  F2->resize(   nf, 1 );

  sqreps=sqrt( PSOPT_extras::GetEPS() );
  
  int tcount=0;
  
  for (int jj=0; jj<nvar; jj++){

        if(  x(jj) <= (xub(jj)-sqreps)||   x(jj)>=(xlb(jj)+sqreps) ) 
           tcount++; 
  }

  if ( tcount )
  {
     F3->resize( nf, 1);
     fun(x, F3, workspace);
  }

  j = jCol;
      delj = sqreps*(1.+fabs(x(j)));
//      delj = sqreps;
      xs   = x(j);
      if ((xs < xub(j)-delj && xs>xlb(j)+delj) || (xub(j)==xlb(j))) {
        // Use central difference formula
      	x(j)+= delj;
      	fun( x, F1, workspace );
      	x(j) = xs-delj;
      	fun( x, F2, workspace );
      	*dfdx_j=( (*F1)-(*F2) )/(2*delj);
      	x(j)=xs;
      }
      else if (xs >= xub(j)-delj) {
        // Variable at upper bound, use backward difference formula
      	x(j)= xs - delj;
      	fun( x, F1, workspace );
        x(j)= xs;
      	*dfdx_j=( (*F1)-(*F3) )/(-delj);
      }
      else if (xs <= xlb(j)+delj) {
        // Variable at lower bound, use forward difference formula
      	x(j)= xs + delj;
      	fun( x, F1, workspace );
        x(j)= xs;
      	*dfdx_j=( (*F1)-(*F3) )/(delj);
      }
//      for(k=0;k<nf;k++) { (*JacColumn)(k,0)=(*dfdx_j)(k,0); } //EIGEN_UPDATE
      *JacColumn = *dfdx_j;
}




void JacobianRow( void fun(MatrixXd& x, MatrixXd* f, Workspace* ), MatrixXd& x, int iRow, int nf,
                  MatrixXd* JacRow, GRWORK* grw, Workspace* workspace )
{
  // Computes only one row of the Jacobian matrix

  int  j;
  double delj;
  MatrixXd sqreps;
  double xs;
  long nvar= x.rows();

  MatrixXd *dfdx_j = grw->dfdx_j;
  MatrixXd *F1   = grw->F1;
  MatrixXd *F2   = grw->F2;
  MatrixXd *F3   = grw->F3;

  MatrixXd& xlb    = *workspace->xlb;
  MatrixXd& xub    = *workspace->xub;


  dfdx_j->resize( nf, 1 );
  F1->resize(   nf, 1 );
  F2->resize(   nf, 1 );

  sqreps=sqrt( PSOPT_extras::GetEPS() )*ones(nvar,1);
  
  int tcount=0;
  
  for (int jj=0; jj<nvar; jj++){

        if(  x(jj) <= (xub(jj)-sqreps(jj))||   x(jj)>=(xlb(jj)+sqreps(jj)) ) 
           tcount++; 
  }

  if ( tcount )
  {
     F3->resize( nf, 1);
     fun(x, F3, workspace);
  }

  for(j=0;j<nvar;j++) { // EIGEN_UPDATE
      delj = sqreps(0);
      xs   = x(j);
      if ((xs < xub(j)-delj && xs>xlb(j)+delj) || (xub(j)==xlb(j))) {
        // Use central difference formula
      	x(j)+= delj;
      	fun( x, F1, workspace );
      	x(j) = xs-delj;
      	fun( x, F2, workspace );
      	*dfdx_j=( (*F1)-(*F2) )/(2*delj);
      	x(j)=xs;
      }
      else if (xs >= xub(j)-delj) {
        // Variable at upper bound, use backward difference formula
      	x(j)= xs - delj;
      	fun( x, F1, workspace );
        x(j)= xs;
      	*dfdx_j=( (*F1)-(*F3) )/(-delj);
      }
      else if (xs <= xlb(j)+delj) {
        // Variable at lower bound, use forward difference formula
      	x(j)= xs + delj;
      	fun( x, F1, workspace );
        x(j)= xs;
      	*dfdx_j=( (*F1)-(*F3) )/(delj);
      }
      (*JacRow)(1,j)=(*dfdx_j)(iRow,1);
  }
}



void ComputeJacobianNonZeros( void fun(MatrixXd& x, MatrixXd* f ), MatrixXd& x,
                int nf, double *nzvalue, int nnz, int* iArow, int* jAcol, GRWORK* grw, Workspace* workspace )
{

  int nvar, I, i;
  double delj;
  MatrixXd sqreps;
  double xs = 0.0;
  nvar= x.rows(); // EIGEN_UPDATE
  int iflag=1;



  MatrixXd *F1   = grw->F1;
  MatrixXd *F2   = grw->F2;
  MatrixXd *F3   = grw->F3;

  MatrixXd& xlb    = *workspace->xlb;
  MatrixXd& xub    = *workspace->xub;

  F1->resize(   nf, 1 );
  F2->resize(   nf, 1 );

  sqreps=sqrt( PSOPT_extras::GetEPS() )*ones(nvar,1);
  
  int tcount=0;
  
  for (int jj=0; jj<nvar; jj++){

        if(  x(jj) <= (xub(jj)-sqreps(jj))||   x(jj)>=(xlb(jj)+sqreps(jj)) ) 
           tcount++; 
  }



  if ( tcount )
  {
     F3->resize( nf, 1);
     fun(x, F3);
  }


  for (i=0;i<nvar;i++)  // EIGEN_UPDATE: Index i shifted by -1.
  {
     for(I=0; I<nnz; I++)
     {
         if( jAcol[I]==i ) {
            if (iflag) {
                      delj = sqreps(0);
                      xs   = x(i);
                      if (xs< xub(i)-delj  || (xub(i)==xlb(i)))
                      {
		          x(i) += delj;
      		          fun( x, F1 );
                      }
                      if (xs> xlb(i)+delj || (xub(i)==xlb(i)))
                      {
                         x(i)  = xs-delj;
                         fun( x, F2 );
                      }
                      x(i) = xs;
                      iflag=0;
            }
            if ((xs< (xub(i)-delj) && (xs> xlb(i)+delj))  || (xub(i)==xlb(i))) {
            	// Use central difference formula
            	nzvalue[I] = ((*F1)(iArow[I]) - (*F2)(iArow[I]))/(2*delj);
            }
            else if (xs>= (xub(i)-delj) ) {
            	// Use backward difference formula
            	nzvalue[I] = ((*F2)(iArow[I]) - (*F3)(iArow[I]))/(-delj);
            }
            else if (xs<= (xlb(i)+delj) ) {
            	// Use forward difference formula
            	nzvalue[I] = ((*F1)(iArow[I]) - (*F3)(iArow[I]))/(delj);
            }

         }
     }
     iflag=1;

  }
}

void deleteIndexGroups(IGroup* igroup, int ncols )
{
   int i;

   for(i=0;i< ncols;i++)
   {
         delete[] igroup->colindex[i];
   }

   delete[] igroup->size;

   delete[] igroup->colindex;
}


void getIndexGroups( IGroup* igroup, int nrows, int ncols, int nnz, int* iArow, int* jAcol, Workspace* workspace)
{
/* This function uses the method of Curtis, Powell and Reid (1974) to find groups of variables
 * to evaluate efficiently the sparse Jacobian by perturbing simultaneously groups of variables.
 * Reference:
 * A. R. Curtis, M.J.D. Powell and J.K. Reid
 * "On the estimation of Sparse Jacobian Matrices"
 * J Inst Maths Applics (1974) 13, 117-119
 *
 */

   int i, j, l, q, r;
   int group_index;
   MatrixXd& C1 = *workspace->JacCol1;
   MatrixXd& C2 = *workspace->JacCol2;

   double dotCols;

   // Form dummy Jacobian matrix with ones at the non-zero elements

   double* ones_pr = workspace->jac_Gij;


   for(i=0;i<nnz;i++)  ones_pr[i] = 1.0;

   TripletSparseMatrix J(ones_pr, nrows, ncols, nnz, iArow, jAcol );

//   J.Save("J.dat");

//   J.SaveSparsityPattern("pattern.txt");

   // Now allocate the pointer to the groups

   igroup->colindex = new int*[ncols];
   int* col_done  = new int[ncols];


   for(i=0;i<ncols;i++) col_done[i]=0;

   igroup->size = new int[ncols];

   for(i=0;i< ncols;i++) {
      igroup->colindex[i] = new int[ncols];
      igroup->size[i]     = 0;
   }


// To form the first group we inspect the columns
//  in turn and include each that has no unknowns in common with those columns already
//  included.
   // Add the first column to the first group
//   igroup->colindex[0][0] = 1;
   igroup->colindex[0][0] = 0; // EIGEN_UPDATE: The first group has index 0.
   col_done[0]=1;
   int gcount   = 1;
   int colcount = 1;


   // Now form the first group
   bool ok;
   for(j=1;j<ncols;j++) {  // EIGEN_UPDATE index j shifted by -1.
          ok = true;
	  for(l=0;l<gcount;l++) {
          	if (j== igroup->colindex[0][l]) {
                        ok=false;
               		break;
                }
                C1 = J.col(igroup->colindex[0][l]);
                C2 = J.col(j);
                  dotCols = (C1.transpose()*C2)(0);

          	if  ( dotCols>0.0 ) {
                        ok=false;
          		break;
		}
          }
  	  if (ok) {
                	igroup->colindex[0][gcount]=j;
                	gcount++;
                	colcount++;
//                        col_done[j-1]=1; //EIGEN_UPDATE
                          col_done[j]=1;
          }

   }


   igroup->size[0]= gcount;

// Now form the remaining groups.
//The other groups are formed successively by applying the same procedure to
//those columns not already included in a group.


   group_index=0;
   while (colcount<ncols) {
        group_index++;
        gcount = 0;
   	for(j=1;j<ncols;j++) { // EIGEN_UPDATE: index j shifted by -1.
                  ok=true;
	 	  for(q = 0; q< group_index; q++)
                  {
			for(r=0; r< igroup->size[q]; r++) {

				if(j==igroup->colindex[q][r]) {
                                	ok=false;
                                        break;
                                }
			}
			if (ok==false)
				break;
                  }
		  for(l=0;l<gcount;l++) {

          		if (j== igroup->colindex[group_index][l]) {
                                ok=false;
               			break;
                        }
                        C1 = J.col(igroup->colindex[group_index][l]);
                        C2 = J.col(j);

                        dotCols = (C1.adjoint()*C2)(0);
          	  	if ( dotCols>0.0 ) {
                                ok = false;
                                break;
                  	}
                  }
                  if (ok) {
                	igroup->colindex[group_index][gcount]=j;
                	gcount++;
                        colcount++;
                          col_done[j] = 1;
          	  }
   	}



        igroup->size[group_index] = gcount;

    }



    for(i=0;i<ncols;i++)
    {
         if (igroup->size[i]==0)
         {
		igroup->number=i;
                break;
         }
    }



    sprintf(workspace->text,"\nNumber of index sets for sparse finite differences = %i\n", igroup->number);
    psopt_print(workspace,workspace->text);



}

void EfficientlyComputeJacobianNonZeros( void fun(MatrixXd& x, MatrixXd* f, Workspace* ), MatrixXd& x, int nf,
            double *nzvalue, int nnz, int* iArow, int* jAcol, IGroup* igroup, GRWORK* grw, Workspace* workspace )
{
/* This function uses the method of Curtis, Powell and Reid (1974) to
 * evaluate efficiently the sparse Jacobian by perturbing simultaneously groups of variables.
 * Reference:
 * A. R. Curtis, M.J.D. Powell and J.K. Reid
 * "On the estimation of Sparse Jacobian Matrices"
 * J Inst Maths Applics (1974) 13, 117-119
 *
 */

  int  j, k, i;
  double delj;
  double sqreps;
  long nvar= x.rows();

  MatrixXd *F1   = grw->F1;
  MatrixXd *F2   = grw->F2;


  MatrixXd xp(nvar,1);

  F1->resize(   nf, 1 );
  F2->resize(   nf, 1 );

  sqreps = sqrt( PSOPT_extras::GetEPS() );

  for (i=0;i<igroup->number;i++)
  {
        xp=x;
	for(j=0; j<igroup->size[i]; j++) {
                      delj = sqreps;
		      xp(igroup->colindex[i][j]) += delj;
        }
        fun( xp, F1, workspace );
        for(j=0; j<igroup->size[i]; j++) {
                      delj=sqreps;
                      xp(igroup->colindex[i][j]) -= 2*delj;
        }
        fun( xp, F2, workspace );
        for(j=0; j<igroup->size[i]; j++) {
             for(k=0;k<nnz;k++) {
              	if (jAcol[k] == igroup->colindex[i][j] )
                {
                     nzvalue[k] = ((*F1)(iArow[k]) - (*F2)(iArow[k]))/(2*delj);
                }
             }
        }
   }

}


void DetectJacobianSparsity(void fun(MatrixXd& x, MatrixXd* f, Workspace* ), MatrixXd& x, int nf,
                           int* nnzA, int* iArow, int* jAcol, double* Aij,
                           int* nnzG, int* jGrow, int* jGcol,
                           GRWORK* grw, Workspace* workspace)
{



  long nvars = x.rows();
  long i,j;
  int nzcount_A=0;
  int nzcount_G=0;

  double s = 1.0e6*sqrt(PSOPT_extras::GetEPS());
  double tol  = 1.e-16*pow( PSOPT_extras::GetEPS(), 0.8)* std::max( 1.0, x.norm() );




  MatrixXd& JacCol1 = *workspace->JacCol1;
  MatrixXd& JacCol2 = *workspace->JacCol2;
  MatrixXd& JacCol3 = *workspace->JacCol3;
  MatrixXd& xp      = *workspace->xp;
  MatrixXd& xlb     = *workspace->xlb;
  MatrixXd& xub     = *workspace->xub;

  for(j=0;j<nvars;j++) {   // EIGEN_UPDATE: j index shifted by -1.

     xp = x;

     clip_vector_given_bounds( xp, xlb, xub);

     JacobianColumn( fun, xp, xlb, xub, j, &JacCol1,  grw, workspace);
     xp = x + 0.1*x.cwiseAbs() + s*ones(nvars,1);

     clip_vector_given_bounds( xp, xlb, xub);

     JacobianColumn( fun, xp, xlb, xub, j, &JacCol2,  grw, workspace);
     xp = x - 0.15*x.cwiseAbs() - 1.1*s*ones(nvars,1);

     clip_vector_given_bounds( xp, xlb, xub);

     JacobianColumn( fun, xp, xlb, xub,j, &JacCol3, grw, workspace);



      for(i=0; i<nf; i++) { // EIGEN_UPDATE: index i shifted by -1
            if ( ( fabs(JacCol1(i,0)) +  fabs(JacCol2(i,0)) + fabs(JacCol3(i,0)) )>=tol ) {
         	if ( fabs(JacCol1(i,0)-JacCol2(i,0))==0.0 && fabs(JacCol1(i,0)-JacCol3(i,0))==0.0 ) {
                        // Constant Jacobian element detected
              		      iArow[nzcount_A]=i;
              		      jAcol[nzcount_A]=j;
                        Aij[nzcount_A]    = JacCol1(i,0);

                        nzcount_A++;
              }
              else {
		       // Non-constant Jacobian element
			               jGrow[nzcount_G]= i;
                        jGcol[nzcount_G]= j;
                        nzcount_G++;
              }
           }
      }

  }
  *nnzA=nzcount_A;
  *nnzG=nzcount_G;

  workspace->jac_nnz  = nzcount_A + nzcount_G;
  workspace->jac_nnzA = nzcount_A;
  workspace->jac_nnzG = nzcount_G;

}


void DetectJacobianSparsityAD(void fun(MatrixXd& x, MatrixXd* f, Workspace* ), MatrixXd& x, int nf,
                           int* nnzA, int* iArow, int* jAcol, double* Aij,
                           int* nnzG, int* jGrow, int* jGcol,
                           GRWORK* grw, Workspace* workspace)
{


  long i,j;
  int nzcount_A=0;
  int nzcount_G=0;
  int n     =  length(x);
  int neF   =  nf;  
  
  double s = 1.0e6*sqrt(PSOPT_extras::GetEPS());


  MatrixXd& xp      = *workspace->xp;
  MatrixXd& xlb     = *workspace->xlb;
  MatrixXd& xub     = *workspace->xub;


        int nF = neF;

        int nvars = n;

        int k;

       xp = x;

       clip_vector_given_bounds( xp, xlb, xub);

       // Compute the full Jacobian using ADOL-C:  J = G(x)+ A
       int options[4];
       int repeat = 1; // To use previously determined sparsity pattern.
       options[0]=0; options[1]=0; options[2]=0;options[3]=0; 
	    sparse_jac(workspace->tag_fg, nF, nvars, repeat, &xp(0), &workspace->F_nnz, &workspace->iGfun2, &workspace->jGvar2, &workspace->G2, options);

       xp = x + 0.05*x.cwiseAbs() + s*ones(nvars,1);
       clip_vector_given_bounds( xp, xlb, xub);

       // Compute the full Jacobian using ADOL-C:  J = G(x)+ A

	    sparse_jac(workspace->tag_fg, nF, nvars, repeat, &xp(0), &workspace->F_nnz, &workspace->iGfun2, &workspace->jGvar2, &workspace->G3, options);

       xp = x - 0.06*x.cwiseAbs() - 0.95*s*ones(nvars,1);

       clip_vector_given_bounds( xp, xlb, xub);
       // Compute the full Jacobian using ADOL-C:  J = G(x)+ A

	    sparse_jac(workspace->tag_fg, nF, nvars, repeat, &xp(0), &workspace->F_nnz, &workspace->iGfun2, &workspace->jGvar2, &workspace->G4, options);


        if (workspace->enable_nlp_counters) {
             workspace->solution->mesh_stats[ workspace->current_mesh_refinement_iteration-1 ].n_jacobian_evals+=3;
        }


  

    for(i=0; i<workspace->F_nnz; i++) { // EIGEN_UPDATE: index i shifted by -1
//            if ( ( fabs(workspace->G2[i]) +  fabs(workspace->G3[i]) + fabs(workspace->G4[i])>=tol ) ) {
//              if ( fabs(workspace->G2[i] - workspace->G3[i])<=tol && fabs(workspace->G2[i] - workspace->G4[i])<=tol ) {
             if ( (workspace->G2[i] == workspace->G3[i] && workspace->G2[i] == workspace->G4[i]) ) {
                        // Constant Jacobian element detected
              		      iArow[nzcount_A]=workspace->iGfun2[i];
              		      jAcol[nzcount_A]=workspace->jGvar2[i];
                        Aij[nzcount_A]    = workspace->G2[i];
                        nzcount_A++;
              }
              else {
		       // Non-constant Jacobian element
			               jGrow[nzcount_G]= workspace->iGfun2[i];
                        jGcol[nzcount_G]= workspace->jGvar2[i];
                        nzcount_G++;
              }
//           }
      }

  
  *nnzA=nzcount_A;
  *nnzG=nzcount_G;

  workspace->jac_nnz  = nzcount_A + nzcount_G;
  workspace->jac_nnzA = nzcount_A;
  workspace->jac_nnzG = nzcount_G;
  
    
//  cout << "\nneA = " << nzcount_A;
//  cout << "\n[iArow , jAcol]";
//  for(int iA=0; iA< nzcount_A; iA++ ) {
//    cout << "\n[" << iArow[iA] << " , " << jAcol[iA] << "]";    
// }

}




void ScalarGradient( double (*fun)(MatrixXd& x, Workspace* workspace), MatrixXd& x,
                MatrixXd* grad, GRWORK* grw, Workspace* workspace )
{

  int j = 0, nf;
  double delj;
  MatrixXd sqreps;
  double xs = 0.0;
  long nvar= x.rows();
  nf  = 1;
  double F1 = 0.0;
  double F2 = 0.0;
  double F3 = 0.0;
  double dfdx = 0.0;;

  MatrixXd& xlb    = *workspace->xlb;
  MatrixXd& xub    = *workspace->xub;

  sqreps=sqrt( PSOPT_extras::GetEPS() )*ones(nvar,1);
  
  int tcount=0;
  
  for (int jj=0; jj<nvar; jj++){

        if(  x(jj) <= (xub(jj)-sqreps(jj))||   x(jj)>=(xlb(jj)+sqreps(jj)) ) 
           tcount++; 
  }

  if (tcount)
  {
     F3 = fun(x, workspace);
  }

  for(j=0;j<nvar;j++) { // EIGEN_UPDATE: index j shifted by -1
      delj = sqreps(0)*(1.0+fabs(x(j)));
      xs   = x(j);
      if (xs< xub(j)-delj || (xub(j)==xlb(j)))
      {
      	x(j) += delj;
      	F1 = fun( x, workspace);
      }
      if (xs> xlb(j)+delj || (xub(j)==xlb(j)))
      {
      	x(j) = xs-delj;
      	F2 = fun( x, workspace );
      }
      if (( (xs< (xub(j)-delj)) && (xs> (xlb(j)+delj)) ) || (xub(j)==xlb(j)) ) {
        // Use central difference formula
      	dfdx = ( F1 - F2 )/(2*delj);
      }
      else if (xs>= (xub(j)-delj)) {
        // Variable at upper bound, use backward difference formula
      	dfdx = ( F2 - F3 )/(-delj);
      }
      else if (xs<= (xlb(j)+delj)) {
        // Variable at lower bound, use forward difference formula
      	dfdx = ( F1 - F3 )/(delj);
      }
      x(j) = xs;
      (*grad)(j) = dfdx;
  }


}

void ScalarGradientAD( adouble (*fun)(adouble *, Workspace*), MatrixXd& x, MatrixXd* grad, bool* trace_done, int itag, Workspace* workspace )
{
    // Compute the gradient of a scalar function using automatic differentiation
    int      n = x.rows();
    int i;
    double  yp = 0.0;
    adouble *xad = workspace->xad;
    adouble  yad;

    if( !(*trace_done) ) {
    	trace_on(itag);
    	for(i=0;i<n;i++) {
            xad[i] <<= (&x(0))[i];
	}
    	yad = (*fun)(xad, workspace);
    	yad >>= yp;
    	trace_off();
        *trace_done = true;
    }

    gradient(itag,n,&x(0),&(*grad)(0));

}


void compute_jacobian_of_constraints_with_respect_to_variables(MatrixXd& Jc, MatrixXd& X, MatrixXd& XL, MatrixXd& XU, Workspace* workspace)
{

    Alg& algorithm = *workspace->algorithm;
    Prob& problem  = *workspace->problem;

    int i, j, k;

    int nvars = get_number_nlp_vars(problem, workspace);

    int ncons = get_number_nlp_constraints(problem, workspace);

    MatrixXd Jctmp(ncons,nvars);

    Jc.resize(ncons,nvars);

    MatrixXd& JacCol1 = *workspace->JacCol1;
    MatrixXd& xlb     = *workspace->xlb;
    MatrixXd& xub     = *workspace->xub;
    MatrixXd& xp      = *workspace->xp;
    MatrixXd jtemp;

    workspace->use_constraint_scaling = 0;

     xp = X;
//     clip_vector_given_bounds( xp, xlb, xub);

     if ( useAutomaticDifferentiation(algorithm) ) {

    unsigned int *jac_rind  = NULL;
	unsigned int *jac_cind  = NULL;
	double       *jac_values = NULL;
	int           nnz;

	adouble *xad = workspace->xad;
	adouble *gad = workspace->gad;
	double  *g   = workspace->fg;

    double  *x   = &xp(0);

	/* Tracing of function gg() */
	trace_on(workspace->tag_gc);
	for(i=0;i<nvars;i++)
		xad[i] <<= x[i];

	gg_ad(xad, gad, workspace);

	for(i=0;i<ncons;i++)
		gad[i] >>= g[i];
	trace_off();


    int options[4];
    options[0]=0; options[1]=0; options[2]=0;options[3]=0;
	sparse_jac(workspace->tag_gc, ncons, nvars, 0, x, &nnz, &jac_rind, &jac_cind, &jac_values, options);


    for(j=0;j<nvars;j++)
    {
         Jctmp( jac_rind[j], jac_cind[j]) = jac_values[j];
    }


   }
   else {

    	MatrixXd& xlb = *(workspace->xlb);
	    MatrixXd& xub = *(workspace->xub);
	    for(j=0;j<nvars;j++) { // EIGEN_UPDATE: index j shifted by -1
	      JacobianColumn( gg_num, xp, xlb, xub,j, &JacCol1, workspace->grw, workspace);
          long nrows = JacCol1.rows();
          Jctmp.block(0,j,nrows,1)= JacCol1;
	    }

   } // end if-else


   MatrixXd& lambda = *workspace->lambda;

   int icount=0;  // EIGEN_UPDATE: icount starts from 0.

   int lam_phase_offset=0;

   int iphase;

   for(iphase=1;iphase<=problem.nphases;iphase++) {
       int ncons_phase_i =  get_ncons_phase_i(problem,iphase-1, workspace);
       for(j=0;j<ncons_phase_i;j++) {  // EIGEN_UPDATE: j index shifted by -1
           int nstates = problem.phases(iphase).nstates;
           int norder  = problem.phases(iphase).current_number_of_intervals;
           i = lam_phase_offset + j;
           if (j<= nstates*(norder+1)-1) {
           // Copy rows corresponding to differential defect constraints.
              for (k=0;k<nvars;k++) {  // EIGEN_UPDATE: j index shifted by -1
                    Jc(icount,k) = Jctmp(i,k);
              }
              icount++;
           }
           else if( j< ncons_phase_i-1) { // discard constraint tf>=t0 for each phase
              if ( lambda(i)!=0.0 ) {
                // Only copy rows corresponding to active inequality constraints.
                 for (k=0;k<nvars;k++) {  // EIGEN_UPDATE: index k shifted by -1
                    Jc(icount,k) = Jctmp(i,k);
                 }
                 icount++;
              }
           }
       }
       lam_phase_offset+= ncons_phase_i;
   }

     Jc = Jc.block(0,0, icount-1, Jc.cols() );

   workspace->use_constraint_scaling = 1;

}


void compute_jacobian_of_residual_vector_with_respect_to_variables(MatrixXd& Jr, MatrixXd& X, MatrixXd& XL, MatrixXd& XU, Workspace* workspace)
{
    MatrixXd Jcol;
	Prob & problem = *(workspace->problem);
    int nvar, nr, iphase, j;

	nvar = get_number_nlp_vars(problem, workspace);

	nr    = 0;

	for(iphase=1; iphase<=problem.nphases;iphase++)
	{
	  nr    += problem.phases(iphase).nobserved*problem.phases(iphase).nsamples;
	}

	Jr.resize(nr,nvar);
	Jcol.resize(nr,1);


	for(j=0;j<nvar;j++) {  // EIGEN_UPDATE: index j shifted by -1
	    JacobianColumn( rr_num, X, XL, XU, j, &Jcol, workspace->grw, workspace);
        Jr.block(0,j,Jcol.rows(), 1) = Jcol;
	}
}

