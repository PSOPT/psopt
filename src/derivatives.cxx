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

// Bring std names into this translation unit (formerly leaked via psopt.h).
using namespace std;

// Numerical Gradient Functions


void JacobianColumn( void fun(MatrixXd& x, MatrixXd* f, Workspace* ), MatrixXd& x, MatrixXd& xlb, MatrixXd& xub, int jCol,
                MatrixXd* JacColumn, GRWORK* grw, Workspace* workspace )
{
  // Computes only one column of the Jacobian matrix

  int  j;
  double delj;
  double sqreps;
  double xs;
  long nf  = JacColumn->rows();

  MatrixXd *dfdx_j = grw->dfdx_j.get();
  MatrixXd *F1   = grw->F1.get();
  MatrixXd *F2   = grw->F2.get();
  MatrixXd *F3   = grw->F3.get();
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

  MatrixXd *dfdx_j = grw->dfdx_j.get();
  MatrixXd *F1   = grw->F1.get();
  MatrixXd *F2   = grw->F2.get();
  MatrixXd *F3   = grw->F3.get();

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



  MatrixXd *F1   = grw->F1.get();
  MatrixXd *F2   = grw->F2.get();
  MatrixXd *F3   = grw->F3.get();

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

   double* ones_pr = workspace->jac_Gij.get();


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



    snprintf(workspace->text,sizeof(workspace->text),"\nNumber of index sets for sparse finite differences = %i\n", igroup->number);
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

  MatrixXd *F1   = grw->F1.get();
  MatrixXd *F2   = grw->F2.get();


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
                // Guard the numerical-path writes against the allocated buffer.
                // jac_nnz_capacity is set (>0) only for the workspace-owned buffers;
                // where a caller passes its own it is 0 and the guard is inactive.
                // error_message throws, aborting cleanly before any out-of-bounds
                // write.
                if ( workspace->jac_nnz_capacity > 0 &&
                     ( nzcount_A >= workspace->jac_nnz_capacity ||
                       nzcount_G >= workspace->jac_nnz_capacity ) ) {
                    // Grow the workspace buffers to the detected count and refresh
                    // the local pointers, so the buffer tracks the non-zero count
                    // instead of aborting. These buffers are workspace-owned; a
                    // caller-owned buffer has jac_nnz_capacity == 0 and is skipped.
                    long need = ((nzcount_A > nzcount_G) ? nzcount_A : nzcount_G) + 1;
                    psopt_grow_jacobian_buffers(workspace, need);
                    iArow = workspace->iArow.get();
                    jAcol = workspace->jAcol.get();
                    Aij   = workspace->jac_Aij.get();
                    jGrow = workspace->iGrow.get();
                    jGcol = workspace->jGcol.get();
                }
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


  long i;
  int nzcount_A=0;
  int nzcount_G=0;
  int n     =  length(x);
  
  double s = 1.0e6*sqrt(PSOPT_extras::GetEPS());


  MatrixXd& xp      = *workspace->xp;
  MatrixXd& xlb     = *workspace->xlb;
  MatrixXd& xub     = *workspace->xub;



        int nvars = n;


       xp = x;

       clip_vector_given_bounds( xp, xlb, xub);

       // Compute the full Jacobian using ADOL-C:  J = G(x)+ A
	    { psopt_ad::SparseTriplet Jg = psopt_ad::ad_sparse_jacobian(workspace->ad_fg, &xp(0), /*reuse=*/true);
	      workspace->iGfun2.assign(Jg.row.begin(), Jg.row.end()); workspace->jGvar2.assign(Jg.col.begin(), Jg.col.end());
	      workspace->G2.assign(Jg.val.begin(), Jg.val.end()); }

       xp = x + 0.05*x.cwiseAbs() + s*ones(nvars,1);
       clip_vector_given_bounds( xp, xlb, xub);

       // Compute the full Jacobian using ADOL-C:  J = G(x)+ A

	    { psopt_ad::SparseTriplet Jg = psopt_ad::ad_sparse_jacobian(workspace->ad_fg, &xp(0), /*reuse=*/true);
	      workspace->G3.assign(Jg.val.begin(), Jg.val.end()); }

       xp = x - 0.06*x.cwiseAbs() - 0.95*s*ones(nvars,1);

       clip_vector_given_bounds( xp, xlb, xub);
       // Compute the full Jacobian using ADOL-C:  J = G(x)+ A

	    { psopt_ad::SparseTriplet Jg = psopt_ad::ad_sparse_jacobian(workspace->ad_fg, &xp(0), /*reuse=*/true);
	      workspace->G4.assign(Jg.val.begin(), Jg.val.end()); }


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

  int j = 0;
  double delj;
  MatrixXd sqreps;
  double xs = 0.0;
  long nvar= x.rows();
  double F1 = 0.0;
  double F2 = 0.0;
  double F3 = 0.0;
  double dfdx = 0.0;;

  MatrixXd& xlb    = *workspace->xlb;
  MatrixXd& xub    = *workspace->xub;

  sqreps=sqrt( PSOPT_extras::GetEPS() )*ones(nvar,1);
  
  // The base point, needed by the one-sided formulas below. The guard this replaces
  // counted variables satisfying "x <= xub - sqrt(eps) OR x >= xlb + sqrt(eps)", which
  // is true unless a variable is somehow close to both ends of its range at once, so it
  // always fired and F3 was always computed. Reading as though the AND was meant, it was
  // one edit away from leaving the one-sided branches differencing against zero.
  F3 = fun(x, workspace);

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

  // A derivative checker, for when a run behaves as though its gradient were wrong.
  // Setting PSOPT_GRAD_VERIFY recomputes the same gradient by automatic
  // differentiation and reports the worst disagreement, relative to the largest
  // component. Finite differences should agree to something like 1e-8; a
  // disagreement far larger than that says the numerical gradient is not to be
  // trusted on this problem, which is worth knowing before drawing conclusions from
  // the solve. Costs nothing when the variable is unset.
  if (getenv("PSOPT_GRAD_VERIFY")) {
      MatrixXd gad(nvar,1);
      bool done = false;
      ScalarGradientAD( ff_ad, x, &gad, &done, workspace->ad_f, workspace );
      double worst = 0.0; long jw = -1; double gs = 0.0;
      for (long t = 0; t < nvar; t++) gs = max(gs, fabs(gad(t)));
      for (long t = 0; t < nvar; t++) {
          const double e = fabs((*grad)(t) - gad(t))/max(1.0e-30, gs);
          if (e > worst) { worst = e; jw = t; }
      }
      fprintf(stderr, "[grad-verify] n=%ld worst relative difference %.3e at j=%ld "
                      "(fd %.8e, ad %.8e)\n", (long) nvar, worst, jw,
                      jw >= 0 ? (*grad)(jw) : 0.0, jw >= 0 ? gad(jw) : 0.0);
  }

}

void ScalarGradientAD( adouble (*fun)(adouble *, Workspace*), MatrixXd& x, MatrixXd* grad, bool* trace_done, psopt_ad::ADHandle& adh, Workspace* workspace )
{
    int n = x.rows();
    if( !(*trace_done) ) {
        psopt_ad::ad_record(adh, n, 1, &x(0),
            [&](const adouble* xin, adouble* yout){ yout[0] = (*fun)(const_cast<adouble*>(xin), workspace); });
        *trace_done = true;
    }
    std::vector<double> g = psopt_ad::ad_gradient(adh, &x(0));
    for(int t=0;t<n;t++) (*grad)(t) = g[t];
}


void compute_jacobian_of_constraints_with_respect_to_variables(MatrixXd& Jc, MatrixXd& X, MatrixXd& XL, MatrixXd& XU, Workspace* workspace)
{

    Alg& algorithm = *workspace->algorithm;
    Prob& problem  = *workspace->problem;

    int i, j, k;

    int nvars = get_number_nlp_vars(problem, workspace);

    int ncons = get_number_nlp_constraints(problem, workspace);

    // Eigen leaves the storage of a newly constructed matrix uninitialised, and the
    // sparse automatic-differentiation branch below writes only the structurally nonzero
    // entries. Every other entry of Jctmp would then be whatever happened to be in memory.
    MatrixXd Jctmp = MatrixXd::Zero(ncons,nvars);

    // Room for the constraint rows plus one unit row per decision variable sitting at a
    // bound (appended below).
    Jc.resize(ncons+nvars,nvars);

    MatrixXd& JacCol1 = *workspace->JacCol1;
    MatrixXd& xp      = *workspace->xp;
    MatrixXd jtemp;

    workspace->use_constraint_scaling = 0;

     xp = X;
//     clip_vector_given_bounds( xp, xlb, xub);

     if ( useAutomaticDifferentiation(algorithm) ) {

    double  *x   = &xp(0);
	psopt_ad::ad_record(workspace->ad_gc, nvars, ncons, x,
		[&](const adouble* xin, adouble* yout){ gg_ad(const_cast<adouble*>(xin), yout, workspace); });
	psopt_ad::SparseTriplet J = psopt_ad::ad_sparse_jacobian(workspace->ad_gc, x, /*reuse=*/false);
	for(int t=0;t<J.nnz();t++)
		Jctmp( J.row[t], J.col[t]) = J.val[t];


   }
   else {

    	MatrixXd& xlb = *(workspace->xlb);
	    MatrixXd& xub = *(workspace->xub);
	    for(j=0;j<nvars;j++) { // EIGEN_UPDATE: index j shifted by -1
	      JacobianColumn( gg_num, xp, xlb, xub,j, &JacCol1, workspace->grw.get(), workspace);
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

   // icount is post-incremented as each row is copied, so rows 0..icount-1 are valid and
   // the assembled block has icount rows. Truncating it to icount-1 silently dropped one
   // active constraint, enlarging the computed null space by one dimension.

   // Simple bounds that are active at the solution are constraints as much as the
   // equalities are: a variable pinned at a bound is not free to vary and must not
   // contribute a direction to the null space. In a parameter estimation problem this
   // always includes the fixed initial and final times, and any parameter that has
   // converged onto one of its bounds. Each active bound contributes a unit row. Rows that
   // merely repeat an equality already present are harmless, because the null space is
   // taken from a rank-revealing decomposition.

   // A bound counts only if the solution is not free to move away from it. A variable
   // whose two bounds coincide is fixed outright and always counts. Otherwise the bound
   // must be active with a non-zero multiplier: a trajectory that merely grazes a state
   // bound at one node, with a zero multiplier, is still free to move inwards and does not
   // remove a degree of freedom. Where the solver supplied no bound multipliers, only
   // fixed variables are counted, which errs towards a larger null space and so towards
   // wider, not narrower, confidence intervals.

   const double bound_tol = 1.0e-8;

   MatrixXd& zbnd = workspace->bound_multipliers;
   bool have_z = ( zbnd.size() == nvars );
   double zmax = 0.0;
   if ( have_z ) for (j=0;j<nvars;j++) zmax = std::max(zmax, zbnd(j));
   const double z_tol = 1.0e-8*(1.0+zmax);

   for (j=0; j<nvars; j++) {
       double xj = X(j), lo = XL(j), up = XU(j);
       bool at_lower = ( lo > -PSOPT::inf/2.0 ) && ( fabs(xj-lo) <= bound_tol*(1.0+fabs(lo)) );
       bool at_upper = ( up <  PSOPT::inf/2.0 ) && ( fabs(xj-up) <= bound_tol*(1.0+fabs(up)) );
       bool fixed    = ( lo > -PSOPT::inf/2.0 ) && ( up < PSOPT::inf/2.0 )
                       && ( fabs(up-lo) <= bound_tol*(1.0+fabs(lo)) );
       bool binding  = fixed || ( (at_lower || at_upper) && have_z && zbnd(j) > z_tol );
       if ( binding ) {
           for (k=0;k<nvars;k++) Jc(icount,k) = 0.0;
           Jc(icount,j) = 1.0;
           icount++;
       }
   }

   Jc = Jc.block(0,0, icount, Jc.cols() );

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
	    JacobianColumn( rr_num, X, XL, XU, j, &Jcol, workspace->grw.get(), workspace);
        Jr.block(0,j,Jcol.rows(), 1) = Jcol;
	}
}

