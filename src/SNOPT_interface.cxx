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

static Workspace* workspace= NULL; // Temporary measure to get this file to compile.
                                      // Consider changing to C++ interface to see if
                                      // the same thing as IPOPT can be done with workspace.

void snPSOPTusrf_(int    *Status, int *n,    double x[],
	     int    *needF,  int *neF,  double F[],
	     int    *needG,  int *neG,  double G[],
	     char       *cu,     int *lencu,
	     int    iu[],    int *leniu,
	     double ru[],    int *lenru )

{

  workspace = tempsnoptworkspace;

  Alg& algorithm    = *workspace->algorithm;

  MatrixXd&        X = *workspace->Xsnopt;
  MatrixXd&        g = *workspace->gsnopt;
  TripletSparseMatrix&  Ax = *workspace->Ax;
  TripletSparseMatrix&  As = *workspace->As;

  int i;


  memcpy( X.GetPr(), x, (*n)*sizeof(double) );


  if (*needF) {

  	// Compute objective function value and assign to F[0]

        F[0] = ff_num(X, workspace);

  	// Compute constraint functions
        gg_num(X, &g, workspace);


  	//  assign g[i] to F[i]
  	for(i=1;i<(*neF);i++) {
	   	F[i] = g(i);
  	}


  	if (workspace->jac_done)
  	{
  		// Compute the linear part of the constraints:
       			Ax = As*X;
  		// Now subtract linear part of the constraints from F[i]
  		// This is done to return only the nonlinear part of the constraint function, so
  		// that SNOPT can handle separately the linear constraints.

  		for(i=0; i<(*neF);i++) {
         		F[i] -= Ax(i,1); // EIGEN_UPDATE
   		}


  	}


  }


  if( useAutomaticDifferentiation(algorithm) && workspace->jac_done && *needG ) {


        int nF = *neF;

        int nvars = *n;

        int k;

        double *xvars = &X(0);



       int options[4];
       options[0]=0; options[1]=0; options[2]=0;options[3]=0;
	    sparse_jac(workspace->tag_fg, nF, nvars, 0, xvars, &workspace->F_nnz, &workspace->iGfun2, &workspace->jGvar2, &workspace->G2, options);



        TripletSparseMatrix GS2(workspace->G2, nF, nvars, workspace->F_nnz, (int*) workspace->iGfun1, (int*) workspace->jGvar1);


        // Copy the result into G[] using the right indices to send it back to SNOPT.
        for(k=0;k<(*neG);k++) {
               G[k] = GS2(workspace->iGfun[k],workspace->jGvar[k]);
        }

        if (workspace->enable_nlp_counters) {
             workspace->solution->mesh_stats[ workspace->current_mesh_refinement_iteration-1 ].n_jacobian_evals++;
        }

  }



}


void fg_ad( adouble* x, adouble* fg, Workspace* workspace)
{

   int j;

   adouble fval;

   adouble* xad = workspace->xad;

   adouble* gad = workspace->gad;

   for(j=0; j<workspace->nvars; j++)
   {
        xad[j] = x[j];
   }

   fval = ff_ad( xad, workspace );

   gg_ad( xad, gad, workspace );

   fg[0] = fval;

   for(j=1; j<=workspace->ncons; j++)
   {
        fg[j] = gad[j-1];
   }

}



#endif // USE_SNOPT

