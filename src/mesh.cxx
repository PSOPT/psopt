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

using namespace Eigen;

// NOTE: compute_next_mesh_size (the legacy global pseudospectral node-count refinement
// heuristic) was retired when hp-adaptive refinement became the default automatic strategy
// for the global pseudospectral methods (Radau/Gauss/Legendre/Chebyshev). The remaining
// functions in this file serve the local (Betts) refinement path.

bool check_for_equidistributed_error(Prob& problem,Alg& algorithm,Sol& solution)
{
    int nphases = problem.nphases;
    int iphase;

    int ecount = 0;

    bool retval = false;


    for(iphase=1;iphase<=nphases;iphase++) {
        	MatrixXd& epsilon = solution.relative_errors[iphase-1];

			double epsilon_max = Max(epsilon);

         double epsilon_mean = mean(epsilon);

			if ( epsilon_max <= 2*epsilon_mean ) {
	    		ecount++;
			}

    }

    if ( ecount == nphases ) {
			retval = true;
    }

    return retval;


}



void estimate_order_reduction(Prob& problem,Alg& algorithm,Sol& solution, Workspace* workspace)
{
  // This function estimates the local order reduction for the local mesh refinement algorithm.
  // Reference: Betts (2001), page 114
  //
  int nphases = problem.nphases;
  int iphase, k;
  // Default to the trapezoidal order. Local mesh refinement only runs with
  // trapezoidal (p=2) or Hermite-Simpson (p=4) defects, but initialising here
  // avoids using p uninitialised should differential_defects hold another value.
  int p = 2;
  double rhat, eta;
  int r;
  int j;



  if (workspace->differential_defects=="trapezoidal")
    p=2;
  else if (workspace->differential_defects=="Hermite-Simpson")
    p=4;

  for(iphase=1;iphase<=nphases;iphase++) {

      MatrixXd& order_reduction = workspace->order_reduction[iphase-1];

      MatrixXd& epsilon = solution.relative_errors[iphase-1];

      MatrixXd& old_epsilon = workspace->old_relative_errors[iphase-1];

      MatrixXd& snodes = workspace->snodes[iphase-1];

      MatrixXd& old_snodes = workspace->old_snodes[iphase-1];

      j = 0; // EIGEN_UPDATE

      int I = 1;

      for (k=0;k< problem.phase[iphase-1].current_number_of_intervals;k++) { // EIGEN_UPDATE: index k shifted by -1

	    eta   = epsilon(k);

	    if ( snodes(k) > old_snodes(j) ) {
	       j++; I=1;
	    }
	    else {
	       I++;
	    }

	    double theta = old_epsilon(j);

	    rhat = p+1.0 - log(theta/eta)/(log(1.0+I));


	    r = std::min(nint(rhat),(double) p);

	    r = std::max( 0, r );

	    order_reduction(k) = r;

      }

  }

}


void zero_order_reduction(Prob& problem,Alg& algorithm,Sol& solution, Workspace* workspace)
{
  // This function sets the order reduction to zero.
  //
  int nphases = problem.nphases;
  int iphase, k;


  for(iphase=1;iphase<=nphases;iphase++) {

      MatrixXd& order_reduction = workspace->order_reduction[iphase-1];

      for (k=0;k< problem.phase[iphase-1].current_number_of_intervals;k++) { // EIGEN_UPDATE

	      order_reduction(k) = 0;

      }

  }

}


void construct_new_mesh(Prob& problem,Alg& algorithm,Sol& solution, Workspace* workspace)
{
  // This function constructs the new mesh as part of the local mesh refinement algorithm.
  // Reference: Betts (2001), page 118
  //
  int nphases = problem.nphases;
  int iphase;
  // Default to the trapezoidal order (see estimate_order_reduction): avoids
  // using p uninitialised should differential_defects hold another value.
  double p = 2.0;
  long imax,rmax;
  bool terminate_flag = false;
  int M1 = algorithm.mr_M1;
  int M;
  double kappa = algorithm.mr_kappa;
  int Mdash;
  int Icount;
  MatrixXd epsilon;
  MatrixXd I;
  int i,l;

  if (workspace->differential_defects=="trapezoidal")
    p=2.0;
  else if (workspace->differential_defects=="Hermite-Simpson")
    p=4.0;

  for(iphase=1;iphase<=nphases;iphase++) {
        Icount = 0;
	M = problem.phase[iphase-1].current_number_of_intervals;
	Mdash =  min( (double) M1, kappa*M)+1;
	I.resize(1,M);
	I=zeros(1,M);
	terminate_flag = false;
        epsilon = solution.relative_errors[iphase-1];
        MatrixXd& r       = workspace->order_reduction[iphase-1];

	while (!terminate_flag) {
	      // Check which interval has maximum relative error
	      double epsilon_max = epsilon.maxCoeff(&rmax,&imax);

	      if ( (Icount>Mdash) && (epsilon_max <= algorithm.ode_tolerance) && (I(imax)==0) )
				terminate_flag=true;
	      if ( (epsilon_max <= kappa*algorithm.ode_tolerance) && (I(imax)<M1) && (I(imax)>0) )
				terminate_flag = true;
	      if ( Icount >= (algorithm.mr_max_growth_factor)*(M-1) )
				terminate_flag = true;
		   for(int iii=0;iii<I.cols(); iii++) {
	         if (  I(iii) == M1  ) {
		         terminate_flag = true;
		      }
		   }

	      if (!terminate_flag) {
	         // Add a point to interval imax
		  		I(imax) = I(imax)+1; Icount++;
		  		// Update the predicted error for interval imax
		  		epsilon(imax) = epsilon_max*pow( 1.0/(1.0+I(imax)), p-r(imax)+1.0);

	      }
	}

	// Now construct the new snodes array and sort it



	MatrixXd& snodes = workspace->snodes[iphase-1];


	for(i=0;i< M; i++) { // EIGEN_UPDATE: index i shifted by -1.
	    int Ii = (int) I(i);
	    if ( Ii> 0 ) {
	         double delta = snodes(i+1)-snodes(i);
	         for(l=1;l<=Ii;l++) { //EIGEN_UPDATE: NO NEED TO SHIFT THIS INDEX AS IT IS NOT USED BY A MATRIX OR VECTOR.

	             MatrixXd snodes_prev = snodes;
	             snodes.resize(snodes.rows(), snodes.cols()+1);
	             snodes.block(0,0, 1, snodes_prev.cols()) = snodes_prev;  
		          snodes(0, snodes.cols() -1 ) = (snodes(i) + (l)*delta/(Ii+1))  ;
		      }
	    }

	}

	sort(snodes);
   problem.phase[iphase-1].current_number_of_intervals = length(snodes)-1;

   fprintf(stderr,"\n >>> Local mesh refinement added %i new nodes in phase %i", (int) ( I.transpose() ).sum(), iphase );


  }



}


