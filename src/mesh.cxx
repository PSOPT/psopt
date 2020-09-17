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

using namespace Eigen;

void compute_next_mesh_size( Prob& problem, Alg& algorithm, Sol& solution, Workspace* workspace )
{
   int iphase;
	MatrixXd PHI;
	MatrixXd theta(2,1);
	double  yd;
	MatrixXd x;
	MatrixXd evec;
	int i;
	MatrixXd y, yls;
	int xd;
	char msg[100];
	int max_increment;
	double yerror;

  	for(iphase=1;iphase<=problem.nphases;iphase++) {
    	int decrease_count=0;
		int Ncurrent = problem.phase[iphase-1].current_number_of_intervals+1;


		int Nd;
	  	MatrixXd& emax_history = workspace->emax_history[iphase-1];

      x    = emax_history.block(0,0, workspace->current_mesh_refinement_iteration, 1); 

      evec = emax_history.block(0,1, workspace->current_mesh_refinement_iteration, 1); 

      y.resize(evec.rows(),1);
      
    	for(int iii=0;iii<evec.rows();iii++)   y(iii)= log(evec(iii));


		// If the error in this phase is not below the tolerance, then increment the nodes

		if ( evec( evec.rows() -1  ) > algorithm.ode_tolerance )   {   

		    // Check convergence behaviour

		   for (i=length(y)-2; i>=0; i--) {   // EIGEN_UPDATE: index i shifted by -1
				if (y(i+1)<y(i)) {
			   	 decrease_count++;
				}
				else break;
		   }

		   if (decrease_count>=1) {
			  // If the error has decreased calculate the next mesh size by using regression.

			  // Create regression matrix

			  PHI.resize(1,2);

			  int istart = length(y)-decrease_count-1;
			  int iend   = length(y)-1;

           int deltai = iend-istart+1;


           y = y.block(istart,0,deltai,1);

           x = x.block(istart,0,deltai,1);

           evec = (evec.block(istart,0,deltai,1));

			  
			  PHI.resize(length(y),2);
			  
           PHI(0,0) = log( x(0)); PHI(0,1) = 1.0;

			  for (i=1;i<length(y);i++) { // EIGEN_UPDATE: Index i shifted by -1.
             PHI(i,0) = log(x(i)); PHI(i,1) = 1.0;
			  }

//	    	          Solve the least squares problem

			  theta = (PHI.transpose()*PHI).colPivHouseholderQr().solve(PHI.transpose()*y);
			  

			  yls = PHI*theta;

           yerror = (y-yls).norm();

			  sprintf(msg, "\nLeast squares problem solved for global mesh refinement. Norm of residual = %e", yerror);

			  psopt_print(workspace,msg);



           yd = std::max( log(0.25*evec(length(evec)-1)), log(0.99*algorithm.ode_tolerance) );

           xd = (int) exp( (yd- theta(1))/theta(0) );

			  max_increment = std::max( (int) algorithm.mr_max_increment_factor*Ncurrent, 1);

			  Nd = std::min( Ncurrent + max_increment, std::max( xd, (Ncurrent+1)  ) );

			  sprintf(msg,"\nPhase %i: extrapolated number of nodes: %i, accepted number of nodes %d\n", iphase, xd, Nd);

		    }

		    else {
			// Take a cautious step forward and hope for the best
			Nd = Ncurrent + 1;
			sprintf(msg,"\nPhase %i: cautious step forward, next number of nodes %d\n", iphase, Nd);
		    }

		    psopt_print(workspace,msg);

	            problem.phase[iphase-1].current_number_of_intervals = Nd-1;

	       }  // End if



	}


}


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
  int p;
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
  double p;
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
	      if ( Icount >= (algorithm.mr_max_increment_factor)*(M-1) )
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


