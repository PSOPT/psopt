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



void compute_next_mesh_size( Prob& problem, Alg& algorithm, Sol& solution, Workspace* workspace )
{
        int iphase;
        int p=1;
	DMatrix PHI;
	DMatrix theta(2,1);
	double  yd;
	DMatrix x;
	DMatrix evec;
	int i;
	DMatrix y, yls;
	int xd;
	char msg[100];
	int max_increment;
	double yerror;

  	for(iphase=1;iphase<=problem.nphases;iphase++) {
    	        int decrease_count=0;
		int Ncurrent = problem.phase[iphase-1].current_number_of_intervals+1;
     	        bool increase_flag = false;

		int Nd;
	  	DMatrix& emax_history = workspace->emax_history[iphase-1];
		x    = emax_history(colon(1,workspace->current_mesh_refinement_iteration), 1);
		evec = emax_history(colon(1,workspace->current_mesh_refinement_iteration), 2);
		y    = log( evec );

		// If the error in this phase is not below the tolerance, then increment the nodes
		if ( evec("end") > algorithm.ode_tolerance )   {

		    // Check convergence behaviour

		    for (i=length(y)-1; i>=1; i--) {
			if (y(i+1)<y(i)) {
			    decrease_count++;
			}
			else break;
		    }

		    if (decrease_count>=1) {
			  // If the error has decreased calculate the next mesh size by using regression.

			  // Create regression matrix

			  PHI.Resize(1,2);

			  int istart = length(y)-decrease_count;
			  int iend   = length(y);


			  y = y(colon(istart,iend));
			  x = x(colon(istart,iend));
			  evec = evec(colon(istart,iend));

			  double yerr_best = inf;

			  int pbest;

			  PHI =  ( log( x(1) )*ones(1,1) || ones(1,1) );

			  for (i=2;i<=length(y);i++) {
			    PHI = PHI && ( log(x(i))*ones(1,1) || ones(1,1) );
			  }

//	    	          Solve the least squares problem

			  theta = inv(tra(PHI)*PHI)*tra(PHI)*y;

			  yls = PHI*theta;

			  yerror = enorm(y-yls);

			  sprintf(msg, "\nLeast squares problem solved for global mesh refinement. Norm of residual = %e", yerror);

			  psopt_print(workspace,msg);

//			  PHI.Print("PHI");
//			  theta.Print("theta");
//			  yls.Print("yls");
//			  y.Print("y");


			  yd = MAX( log(0.25*evec("end")), log(0.99*algorithm.ode_tolerance) );

			  xd = (int) exp( (yd- theta(2))/theta(1) );

			  max_increment = MAX(algorithm.mr_max_increment_factor*Ncurrent,1);

			  Nd = min( Ncurrent + max_increment, max(xd,Ncurrent+1) );

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
        DMatrix& epsilon = solution.relative_errors[iphase-1];
	double epsilon_max = Max(tra(epsilon));
	double epsilon_mean = (mean(tra(epsilon)))(1);

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
  double rhat, eta, theta;
  int r;
  int j;
  int I;


  if (workspace->differential_defects=="trapezoidal")
    p=2;
  else if (workspace->differential_defects=="Hermite-Simpson")
    p=4;

  for(iphase=1;iphase<=nphases;iphase++) {

      DMatrix& order_reduction = workspace->order_reduction[iphase-1];

      DMatrix& epsilon = solution.relative_errors[iphase-1];

      DMatrix& old_epsilon = workspace->old_relative_errors[iphase-1];

      DMatrix& snodes = workspace->snodes[iphase-1];

      DMatrix& old_snodes = workspace->old_snodes[iphase-1];

      j = 1;

      int I = 1;

      for (k=1;k<= problem.phase[iphase-1].current_number_of_intervals;k++) {

	    eta   = epsilon(k);

	    if ( snodes(k) > old_snodes(j) ) {
	       j++; I=1;
	    }
	    else {
	       I++;
	    }

	    double theta = old_epsilon(j);

	    rhat = p+1.0 - log(theta/eta)/(log(1.0+I));


	    r = MIN(nint(rhat),p);

	    r = MAX( 0, r );

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

      DMatrix& order_reduction = workspace->order_reduction[iphase-1];

      for (k=1;k<= problem.phase[iphase-1].current_number_of_intervals;k++) {

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
  int iphase, k;
  double p;
  int imax;
  bool terminate_flag = false;
  int M1 = algorithm.mr_M1;
  int M;
  double kappa = algorithm.mr_kappa;
  int Mdash;
  int Icount;
  DMatrix epsilon;
  DMatrix I;
  int i,l;

  if (workspace->differential_defects=="trapezoidal")
    p=2.0;
  else if (workspace->differential_defects=="Hermite-Simpson")
    p=4.0;

  for(iphase=1;iphase<=nphases;iphase++) {
        Icount = 0;
	M = problem.phase[iphase-1].current_number_of_intervals;
	Mdash =  MIN( M1, kappa*M)+1;
	I.Resize(1,M);
	I.FillWithZeros();
	terminate_flag = false;
        epsilon = solution.relative_errors[iphase-1];
        DMatrix& r       = workspace->order_reduction[iphase-1];
//	r.Print("order reduction");
	while (!terminate_flag) {
	      // Check which interval has maximum relative error
	      double epsilon_max = Max(tra(epsilon), &imax);
//	      fprintf(stderr,"\n epsilon_max = %f, imax = %i", epsilon_max, imax );
	      if ( (Icount>Mdash) && (epsilon_max <= algorithm.ode_tolerance) && (I(imax)==0) )
		terminate_flag=true;
	      if ( (epsilon_max <= kappa*algorithm.ode_tolerance) && (I(imax)<M1) && (I(imax)>0) )
		terminate_flag = true;
	      if ( Icount >= (algorithm.mr_max_increment_factor)*(M-1) )
		terminate_flag = true;
	      if ( any( I == M1 ) )
		terminate_flag = true;

	      if (!terminate_flag) {
	          // Add a point to interval imax
		  I(imax) = I(imax)+1; Icount++;
		  // Update the predicted error for interval imax
		  epsilon(imax) = epsilon_max*pow( 1.0/(1.0+I(imax)), p-r(imax)+1.0);

//		  fprintf(stderr,"\n iphase: %i, imax: %i, old error=%f, new error=%f, r=%f, I(imax)=%f", iphase, imax, epsilon_max, epsilon(imax), r(imax), I(imax));
	      }
	}

	// Now construct the new snodes array and sort it

//	I.Print(">>> Added nodes per interval");

	DMatrix& snodes = workspace->snodes[iphase-1];


	for(i=1;i<= M; i++) {
	    int Ii = (int) I(i);
	    if ( Ii> 0 ) {
	         double delta = snodes(i+1)-snodes(i);
	         for(l=1;l<=Ii;l++) {
		     snodes = (snodes || (snodes(i) + (l)*delta/(Ii+1))*ones(1,1) ) ;
		 }
	    }

	}

	sort(snodes);
        problem.phase[iphase-1].current_number_of_intervals = length(snodes)-1;
	fprintf(stderr,"\n >>> Local mesh refinement added %i new nodes in phase %i", (int) sum(tra(I)).elem(1,1), iphase );


  }



}


