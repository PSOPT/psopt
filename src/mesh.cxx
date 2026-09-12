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
e-mail:    vmbecerra@vmb1.com

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

      const int M_new = problem.phase[iphase-1].current_number_of_intervals;
      const int M_old = length(old_snodes) - 1;

      // Map each interval of the new mesh to the interval of the old mesh that
      // contains it, and count how many pieces each old interval was cut into.
      //
      // The lookup has to advance on the *right* endpoint of the old interval:
      // a new sub-interval strictly inside old interval j has a left endpoint
      // greater than old_snodes(j), so comparing against the left endpoint
      // advances j once per new interval and attributes every sub-interval but
      // the first to the wrong parent. On a uniformly refined mesh the two agree;
      // on the non-uniform meshes this algorithm actually produces they do not,
      // and theta was then read from an unrelated interval.
      std::vector<int> parent(M_new, 0);
      std::vector<int> pieces(M_old > 0 ? M_old : 1, 0);

      j = 0;
      for (k=0; k<M_new; k++) {
	    const double sk = snodes(k);
	    while ( j < M_old-1 && sk >= old_snodes(j+1) - PSOPT_extras::GetEPS() ) j++;
	    parent[k] = j;
	    pieces[j]++;
      }

      for (k=0; k<M_new; k++) {

	    eta   = epsilon(k);

	    j = parent[k];

	    // Points added to the parent interval: one fewer than the number of
	    // pieces it was cut into. Betts's estimate divides by log(1+I_k), so an
	    // interval that was not subdivided carries no information about the
	    // order and is left at zero rather than dividing by log(1).
	    const int Ik = pieces[j] - 1;

	    if ( Ik <= 0 ) { order_reduction(k) = 0; continue; }

	    double theta = old_epsilon(j);

	    if ( !(theta > 0.0) || !(eta > 0.0) ) { order_reduction(k) = 0; continue; }

	    rhat = p+1.0 - log(theta/eta)/(log(1.0+ (double) Ik));

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



// ===========================================================================================
// Automatic mesh refinement for the integrated-residual transcription, which refines
// ELEMENTS rather than inserting nodes.
//
// Betts refinement puts new nodes at interval midpoints, and that is precisely what an
// element basis cannot take. The nodes of a Nie-Kerrigan element sit at its own LGL
// abscissae, so a node inserted between them belongs to no element, and the divisibility
// rule norder % d == 0 stops holding the moment the count changes by anything but a
// multiple of d. PSOPT refused the combination rather than run it and return an answer to a
// discretisation nobody had defined.
//
// What was wanted instead is refinement in the currency the transcription is written in. An
// element is either kept or split into k equal sub-elements, each of which is a proper
// element carrying its own abscissae, and the partition that results is one the
// transcription can state exactly. The node count moves in multiples of the stride by
// construction, so the divisibility rule is not something to check afterwards -- it cannot
// be broken.
//
// It composes with the flexible mesh, which is the reason it takes this form. The solved
// widths are already in snodes by the time this runs (ir_write_back_snodes), so the split
// works on the partition the previous solve CHOSE rather than on the uniform one it started
// from: a boundary the optimiser moved onto a switch is a boundary of the new partition too,
// and the elements on either side of it are then refined independently of each other. The
// alternative -- rebuilding a uniform partition at a larger node count -- throws away the one
// thing the flexible mesh found, every iteration, and asks the next solve to find it again.
//
// Called after the solve, like hp_refine_driver and unlike construct_new_mesh, because it
// needs the error estimate of the solve just finished and because old_snodes must already
// hold the solved mesh for the hot start to interpolate from.
// ===========================================================================================
void ir_refine_driver(Prob& problem, Alg& algorithm, Sol& solution, Workspace* workspace)
{
    // How many pieces one element may be cut into in a single iteration. A large error is a
    // reason to subdivide repeatedly across iterations, not to shatter an element in one:
    // the error estimate is computed on a trajectory the previous mesh could not represent,
    // so it says where the trouble is far more reliably than how much of it there is.
    const int    IR_MAX_SPLIT = 4;
    const double tol          = algorithm.ode_tolerance;

    const int stride = ir_element_stride(algorithm);

    for (int i = 0; i < problem.nphases; i++) {

        const int norder = problem.phase[i].current_number_of_intervals;
        const int M      = ir_num_elements(norder, algorithm);
        if ( M <= 0 ) continue;

        MatrixXd& sn  = workspace->snodes[i];
        MatrixXd& eps = solution.relative_errors[i];        // 1 x norder, per interval

        if ( (int) eps.size() < norder ) continue;          // no estimate to refine on

        // The error of an element is the worst of the intervals inside it, and its width is
        // read from the mesh that was actually solved on.
        std::vector<double> err(M, 0.0), h(M, 0.0);
        for (int e = 0; e < M; e++) {
            for (int r = e*stride; r < (e+1)*stride; r++)
                err[e] = std::max( err[e], eps(0,r) );
            h[e] = sn((e+1)*stride) - sn(e*stride);
        }

        // The order the local error converges at: degree d for a Nie-Kerrigan element, and
        // the Hermite-Simpson order for the cubic-Hermite one. The exponent only grades how
        // aggressively an element is cut; the decision to cut it at all is the tolerance.
        const int q = ( stride > 1 ) ? stride : 3;

        double emax = 0.0;
        for (int e = 0; e < M; e++) emax = std::max( emax, err[e] );

        std::vector<int> k(M, 1);
        int M_new = M;

        if ( !algorithm.ir_flexible_mesh ) {
            // A fixed mesh: the refinement decides both how many elements and where, which is
            // the classic arrangement, so split the elements whose error exceeds the tolerance
            // and grade the cut by how far it exceeds it.
            for (int e = 0; e < M; e++) {
                if ( err[e] <= tol ) continue;
                const double ratio = pow( err[e]/tol, 1.0/((double) q + 1.0) );
                int ke = (int) ceil(ratio);
                if ( ke < 2 )            ke = 2;
                if ( ke > IR_MAX_SPLIT ) ke = IR_MAX_SPLIT;
                k[e] = ke;
            }
            M_new = 0;
            for (int e = 0; e < M; e++) M_new += k[e];
        }
        else {
            // A flexible mesh: the two mechanisms are given disjoint jobs, because when they
            // are given the same one they fight.
            //
            // Measured, splitting by error with the flexible mesh on makes the answer WORSE.
            // On the minimum-time problem with a switch at tf/3 the flexible mesh alone
            // reaches 9.1e-8 on nine nodes; error-directed refinement to thirty-seven nodes
            // gives 1.9e-6, and the cubic-Hermite form degrades to 6.2e-5 with thirty tiny
            // elements packed around the switch. The reason is visible in the partition and
            // is not a tuning problem. Once the flexible mesh has isolated a discontinuity
            // inside a thin element, the local error of that element stays large however thin
            // it is -- the error is the jump, not the resolution -- so an estimator built for
            // smooth solutions flags it every iteration, the refinement splits it every
            // iteration, and the mesh starves the rest of the trajectory to feed a point that
            // was already handled.
            //
            // So the estimator is asked the only question it can answer well here, which is
            // HOW MANY elements the phase needs; where they go is the flexible mesh's job,
            // and it will re-place every boundary in the next solve in any case. The new
            // elements are seeded by splitting the WIDEST elements, which is where resolution
            // is cheap and which leaves a partition the next solve can start from without a
            // degenerate element in it.
            if ( emax <= tol ) continue;
            const double ratio  = pow( emax/tol, 1.0/((double) q + 1.0) );
            int M_target = (int) ceil( M*std::min( ratio, 1.0 + algorithm.mr_max_growth_factor ) );
            if ( M_target < M+1 ) M_target = M+1;

            while ( M_new < M_target ) {
                int widest = -1; double best = 0.0;
                for (int e = 0; e < M; e++) {
                    if ( k[e] >= IR_MAX_SPLIT ) continue;
                    const double sub = h[e]/((double) k[e]);
                    if ( widest < 0 || sub > best ) { widest = e; best = sub; }
                }
                if ( widest < 0 ) break;
                k[widest]++; M_new++;
            }
        }

        if ( M_new == M ) continue;                          // nothing exceeded the tolerance

        // Two ceilings. The growth factor is the user's limit on how fast the mesh may grow,
        // and the workspace ceiling is not negotiable at all: max_nodes sized xad, and a mesh
        // past it writes off the end of the tape.
        const int max_nodes  = get_max_nodes(problem, i+1, &algorithm);
        int cap = M + (int) floor( M*algorithm.mr_max_growth_factor );
        if ( cap < M+1 )              cap = M+1;
        if ( cap > max_nodes/stride ) cap = max_nodes/stride;
        if ( cap < M )                cap = M;

        // Over budget: give up the splits of the least troublesome elements first, so that
        // the budget is spent where the error is.
        while ( M_new > cap ) {
            int worst = -1; double smallest = 0.0;
            for (int e = 0; e < M; e++) {
                if ( k[e] <= 1 ) continue;
                if ( worst < 0 || err[e] < smallest ) { worst = e; smallest = err[e]; }
            }
            if ( worst < 0 ) break;
            k[worst]--; M_new--;
        }
        if ( M_new <= M ) continue;

        // A split may not push an element below the floor the flexible mesh will impose on
        // it, or the mesh just built would be outside its own bounds. The floor falls as the
        // partition grows, so removing a split can make another one legal again; the loop
        // settles because M_new only decreases.
        if ( algorithm.ir_flexible_mesh ) {
            bool changed = true;
            while ( changed && M_new > M ) {
                changed = false;
                const double hlo = algorithm.ir_min_element_fraction * 2.0/((double) M_new);
                for (int e = 0; e < M; e++) {
                    while ( k[e] > 1 && h[e]/((double) k[e]) < hlo ) { k[e]--; M_new--; changed = true; }
                }
            }
            if ( M_new <= M ) continue;
        }

        // Build the new partition. Each split element contributes k equal pieces of its own
        // width, so the boundaries of the old partition all survive into the new one.
        std::vector<double> hn;
        hn.reserve(M_new);
        for (int e = 0; e < M; e++)
            for (int j = 0; j < k[e]; j++) hn.push_back( h[e]/((double) k[e]) );

        const int norder_new = M_new*stride;
        sn.resize(1, norder_new+1);
        MatrixXd& lgl01 = workspace->ir_lgl01;

        double a = -1.0;
        for (int e = 0; e < M_new; e++) {
            if ( stride == 1 ) sn(e) = a;
            else for (int r = 0; r < stride; r++) sn(e*stride + r) = a + lgl01(r)*hn[e];
            a += hn[e];
        }
        sn(norder_new) = 1.0;      // pinned, for the reason ir_write_back_snodes pins it

        problem.phase[i].current_number_of_intervals = norder_new;

        snprintf(workspace->text, sizeof(workspace->text),
                 "\n>>> Phase %d: integrated-residual element refinement, %d -> %d elements "
                 "(%d -> %d intervals), worst element error %e\n",
                 i+1, M, M_new, norder, norder_new,
                 *std::max_element(err.begin(), err.end()));
        psopt_print(workspace, workspace->text);
    }
}
