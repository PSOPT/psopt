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
//
// The splitting rule itself is shared with multiple shooting's segment refinement, because
// the rule is the same and only the indicator differs: split what exceeds the tolerance and
// grade the cut by how far it exceeds it, unless a flexible partition is in force, in which
// case decide only HOW MANY pieces the phase needs and split the widest. Two copies of that
// would have drifted the first time either was tuned.
// ===========================================================================================

// How many pieces one piece may be cut into in a single iteration. A large error is a reason
// to subdivide repeatedly across iterations, not to shatter a piece in one: the indicator is
// computed on a trajectory the previous mesh could not represent, so it says where the trouble
// is far more reliably than how much of it there is.
#define PSOPT_MAX_SPLIT 4

// Decides the new partition of one phase from a per-piece indicator and the widths that were
// solved on, and writes snodes and current_number_of_intervals. Returns true when the mesh
// changed.
//
//   err     one entry per piece: the indicator the refinement acts on
//   h       one entry per piece: its width on the normalised interval, as solved
//   stride  intervals per piece (1 for a shooting segment or a cubic-Hermite element)
//   q       the order the piece's local error converges at; grades the cut, nothing else
//   tol     the tolerance err is compared against
//   label   what to call the thing in the message
static bool psopt_refine_partition(Prob& problem, Alg& algorithm, Workspace* workspace,
                                   int i, const std::vector<double>& err,
                                   const std::vector<double>& h,
                                   int stride, int q, double tol, const char* label)
{
    const int M      = (int) err.size();
    const int norder = problem.phase[i].current_number_of_intervals;
    MatrixXd& sn     = workspace->snodes[i];

    double emax = 0.0;
    for (int e = 0; e < M; e++) emax = std::max( emax, err[e] );

    std::vector<int> k(M, 1);
    int M_new = M;

    // The budget, computed BEFORE anything is split rather than after. Two ceilings: the
    // growth factor is the user's limit on how fast the mesh may grow, and the workspace
    // ceiling is not negotiable at all -- max_nodes sized xad, and a mesh past it writes off
    // the end of the tape.
    const int max_nodes  = get_max_nodes(problem, i+1, &algorithm);
    int cap = M + (int) floor( M*algorithm.mr_max_growth_factor );
    if ( cap < M+1 )              cap = M+1;
    if ( cap > max_nodes/stride ) cap = max_nodes/stride;
    if ( cap < M )                cap = M;

    if ( !flexible_partition_active(algorithm) ) {
        // A fixed mesh: the refinement decides both how many pieces and where, which is the
        // classic arrangement. Ask of each piece what its own error alone suggests -- split
        // it, graded by how far it exceeds the tolerance -- and let the budget below decide
        // how much of that can be afforded.
        for (int e = 0; e < M; e++) {
            if ( err[e] <= tol ) continue;
            const double ratio = pow( err[e]/tol, 1.0/((double) q + 1.0) );
            int ke = (int) ceil(ratio);
            if ( ke < 2 )               ke = 2;
            if ( ke > PSOPT_MAX_SPLIT ) ke = PSOPT_MAX_SPLIT;
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
        // HOW MANY pieces the phase needs; where they go is the flexible partition's job,
        // and it will re-place every boundary in the next solve in any case. The new
        // pieces are seeded by splitting the WIDEST, which is where resolution is cheap and
        // which leaves a partition the next solve can start from without a degenerate
        // piece in it.
        if ( emax <= tol ) return false;
        const double ratio  = pow( emax/tol, 1.0/((double) q + 1.0) );
        int M_target = (int) ceil( M*std::min( ratio, 1.0 + algorithm.mr_max_growth_factor ) );
        if ( M_target < M+1 ) M_target = M+1;

        while ( M_new < M_target ) {
            int widest = -1; double best = 0.0;
            for (int e = 0; e < M; e++) {
                if ( k[e] >= PSOPT_MAX_SPLIT ) continue;
                const double sub = h[e]/((double) k[e]);
                if ( widest < 0 || sub > best ) { widest = e; best = sub; }
            }
            if ( widest < 0 ) break;
            k[widest]++; M_new++;
        }
    }

    if ( M_new == M ) return false;                      // nothing exceeded the tolerance

    // Over budget, which on a fixed mesh is the normal case rather than the exception: what
    // each piece asks for is decided by its own error alone and the budget is a limit on the
    // total, so the interesting question is not what to ask for but what to give up.
    //
    // Give up the split whose loss costs least, measured by the PREDICTED error the piece
    // would be left with -- err/(k-1), one fewer piece than it currently holds. That is the
    // change from what this did before, which was to give up the splits of the piece with the
    // smallest raw error first. The two agree when one piece is far worse than the rest, and
    // they disagree completely when every piece is about as bad as every other: the old rule
    // drove the small-error pieces to a single piece each and then cut the worst back, which
    // on a smooth problem meant the entire budget went into quartering a seventh of the mesh
    // while the widest pieces -- the ones that set the error -- were left untouched. Measured
    // on the oscillator with a held control: six refinements took five segments to thirty and
    // the cost error fell by a factor of 1.8, where a UNIFORM mesh of thirty segments gives a
    // factor of 50. Equalising the predicted error instead spreads the budget and the mesh
    // stays near-uniform where the difficulty is.
    while ( M_new > cap ) {
        int give = -1; double smallest = 0.0;
        for (int e = 0; e < M; e++) {
            if ( k[e] <= 1 ) continue;
            const double left = err[e]/((double) k[e] - 1.0);
            if ( give < 0 || left < smallest ) { give = e; smallest = left; }
        }
        if ( give < 0 ) break;
        k[give]--; M_new--;
    }
    if ( M_new <= M ) return false;

    // A split may not push a piece below the floor the flexible partition will impose on
    // it, or the mesh just built would be outside its own bounds. The floor falls as the
    // partition grows, so removing a split can make another one legal again; the loop
    // settles because M_new only decreases.
    if ( flexible_partition_active(algorithm) ) {
        bool changed = true;
        while ( changed && M_new > M ) {
            changed = false;
            const double hlo = min_partition_fraction(algorithm) * 2.0/((double) M_new);
            for (int e = 0; e < M; e++) {
                while ( k[e] > 1 && h[e]/((double) k[e]) < hlo ) { k[e]--; M_new--; changed = true; }
            }
        }
        if ( M_new <= M ) return false;
    }

    // Build the new partition. Each split piece contributes k equal parts of its own
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
             "\n>>> Phase %d: %s, %d -> %d (%d -> %d intervals), worst indicator %e\n",
             i+1, label, M, M_new, norder, norder_new, emax);
    psopt_print(workspace, workspace->text);
    return true;
}


void ir_refine_driver(Prob& problem, Alg& algorithm, Sol& solution, Workspace* workspace)
{
    const double tol    = algorithm.ode_tolerance;
    const int    stride = ir_element_stride(algorithm);

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

        (void) psopt_refine_partition(problem, algorithm, workspace, i, err, h, stride, q, tol,
                                      "integrated-residual element refinement");
    }
}


// ===========================================================================================
// Multiple shooting's automatic segment refinement.
//
// The question is not the one the other drivers answer, and mistaking it for that question is
// the whole risk here. On a collocation mesh more nodes means a better approximation of the
// DYNAMICS. On a shooting mesh the dynamics are integrated to whatever ms_steps_per_segment
// and ms_integrator buy, however many segments there are, so refining segments to drive the
// reported ODE error down would be spending decision variables on an error that a step count
// buys for none. What the segment count does control is two other things:
//
//   the resolution of the CONTROL PARAMETERISATION, which is what caps the accuracy of the
//   answer once the integrator is adequate; and
//
//   the coverage of the PATH CONSTRAINTS, which between two boundaries are enforced at
//   ms_path_samples interior points and nowhere else.
//
// The indicator measures both and takes the worse of them.
//
// THE CONTROL PART, AND THE TRAP IN IT. The natural estimate of a representation's error is
// its departure from a one-degree-higher reconstruction of the same data. Written naively
// that flags a corner of the optimal control sitting exactly ON a segment boundary -- where
// the representation is not in error at all, a piecewise form being free to have a corner at
// a boundary -- and it flags it just as loudly however fine the mesh becomes, because the
// jump does not shrink with the mesh. That is the patch-177 trap in a new place: an estimator
// built for smooth solutions, pointed at a discontinuity the mesh has already resolved.
//
// The remedy is to form the departure twice, from the window extended to the LEFT and from
// the window extended to the RIGHT, and take the SMALLER of the two. A corner at a boundary
// spoils exactly one of the two windows, so the minimum ignores it; a corner INSIDE a segment
// spoils both, so the minimum sees it; on a smooth arc the two agree and the minimum is the
// estimate. What the indicator then measures is "how badly is this segment's own interior
// represented", which is the question the segment count answers, and not "how fast is the
// control changing", which it does not.
// ===========================================================================================

// How many points inside a segment the indicator looks at. They are for measuring, not for
// constraining, so they are independent of ms_path_samples -- the point of looking there is
// precisely that nothing is being enforced at those places.
#define MS_INDICATOR_SAMPLES 3

// Lagrange interpolation through n samples, evaluated at t. n = 1 returns the sample.
static double ms_lagrange_at(const double* ts, const double* us, int n, double t)
{
    double v = 0.0;
    for (int j = 0; j < n; j++) {
        double L = 1.0;
        for (int m = 0; m < n; m++) {
            if ( m == j ) continue;
            const double den = ts[j] - ts[m];
            if ( fabs(den) < 1.0e-300 ) return us[j];     // coincident samples: degenerate
            L *= (t - ts[m])/den;
        }
        v += L*us[j];
    }
    return v;
}


void ms_segment_indicators(Prob& problem, Alg& algorithm, Sol& solution, Workspace* workspace,
                           int iphase, std::vector<double>& indicator)
{
    const int i         = iphase - 1;
    const int M         = problem.phase[i].current_number_of_intervals;
    const int ncontrols = problem.phase[i].ncontrols;
    const int nstates   = problem.phase[i].nstates;
    const int npath     = problem.phase[i].npath;
    const int nparam    = problem.phase[i].nparameters;
    const int d         = ms_control_degree(algorithm);

    indicator.assign( (M > 0) ? M : 1, 0.0 );
    if ( M <= 0 ) return;

    // ---- the control parameterisation --------------------------------------------------
    //
    // The sample sequence the representation is built from, which is not the same array for
    // the three forms and must not be taken to be.
    //
    //   held      the M segment values, placed at the segment MIDPOINTS. A constant best
    //             approximates a function at the middle of its span, so that is where the
    //             sample is; and the sequence has M entries, not M+1, because the terminal
    //             control slot is not a segment's value at all -- it belongs to no segment
    //             and is pinned to its neighbour. Reading it as one makes the last segment's
    //             difference identically zero, which is a segment that can never be refined
    //             however badly it needs to be. Measured before this was seen: the last
    //             segment kept its original width while the mesh grew around it, ending six
    //             times the width of its neighbours and carrying the whole error.
    //   ramp      the M+1 nodal values.
    //   parabola  the 2M+1 node-and-midpoint values, which is the whole control history there.
    //
    // The reference carries TWO more degrees than the representation, not one. One more is
    // what estimates the leading error term, and the leading term vanishes where it is least
    // safe for it to: a held control at an extremum of the optimal control has no first
    // difference at all, and an error of O(h^2 u'') that a first difference cannot see. The
    // extra degree costs two samples of window and sees it.
    if ( ncontrols > 0 ) {
        const bool use_hs = ( d == 2 && solution.controls_hs != NULL
                              && solution.controls_hs[i].cols() == 2*M + 1 );
        const MatrixXd& Usrc = use_hs ? solution.controls_hs[i] : solution.controls[i];
        const MatrixXd& Tsrc = use_hs ? solution.nodes_hs[i]    : solution.nodes[i];
        const MatrixXd& Tseg = solution.nodes[i];

        const int per = use_hs ? 2 : 1;
        const int own = d + 1;
        const int ref = own + 2;

        // N: how many samples the sequence has, and where sample j sits.
        const int N = ( d == 0 ) ? M : (int) Usrc.cols();

        if ( N >= 2 && (int) Tsrc.cols() >= ((d == 0) ? M+1 : N)
             && (int) Tseg.cols() >= M+1 && (int) Usrc.rows() >= ncontrols ) {

            std::vector<double> ts(N);
            for (int j = 0; j < N; j++)
                ts[j] = ( d == 0 ) ? 0.5*( Tseg(0,j) + Tseg(0,j+1) ) : Tsrc(0,j);

            std::vector<double> tw(ref), uw(ref), to(own), uo(own);
            for (int c = 0; c < ncontrols; c++) {
                double uscale = 0.0;
                for (int j = 0; j < N; j++) uscale = std::max( uscale, fabs(Usrc(c,j)) );
                uscale += 1.0;

                for (int k = 0; k < M; k++) {
                    const int b = k*per;
                    if ( b + own - 1 > N-1 ) continue;
                    const double ta = Tseg(0,k), tb = Tseg(0,k+1);
                    if ( !(tb > ta) ) continue;

                    for (int j = 0; j < own; j++) { to[j] = ts[b+j]; uo[j] = Usrc(c,b+j); }

                    // The two windows. One extends to the left of the segment and one to the
                    // right, and neither reaches across it, which is what lets the smaller of
                    // the two ignore a corner sitting ON a boundary: such a corner spoils
                    // exactly one window. A corner INSIDE the segment spoils both, and is
                    // seen. On a smooth arc the two agree.
                    double dep[2] = { -1.0, -1.0 };
                    for (int side = 0; side < 2; side++) {
                        const int start = ( side == 0 ) ? b + own - ref : b;
                        if ( start < 0 || start + ref - 1 > N-1 ) continue;
                        for (int j = 0; j < ref; j++) { tw[j] = ts[start+j]; uw[j] = Usrc(c,start+j); }
                        double bad = 0.0;
                        for (int q = 0; q <= 8; q++) {
                            const double t = ta + (((double) q)/8.0)*(tb - ta);
                            const double rep = ms_lagrange_at(to.data(), uo.data(), own, t);
                            const double rfn = ms_lagrange_at(tw.data(), uw.data(), ref, t);
                            bad = std::max( bad, fabs(rep - rfn) );
                        }
                        dep[side] = bad/uscale;
                    }
                    double e = 0.0;
                    if      ( dep[0] >= 0.0 && dep[1] >= 0.0 ) e = std::min(dep[0], dep[1]);
                    else if ( dep[0] >= 0.0 )                  e = dep[0];
                    else if ( dep[1] >= 0.0 )                  e = dep[1];
                    indicator[k] = std::max( indicator[k], e );
                }
            }
        }
    }

    // ---- the path constraints, where nothing is enforcing them --------------------------
    //
    // The state between two boundaries is not a decision variable, so a path constraint
    // imposed only at boundaries can be violated freely in between -- and the violation is
    // not visible in anything the NLP reports, because from the NLP's point of view every
    // constraint it was given is satisfied. It has to be looked for, by propagating the
    // segment again and asking the user's own path function what happened inside it.
    if ( npath > 0 && nstates > 0 ) {
        adouble* xad = workspace->xad.get();
        adouble t0, tf;
        get_times(&t0, &tf, xad, iphase, workspace);

        const int NS = MS_INDICATOR_SAMPLES;
        int nsteps = algorithm.ms_steps_per_segment;
        if ( nsteps < NS + 1 ) nsteps = NS + 1;          // room for the samples to sit in

        std::vector<adouble> par( (nparam > 0) ? nparam : 1 );
        get_parameters(par.data(), xad, iphase, workspace);

        std::vector<adouble> xend(nstates), xs(NS*nstates), ts(NS);
        std::vector<adouble> us( NS*((ncontrols > 0) ? ncontrols : 1) );
        std::vector<adouble> deriv(nstates), pth(npath);

        for (int k = 0; k < M; k++) {
            ms_propagate_segment(xend.data(), NULL, k, xad, iphase, t0, tf, par.data(),
                                 workspace, nsteps, xs.data(),
                                 (ncontrols > 0) ? us.data() : NULL, ts.data(), NULL, NS);
            for (int q = 0; q < NS; q++) {
                problem.dae(deriv.data(), pth.data(), xs.data() + q*nstates,
                            (ncontrols > 0) ? (us.data() + q*ncontrols) : us.data(),
                            par.data(), ts[q], xad, iphase, workspace);
                for (int pj = 0; pj < npath; pj++) {
                    const double lo = problem.phase[i].bounds.lower.path(pj);
                    const double up = problem.phase[i].bounds.upper.path(pj);
                    const double v  = pth[pj].value();
                    const double sc = std::max( fabs(lo), fabs(up) ) + 1.0;
                    double viol = 0.0;
                    if ( v > up ) viol = (v - up)/sc;
                    if ( v < lo ) viol = (lo - v)/sc;
                    indicator[k] = std::max( indicator[k], viol );
                }
            }
        }
    }
}


double ms_refine_driver(Prob& problem, Alg& algorithm, Sol& solution, Workspace* workspace,
                        bool do_refine)
{
    const double tol   = algorithm.ms_refine_tolerance;
    double       worst = 0.0;

    for (int i = 0; i < problem.nphases; i++) {

        const int M = problem.phase[i].current_number_of_intervals;
        if ( M <= 0 ) continue;

        std::vector<double> ind;
        ms_segment_indicators(problem, algorithm, solution, workspace, i+1, ind);
        if ( (int) ind.size() < M ) continue;

        for (int k = 0; k < M; k++) worst = std::max( worst, ind[k] );

        if ( !do_refine ) continue;

        // The widths solved on. With ms_flexible_segments they are not uniform, and they are
        // what the split has to work from: a boundary the optimiser moved onto a switch is a
        // boundary of the new partition too.
        MatrixXd& sn = workspace->snodes[i];
        std::vector<double> h(M, 0.0);
        for (int k = 0; k < M; k++) h[k] = sn(k+1) - sn(k);

        // The representation's local error is O(h^(degree+1)), so the exponent that grades a
        // split is 1/(degree+1). A held control is first order and is graded most aggressively
        // of the three, which is right: it is the form that needs the most segments.
        (void) psopt_refine_partition(problem, algorithm, workspace, i, ind, h, 1,
                                      ms_control_degree(algorithm), tol,
                                      "multiple-shooting segment refinement");
    }
    return worst;
}
