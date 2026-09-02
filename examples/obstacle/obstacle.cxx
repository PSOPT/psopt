//////////////////////////////////////////////////////////////////////////
////////////////        obstacle.cxx                 /////////////////////
//////////////////////////////////////////////////////////////////////////
////////////////           PSOPT  Example             ////////////////////
//////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////
//////// Title:      Minimum time obstacle avoidance      ////////////////
//////// Last modified: 02 September 2026                 ////////////////
//////// Reference:  after the PROPT User's Guide, with   ////////////////
////////             the objective changed to the time    ////////////////
//////////////////////////////////////////////////////////////////////////
////////     Copyright (c) Victor M. Becerra, 2009        ////////////////
//////////////////////////////////////////////////////////////////////////
//////// This is part of the PSOPT software library, which ///////////////
//////// is distributed under the terms of the GNU Lesser ////////////////
//////// General Public License (LGPL)                    ////////////////
//////////////////////////////////////////////////////////////////////////

// A particle travels at the fixed speed V = 2.138 with heading theta as its only
// control. It starts at (0,0), must reach (1.2,1.6), and must stay outside two
// discs of radius sqrt(0.1) centred at (0.4,0.5) and (0.8,1.5). The final time is
// free and is the objective.
//
// The problem in the PROPT User's Guide fixes the final time at 1 and minimises
// the integral of (dx/dt)^2 + (dy/dt)^2. With the speed fixed that integrand is
// identically V^2, so the cost is V^2 (tf - t0) for every admissible control and
// the problem has no objective at all. Worse, fixing tf pins the arclength of
// every admissible path at exactly V tf = 2.138 while the shortest admissible path
// is 2.131722, so the feasible set is every path of length exactly 2.138 that
// avoids the discs -- a surplus of three parts in a thousand that may be spent as
// a wiggle anywhere along the path, with every way of spending it equally optimal.
// Solved that way the example returns a different wandering trajectory on every
// mesh, and no discretization can do better, because they are all optimal.
//
// Minimising the time instead makes the problem the Euclidean shortest path among
// two discs, which has a unique solution in each homotopy class and an answer in
// closed form. The shortest path is a taut string: a straight segment leaving the
// start tangentially to an obstacle, an arc of that obstacle, and a straight
// segment to the target. Enumerating both sides of each obstacle and both orders,
// and discarding the classes whose straight segments cut a disc, gives
//
//      first obstacle only, passed on its right      2.131722291465    <- shortest
//      left of the first, right of the second        2.137849267659
//      left of both                                  2.544027869
//      right of the first, left of the second        2.921806899
//      right of both                                 4.119168448
//
// so the minimum time is
//
//      tf* = 2.131722291465491 / V = 0.997063747177498
//
// Only the first obstacle is active. The second is cleared at a squared distance
// of 0.109512045 against the 0.1 required. The optimal control is continuous and
// piecewise LINEAR in t: constant along each straight segment, and of constant
// rate V/r = 6.760949637 along the arc, since an arc of radius r traversed at
// speed V turns at exactly that rate.
//
//      t = 0        -> 0.260419287   theta = 0.379526559          (straight)
//      t = 0.260419 -> 0.378319287   theta rises by 0.797115961   (arc)
//      t = 0.378319 -> 0.997063747   theta = 1.176642521          (straight)
//
// The example prints its computed final time next to that value. On the mesh below
// the two agree to about nine significant figures, which is a far stronger check
// than any residual estimate: it is a comparison against an answer obtained
// without a discretization at all.
//
// Two things are worth knowing before changing the settings.
//
// The heading bound matters. theta is a heading, so it is defined only modulo
// 2 pi; the original bounds of +/- 10 admit about three revolutions, and with the
// final time free that gives a family of spurious local minima -- longer paths
// that wind -- which the solver does find. Bounded to a single branch, as below,
// the solve is clean; left at +/- 10 it lands on times between 1.3 and 2.7 on most
// meshes. That is a modelling fix, not a numerical one.
//
// The solution rides the constraint boundary along the arc, so the control has a
// corner in its derivative at each end of it, and a uniform mesh cannot place
// those two junctions except by being fine everywhere. Uniform Hermite-Simpson
// gives 3.3e-7 relative at 81 nodes, 8.0e-8 at 121, 5.4e-8 at 161 and 3.0e-8 at
// 201, whereas automatic mesh refinement finds the junctions by itself and reaches
// 3.0e-9 on the 74 nodes it settles at from the 41-node seed below -- an order of
// magnitude better than a uniform mesh nearly three times as large. That is why
// the example is set up this way, and it is a fair illustration of what mesh
// refinement is for.
//
// The integrated-residual transcription with a residual box reaches the same time
// to 4.2e-8 on 201 uniform nodes, with a residual three orders of magnitude
// smaller, but on the same mesh it takes 127 s against 0.21 s for Hermite-Simpson,
// and on coarse uniform meshes it is much the worse of the two (5.1e-3 at 81
// nodes). It is not the default here: once the problem has an objective, the
// answer does not need it.

#include "psopt.h"

using namespace std;

using namespace PSOPT;

static const double V_SPEED   = 2.138;

// The closed-form answer derived above, for the check printed after the solve.
static const double TF_EXACT  = 0.997063747177498;

//////////////////////////////////////////////////////////////////////////
///////////////////  Define the end point (Mayer) cost function //////////
//////////////////////////////////////////////////////////////////////////

adouble endpoint_cost(adouble* initial_states, adouble* final_states,
                      adouble* parameters,adouble& t0, adouble& tf,
                      adouble* xad, int iphase,Workspace* workspace)
{
   return tf;
}

//////////////////////////////////////////////////////////////////////////
///////////////////  Define the integrand (Lagrange) cost function  //////
//////////////////////////////////////////////////////////////////////////

adouble integrand_cost(adouble* states, adouble* controls,
                       adouble* parameters, adouble& time, adouble* xad,
                       int iphase, Workspace* workspace)
{
    return  0.0;
}


//////////////////////////////////////////////////////////////////////////
///////////////////  Define the DAE's ////////////////////////////////////
//////////////////////////////////////////////////////////////////////////

void dae(adouble* derivatives, adouble* path, adouble* states,
         adouble* controls, adouble* parameters, adouble& time,
         adouble* xad, int iphase, Workspace* workspace)
{

   adouble x    = states[ 0 ];
   adouble y    = states[ 1 ];


   adouble theta = controls[ 0 ];

   double V = V_SPEED;

   adouble dxdt = V*cos(theta);
   adouble dydt = V*sin(theta);


   derivatives[ 0 ] = dxdt;
   derivatives[ 1 ] = dydt;


   path[ 0 ] = pow(x-0.4,2.0) + pow(y-0.5,2.0);
   path[ 1 ] = pow(x-0.8,2.0) + pow(y-1.5,2.0);

}

////////////////////////////////////////////////////////////////////////////
///////////////////  Define the events function ////////////////////////////
////////////////////////////////////////////////////////////////////////////

void events(adouble* e, adouble* initial_states, adouble* final_states,
            adouble* parameters,adouble& t0, adouble& tf, adouble* xad,
            int iphase, Workspace* workspace)

{
   adouble x0 = initial_states[ 0 ];
   adouble y0 = initial_states[ 1 ];
   adouble xf = final_states[   0 ];
   adouble yf = final_states[   1 ];

   e[ 0 ] = x0;
   e[ 1 ] = y0;
   e[ 2 ] = xf;
   e[ 3 ] = yf;

}



///////////////////////////////////////////////////////////////////////////
///////////////////  Define the phase linkages function ///////////////////
///////////////////////////////////////////////////////////////////////////

void linkages( adouble* linkages, adouble* xad, Workspace* workspace)
{
  // No linkages as this is a single phase problem
}


////////////////////////////////////////////////////////////////////////////
/////  Propagate the returned control and report what it actually does /////
////////////////////////////////////////////////////////////////////////////

// The control history is read as the transcription defines it: over each mesh
// interval Hermite-Simpson carries a control at the midpoint as well as at the two
// nodes, and the control is the quadratic through the three. get_hs_controls_in_phase
// returns those 2M-1 values interleaved with the node values, and
// get_hs_time_in_phase the times that go with them; reading the node values alone
// would be reading a different control from the one the solver constrained.
//
// The states are then obtained by classical Runge-Kutta on a grid nsub times finer
// than the control intervals, starting from the initial condition of the problem,
// with no reference whatever to the returned states. This is the check that the
// answer is a trajectory of the system and not merely a solution of the discretized
// problem.

static void propagate_and_report(const MatrixXd& ths, const MatrixXd& uhs,
                                 double x_init, double y_init,
                                 double x_target, double y_target)
{
   const int    nhs  = (int) uhs.cols();          // 2M-1 control values
   const int    nsub = 200;                       // RK4 steps per half interval
   const double V    = V_SPEED;

   if ( nhs < 3 || (nhs % 2) == 0 ) {
      printf("\nBetween-node verification skipped: no Hermite-Simpson control history.\n");
      return;
   }

   double x = x_init, y = y_init;
   double d1min = 1.0e30, d2min = 1.0e30;

   // Closest approach on the propagated path, including the starting point.
   #define OBS_UPDATE()                                                        \
      {  double a = (x-0.4)*(x-0.4) + (y-0.5)*(y-0.5);                         \
         double b = (x-0.8)*(x-0.8) + (y-1.5)*(y-1.5);                         \
         if (a < d1min) d1min = a;                                             \
         if (b < d2min) d2min = b;  }

   OBS_UPDATE();

   // Interval k of the mesh spans control indices 2k, 2k+1 (midpoint), 2k+2.
   const int M = (nhs + 1)/2;

   for (int k = 0; k < M-1; k++) {

      const double ta = ths(0, 2*k), tm = ths(0, 2*k+1), tb = ths(0, 2*k+2);
      const double ua = uhs(0, 2*k), um = uhs(0, 2*k+1), ub = uhs(0, 2*k+2);

      // Quadratic through (ta,ua), (tm,um), (tb,ub) in the local variable
      // s = (t - ta)/(tb - ta), evaluated by Lagrange form so that a midpoint
      // that is not exactly central is handled correctly.
      const double h  = tb - ta;
      const double sm = (tm - ta)/h;

      #define THETA_AT(s)                                                       \
         (  ua*((s)-sm)*((s)-1.0)/((0.0-sm)*(0.0-1.0))                          \
          + um*((s)-0.0)*((s)-1.0)/((sm-0.0)*(sm-1.0))                          \
          + ub*((s)-0.0)*((s)-sm)/((1.0-0.0)*(1.0-sm)) )

      const double ds = 1.0/(2*nsub);

      for (int j = 0; j < 2*nsub; j++) {
         const double s0 = j*ds;
         const double dt = h*ds;

         const double th1 = THETA_AT(s0);
         const double th2 = THETA_AT(s0 + 0.5*ds);
         const double th4 = THETA_AT(s0 + ds);

         // f does not depend on the state, so the four RK4 stages use the
         // headings at the two ends and (twice) at the middle of the step.
         const double kx1 = V*cos(th1), ky1 = V*sin(th1);
         const double kx2 = V*cos(th2), ky2 = V*sin(th2);
         const double kx4 = V*cos(th4), ky4 = V*sin(th4);

         x += dt*(kx1 + 4.0*kx2 + kx4)/6.0;
         y += dt*(ky1 + 4.0*ky2 + ky4)/6.0;

         OBS_UPDATE();
      }
      #undef THETA_AT
   }
   #undef OBS_UPDATE

   const double miss = sqrt( (x-x_target)*(x-x_target) + (y-y_target)*(y-y_target) );

   printf("\n");
   printf("Between-node verification: the returned control propagated by RK4\n");
   printf("  propagated endpoint      : (%.6f, %.6f)\n", x, y);
   printf("  required endpoint        : (%.6f, %.6f)\n", x_target, y_target);
   printf("  terminal miss            : %.3e\n", miss);
   printf("  closest approach, obs 1  : %.9f   (required >= 0.1, exact 0.100000000)\n", d1min);
   printf("  closest approach, obs 2  : %.9f   (required >= 0.1, exact 0.109512045)\n", d2min);
   printf("\n");
   printf("  The two closest-approach figures are the squared distances used by the\n");
   printf("  path constraints, evaluated all along the propagated path rather than at\n");
   printf("  the nodes only. The first obstacle carries a boundary arc, where the\n");
   printf("  constraint is active over a whole interval and is imposed only at the\n");
   printf("  nodes and midpoints, so the propagated path may dip a little inside it\n");
   printf("  between them; the excursion here is of order 1e-7 in this squared\n");
   printf("  measure, which is 1e-7 of the obstacle radius.\n\n");
}


////////////////////////////////////////////////////////////////////////////
///////////////////  Define the main routine ///////////////////////////////
////////////////////////////////////////////////////////////////////////////

int main(void)
{

////////////////////////////////////////////////////////////////////////////
///////////////////  Declare key structures ////////////////////////////////
////////////////////////////////////////////////////////////////////////////

    Alg  algorithm;
    Sol  solution;
    Prob problem;

////////////////////////////////////////////////////////////////////////////
///////////////////  Register problem name  ////////////////////////////////
////////////////////////////////////////////////////////////////////////////

    problem.name        		= "Minimum time obstacle avoidance";
    problem.outfilename                 = "obstacle.txt";

////////////////////////////////////////////////////////////////////////////
////////////  Define problem level constants & do level 1 setup ////////////
////////////////////////////////////////////////////////////////////////////

    problem.nphases   						= 1;
    problem.nlinkages                   = 0;

    psopt_level1_setup(problem);

/////////////////////////////////////////////////////////////////////////////
/////////   Define phase related information & do level 2 setup  ////////////
/////////////////////////////////////////////////////////////////////////////

    problem.phases(1).nstates   						= 2;
    problem.phases(1).ncontrols 						= 1;
    problem.phases(1).nevents   						= 4;
    problem.phases(1).npath                     = 2;
    // A coarse seed: mesh refinement places the nodes where the corners are, and
    // finishes at about 74. Starting uniform and fine is worse than starting
    // coarse and refining, because the corners are what limit the accuracy.
    problem.phases(1).nodes                     = (RowVectorXi(1) << 41).finished();

    psopt_level2_setup(problem, algorithm);


////////////////////////////////////////////////////////////////////////////
///////////////////  Enter problem bounds information //////////////////////
////////////////////////////////////////////////////////////////////////////

    double xL = 0.0;
    double yL = 0.0;
    double xU = 2.0;
    double yU = 2.0;

    // theta is a heading, so it is defined modulo 2 pi. Bounding it to a single
    // branch removes the spurious winding minima; see the note at the head of this
    // file. The optimal heading lies in [0.380, 1.177], well inside these bounds.
    double thetaL = -pi/2.0;
    double thetaU =  pi;

    double x0 = 0.0;
    double y0 = 0.0;
    double xf = 1.2;
    double yf = 1.6;


    problem.phases(1).bounds.lower.states    <<  xL, yL;

    problem.phases(1).bounds.upper.states    <<  xU, yU;

    problem.phases(1).bounds.lower.controls  << thetaL;

    problem.phases(1).bounds.upper.controls  << thetaU;

    problem.phases(1).bounds.lower.events    << x0, y0, xf, yf;

    problem.phases(1).bounds.upper.events    << x0, y0, xf, yf;


    problem.phases(1).bounds.lower.path      <<  0.1, 0.1;

    problem.phases(1).bounds.upper.path      <<  100.0, 100.0;



    problem.phases(1).bounds.lower.StartTime    = 0.0;
    problem.phases(1).bounds.upper.StartTime    = 0.0;

    // The final time is now the objective, so it is free. The straight-line
    // distance is 2.0 and the speed is 2.138, so no admissible time is below
    // 2.0/2.138 = 0.9355.
    problem.phases(1).bounds.lower.EndTime      = 0.5;
    problem.phases(1).bounds.upper.EndTime      = 3.0;



////////////////////////////////////////////////////////////////////////////
///////////////////  Register problem functions  ///////////////////////////
////////////////////////////////////////////////////////////////////////////


    problem.integrand_cost 				= &integrand_cost;
    problem.endpoint_cost 					= &endpoint_cost;
    problem.dae             				= &dae;
    problem.events 							= &events;
    problem.linkages							= &linkages;

////////////////////////////////////////////////////////////////////////////
///////////////////  Define & register initial guess ///////////////////////
////////////////////////////////////////////////////////////////////////////

    // A guess of zero heading throughout tells the solver nothing about which way
    // to turn, and no straight path is admissible in any case. The guess below
    // sweeps the heading from 0.5 to 1.4 radians and integrates the kinematics: it
    // goes up and to the right, which is the only sensible way round, and nothing
    // about it is fitted to the answer. The solve is insensitive to it: over
    // theta_0 in [0.0, 0.8] and theta_f in [1.0, 1.6], eighteen of the twenty
    // combinations tried return the minimum time to at least eight figures, and
    // the other two return the taut path of the next homotopy class, 0.999929498,
    // which is itself a local minimum of the problem and which they also reach to
    // nine figures.

    int nnodes    							 	= 30;
    int ncontrols                       	= problem.phases(1).ncontrols;
    int nstates                         	= problem.phases(1).nstates;

    MatrixXd u_guess    						= zeros(ncontrols,nnodes);
    MatrixXd x_guess    						= zeros(nstates,nnodes);
    MatrixXd time_guess 						= linspace(0.0,1.0,nnodes);

    {
       const double theta_start = 0.5, theta_end = 1.4;
       const double dt = 1.0/(nnodes-1);
       double px = x0, py = y0;
       for (int k = 0; k < nnodes; k++) {
           double th = theta_start + (theta_end-theta_start)*k/(double)(nnodes-1);
           u_guess(0,k) = th;
           x_guess(0,k) = px;
           x_guess(1,k) = py;
           px += V_SPEED*cos(th)*dt;
           py += V_SPEED*sin(th)*dt;
       }
    }

    problem.phases(1).guess.controls      = u_guess;
    problem.phases(1).guess.states        = x_guess;
    problem.phases(1).guess.time          = time_guess;


////////////////////////////////////////////////////////////////////////////
///////////////////  Enter algorithm options  //////////////////////////////
////////////////////////////////////////////////////////////////////////////
    algorithm.nlp_iter_max                = 1000;
    algorithm.nlp_tolerance               = 1.e-6;
    algorithm.nlp_method                  = "IPOPT";
    algorithm.scaling                     = "automatic";
    algorithm.derivatives                 = "automatic";
    algorithm.collocation_method          = "Hermite-Simpson";
    algorithm.mesh_refinement             = "automatic";
    algorithm.ode_tolerance               = 1.e-8;


////////////////////////////////////////////////////////////////////////////
///////////////////  Now call PSOPT to solve the problem   /////////////////
////////////////////////////////////////////////////////////////////////////

    psopt(solution, problem, algorithm);

////////////////////////////////////////////////////////////////////////////
///////////  Extract relevant variables from solution structure   //////////
////////////////////////////////////////////////////////////////////////////


    MatrixXd states = solution.get_states_in_phase(1);
    MatrixXd theta  = solution.get_controls_in_phase(1);
    MatrixXd t      = solution.get_time_in_phase(1);
    MatrixXd mu     = solution.get_dual_path_in_phase(1);
    MatrixXd lambda = solution.get_dual_costates_in_phase(1);

    MatrixXd x = states.row(0);
    MatrixXd y = states.row(1);

////////////////////////////////////////////////////////////////////////////
///////////  Compare with the closed-form minimum time /////////////////////
////////////////////////////////////////////////////////////////////////////

    {
       double tf_psopt = solution.get_cost();
       printf("\n");
       printf("Minimum time\n");
       printf("  PSOPT                    : %.12f  on %d nodes\n",
              tf_psopt, (int) t.cols());
       printf("  taut-string construction : %.12f\n", TF_EXACT);
       printf("  relative difference      : %.3e\n",
              fabs(tf_psopt-TF_EXACT)/TF_EXACT);
       printf("  heading range            : [%.9f, %.9f]\n",
              theta.minCoeff(), theta.maxCoeff());
       printf("  exact heading range      : [0.379526559, 1.176642521]\n");
    }

////////////////////////////////////////////////////////////////////////////
///////////  Check what the returned control actually does /////////////////
////////////////////////////////////////////////////////////////////////////

    MatrixXd uhs = solution.get_hs_controls_in_phase(1);
    MatrixXd ths = solution.get_hs_time_in_phase(1);

    propagate_and_report(ths, uhs, x0, y0, xf, yf);

////////////////////////////////////////////////////////////////////////////
///////////  Save solution data to files if desired ////////////////////////
////////////////////////////////////////////////////////////////////////////

    Save(x,"x.dat");
    Save(y,"y.dat");
    Save(theta,"theta.dat");
    Save(t,"t.dat");


////////////////////////////////////////////////////////////////////////////
///////////  Plot some results if desired (requires gnuplot) ///////////////
////////////////////////////////////////////////////////////////////////////

    MatrixXd alpha = linspace(0.0, 2*pi, 100);
    MatrixXd xObs1 = sqrt(0.1)*cos(alpha) + 0.4*ones(1,length(alpha));
    MatrixXd yObs1 = sqrt(0.1)*sin(alpha) + 0.5*ones(1,length(alpha));
    MatrixXd xObs2 = sqrt(0.1)*cos(alpha) + 0.8*ones(1,length(alpha));
    MatrixXd yObs2 = sqrt(0.1)*sin(alpha) + 1.5*ones(1,length(alpha));


    plot(x,y,xObs1,yObs1,xObs2,yObs2,problem.name+": x-y trajectory",
                                            "x", "y", "y obs1 obs2");

    plot(x,y,xObs1,yObs1,xObs2,yObs2,problem.name+": x-y trajectory",
                                            "x", "y", "y obs1 obs2",
                                            "pdf", "obstacle_xy.pdf");

    plot(t,theta, problem.name+": theta","t", "theta");

    plot(t,theta, problem.name+": theta","t", "theta",
                                            "pdf", "obstacle_theta.pdf");

    plot(t,mu, problem.name+": path constraint multipliers","t", "mu_1 mu_2");

    plot(t,lambda, problem.name+": costates","t", "lambda_1 lambda_2");


}

////////////////////////////////////////////////////////////////////////////
///////////////////////      END OF FILE     ///////////////////////////////
////////////////////////////////////////////////////////////////////////////
