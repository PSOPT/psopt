//////////////////////////////////////////////////////////////////////////
////////////////        obstacle.cxx                 /////////////////////
//////////////////////////////////////////////////////////////////////////
////////////////           PSOPT  Example             ////////////////////
//////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////
//////// Title:         Obstacle avoidance problem        ////////////////
//////// Last modified: 01 September 2026                 ////////////////
//////// Reference:  PROPT User's Guide                   ////////////////
//////// (See PSOPT handbook for full reference)          ////////////////
//////////////////////////////////////////////////////////////////////////
////////     Copyright (c) Victor M. Becerra, 2009        ////////////////
//////////////////////////////////////////////////////////////////////////
//////// This is part of the PSOPT software library, which ///////////////
//////// is distributed under the terms of the GNU Lesser ////////////////
//////// General Public License (LGPL)                    ////////////////
//////////////////////////////////////////////////////////////////////////

// A particle travels at the fixed speed V = 2.138 with heading theta as the only
// control, from (0,0) to (1.2,1.6) in one unit of time, and must stay outside two
// discs of radius sqrt(0.1) centred at (0.4,0.5) and (0.8,1.5).
//
// Two properties of the problem as posed deserve a note, because they explain why
// this example is discretized differently from most of the others in this
// directory.
//
// First, the problem has no objective. The integrand is
//
//        (dx/dt)^2 + (dy/dt)^2 = V^2 cos^2(theta) + V^2 sin^2(theta) = V^2,
//
// so the cost is identically V^2 (tf - t0) = 4.571044 for every admissible control.
// The problem is a pure feasibility problem, and the reported cost carries no
// information about the quality of the answer: two solutions can be equally optimal
// and yet one of them can be nonsense.
//
// Second, theta is a heading, so it is defined only modulo 2 pi, and the bounds
// below admit about three revolutions. Nothing in the formulation prefers one
// branch over another, so an optimizer is free to wander between them.
//
// Together these two facts make the problem an unusually severe test of both the
// initial guess and the transcription. With no objective to pull the iterate
// anywhere and no penalty for crossing branches, a discretization is free to return
// node values that satisfy every collocation and path condition at the nodes while
// describing a control that is nonsense between them.
//
// The four columns below are the same problem under four settings, all measured with
// PSOPT's own local error estimate and by propagating the returned control
// accurately from (0,0) and comparing where it arrives with the required (1.2,1.6).
// The last column is the closest the propagated path comes to the centre of the
// first obstacle, in the squared measure the path constraint uses, so 0.1 is the
// requirement and anything below it means the path goes inside the disc.
//
//                                          error est.  heading at nodes   miss   obs 1
//   Legendre, 60 then 80 nodes, zero guess   1.77e-2   [-5.72,  7.65]    0.819  0.0280
//   Legendre, 81 nodes, guess below          2.60e-6   [ 0.348, 1.515]   1.4e-7 0.1012
//   Hermite-Simpson, 81 nodes, guess below   4.99e-3   [-5.99,  1.42]    1.5e-2 0.0942
//   integrated residual, box 1e-6            6.08e-9   [ 0.341, 1.273]   3.2e-10 0.1010
//
// The first line is how this example was set up until now: the heading crosses
// branches, the error estimate stalls near 2e-2, the propagated path ends 0.8 away
// from the target and passes well inside the first obstacle. The second line shows
// that most of that is the fault of the initial guess rather than of the pseudo-
// spectral method. The third shows that a good guess is not by itself enough, since
// Hermite-Simpson still crosses a branch and its propagated path enters the first
// obstacle even though every node and midpoint satisfies the path constraint.
//
// The integrated residual transcription with a residual box is the setting that is
// both accurate and reliable here, because it bounds the residual of the differential
// equations at points interior to each interval rather than only at the collocation
// points: a control that is nonsense between the nodes is then inadmissible rather
// than merely undetected.
//
// The verification block after the solve is part of the example on purpose. It
// propagates the control that the transcription itself returns and reports the
// terminal miss and the closest approach to each obstacle, which is the check that
// separates the four settings above; the objective cannot, since it is constant.
//
// If a formulation change is acceptable in an application of this kind, bounding
// theta to a single branch (for instance [-pi/2, pi] here) is the underlying
// modelling fix and helps under any discretization. The bounds are left as posed
// here so that the example matches its reference.

#include "psopt.h"

using namespace std;

using namespace PSOPT;

static const double V_SPEED = 2.138;

//////////////////////////////////////////////////////////////////////////
///////////////////  Define the end point (Mayer) cost function //////////
//////////////////////////////////////////////////////////////////////////

adouble endpoint_cost(adouble* initial_states, adouble* final_states,
                      adouble* parameters,adouble& t0, adouble& tf,
                      adouble* xad, int iphase,Workspace* workspace)
{
   return 0;
}

//////////////////////////////////////////////////////////////////////////
///////////////////  Define the integrand (Lagrange) cost function  //////
//////////////////////////////////////////////////////////////////////////

adouble integrand_cost(adouble* states, adouble* controls,
                       adouble* parameters, adouble& time, adouble* xad,
                       int iphase, Workspace* workspace)
{
    double V = V_SPEED;
    adouble theta = controls[ 0 ];

    adouble dxdt = V*cos(theta);
    adouble dydt = V*sin(theta);

    // Identically V^2: see the note at the head of this file.
    adouble L =  pow(dxdt,2.0) + pow(dydt,2.0);

    return  L;
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
// with no reference whatever to the returned states.

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
   printf("  closest approach, obs 1  : %.6f   (required >= 0.100000)\n", d1min);
   printf("  closest approach, obs 2  : %.6f   (required >= 0.100000)\n", d2min);
   printf("\n");
   printf("  The two closest-approach figures are the squared distances used by the\n");
   printf("  path constraints, evaluated all along the propagated path rather than at\n");
   printf("  the nodes only. The objective cannot be used as a check here, since it is\n");
   printf("  the same constant for every admissible control.\n\n");
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

    problem.name        		= "Obstacle avoidance problem";
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
    // A single fixed mesh: the residual box, not mesh refinement, is what controls
    // the accuracy here, and the error estimate is already about 3e-8 on this mesh.
    problem.phases(1).nodes                     = (RowVectorXi(1) << 81).finished();

    psopt_level2_setup(problem, algorithm);


////////////////////////////////////////////////////////////////////////////
///////////////////  Enter problem bounds information //////////////////////
////////////////////////////////////////////////////////////////////////////

    double xL = 0.0;
    double yL = 0.0;
    double xU = 2.0;
    double yU = 2.0;

    // theta is a heading, so these bounds admit about three revolutions and do not
    // single out a branch: see the note at the head of this file.
    double thetaL = -10.0;
    double thetaU = 10.0;

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

    problem.phases(1).bounds.lower.EndTime      = 1.0;
    problem.phases(1).bounds.upper.EndTime      = 1.0;



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

    // The straight line from (0,0) to (1.2,1.6) is 2.0 long and the particle covers
    // 2.138 in the available time, so no straight path is admissible and a guess of
    // zero heading throughout tells the solver nothing about which way to turn. The
    // guess below instead sweeps the heading from 0.5 to 1.4 radians and integrates
    // the kinematics: it goes up and to the right, which is the only sensible way
    // round, and nothing about it is fitted to the answer. Started from the zero
    // heading the solve here fails, which is a fair reflection of a problem with no
    // objective to guide the iterate.

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
    algorithm.nlp_iter_max                = 3000;
    algorithm.nlp_tolerance               = 1.e-6;
    algorithm.nlp_method                  = "IPOPT";
    algorithm.scaling                     = "automatic";
    algorithm.derivatives                 = "automatic";

    // Integrated residual transcription with a residual box: minimise the cost
    // subject to the residual of the differential equations being no larger than
    // ir_residual_bound at ir_residual_nodes Gauss-Legendre points inside each
    // interval. Bounding the residual between the nodes is what rules out the
    // control that local collocation returns here.
    algorithm.collocation_method          = "Hermite-Simpson";
    algorithm.transcription_method        = "integrated-residual";
    algorithm.ir_objective                = "cost";
    algorithm.ir_residual_bound           = 1.e-6;
    algorithm.ir_residual_nodes           = 4;
    algorithm.mesh_refinement             = "manual";


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

    plot(t,mu, problem.name+": path constraint multipliers","t", "mu_1 mu_2");

    plot(t,lambda, problem.name+": costates","t", "lambda_1 lambda_2");


}

////////////////////////////////////////////////////////////////////////////
///////////////////////      END OF FILE     ///////////////////////////////
////////////////////////////////////////////////////////////////////////////
