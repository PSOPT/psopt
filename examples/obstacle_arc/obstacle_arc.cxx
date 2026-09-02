//////////////////////////////////////////////////////////////////////////
////////////////      obstacle_arc.cxx               /////////////////////
//////////////////////////////////////////////////////////////////////////
////////////////           PSOPT  Example             ////////////////////
//////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////
//////// Title:   Minimum time obstacle avoidance, as three ///////////////
////////          phases with the boundary arc in its own   ///////////////
//////// Last modified: 02 September 2026                 ////////////////
//////////////////////////////////////////////////////////////////////////
////////     Copyright (c) Victor M. Becerra, 2026        ////////////////
//////////////////////////////////////////////////////////////////////////
//////// This is part of the PSOPT software library, which ///////////////
//////// is distributed under the terms of the GNU Lesser ////////////////
//////// General Public License (LGPL)                    ////////////////
//////////////////////////////////////////////////////////////////////////

// The same problem as examples/obstacle: a particle at fixed speed V = 2.138 with
// heading theta as its only control, from (0,0) to (1.2,1.6) in minimum time,
// keeping outside two discs of radius sqrt(0.1) centred at (0.4,0.5) and (0.8,1.5).
//
// examples/obstacle solves it in one phase and knows nothing about the structure of
// the answer. This file solves the same problem knowing that structure, which is
// how a state-constrained problem is normally handled once its structure is known:
// the solution leaves the start on a straight line, runs along the boundary of the
// first obstacle for a while, and leaves on another straight line. The constrained
// arc is given a phase of its own, in which the path constraint is imposed as an
// EQUALITY rather than an inequality, and the two junction times are free variables
// determined by the solve.
//
// The point of doing so is that the difficulty in the one-phase version is entirely
// at the two junctions. There the optimal heading has a corner -- it is constant,
// then rises at the constant rate V/r, then is constant again -- and a mesh that
// does not have a node at the corner smooths it and pays for that in the objective.
// One phase with automatic mesh refinement reaches 3.0e-9 relative on 74 nodes;
// uniform meshes do far worse, 3.0e-8 at 201 nodes. Splitting at the junctions
// removes the corners from the interior of every phase: phases 1 and 3 have a
// constant heading and a straight trajectory, which Hermite-Simpson represents
// exactly, and phase 2 has a heading linear in t. Only the states of phase 2 are
// not polynomial, being a circular arc, and that is the whole of the remaining
// discretization error. It converges at the fourth order of the method, with the
// arc phase carrying 11, 21, 41, 81 and 161 nodes and the two straight phases 11
// each:
//
//     arc nodes    tf, rel.     first junction   turn rate on the arc
//        11        6.66e-09        3.00e-05           1.32e-04
//        21        4.16e-10        7.48e-06           3.31e-05
//        41        2.60e-11        1.86e-06           8.27e-06
//        81        1.63e-12        4.58e-07           2.07e-06
//       161        9.02e-14        1.31e-07           6.03e-07
//
// The objective falls by a factor of sixteen per halving, which is the order of
// Hermite-Simpson; the junction times and the turn rate, being read off the ends of
// the arc rather than integrated over it, converge at second order. The setting
// below is the third line: 63 nodes in all, a third of a second, and a minimum time
// correct to eleven significant figures against a value obtained without a
// discretization at all. The one-phase version of the same problem manages nine on
// 74 nodes. An optional argument changes the arc node count so the rest of the table
// can be reproduced without editing this file.
//
// The junction conditions are the interesting part. Continuity of position at each
// junction is the obvious requirement. Continuity of the HEADING is the other one,
// and it is what makes the arc a boundary arc rather than a corner: for a state
// constraint of first order -- one differentiation of
// g = (x-0.4)^2 + (y-0.5)^2 - 0.1 exposes the control -- the trajectory must enter
// and leave the boundary tangentially, and here that is exactly theta being
// continuous. Imposing it as a linkage is what pins the junction times; without it
// the phase times are free to slide and the problem is degenerate.
//
// The answer is known in closed form. The shortest path avoiding the two discs is
// a taut string, and the enumeration in the header of examples/obstacle gives
//
//      L*  = 2.131722291465491    (only the first obstacle active)
//      tf* = L*/V = 0.997063747177498
//
// with the structure
//
//      straight   theta = 0.379526559413    t = 0        -> 0.260419287317
//      arc        theta rises at V/r        t = 0.260419 -> 0.378319287153
//                 = 6.760949637440, by 0.797115961148
//      straight   theta = 1.176642520561    t = 0.378319 -> 0.997063747177
//
// and the tangency points (0.517156303191, 0.206274957447) and
// (0.691980007889, 0.378559994263). The example prints its junction times, its two
// straight-line headings and its turn rate next to all of these, so that every
// element of the structure is checked and not only the objective.
//
// The initial guess is deliberately crude: junction times of 0.3 and 0.5, a final
// time of 1, and a heading swept linearly from 0.5 to 1.4 with the kinematics
// integrated. None of the numbers above is used to start the solve.

#include "psopt.h"

using namespace std;

using namespace PSOPT;

static const double V_SPEED = 2.138;
static const double R2_OBS    = 0.1;                // squared obstacle radius
static const double THETA_TOL = 1.0e-6;             // tolerance on "the heading is constant"

// The closed-form answer, for the comparison printed after the solve.
static const double TF_EXACT     = 0.997063747177498;
static const double T1_EXACT     = 0.260419287317;
static const double T2_EXACT     = 0.378319287153;
static const double THETA1_EXACT = 0.379526559413;
static const double THETA2_EXACT = 1.176642520561;
static const double RATE_EXACT   = 6.760949637440;  // V/r on the arc

//////////////////////////////////////////////////////////////////////////
///////////////////  Define the end point (Mayer) cost function //////////
//////////////////////////////////////////////////////////////////////////

adouble endpoint_cost(adouble* initial_states, adouble* final_states,
                      adouble* parameters,adouble& t0, adouble& tf,
                      adouble* xad, int iphase,Workspace* workspace)
{
   // The objective is the time at which the last phase ends.
   if (iphase == 3) return tf;
   return 0.0;
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

// The dynamics are the same in all three phases. What differs is how the path
// constraints are used: on phase 2 the first one is an equality, which is what puts the
// trajectory on the boundary of the first obstacle, and phases 1 and 3 carry a third
// one that holds the heading constant.
//
// That third constraint is not an assumption about the answer, it is a consequence of
// the maximum principle for this problem. Off the obstacle boundary the Hamiltonian
// V(lambda_1 cos theta + lambda_2 sin theta) does not depend on the state, so the
// costates are constant there, and the heading that minimises it is constant with them:
// the unconstrained arcs of the optimal trajectory are straight lines. Saying so in the
// model does two things. It makes phases 1 and 3 exactly representable by the
// discretization, since a straight line at constant speed is linear in t. And it removes
// a degeneracy: without it, phase 1 is free to ride the obstacle boundary as well --
// its inequality constraint permits that -- so the first junction can sit anywhere along
// the arc, and although the total time comes out right the junction time does not.

void dae(adouble* derivatives, adouble* path, adouble* states,
         adouble* controls, adouble* parameters, adouble& time,
         adouble* xad, int iphase, Workspace* workspace)
{
   adouble x     = states[ 0 ];
   adouble y     = states[ 1 ];
   adouble theta = controls[ 0 ];

   derivatives[ 0 ] = V_SPEED*cos(theta);
   derivatives[ 1 ] = V_SPEED*sin(theta);

   path[ 0 ] = pow(x-0.4,2.0) + pow(y-0.5,2.0);
   path[ 1 ] = pow(x-0.8,2.0) + pow(y-1.5,2.0);

   if ( iphase != 2 ) path[ 2 ] = theta - parameters[ 0 ];
}

////////////////////////////////////////////////////////////////////////////
///////////////////  Define the events function ////////////////////////////
////////////////////////////////////////////////////////////////////////////

void events(adouble* e, adouble* initial_states, adouble* final_states,
            adouble* parameters,adouble& t0, adouble& tf, adouble* xad,
            int iphase, Workspace* workspace)
{
   if (iphase == 1) {
      e[ 0 ] = initial_states[ 0 ];
      e[ 1 ] = initial_states[ 1 ];
   }
   else if (iphase == 3) {
      e[ 0 ] = final_states[ 0 ];
      e[ 1 ] = final_states[ 1 ];
   }
}


///////////////////////////////////////////////////////////////////////////
///////////////////  Define the phase linkages function ///////////////////
///////////////////////////////////////////////////////////////////////////

// Two junctions, each carrying continuity of x, y, t and theta: eight linkages, plus a
// ninth that keeps the arc phase from collapsing.
//
// The heading linkages are the junction conditions of the boundary arc. The path
// constraint is of first order, so the trajectory must join the boundary
// tangentially; theta continuous is that condition written in the variables of this
// problem. They are also what determines the junction times: with position and time
// continuity alone the pair (t1, t2) could slide along the arc, and the solve is
// degenerate.

void linkages( adouble* linkages, adouble* xad, Workspace* workspace)
{
   adouble xf[2], xi[2], uf[1], ui[1], tfa, t0b;

   // Phases 1 and 2
   get_final_states  ( xf, xad, 1, workspace );
   get_initial_states( xi, xad, 2, workspace );
   get_final_controls  ( uf, xad, 1, workspace );
   get_initial_controls( ui, xad, 2, workspace );
   tfa = get_final_time  ( xad, 1, workspace );
   t0b = get_initial_time( xad, 2, workspace );

   linkages[ 0 ] = xf[ 0 ] - xi[ 0 ];
   linkages[ 1 ] = xf[ 1 ] - xi[ 1 ];
   linkages[ 2 ] = tfa     - t0b;
   linkages[ 3 ] = uf[ 0 ] - ui[ 0 ];

   // Phases 2 and 3
   get_final_states  ( xf, xad, 2, workspace );
   get_initial_states( xi, xad, 3, workspace );
   get_final_controls  ( uf, xad, 2, workspace );
   get_initial_controls( ui, xad, 3, workspace );
   tfa = get_final_time  ( xad, 2, workspace );
   t0b = get_initial_time( xad, 3, workspace );

   linkages[ 4 ] = xf[ 0 ] - xi[ 0 ];
   linkages[ 5 ] = xf[ 1 ] - xi[ 1 ];
   linkages[ 6 ] = tfa     - t0b;
   linkages[ 7 ] = uf[ 0 ] - ui[ 0 ];

   // The duration of the arc phase. Nothing above prevents t2 = t1, and a phase of zero
   // length satisfies every continuity condition trivially, which is what the solver
   // does if allowed: the arc disappears, the two straight segments meet at a corner,
   // and the problem becomes infeasible. Bounding this below asserts what the three-phase
   // formulation assumes in the first place, that there is an arc.
   linkages[ 8 ] = get_final_time( xad, 2, workspace ) - get_initial_time( xad, 2, workspace );
}


////////////////////////////////////////////////////////////////////////////
///////////////////  Define the main routine ///////////////////////////////
////////////////////////////////////////////////////////////////////////////

int main(int argc, char* argv[])
{
// An optional argument sets the number of nodes on the arc phase, so that the
// convergence of the junction times and the turn rate can be watched without editing
// the source:   ./obstacle_arc [nodes]
    int arc_nodes = (argc > 1) ? atoi(argv[1]) : 41;


////////////////////////////////////////////////////////////////////////////
///////////////////  Declare key structures ////////////////////////////////
////////////////////////////////////////////////////////////////////////////

    Alg  algorithm;
    Sol  solution;
    Prob problem;

////////////////////////////////////////////////////////////////////////////
///////////////////  Register problem name  ////////////////////////////////
////////////////////////////////////////////////////////////////////////////

    problem.name                        = "Minimum time obstacle avoidance with a boundary arc";
    problem.outfilename                 = "obstacle_arc.txt";

////////////////////////////////////////////////////////////////////////////
////////////  Define problem level constants & do level 1 setup ////////////
////////////////////////////////////////////////////////////////////////////

    problem.nphases                     = 3;
    problem.nlinkages                   = 9;

    psopt_level1_setup(problem);

/////////////////////////////////////////////////////////////////////////////
/////////   Define phase related information & do level 2 setup  ////////////
/////////////////////////////////////////////////////////////////////////////

    // Phases 1 and 3 carry a straight line at constant heading, which
    // Hermite-Simpson represents exactly, so they need almost no nodes. Phase 2
    // carries a circular arc, whose states are not polynomial, and is the only
    // source of discretization error left; it is given the most nodes for that
    // reason. There is no mesh refinement: there is nothing left for it to find.
    for (int i = 1; i <= 3; i++) {
        problem.phases(i).nstates   = 2;
        problem.phases(i).ncontrols = 1;
    }
    // Phases 1 and 3 carry a third path constraint, theta - p = 0, with p a static
    // parameter of the phase: the constant heading of the straight arc.
    problem.phases(1).npath      = 3;
    problem.phases(2).npath      = 2;
    problem.phases(3).npath      = 3;
    problem.phases(1).nparameters = 1;
    problem.phases(2).nparameters = 0;
    problem.phases(3).nparameters = 1;
    problem.phases(1).nevents = 2;
    problem.phases(2).nevents = 0;
    problem.phases(3).nevents = 2;

    problem.phases(1).nodes    << 11;
    problem.phases(2).nodes    << arc_nodes;
    problem.phases(3).nodes    << 11;

    psopt_level2_setup(problem, algorithm);


////////////////////////////////////////////////////////////////////////////
///////////////////  Enter problem bounds information //////////////////////
////////////////////////////////////////////////////////////////////////////

    double x0 = 0.0, y0 = 0.0, xf = 1.2, yf = 1.6;

    for (int i = 1; i <= 3; i++) {
        problem.phases(i).bounds.lower.states   <<  0.0, 0.0;
        problem.phases(i).bounds.upper.states   <<  2.0, 2.0;
        // theta is a heading, defined modulo 2 pi; a single branch, as in
        // examples/obstacle.
        problem.phases(i).bounds.lower.controls << -pi/2.0;
        problem.phases(i).bounds.upper.controls <<  pi;
    }

    // Phases 1 and 3: stay outside both obstacles, and hold the heading constant.
    // The heading is required constant to within THETA_TOL rather than exactly. With an
    // exact equality at every node and midpoint the NLP is square -- the three-phase
    // structure, its tangency conditions and the boundary conditions determine the answer
    // without reference to the objective at all -- and IPOPT declines a problem with no
    // degrees of freedom. A band of 1e-6 radians leaves the rows as inequalities. It costs
    // nothing: a heading wrong by 1e-6 over a segment of length 1 moves the arrival by
    // 1e-6 and the time by parts in 1e13.
    problem.phases(1).bounds.lower.path <<  R2_OBS, R2_OBS, -THETA_TOL;
    problem.phases(1).bounds.upper.path <<  100.0,  100.0,   THETA_TOL;
    problem.phases(3).bounds.lower.path <<  R2_OBS, R2_OBS, -THETA_TOL;
    problem.phases(3).bounds.upper.path <<  100.0,  100.0,   THETA_TOL;

    problem.phases(1).bounds.lower.parameters << -pi/2.0;
    problem.phases(1).bounds.upper.parameters <<  pi;
    problem.phases(3).bounds.lower.parameters << -pi/2.0;
    problem.phases(3).bounds.upper.parameters <<  pi;

    // Phase 2: ON the boundary of the first obstacle, and outside the second.
    problem.phases(2).bounds.lower.path <<  R2_OBS, R2_OBS;
    problem.phases(2).bounds.upper.path <<  R2_OBS, 100.0;

    problem.phases(1).bounds.lower.events << x0, y0;
    problem.phases(1).bounds.upper.events << x0, y0;
    problem.phases(3).bounds.lower.events << xf, yf;
    problem.phases(3).bounds.upper.events << xf, yf;

    // t0 is fixed at zero; the two junction times and the final time are free, with
    // wide brackets that do not encode the answer.
    problem.phases(1).bounds.lower.StartTime = 0.0;
    problem.phases(1).bounds.upper.StartTime = 0.0;
    problem.phases(1).bounds.lower.EndTime   = 0.02;
    problem.phases(1).bounds.upper.EndTime   = 2.0;

    problem.phases(2).bounds.lower.StartTime = 0.02;
    problem.phases(2).bounds.upper.StartTime = 2.0;
    problem.phases(2).bounds.lower.EndTime   = 0.04;
    problem.phases(2).bounds.upper.EndTime   = 2.5;

    problem.phases(3).bounds.lower.StartTime = 0.04;
    problem.phases(3).bounds.upper.StartTime = 2.5;
    problem.phases(3).bounds.lower.EndTime   = 0.5;
    problem.phases(3).bounds.upper.EndTime   = 3.0;

    // The first eight linkages are continuity conditions; the ninth is the duration of
    // the arc phase, which must be positive and is otherwise unconstrained.
    problem.bounds.lower.linkage = zeros(9,1);
    problem.bounds.upper.linkage = zeros(9,1);
    problem.bounds.lower.linkage(8) = 0.02;
    problem.bounds.upper.linkage(8) = 2.0;


////////////////////////////////////////////////////////////////////////////
///////////////////  Register problem functions  ///////////////////////////
////////////////////////////////////////////////////////////////////////////

    problem.integrand_cost      = &integrand_cost;
    problem.endpoint_cost       = &endpoint_cost;
    problem.dae                 = &dae;
    problem.events              = &events;
    problem.linkages            = &linkages;

////////////////////////////////////////////////////////////////////////////
///////////////////  Define & register initial guess ///////////////////////
////////////////////////////////////////////////////////////////////////////

    // A crude guess built from the structural assumption and nothing else: a straight
    // run at a heading of 0.4 rad for 0.3 s, an arc on the first obstacle for 0.15 s,
    // and a straight run to the target ending at t = 1. The arc is started by projecting
    // the end of the first segment radially onto the obstacle boundary, which is the
    // generic way to put a guess on a constraint boundary; it is swept anticlockwise,
    // which is the sense in which the obstacle is passed. None of the closed-form
    // numbers in the header is used.
    {
       const double cx = 0.4, cy = 0.5, r = sqrt(R2_OBS);
       const double tA = 0.3, tB = 0.45, tC = 1.0;
       const double th_a = 0.4;
       const int    ng = 20;

       // end of phase 1, and the radial projection of it onto the circle
       const double ax = V_SPEED*cos(th_a)*tA, ay = V_SPEED*sin(th_a)*tA;
       const double dn = sqrt((ax-cx)*(ax-cx) + (ay-cy)*(ay-cy));
       const double phi0 = atan2(ay-cy, ax-cx);
       const double dphi = (V_SPEED/r)*(tB-tA);          // swept anticlockwise
       const double bx = cx + r*cos(phi0+dphi), by = cy + r*sin(phi0+dphi);
       const double th_c = atan2(yf-by, xf-bx);
       (void) dn;

       MatrixXd xg(2,ng), ug(1,ng), tt(1,ng);

       for (int k = 0; k < ng; k++) {                    // phase 1: straight
           double tk = tA*k/(double)(ng-1);
           tt(0,k) = tk; ug(0,k) = th_a;
           xg(0,k) = V_SPEED*cos(th_a)*tk;
           xg(1,k) = V_SPEED*sin(th_a)*tk;
       }
       problem.phases(1).guess.time = tt;
       problem.phases(1).guess.controls = ug;
       problem.phases(1).guess.states = xg;
       problem.phases(1).guess.parameters = th_a*ones(1,1);

       for (int k = 0; k < ng; k++) {                    // phase 2: on the boundary
           double f  = k/(double)(ng-1);
           double tk = tA + (tB-tA)*f;
           double ph = phi0 + dphi*f;
           tt(0,k) = tk;
           xg(0,k) = cx + r*cos(ph);
           xg(1,k) = cy + r*sin(ph);
           ug(0,k) = atan2(cos(ph), -sin(ph));           // anticlockwise tangent
       }
       problem.phases(2).guess.time = tt;
       problem.phases(2).guess.controls = ug;
       problem.phases(2).guess.states = xg;

       for (int k = 0; k < ng; k++) {                    // phase 3: straight to the target
           double f  = k/(double)(ng-1);
           tt(0,k) = tB + (tC-tB)*f;
           xg(0,k) = bx + (xf-bx)*f;
           xg(1,k) = by + (yf-by)*f;
           ug(0,k) = th_c;
       }
       problem.phases(3).guess.time = tt;
       problem.phases(3).guess.controls = ug;
       problem.phases(3).guess.states = xg;
       problem.phases(3).guess.parameters = th_c*ones(1,1);
    }


////////////////////////////////////////////////////////////////////////////
///////////////////  Enter algorithm options  //////////////////////////////
////////////////////////////////////////////////////////////////////////////

    algorithm.nlp_iter_max                = 1000;
    algorithm.nlp_tolerance               = 1.e-8;
    algorithm.nlp_method                  = "IPOPT";
    algorithm.scaling                     = "automatic";
    algorithm.derivatives                 = "automatic";
    algorithm.collocation_method          = "Hermite-Simpson";


////////////////////////////////////////////////////////////////////////////
///////////////////  Now call PSOPT to solve the problem   /////////////////
////////////////////////////////////////////////////////////////////////////

    psopt(solution, problem, algorithm);

////////////////////////////////////////////////////////////////////////////
///////////  Extract relevant variables from solution structure   //////////
////////////////////////////////////////////////////////////////////////////

    MatrixXd x1 = solution.get_states_in_phase(1);
    MatrixXd x2 = solution.get_states_in_phase(2);
    MatrixXd x3 = solution.get_states_in_phase(3);
    MatrixXd u1 = solution.get_controls_in_phase(1);
    MatrixXd u2 = solution.get_controls_in_phase(2);
    MatrixXd u3 = solution.get_controls_in_phase(3);
    MatrixXd t1 = solution.get_time_in_phase(1);
    MatrixXd t2 = solution.get_time_in_phase(2);
    MatrixXd t3 = solution.get_time_in_phase(3);

    const int n1 = (int) t1.cols(), n2 = (int) t2.cols(), n3 = (int) t3.cols();
    MatrixXd x(2, n1+n2+n3), u(1, n1+n2+n3), t(1, n1+n2+n3);
    x << x1, x2, x3;
    u << u1, u2, u3;
    t << t1, t2, t3;

////////////////////////////////////////////////////////////////////////////
///////////  Compare the whole structure with the closed form //////////////
////////////////////////////////////////////////////////////////////////////

    {
       const double tj1  = t1(0, t1.cols()-1);          // first junction
       const double tj2  = t2(0, t2.cols()-1);          // second junction
       const double tend = t3(0, t3.cols()-1);
       // The heading is constant on phases 1 and 3; report its mean and its spread.
       const double th1  = u1.row(0).mean();
       const double th3  = u3.row(0).mean();
       const double sp1  = u1.maxCoeff() - u1.minCoeff();
       const double sp3  = u3.maxCoeff() - u3.minCoeff();
       // On the arc the heading is linear in t; the turn rate is its slope.
       const double rate = (u2(0,u2.cols()-1) - u2(0,0)) / (tj2 - tj1);

       printf("\n");
       printf("                          PSOPT                closed form           rel. diff\n");
       printf("  final time         %18.12f   %18.12f   %9.2e\n",
              tend, TF_EXACT, fabs(tend-TF_EXACT)/TF_EXACT);
       printf("  first junction     %18.12f   %18.12f   %9.2e\n",
              tj1, T1_EXACT, fabs(tj1-T1_EXACT)/T1_EXACT);
       printf("  second junction    %18.12f   %18.12f   %9.2e\n",
              tj2, T2_EXACT, fabs(tj2-T2_EXACT)/T2_EXACT);
       printf("  heading, phase 1   %18.12f   %18.12f   %9.2e\n",
              th1, THETA1_EXACT, fabs(th1-THETA1_EXACT)/THETA1_EXACT);
       printf("  heading, phase 3   %18.12f   %18.12f   %9.2e\n",
              th3, THETA2_EXACT, fabs(th3-THETA2_EXACT)/THETA2_EXACT);
       printf("  turn rate on arc   %18.12f   %18.12f   %9.2e\n",
              rate, RATE_EXACT, fabs(rate-RATE_EXACT)/RATE_EXACT);
       printf("\n");
       printf("  spread of the heading over phase 1: %.3e, over phase 3: %.3e\n", sp1, sp3);
       printf("  (both should be zero: the straight arcs carry a constant heading)\n\n");
    }

////////////////////////////////////////////////////////////////////////////
///////////  Save solution data to files if desired ////////////////////////
////////////////////////////////////////////////////////////////////////////

    Save(x,"x.dat");
    Save(u,"u.dat");
    Save(t,"t.dat");

////////////////////////////////////////////////////////////////////////////
///////////  Plot some results if desired (requires gnuplot) ///////////////
////////////////////////////////////////////////////////////////////////////

    MatrixXd alpha = linspace(0.0, 2*pi, 100);
    MatrixXd xObs1 = sqrt(R2_OBS)*cos(alpha) + 0.4*ones(1,length(alpha));
    MatrixXd yObs1 = sqrt(R2_OBS)*sin(alpha) + 0.5*ones(1,length(alpha));
    MatrixXd xObs2 = sqrt(R2_OBS)*cos(alpha) + 0.8*ones(1,length(alpha));
    MatrixXd yObs2 = sqrt(R2_OBS)*sin(alpha) + 1.5*ones(1,length(alpha));

    MatrixXd xx = x.row(0);
    MatrixXd yy = x.row(1);

    plot(xx,yy,xObs1,yObs1,xObs2,yObs2,problem.name+": x-y trajectory",
                                            "x", "y", "y obs1 obs2");

    plot(xx,yy,xObs1,yObs1,xObs2,yObs2,problem.name+": x-y trajectory",
                                            "x", "y", "y obs1 obs2",
                                            "pdf", "obstacle_arc_xy.pdf");

    plot(t,u, problem.name+": theta","t", "theta");

    plot(t,u, problem.name+": theta","t", "theta",
                                            "pdf", "obstacle_arc_theta.pdf");

}

////////////////////////////////////////////////////////////////////////////
///////////////////////      END OF FILE     ///////////////////////////////
////////////////////////////////////////////////////////////////////////////
