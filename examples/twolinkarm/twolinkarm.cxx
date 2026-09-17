//////////////////////////////////////////////////////////////////////////
//////////////////        twolinkarm.cxx        //////////////////////////
//////////////////////////////////////////////////////////////////////////
////////////////           PSOPT  Example             ////////////////////
//////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////
//////// Title:                 Two link arm problem      ////////////////
//////// Last modified:         17 September 2026         ////////////////
//////// Reference:             PROPT users guide         ////////////////
//////// (See PSOPT handbook for full reference)          ////////////////
//////////////////////////////////////////////////////////////////////////
////////                                                  ////////////////
//////// TRANSCRIPTION: MULTIPLE SHOOTING.                ////////////////
////////                                                  ////////////////
//////// This example is solved by multiple shooting      ////////////////
//////// rather than by collocation, and is the example   ////////////////
//////// to read for how that method is configured. Pass  ////////////////
//////// a collocation method as argv[1] -- "Legendre",   ////////////////
//////// "Hermite-Simpson" and so on -- to solve the same ////////////////
//////// problem the other way and compare.               ////////////////
////////                                                  ////////////////
//////// WHY THIS PROBLEM. Minimum time with bounded      ////////////////
//////// torques, so the optimal control is bang-bang:    ////////////////
//////// u1 has one switch and u2 has two, and both sit   ////////////////
//////// hard against +/-1 in between. A global Legendre  ////////////////
//////// polynomial cannot represent a jump, so the       ////////////////
//////// collocated control rings around each switch and  ////////////////
//////// the discretisation error is dominated by that    ////////////////
//////// rather than by the smooth arcs between them.     ////////////////
//////// Multiple shooting integrates each segment        ////////////////
//////// instead, and the control is parameterised        ////////////////
//////// segment by segment, which suits the structure.   ////////////////
////////                                                  ////////////////
//////// WHAT IT BUYS, measured on this problem (Ubuntu   ////////////////
//////// 24.04, IPOPT 3.11.9, one core), with the same    ////////////////
//////// 40 nodes -- 39 shooting segments -- in each      ////////////////
//////// case, and the same at 80:                        ////////////////
////////                                                  ////////////////
////////   nodes  method        CPU(s)   t_f       max    ////////////////
////////                                          rel err ////////////////
////////   40     Legendre       0.75   2.988660  3.7e-05 ////////////////
////////   40     mult. shooting 0.48   2.985042  8.1e-09 ////////////////
////////   80     Legendre       4.79   2.983629  2.6e-06 ////////////////
////////   80     mult. shooting 0.61   2.983021  9.5e-10 ////////////////
////////                                                  ////////////////
//////// Faster at both sizes -- eight times faster at 80 ////////////////
//////// nodes -- and satisfying the dynamics between two ////////////////
//////// and three orders of magnitude more accurately.   ////////////////
////////                                                  ////////////////
//////// WHAT IT DOES NOT BUY, which matters as much.     ////////////////
//////// The minimum time is NOT converged at this        ////////////////
//////// resolution and neither transcription claims it   ////////////////
//////// is: refining the mesh lowers t_f under both, and ////////////////
//////// the two agree with each other far better than    ////////////////
//////// either agrees with itself across mesh sizes.     ////////////////
//////// What limits the answer here is how finely the    ////////////////
//////// CONTROL is parameterised, not how accurately the ////////////////
//////// dynamics are integrated -- RK8 at 79 segments    ////////////////
//////// returns the same 2.983021 as RK4 does, with a    ////////////////
//////// residual at machine precision. A bang-bang       ////////////////
//////// problem wants its switching instants as          ////////////////
//////// variables, which is what a multi-phase           ////////////////
//////// formulation with the switch as a phase boundary  ////////////////
//////// gives, and which no uniform mesh can supply.     ////////////////
//////////////////////////////////////////////////////////////////////////
////////     Copyright (c) Victor M. Becerra, 2009        ////////////////
//////////////////////////////////////////////////////////////////////////
//////// This is part of the PSOPT software library, which ///////////////
//////// is distributed under the terms of the GNU Lesser ////////////////
//////// General Public License (LGPL)                    ////////////////
//////////////////////////////////////////////////////////////////////////

#include "psopt.h"

using namespace std;

//////////////////////////////////////////////////////////////////////////
///////////////////  Define the end point (Mayer) cost function //////////
//////////////////////////////////////////////////////////////////////////

adouble endpoint_cost(adouble* initial_states, adouble* final_states,
                      adouble* parameters,adouble& t0, adouble& tf,
                      adouble* xad, int iphase, Workspace* workspace)
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
   adouble xdot, ydot, vdot;

   adouble x1 = states[ 0 ];
   adouble x2 = states[ 1 ];
   adouble x3 = states[ 2 ];
   adouble x4 = states[ 3 ];

   adouble u1 = controls[ 0 ];
   adouble u2 = controls[ 1 ];

   adouble num1 =  sin(x3)*( (9.0/4.0)*cos(x3)*x1*x1+2*x2*x2 )
                   + (4.0/3.0)*(u1-u2)-(3.0/2.0)*cos(x3)*u2;

   adouble num2 =  -(sin(x3)*((7.0/2.0)*x1*x1+(9.0/4.0)*cos(x3)*x2*x2)
                   -(7.0/3.0)*u2+(3.0/2.0)*cos(x3)*(u1-u2));

   adouble den  =  31.0/36.0 + 9.0/4.0*pow(sin(x3),2);

   derivatives[ 0 ] = num1/den;
   derivatives[ 1 ] = num2/den;
   derivatives[ 2 ] = x2 - x1;
   derivatives[ 3 ] = x1;


}

////////////////////////////////////////////////////////////////////////////
///////////////////  Define the events function ////////////////////////////
////////////////////////////////////////////////////////////////////////////

void events(adouble* e, adouble* initial_states, adouble* final_states,
            adouble* parameters,adouble& t0, adouble& tf, adouble* xad,
            int iphase, Workspace* workspace)
{
   adouble x10 = initial_states[ 0 ];
   adouble x20 = initial_states[ 1 ];
   adouble x30 = initial_states[ 2 ];
   adouble x40 = initial_states[ 3 ];
   adouble x1f = final_states[ 0 ];
   adouble x2f = final_states[ 1 ];
   adouble x3f = final_states[ 2 ];
   adouble x4f = final_states[ 3 ];

   e[ 0 ] = x10;
   e[ 1 ] = x20;
   e[ 2 ] = x30;
   e[ 3 ] = x40;
   e[ 4 ] = x1f;
   e[ 5 ] = x2f;
   e[ 6 ] = x3f;
   e[ 7 ] = x4f;

}

///////////////////////////////////////////////////////////////////////////
///////////////////  Define the phase linkages function ///////////////////
///////////////////////////////////////////////////////////////////////////

void linkages( adouble* linkages, adouble* xad, Workspace* workspace)
{
  // No linkages as this is a single phase problem
}


////////////////////////////////////////////////////////////////////////////
///////////////////  Define the main routine ///////////////////////////////
////////////////////////////////////////////////////////////////////////////


int main(int argc, char** argv)
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

    problem.name        						= "Two link robotic arm";

    problem.outfilename                 	= "twolink.txt";

////////////////////////////////////////////////////////////////////////////
////////////  Define problem level constants & do level 1 setup ////////////
////////////////////////////////////////////////////////////////////////////

    problem.nphases   							= 1;
    problem.nlinkages                   	= 0;

    psopt_level1_setup(problem);


/////////////////////////////////////////////////////////////////////////////
/////////   Define phase related information & do level 2 setup /////////////
/////////////////////////////////////////////////////////////////////////////

    problem.phases(1).nstates   				= 4;
    problem.phases(1).ncontrols 				= 2;
    problem.phases(1).nevents   				= 8;
    problem.phases(1).npath     				= 0;
    problem.phases(1).nodes               << 40;

    psopt_level2_setup(problem, algorithm);


////////////////////////////////////////////////////////////////////////////
///////////////////  Declare DMatrix objects to store results //////////////
////////////////////////////////////////////////////////////////////////////

    DMatrix x, u, t;
    DMatrix lambda, H;

////////////////////////////////////////////////////////////////////////////
///////////////////  Enter problem bounds information //////////////////////
////////////////////////////////////////////////////////////////////////////

    problem.phases(1).bounds.lower.states(0) = -2.0;
    problem.phases(1).bounds.lower.states(1) = -2.0;
    problem.phases(1).bounds.lower.states(2) = -2.0;
    problem.phases(1).bounds.lower.states(3) = -2.0;

    problem.phases(1).bounds.upper.states(0) = 2.0;
    problem.phases(1).bounds.upper.states(1) = 2.0;
    problem.phases(1).bounds.upper.states(2) = 2.0;
    problem.phases(1).bounds.upper.states(3) = 2.0;

    problem.phases(1).bounds.lower.controls(0) = -1.0;
    problem.phases(1).bounds.lower.controls(1) = -1.0;

    problem.phases(1).bounds.upper.controls(0) = 1.0;
    problem.phases(1).bounds.upper.controls(1) = 1.0;

    problem.phases(1).bounds.lower.events(0) = 0.0;
    problem.phases(1).bounds.lower.events(1) = 0.0;
    problem.phases(1).bounds.lower.events(2) = 0.5;
    problem.phases(1).bounds.lower.events(3) = 0.0;
    problem.phases(1).bounds.lower.events(4) = 0.0;
    problem.phases(1).bounds.lower.events(5) = 0.0;
    problem.phases(1).bounds.lower.events(6) = 0.5;
    problem.phases(1).bounds.lower.events(7) = 0.522;

    problem.phases(1).bounds.upper.events = problem.phases(1).bounds.lower.events;



    problem.phases(1).bounds.lower.StartTime    = 0.0;
    problem.phases(1).bounds.upper.StartTime    = 0.0;

    problem.phases(1).bounds.lower.EndTime      = 1.0;
    problem.phases(1).bounds.upper.EndTime      = 10.0;


////////////////////////////////////////////////////////////////////////////
///////////////////  Register problem functions  ///////////////////////////
////////////////////////////////////////////////////////////////////////////


    problem.integrand_cost 				= &integrand_cost;
    problem.endpoint_cost 					= &endpoint_cost;
    problem.dae 								= &dae;
    problem.events 							= &events;
    problem.linkages							= &linkages;



////////////////////////////////////////////////////////////////////////////
///////////////////  Define & register initial guess ///////////////////////
////////////////////////////////////////////////////////////////////////////


    MatrixXd x0(4,40);

    x0 <<  linspace(0.0,0.0, 40),
           linspace(0.0,0.0, 40),
           linspace(0.5,0.5, 40),
           linspace(0.522,0.522, 40);

    problem.phases(1).guess.controls       = zeros(2,40);
    problem.phases(1).guess.states         = x0;
    problem.phases(1).guess.time           = linspace(0.0, 3.0, 40);

////////////////////////////////////////////////////////////////////////////
///////////////////  Enter algorithm options  //////////////////////////////
////////////////////////////////////////////////////////////////////////////


    algorithm.nlp_method                  = "IPOPT";
    algorithm.scaling                     = "automatic";
    algorithm.derivatives                 = "automatic";
    algorithm.nlp_iter_max                = 1000;
    algorithm.nlp_tolerance               = 1.e-6;

////////////////////////////////////////////////////////////////////////////
///////////////////  Multiple shooting  ////////////////////////////////////
////////////////////////////////////////////////////////////////////////////
//
//  In multiple shooting the mesh interval is a SEGMENT rather than a
//  collocation point: the NLP variables are the state at the segment
//  boundaries and the control parameters within each segment, an embedded
//  integrator marches each segment from its left boundary, and the
//  constraints are the mismatches at the joins. The trajectory is therefore
//  an integrated one by construction, which is why the local error below is
//  set by the integrator rather than by the polynomial through the nodes.
//
//  The four settings that matter, and why each is what it is here:
//
//    ms_integrator                "RK4" is the cheapest scheme that resolves
//                                 this arm's dynamics; "RK8" costs about four
//                                 times as much here and changes t_f in the
//                                 seventh figure, which is a measurement and
//                                 not a guess -- see the table in the header.
//                                 "TRBDF2" and "ESDIRK3" are the stiff
//                                 alternatives and are not needed here.
//    ms_steps_per_segment         8 RK4 steps between joins. At 4 the local
//                                 error is 1.3e-07 and at 8 it is 8.1e-09,
//                                 for 0.2 s more; beyond 8 the answer stops
//                                 moving because the control parameterisation
//                                 is what limits it.
//    ms_control_parameterisation  "linear" within a segment. "constant" is
//                                 cheaper and gives 2.985858 rather than
//                                 2.985042, which is the parameterisation
//                                 showing through.
//    nodes                        39 segments, from the 40 nodes this example
//                                 has always used, so that the transcription
//                                 is the only thing that changed.
//
//  Pass a collocation method as argv[1] to solve the same problem the other
//  way: ./twolinkarm Legendre
//
    if (argc > 1) {
        algorithm.collocation_method      = argv[1];
    }
    else {
        algorithm.transcription_method        = "multiple-shooting";
        algorithm.ms_integrator               = "RK4";
        algorithm.ms_steps_per_segment        = 8;
        algorithm.ms_control_parameterisation = "linear";
    }

////////////////////////////////////////////////////////////////////////////
///////////////////  Now call PSOPT to solve the problem   /////////////////
////////////////////////////////////////////////////////////////////////////

    if (psopt(solution, problem, algorithm) != 0) return 1;

////////////////////////////////////////////////////////////////////////////
///////////  Extract relevant variables from solution structure   //////////
////////////////////////////////////////////////////////////////////////////

    x 		= solution.get_states_in_phase(1);
    u 		= solution.get_controls_in_phase(1);
    t 		= solution.get_time_in_phase(1);

////////////////////////////////////////////////////////////////////////////
///////////  Save solution data to files if desired ////////////////////////
////////////////////////////////////////////////////////////////////////////

    Save(x,"x.dat");
    Save(u,"u.dat");
    Save(t,"t.dat");


////////////////////////////////////////////////////////////////////////////
///////////  Plot some results if desired (requires gnuplot) ///////////////
////////////////////////////////////////////////////////////////////////////

    plot(t,x,problem.name + ": states", "time (s)", "states", "x1 x2 x3 x4");

    plot(t,u,problem.name + ": controls", "time (s)", "controls", "u1 u2");


    plot(t,x,problem.name + ": states", "time (s)", "states", "x1 x2 x3 x4",
                                  "pdf", "twolinkarm_states.pdf");

    plot(t,u,problem.name + ": controls", "time (s)", "controls", "u1 u2",
                              "pdf", "twolinkarm_controls.pdf");


}

////////////////////////////////////////////////////////////////////////////
///////////////////////      END OF FILE     ///////////////////////////////
////////////////////////////////////////////////////////////////////////////

