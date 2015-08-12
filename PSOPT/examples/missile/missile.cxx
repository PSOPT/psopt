//////////////////////////////////////////////////////////////////////////
////////////////           missile_min_time.cxx      /////////////////////
//////////////////////////////////////////////////////////////////////////
////////////////           PSOPT Example             /////////////////////
//////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////
//////// Title: 	Missile terminal burn maneouvre   ////////////////
////////                (minimum time version)            ////////////////
//////// Last modified: 12 July 2009                      ////////////////
//////// Reference:     Subchan and Zbikowski (2009)      ////////////////
////////                                                  ////////////////
//////////////////////////////////////////////////////////////////////////
////////     Copyright (c) Victor M. Becerra, 2009        ////////////////
//////////////////////////////////////////////////////////////////////////
//////// This is part of the PSOPT software library, which ///////////////
//////// is distributed under the terms of the GNU Lesser ////////////////
//////// General Public License (LGPL)                    ////////////////
//////////////////////////////////////////////////////////////////////////

#include "psopt.h"


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

adouble integrand_cost(adouble* states, adouble* controls, adouble* parameters,
                     adouble& time, adouble* xad, int iphase, Workspace* workspace)
{
   return 0;
}


//////////////////////////////////////////////////////////////////////////
///////////////////  Define the DAE's ////////////////////////////////////
//////////////////////////////////////////////////////////////////////////


void dae(adouble* derivatives, adouble* path, adouble* states,
         adouble* controls, adouble* parameters, adouble& time,
         adouble* xad, int iphase, Workspace* workspace)
{

    adouble T 		= controls[ CINDEX(1) ]; // Thrust
    adouble alpha       = controls[ CINDEX(2) ]; // Angle of attack

    adouble gamma 	= states[ CINDEX(1) ]; // flight path angle
    adouble V 		= states[ CINDEX(2) ]; // speed
    adouble x		= states[ CINDEX(3) ]; // horizontal position
    adouble h 		= states[ CINDEX(4) ]; // altitude

    double  m		= 1005.0;  // kg
    double  g  		= 9.81;    // m/s^2
    double  Sref  	= 0.3376;  // m^2
    double  A1		= -1.9431;
    double  A2 		= -0.1499;
    double  A3 		= 0.2359;
    double  B1 		= 21.9;
    double  B2		= 0.0;
    double  C1  	= 3.312e-9;
    double  C2         	= -1.142e-4;
    double  C3  	= 1.224;

    adouble sina 	= sin(alpha);
    adouble cosa 	= cos(alpha);
    adouble sing	= sin(gamma);
    adouble cosg        = cos(gamma);

    adouble Cd, D, L, Cl, rho;

    rho = C1*h*h + C2*h + C3;  // Air density

    Cd = A1*alpha*alpha + A2*alpha + A3;

    Cl = B1*alpha + B2;

    D = 0.5*Cd*rho*V*V*Sref; // Axial aerodynamic force

    L = 0.5*Cl*rho*V*V*Sref; // Normal aerodynamic force

    derivatives[ CINDEX(1) ]  = (T-D)/(m*V)*sina + L/(m*V)*cosa - g*cosg/V;
    derivatives[ CINDEX(2) ]  = (T-D)/m*cosa - L/m*sina - g*sing;
    derivatives[ CINDEX(3) ]  = V*cosg;
    derivatives[ CINDEX(4) ]  = V*sing;

    path[ CINDEX(1) ]  	      = L/(m*g);

    adouble H = smooth_heaviside(7500.0-x, 10.0 );
    path[ CINDEX(2) ]         =  H*(h-30.0)  +  (1.0-H)*h;


}

////////////////////////////////////////////////////////////////////////////
///////////////////  Define the events function ////////////////////////////
////////////////////////////////////////////////////////////////////////////

void events(adouble* e, adouble* initial_states, adouble* final_states,
            adouble* parameters,adouble& t0, adouble& tf, adouble* xad,
            int iphase, Workspace* workspace)

{

    adouble x1_i 	= initial_states[ CINDEX(1) ];
    adouble x2_i 	= initial_states[ CINDEX(2) ];
    adouble x3_i 	= initial_states[ CINDEX(3) ];
    adouble x4_i 	= initial_states[ CINDEX(4) ];

    adouble x1_f 	= final_states[ CINDEX(1) ];
    adouble x2_f	= final_states[ CINDEX(2) ];
    adouble x3_f	= final_states[ CINDEX(3) ];
    adouble x4_f 	= final_states[ CINDEX(4) ];



    e[ CINDEX(1) ] 	= x1_i;
    e[ CINDEX(2) ]      = x2_i;
    e[ CINDEX(3) ]      = x3_i;
    e[ CINDEX(4) ]    	= x4_i;


    e[ CINDEX(5) ]   	= x1_f;
    e[ CINDEX(6) ]      = x2_f;
    e[ CINDEX(7) ]  	= x3_f;
    e[ CINDEX(8) ]  	= x4_f;


}



///////////////////////////////////////////////////////////////////////////
///////////////////  Define the phase linkages function ///////////////////
///////////////////////////////////////////////////////////////////////////

void linkages( adouble* linkages, adouble* xad, Workspace* workspace)
{



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

    problem.name        		= "Missile problem";
    problem.outfilename                 = "missile.txt";

////////////////////////////////////////////////////////////////////////////
////////////  Define problem level constants & do level 1 setup ////////////
////////////////////////////////////////////////////////////////////////////

    problem.nphases   			= 1;
    problem.nlinkages                   = 0;

    psopt_level1_setup(problem);

/////////////////////////////////////////////////////////////////////////////
/////////   Define phase related information & do level 2 setup  ////////////
/////////////////////////////////////////////////////////////////////////////

    problem.phases(1).nstates   		= 4;
    problem.phases(1).ncontrols 		= 2;
    problem.phases(1).nevents   		= 8;
    problem.phases(1).npath     		= 2;
    problem.phases(1).nodes                     = "[20 80]";

    psopt_level2_setup(problem, algorithm);

////////////////////////////////////////////////////////////////////////////
///////////////////  Enter problem bounds information //////////////////////
////////////////////////////////////////////////////////////////////////////


    problem.phases(1).bounds.lower.controls(1) = 1000.0;
    problem.phases(1).bounds.lower.controls(2) = -0.3;

    problem.phases(1).bounds.upper.controls(1) = 6000.0;
    problem.phases(1).bounds.upper.controls(2) =  0.3;

    problem.phases(1).bounds.lower.states(1) = -pi/2.0;
    problem.phases(1).bounds.lower.states(2) = 200.0;
    problem.phases(1).bounds.lower.states(3) = 0.0;
    problem.phases(1).bounds.lower.states(4) = 0.0;

    problem.phases(1).bounds.upper.states(1) = pi/2.0;
    problem.phases(1).bounds.upper.states(2) = 310.0;
    problem.phases(1).bounds.upper.states(3) = 10000.0;
    problem.phases(1).bounds.upper.states(4) = 1900.0;

    problem.phases(1).bounds.lower.events(1) = 0.0;
    problem.phases(1).bounds.lower.events(2) = 272.0;
    problem.phases(1).bounds.lower.events(3) = 0.0;
    problem.phases(1).bounds.lower.events(4) = 30.0;

    problem.phases(1).bounds.lower.events(5) = -pi/2.0;
    problem.phases(1).bounds.lower.events(6) = 310.0;
    problem.phases(1).bounds.lower.events(7) = 10000.0;
    problem.phases(1).bounds.lower.events(8) = 0.0;

    problem.phases(1).bounds.upper.events=    problem.phases(1).bounds.lower.events;


    problem.phases(1).bounds.lower.path(1) = -4.0;
    problem.phases(1).bounds.upper.path(1) =  4.0;

    problem.phases(1).bounds.lower.path(2) = 0.0;
    problem.phases(1).bounds.upper.path(2) = 1900.0;

    problem.phases(1).bounds.lower.StartTime    = 0.0;
    problem.phases(1).bounds.upper.StartTime    = 0.0;

    problem.phases(1).bounds.lower.EndTime      = 10.0;
    problem.phases(1).bounds.upper.EndTime      = 100.0;


////////////////////////////////////////////////////////////////////////////
///////////////////  Register problem functions  ///////////////////////////
////////////////////////////////////////////////////////////////////////////


    problem.integrand_cost 	= &integrand_cost;
    problem.endpoint_cost 	= &endpoint_cost;
    problem.dae             	= &dae;
    problem.events 		= &events;
    problem.linkages		= &linkages;

////////////////////////////////////////////////////////////////////////////
///////////////////  Define & register initial guess ///////////////////////
////////////////////////////////////////////////////////////////////////////

    int nnodes    			= problem.phases(1).nodes(1);
    int ncontrols                       = problem.phases(1).ncontrols;
    int nstates                         = problem.phases(1).nstates;

    DMatrix state_guess    =  zeros(nstates,nnodes);
    DMatrix control_guess  =  zeros(ncontrols,nnodes);
    DMatrix time_guess     =  linspace(0.0,50.0,nnodes);



    state_guess(1,colon()) = linspace(0.0, -pi/2, nnodes);
    state_guess(2,colon()) = linspace(272.0, 310.0, nnodes);
    state_guess(3,colon()) = linspace(0.0, 10000.0, nnodes);
    state_guess(4,colon()) = linspace(30.0, 0.0, nnodes);

    control_guess(1,colon()) = 6000.0*ones(1,nnodes);

    problem.phases(1).guess.controls   = control_guess;
    problem.phases(1).guess.states     = state_guess;
    problem.phases(1).guess.time       = time_guess;

////////////////////////////////////////////////////////////////////////////
///////////////////  Enter algorithm options  //////////////////////////////
////////////////////////////////////////////////////////////////////////////


    algorithm.nlp_iter_max                = 1000;
    algorithm.nlp_tolerance               = 1.e-6;
    algorithm.nlp_method                  = "IPOPT";
    algorithm.scaling                     = "automatic";
    algorithm.derivatives                 = "automatic";
    algorithm.defect_scaling              = "jacobian-based";
    algorithm.collocation_method          = "trapezoidal";
//    algorithm.mesh_refinement             = "automatic";
    algorithm.hessian                     = "exact";
    algorithm.ode_tolerance               = 1.e-5;
    algorithm.mr_max_iterations           = 10;



////////////////////////////////////////////////////////////////////////////
///////////////////  Now call PSOPT to solve the problem   //////////////////
////////////////////////////////////////////////////////////////////////////

    psopt(solution, problem, algorithm);

////////////////////////////////////////////////////////////////////////////
///////////  Extract relevant variables from solution structure   //////////
////////////////////////////////////////////////////////////////////////////


    DMatrix x, u, p, t;

    x      = solution.get_states_in_phase(1);
    u      = solution.get_controls_in_phase(1);
    t      = solution.get_time_in_phase(1);
    p      = solution.get_dual_costates_in_phase(1);


////////////////////////////////////////////////////////////////////////////
///////////  Save solution data to files if desired ////////////////////////
////////////////////////////////////////////////////////////////////////////

    x.Save("x.dat");
    u.Save("u.dat");
    t.Save("t.dat");
    p.Save("p.dat");

////////////////////////////////////////////////////////////////////////////
///////////  Plot some results if desired (requires gnuplot) ///////////////
////////////////////////////////////////////////////////////////////////////

    plot(t,x(1,colon()),problem.name+": flight path angle","time (s)", "gamma", "gamma");

    plot(t,x(1,colon()),problem.name+": flight path angle","time (s)", "{/Symbol g}",
                                        "{/Symbol g}", "postscript eps enhanced color", "missile_gamma.eps");

    plot(t,x(1,colon()),problem.name+": flight path angle","time (s)", "{/Symbol g}",
                                        "{/Symbol g}", "pdf", "missile_gamma.pdf");

    plot(t,x(2,colon()),problem.name+": speed","time (s)", "V", "V");

    plot(t,x(2,colon()),problem.name+": speed (m/s)","time (s)", "V (m/s)",
                                        "V", "pdf", "missile_speed.pdf");

    plot(t,x(3,colon()),problem.name+": horizontal position","time (s)", "x", "x");

    plot(t,x(4,colon()),problem.name+": altitude","time (s)", "h", "h");

    plot(t,x(4,colon()),problem.name+": altitude (m)","time (s)", "h (m)",
                                        "h", "postscript eps enhanced color", "missile_alt.eps");

    plot(t,x(4,colon()),problem.name+": altitude (m)","time (s)", "h (m)",
                                        "h", "pdf", "missile_alt.pdf");

    plot(t,u(1,colon()),problem.name+": Thrust","time (s)", "T", "T");

    plot(t,u(2,colon()),problem.name+": Angle of attack","time (s)", "alpha", "alpha");

    plot(t,u(2,colon()),problem.name+": angle of attack (rad)","time (s)", "{/Symbol a} (rad)",
                                        "{/Symbol a}", "postscript eps enhanced color", "missile_alpha.eps");

    plot(t,u(2,colon()),problem.name+": angle of attack (rad)","time (s)", "alpha (rad)",
                                        "alpha", "pdf", "missile_alpha.pdf");

    plot(x(3,colon()),x(4,colon()),problem.name+": altitude (m)","x (m)", "h (m)",
                                        "h", "postscript eps enhanced color", "missile_alt_vs_x.eps");

    plot(x(3,colon()),x(4,colon()),problem.name+": altitude (m)","x (m)", "h (m)",
                                        "h", "pdf", "missile_alt_vs_x.pdf");


}

////////////////////////////////////////////////////////////////////////////
///////////////////////      END OF FILE     ///////////////////////////////
////////////////////////////////////////////////////////////////////////////
