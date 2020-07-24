//////////////////////////////////////////////////////////////////////////
//////////////////        crane.cxx             //////////////////////////
//////////////////////////////////////////////////////////////////////////
////////////////           PSOPT  Example             ////////////////////
//////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////
//////// Title: Minimum Swing Control for Container Crane ////////////////
//////// Last modified:         06 February 2009          ////////////////
//////// Reference:             Teo and Goh               ////////////////
//////// (See PSOPT handbook for full reference)           ///////////////
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
    return 0;
}

//////////////////////////////////////////////////////////////////////////
///////////////////  Define the integrand (Lagrange) cost function  //////
//////////////////////////////////////////////////////////////////////////

adouble integrand_cost(adouble* states, adouble* controls,
                       adouble* parameters, adouble& time, adouble* xad,
                       int iphase, Workspace* workspace)
{
    adouble x3 = states[CINDEX(3)];
    adouble x6 = states[CINDEX(6)];
    return  4.5*( pow(x3,2) + pow(x6,2) );
}

//////////////////////////////////////////////////////////////////////////
///////////////////  Define the DAE's ////////////////////////////////////
//////////////////////////////////////////////////////////////////////////

void dae(adouble* derivatives, adouble* path, adouble* states,
         adouble* controls, adouble* parameters, adouble& time,
         adouble* xad, int iphase, Workspace* workspace)
{
   adouble xdot, ydot, vdot;

   adouble x1 = states[ CINDEX(1) ];
   adouble x2 = states[ CINDEX(2) ];
   adouble x3 = states[ CINDEX(3) ];
   adouble x4 = states[ CINDEX(4) ];
   adouble x5 = states[ CINDEX(5) ];
   adouble x6 = states[ CINDEX(6) ];

   adouble u1 = controls[ CINDEX(1) ];
   adouble u2 = controls[ CINDEX(2) ];



   derivatives[ CINDEX(1) ] = 9*x4;
   derivatives[ CINDEX(2) ] = 9*x5;
   derivatives[ CINDEX(3) ] = 9*x6;
   derivatives[ CINDEX(4) ] = 9*(u1 + 17.2656*x3);
   derivatives[ CINDEX(5) ] = 9*u2;
   derivatives[ CINDEX(6) ] = -(9/x2)*(u1 + 27.0756*x3 + 2*x5*x6);



}

////////////////////////////////////////////////////////////////////////////
///////////////////  Define the events function ////////////////////////////
////////////////////////////////////////////////////////////////////////////

void events(adouble* e, adouble* initial_states, adouble* final_states,
            adouble* parameters,adouble& t0, adouble& tf, adouble* xad,
            int iphase, Workspace* workspace)
{
   adouble x10 = initial_states[ CINDEX(1) ];
   adouble x20 = initial_states[ CINDEX(2) ];
   adouble x30 = initial_states[ CINDEX(3) ];
   adouble x40 = initial_states[ CINDEX(4) ];
   adouble x50 = initial_states[ CINDEX(5) ];
   adouble x60 = initial_states[ CINDEX(6) ];
   adouble x1f = final_states[ CINDEX(1) ];
   adouble x2f = final_states[ CINDEX(2) ];
   adouble x3f = final_states[ CINDEX(3) ];
   adouble x4f = final_states[ CINDEX(4) ];
   adouble x5f = final_states[ CINDEX(5) ];
   adouble x6f = final_states[ CINDEX(6) ];


   e[ CINDEX(1) ]  = x10;
   e[ CINDEX(2) ]  = x20;
   e[ CINDEX(3) ]  = x30;
   e[ CINDEX(4) ]  = x40;
   e[ CINDEX(5) ]  = x50;
   e[ CINDEX(6) ]  = x60;
   e[ CINDEX(7) ]  = x1f;
   e[ CINDEX(8) ]  = x2f;
   e[ CINDEX(9) ]  = x3f;
   e[ CINDEX(10)]  = x4f;
   e[ CINDEX(11)]  = x5f;
   e[ CINDEX(12)]  = x6f;

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

    problem.name = "Minimum swing control for a container crane";

    problem.outfilename                 = "crane.txt";

////////////////////////////////////////////////////////////////////////////
////////////  Define problem level constants & do level 1 setup ////////////
////////////////////////////////////////////////////////////////////////////

    problem.nphases   			= 1;
    problem.nlinkages                   = 0;

    psopt_level1_setup(problem);


/////////////////////////////////////////////////////////////////////////////
/////////   Define phase related information & do level 2 setup /////////////
/////////////////////////////////////////////////////////////////////////////


    problem.phases(1).nstates   		= 6;
    problem.phases(1).ncontrols 		= 2;
    problem.phases(1).nevents   		= 12;
    problem.phases(1).npath     		= 0;
    problem.phases(1).nodes 	 		= "[40, 60, 80]";

    psopt_level2_setup(problem, algorithm);


////////////////////////////////////////////////////////////////////////////
///////////////////  Enter problem bounds information //////////////////////
////////////////////////////////////////////////////////////////////////////

    problem.phases(1).bounds.lower.states(1) = -5.0;
    problem.phases(1).bounds.lower.states(2) = -5.0;
    problem.phases(1).bounds.lower.states(3) = -5.0;
    problem.phases(1).bounds.lower.states(4) = -2.5;
    problem.phases(1).bounds.lower.states(5) = -1.0;
    problem.phases(1).bounds.lower.states(6) = -5.0;

    problem.phases(1).bounds.upper.states(1) = 15.0;
    problem.phases(1).bounds.upper.states(2) = 25.0;
    problem.phases(1).bounds.upper.states(3) = 15.0;
    problem.phases(1).bounds.upper.states(4) = 2.5;
    problem.phases(1).bounds.upper.states(5) = 1.0;
    problem.phases(1).bounds.upper.states(6) = 15.0;

    problem.phases(1).bounds.lower.controls(1) = -2.83374;
    problem.phases(1).bounds.lower.controls(2) = -0.80865;
    problem.phases(1).bounds.upper.controls(1) = 2.83374;
    problem.phases(1).bounds.upper.controls(2) = 0.71265;

    // Initial states
    problem.phases(1).bounds.lower.events(1) = 0.0;
    problem.phases(1).bounds.lower.events(2) = 22.0;
    problem.phases(1).bounds.lower.events(3) = 0.0;
    problem.phases(1).bounds.lower.events(4) = 0.0;
    problem.phases(1).bounds.lower.events(5) = -1.0;
    problem.phases(1).bounds.lower.events(6) = 0.0;

    // Final states
    problem.phases(1).bounds.lower.events(7) = 10.0;
    problem.phases(1).bounds.lower.events(8) = 14.0;
    problem.phases(1).bounds.lower.events(9) = 0.0;
    problem.phases(1).bounds.lower.events(10)= 2.5;
    problem.phases(1).bounds.lower.events(11)= 0.0;
    problem.phases(1).bounds.lower.events(12)= 0.0;

    problem.phases(1).bounds.upper.events = problem.phases(1).bounds.lower.events;

    problem.phases(1).bounds.lower.StartTime    = 0.0;
    problem.phases(1).bounds.upper.StartTime    = 0.0;

    problem.phases(1).bounds.lower.EndTime      = 1.0;
    problem.phases(1).bounds.upper.EndTime      = 1.0;


////////////////////////////////////////////////////////////////////////////
///////////////////  Register problem functions  ///////////////////////////
////////////////////////////////////////////////////////////////////////////


    problem.integrand_cost 	= &integrand_cost;
    problem.endpoint_cost 	= &endpoint_cost;
    problem.dae 		= &dae;
    problem.events 		= &events;
    problem.linkages		= &linkages;



////////////////////////////////////////////////////////////////////////////
///////////////////  Define & register initial guess ///////////////////////
////////////////////////////////////////////////////////////////////////////

    DMatrix x0(6,20);

    x0(1,colon()) = linspace(0.0,10.0, 20);
    x0(2,colon()) = linspace(22.0,14.0, 20);
    x0(3,colon()) = linspace(0.,0., 20);
    x0(4,colon()) = linspace(0.,2.5, 20);
    x0(5,colon()) = linspace(-1.0,0., 20);
    x0(6,colon()) = linspace(0.,0., 20);

    problem.phases(1).guess.controls       = zeros(2, 20);
    problem.phases(1).guess.states         = x0;
    problem.phases(1).guess.time           = linspace(0.0, 1.0, 20); ;


////////////////////////////////////////////////////////////////////////////
///////////////////  Enter algorithm options  //////////////////////////////
////////////////////////////////////////////////////////////////////////////


    algorithm.nlp_method                  = "IPOPT";
    algorithm.scaling                     = "automatic";
    algorithm.derivatives                 = "automatic";
    algorithm.nlp_iter_max                = 1000;
    algorithm.nlp_tolerance               = 1.e-6;
    algorithm.collocation_method          = "Legendre";
 

////////////////////////////////////////////////////////////////////////////
///////////////////  Now call PSOPT to solve the problem   /////////////////
////////////////////////////////////////////////////////////////////////////

    psopt(solution, problem, algorithm);

////////////////////////////////////////////////////////////////////////////
///////////  Extract relevant variables from solution structure   //////////
////////////////////////////////////////////////////////////////////////////

    DMatrix x 		= solution.get_states_in_phase(1);
    DMatrix u 		= solution.get_controls_in_phase(1);
    DMatrix t 		= solution.get_time_in_phase(1);


////////////////////////////////////////////////////////////////////////////
///////////  Save solution data to files if desired ////////////////////////
////////////////////////////////////////////////////////////////////////////

    x.Save("x.dat");
    u.Save("u.dat");
    t.Save("t.dat");

    DMatrix x13 = x(colon(1,3), colon() );
    DMatrix x46 = x(colon(4,6), colon() );


////////////////////////////////////////////////////////////////////////////
///////////  Plot some results if desired (requires gnuplot) ///////////////
////////////////////////////////////////////////////////////////////////////

    plot(t,x13,problem.name + ": states x1, x2 and x3", "time (s)", "states", "x1 x2 x3");

    plot(t,x46,problem.name + ": states x4, x5 and x6", "time (s)", "states", "x4 x5 x6");

    plot(t,u,problem.name + ": controls", "time", "controls", "u1 u2");


    plot(t,x13,problem.name + ": states x1, x2 and x3", "time (s)", "states", "x1 x2 x3",
                                  "pdf", "crane_states13.pdf");

    plot(t,x46,problem.name + ": states x4, x5 and x6", "time (s)", "states", "x4 x5 x6",
                                  "pdf", "crane_states46.pdf");

    plot(t,u,problem.name + ": controls", "time", "controls", "u1 u2",
                              "pdf", "crane_controls.pdf");


}

////////////////////////////////////////////////////////////////////////////
///////////////////////      END OF FILE     ///////////////////////////////
////////////////////////////////////////////////////////////////////////////

