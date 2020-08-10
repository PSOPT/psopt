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
    adouble x3 = states[2];
    adouble x6 = states[5];
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

   adouble x1 = states[ 0 ];
   adouble x2 = states[ 1 ];
   adouble x3 = states[ 2 ];
   adouble x4 = states[ 3 ];
   adouble x5 = states[ 4 ];
   adouble x6 = states[ 5 ];

   adouble u1 = controls[ 0 ];
   adouble u2 = controls[ 1 ];



   derivatives[ 0 ] = 9*x4;
   derivatives[ 1 ] = 9*x5;
   derivatives[ 2 ] = 9*x6;
   derivatives[ 3 ] = 9*(u1 + 17.2656*x3);
   derivatives[ 4 ] = 9*u2;
   derivatives[ 5 ] = -(9/x2)*(u1 + 27.0756*x3 + 2*x5*x6);



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
   adouble x50 = initial_states[ 4 ];
   adouble x60 = initial_states[ 5 ];
   adouble x1f = final_states[ 0 ];
   adouble x2f = final_states[ 1 ];
   adouble x3f = final_states[ 2 ];
   adouble x4f = final_states[ 3 ];
   adouble x5f = final_states[ 4 ];
   adouble x6f = final_states[ 5 ];


   e[ 0 ]  = x10;
   e[ 1 ]  = x20;
   e[ 2 ]  = x30;
   e[ 3 ]  = x40;
   e[ 4 ]  = x50;
   e[ 5 ]  = x60;
   e[ 6 ]  = x1f;
   e[ 7 ]  = x2f;
   e[ 8 ]  = x3f;
   e[ 9]  = x4f;
   e[ 10]  = x5f;
   e[ 11]  = x6f;

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

    problem.nphases   							= 1;
    problem.nlinkages                  	= 0;

    psopt_level1_setup(problem);


/////////////////////////////////////////////////////////////////////////////
/////////   Define phase related information & do level 2 setup /////////////
/////////////////////////////////////////////////////////////////////////////


    problem.phases(1).nstates   		= 6;
    problem.phases(1).ncontrols 		= 2;
    problem.phases(1).nevents   		= 12;
    problem.phases(1).npath     		= 0;
    problem.phases(1).nodes 	 		= (RowVectorXi(3) << 40, 60, 80).finished(); 

    psopt_level2_setup(problem, algorithm);


////////////////////////////////////////////////////////////////////////////
///////////////////  Enter problem bounds information //////////////////////
////////////////////////////////////////////////////////////////////////////

	 problem.phases(1).bounds.lower.states    << -5.0, -5.0, -5.0, -2.5, -1.0, -5.0;
	 problem.phases(1).bounds.upper.states    << 15.0, 25.0, 15.0,  2.5,  1.0, 15.0;

    problem.phases(1).bounds.lower.controls  << -2.83374, -0.80865;

    problem.phases(1).bounds.upper.controls  <<  2.83374,  0.71265;



    // Initial and final states

    problem.phases(1).bounds.lower.events    <<  0.0, 22.0, 0.0, 0.0, -1.0, 0.0, 10.0, 14.0, 0.0, 2.5, 0.0, 0.0; 

    problem.phases(1).bounds.upper.events    = problem.phases(1).bounds.lower.events;

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

    MatrixXd state_guess(6,20);

    state_guess   << 	linspace(0.0,10.0, 20),
    						 	linspace(22.0,14.0, 20),
    							linspace(0.,0., 20),
    							linspace(0.,2.5, 20),
    							linspace(-1.0,0., 20),
   							linspace(0.,0., 20);

    problem.phases(1).guess.controls       = zeros(2, 20);
    problem.phases(1).guess.states         = state_guess;
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

    MatrixXd x 		= solution.get_states_in_phase(1);
    MatrixXd u 		= solution.get_controls_in_phase(1);
    MatrixXd t 		= solution.get_time_in_phase(1);


////////////////////////////////////////////////////////////////////////////
///////////  Save solution data to files if desired ////////////////////////
////////////////////////////////////////////////////////////////////////////

    Save(x,"x.dat");
    Save(u,"u.dat");
    Save(t,"t.dat");

    MatrixXd x13 = x.block(0,0,3, length(t) ); 
    MatrixXd x46 = x.block(3,0,3, length(t) ); 


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

