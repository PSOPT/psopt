//////////////////////////////////////////////////////////////////////////
//////////////////           coulomb.cxx        //////////////////////////
//////////////////////////////////////////////////////////////////////////
////////////////           PSOPT  Example            /////////////////////
//////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////
//////// Title:                Coulumb friction problem   ////////////////
//////// Last modified:        04 January 2009            ////////////////
//////// Reference:   Driessen and Sadegh (2000)          ////////////////
//////// (See PSOPT handbook for full reference)          ////////////////
//////////////////////////////////////////////////////////////////////////
////////     Copyright (c) Victor M. Becerra, 2009         ///////////////
//////////////////////////////////////////////////////////////////////////
//////// This is part of the PSOPT software library, which ///////////////
//////// is distributed under the terms of the GNU Lesser  ///////////////
//////// General Public License (LGPL)                     ///////////////
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
   adouble q1    = states[ 0 ];
   adouble q1dot = states[ 1 ];
   adouble q2    = states[ 2 ];
   adouble q2dot = states[ 3 ];

   adouble u1 = controls[ 0 ];
   adouble u2 = controls[ 1 ];

   double k1 = 0.95; 
   double k2 = 0.85;
   double mu = 1.0;
   double m1 = 1.1;
   double m2 = 1.2;

   double epsilon = 0.01;

   derivatives[ 0 ] = q1dot;
   derivatives[ 1 ] = ( (-k1-k2)*q1+k2*q2-mu*smooth_sign(q1dot,epsilon)+u1 )/m1;
   derivatives[ 2 ] = q2dot;
   derivatives[ 3 ] = ( k2*q1-k2*q2-mu*smooth_sign(q2dot,epsilon)+u2 )/m2;
}

////////////////////////////////////////////////////////////////////////////
///////////////////  Define the events function ////////////////////////////
////////////////////////////////////////////////////////////////////////////

void events(adouble* e, adouble* initial_states, adouble* final_states, 
            adouble* parameters,adouble& t0, adouble& tf, adouble* xad, 
            int iphase, Workspace* workspace) 
{

   adouble q1_0    = initial_states[ 0 ];
   adouble q1dot_0 = initial_states[ 1 ];
   adouble q2_0    = initial_states[ 2 ];
   adouble q2dot_0 = initial_states[ 3 ];

   adouble q1_f    = final_states[ 0 ];
   adouble q1dot_f = final_states[ 1 ];
   adouble q2_f    = final_states[ 2 ];
   adouble q2dot_f = final_states[ 3 ];


   e[ 0 ] = q1_0;
   e[ 1 ] = q1dot_0;
   e[ 2 ] = q2_0;
   e[ 3 ] = q2dot_0;
   e[ 4 ] = q1_f;
   e[ 5 ] = q1dot_f;
   e[ 6 ] = q2_f;
   e[ 7 ] = q2dot_f;

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

    problem.name        		= "Coulomb friction problem";

    problem.outfilename                 = "coulomb.txt";

////////////////////////////////////////////////////////////////////////////
////////////  Define problem level constants & do level 1 setup ////////////
////////////////////////////////////////////////////////////////////////////

    problem.nphases   			= 1;
    problem.nlinkages                   = 0;

    psopt_level1_setup(problem);


/////////////////////////////////////////////////////////////////////////////
/////////   Define phase related information & do level 2 setup /////////////
/////////////////////////////////////////////////////////////////////////////

    problem.phases(1).nstates   		= 4;
    problem.phases(1).ncontrols 		= 2;
    problem.phases(1).nevents   		= 8;
    problem.phases(1).npath     		= 0;
    problem.phases(1).nodes         << 40; 

    psopt_level2_setup(problem, algorithm);

////////////////////////////////////////////////////////////////////////////
///////////////////  Enter problem bounds information //////////////////////
////////////////////////////////////////////////////////////////////////////

    double q1_0 =           0.0;
    double dotq1_0 =       -1.0;
    double q2_0 =           0.0;
    double dotq2_0 =       -2.0;
    double q1_f =           1.0;
    double dotq1_f =        0.0;
    double q2_f =           2.0;
    double dotq2_f =        0.0;


    problem.phases(1).bounds.lower.states(0) = -2.0;
    problem.phases(1).bounds.lower.states(1) = -20.0;
    problem.phases(1).bounds.lower.states(2) = -2.0;
    problem.phases(1).bounds.lower.states(3) = -20.0;

    problem.phases(1).bounds.upper.states(0) =  2.0;
    problem.phases(1).bounds.upper.states(1) =  20.0;
    problem.phases(1).bounds.upper.states(2) =  2.0;
    problem.phases(1).bounds.upper.states(3) =  20.0;

    problem.phases(1).bounds.lower.controls(0) = -4.0;
    problem.phases(1).bounds.lower.controls(1) = -4.0;
    problem.phases(1).bounds.upper.controls(0) =  4.0;
    problem.phases(1).bounds.upper.controls(1) =  4.0;

    problem.phases(1).bounds.lower.events(0) = q1_0;
    problem.phases(1).bounds.lower.events(1) = dotq1_0;
    problem.phases(1).bounds.lower.events(2) = q2_0;
    problem.phases(1).bounds.lower.events(3) = dotq2_0;
    problem.phases(1).bounds.lower.events(4) = q1_f;
    problem.phases(1).bounds.lower.events(5) = dotq1_f;
    problem.phases(1).bounds.lower.events(6) = q2_f;
    problem.phases(1).bounds.lower.events(7) = dotq2_f;

    problem.phases(1).bounds.upper.events(0) = q1_0;
    problem.phases(1).bounds.upper.events(1) = dotq1_0;
    problem.phases(1).bounds.upper.events(2) = q2_0;
    problem.phases(1).bounds.upper.events(3) = dotq2_0;
    problem.phases(1).bounds.upper.events(4) = q1_f;
    problem.phases(1).bounds.upper.events(5) = dotq1_f;
    problem.phases(1).bounds.upper.events(6) = q2_f;
    problem.phases(1).bounds.upper.events(7) = dotq2_f;


    problem.phases(1).bounds.lower.StartTime    = 0.0;
    problem.phases(1).bounds.upper.StartTime    = 0.0;

    problem.phases(1).bounds.lower.EndTime      = 1.8;
    problem.phases(1).bounds.upper.EndTime      = 4.0;


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


    MatrixXd state_guess(4,40);

    state_guess   << linspace(q1_0,q1_f, 40),
    						linspace(dotq1_0, dotq1_f, 40),
    						linspace(q2_0, q2_f, 40),
							linspace(dotq2_0, dotq2_f, 40);

    problem.phases(1).guess.controls       = zeros(2,40);
    problem.phases(1).guess.states         = state_guess;
    problem.phases(1).guess.time           = linspace(0.0, 4.0,40); 

////////////////////////////////////////////////////////////////////////////
///////////////////  Enter algorithm options  //////////////////////////////
////////////////////////////////////////////////////////////////////////////


    algorithm.nlp_method                  = "IPOPT";
    algorithm.scaling                     = "automatic";
    algorithm.derivatives                 = "automatic";
    algorithm.nlp_iter_max                = 1000;
    algorithm.nlp_tolerance               = 1.e-6;

////////////////////////////////////////////////////////////////////////////
///////////////////  Now call PSOPT to solve the problem   //////////////////
////////////////////////////////////////////////////////////////////////////

    psopt(solution, problem, algorithm);

    if (solution.error_flag) exit(0);

////////////////////////////////////////////////////////////////////////////
///////////  Extract relevant variables from solution structure   //////////
////////////////////////////////////////////////////////////////////////////
    MatrixXd x, u, t;
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

    MatrixXd q12(2, length(t));
    
    q12 = x.row(0), 
          x.row(2); 

    plot(t,q12,problem.name + ": states q1 and q2", 
                        "time (s)", "states", "q1 q2");

    plot(t,u,problem.name + ": controls", "time (s)", "control", "u1 u2");

    plot(t,q12,problem.name + ": states q1 and q2", 
                        "time (s)", "states", "q1 q2", 
                        "pdf", "coulomb_states.pdf");

    plot(t,u,problem.name + ": controls", "time (s)", "controls", "u1 u2", 
                              "pdf", "coulomb_control.pdf");

}

////////////////////////////////////////////////////////////////////////////
///////////////////////      END OF FILE     ///////////////////////////////
////////////////////////////////////////////////////////////////////////////


