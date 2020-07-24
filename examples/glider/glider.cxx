//////////////////////////////////////////////////////////////////////////
////////////////             glider.cxx              /////////////////////
//////////////////////////////////////////////////////////////////////////
////////////////           PSOPT example            /////////////////////
//////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////
//////// Title: 	Hang glider problem               ////////////////
//////// Last modified: 11 July 2009                      ////////////////
//////// Reference:     PROPT User Manual                 ////////////////
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
   adouble xf = final_states[CINDEX(1)];

   return -(xf);
}

//////////////////////////////////////////////////////////////////////////
///////////////////  Define the integrand (Lagrange) cost function  //////
//////////////////////////////////////////////////////////////////////////

adouble integrand_cost(adouble* states, adouble* controls, adouble* parameters,
                     adouble& time, adouble* xad, int iphase, Workspace* workspace)
{
   return 0.0;
}


//////////////////////////////////////////////////////////////////////////
///////////////////  Define the DAE's ////////////////////////////////////
//////////////////////////////////////////////////////////////////////////


void dae(adouble* derivatives, adouble* path, adouble* states,
         adouble* controls, adouble* parameters, adouble& time,
         adouble* xad, int iphase, Workspace* workspace)
{

    adouble CL 		= controls[ CINDEX(1) ];

    adouble x 		= states[ CINDEX(1) ];
    adouble y 	    = states[ CINDEX(2) ];
    adouble vx 		= states[ CINDEX(3) ];
    adouble vy 	    = states[ CINDEX(4) ];

    double m  = 100.0,      g  = 9.80665;
    double uM = 2.5,        R  = 100.0;
    double C0 = 0.034,      k = 0.069662;
    double S  = 14.0,      rho = 1.13;

    adouble sin_eta, cos_eta, D, L, CD, Vy, ua, X, vr, W;


    CD = C0 + k*CL*CL;
    vr = sqrt(vx*vx + vy*vy);
    D = 0.5*CD*rho*S*vr*vr;
    L = 0.5*CL*rho*S*vr*vr;
    X  = pow(x/R - 2.5, 2.0);
    ua = uM*(1.0-X)*exp(-X);
    Vy = vy-ua;
    sin_eta = Vy/vr;
    cos_eta = vx/vr;
    W = m*g;


    derivatives[ CINDEX(1) ] =   vx;
    derivatives[ CINDEX(2) ] =   vy;
    derivatives[ CINDEX(3) ] =   1.0/m*(-L*sin_eta - D*cos_eta    );
    derivatives[ CINDEX(4) ] =   1.0/m*( L*cos_eta - D*sin_eta - W);



}

////////////////////////////////////////////////////////////////////////////
///////////////////  Define the events function ////////////////////////////
////////////////////////////////////////////////////////////////////////////

void events(adouble* e, adouble* initial_states, adouble* final_states,
            adouble* parameters,adouble& t0, adouble& tf, adouble* xad,
            int iphase, Workspace* workspace)

{

    adouble x_i 	= initial_states[ CINDEX(1) ];
    adouble y_i 	= initial_states[ CINDEX(2) ];
    adouble vx_i 	= initial_states[ CINDEX(3) ];
    adouble vy_i 	= initial_states[ CINDEX(4) ];

    adouble x_f 	= final_states[ CINDEX(1) ];
    adouble y_f		= final_states[ CINDEX(2) ];
    adouble vx_f	= final_states[ CINDEX(3) ];
    adouble vy_f 	= final_states[ CINDEX(4) ];

    e[ CINDEX(1) ] 	    = x_i;
    e[ CINDEX(2) ]      = y_i;
    e[ CINDEX(3) ]      = vx_i;
    e[ CINDEX(4) ]    	= vy_i;

    e[ CINDEX(5) ]      = y_f;
    e[ CINDEX(6) ]   	= vx_f;
    e[ CINDEX(7) ]  	= vy_f;


}



///////////////////////////////////////////////////////////////////////////
///////////////////  Define the phase linkages function ///////////////////
///////////////////////////////////////////////////////////////////////////

void linkages( adouble* linkages, adouble* xad, Workspace* workspace)
{

   // Single phase problem

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

    problem.name        		= "Hang glider problem";
    problem.outfilename                 = "hang_glider.txt";

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
    problem.phases(1).ncontrols 		= 1;
    problem.phases(1).nevents   		= 7;
    problem.phases(1).npath     		= 0;
    problem.phases(1).nodes                     = "[30 40 50 80]";


    psopt_level2_setup(problem, algorithm);



////////////////////////////////////////////////////////////////////////////
///////////////////  Enter problem bounds information //////////////////////
////////////////////////////////////////////////////////////////////////////


    problem.phases(1).bounds.lower.states = "[0.0      0.0    0.0  -4.0]";
    problem.phases(1).bounds.upper.states = "[1500.0  1100.0 15.0   4.0]";


    problem.phases(1).bounds.lower.controls(1) = 0.0;

    problem.phases(1).bounds.upper.controls(1) = 1.4;


    problem.phases(1).bounds.lower.events="[0.0,1000.0,13.2275675,-1.28750052,900.00,13.2275675,-1.28750052 ]";
    problem.phases(1).bounds.upper.events="[0.0,1000.0,13.2275675,-1.28750052,900.00,13.2275675,-1.28750052 ]";


    problem.phases(1).bounds.lower.StartTime    = 0.0;
    problem.phases(1).bounds.upper.StartTime    = 0.0;

    problem.phases(1).bounds.lower.EndTime      = 0.1;
    problem.phases(1).bounds.upper.EndTime      = 200.0;



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
    int ncontrols           = problem.phases(1).ncontrols;
    int nstates             = problem.phases(1).nstates;

    DMatrix state_guess    =  zeros(nstates,nnodes);
    DMatrix control_guess  =  1.0*ones(ncontrols,nnodes);
    DMatrix time_guess     =  linspace(0.0,105.0,nnodes);


    state_guess(1,colon()) = linspace(0.0, 1250, nnodes);
    state_guess(2,colon()) = linspace(1000.0, 900.0, nnodes);
    state_guess(3,colon()) = 13.23*ones(1,nnodes);
    state_guess(4,colon()) = -1.288*ones(1,nnodes);

    problem.phases(1).guess.states   = state_guess;
    problem.phases(1).guess.controls = control_guess;
    problem.phases(1).guess.time     = time_guess;



////////////////////////////////////////////////////////////////////////////
///////////////////  Enter algorithm options  //////////////////////////////
////////////////////////////////////////////////////////////////////////////


    algorithm.nlp_iter_max                = 1000;
    algorithm.nlp_tolerance               = 1.e-6;
    algorithm.nlp_method                  = "IPOPT";
    algorithm.scaling                     = "automatic";
    algorithm.derivatives                 = "automatic";


////////////////////////////////////////////////////////////////////////////
///////////////////  Now call PSOPT to solve the problem   //////////////////
////////////////////////////////////////////////////////////////////////////

    psopt(solution, problem, algorithm);

////////////////////////////////////////////////////////////////////////////
///////////  Extract relevant variables from solution structure   //////////
////////////////////////////////////////////////////////////////////////////


    DMatrix states, CL, t, x, y, speeds;

    states      = solution.get_states_in_phase(1);
    CL          = solution.get_controls_in_phase(1);
    t           = solution.get_time_in_phase(1);




////////////////////////////////////////////////////////////////////////////
///////////  Save solution data to files if desired ////////////////////////
////////////////////////////////////////////////////////////////////////////

    states.Save("states.dat");
    CL.Save("cL.dat");
    t.Save("t.dat");


////////////////////////////////////////////////////////////////////////////
///////////  Plot some results if desired (requires gnuplot) ///////////////
////////////////////////////////////////////////////////////////////////////


   x 	= states(1,colon());
   y 	= states(2,colon());
   speeds = (states(3,colon()) && states(4,colon()));

   plot(x,y,problem.name+": trajectory","x [m]", "y [m]", "traj");

   plot(t,speeds,problem.name+": speeds","time (s)", "speeds [m/s]", "dxdt dydt");

   plot(t,CL,problem.name+": control","time (s)", "control", "CL");

   plot(x,y,problem.name+": trajectory","x [m]", "y [m]", "traj", "pdf", "traj.pdf");

   plot(t,speeds,problem.name+": speeds","time (s)", "speeds [m/s]", "dxdt dydt", "pdf", "velocities.pdf");

   plot(t,CL,problem.name+": control - lift coefficient","time (s)", "CL", "CL", "pdf", "control.pdf");


}

////////////////////////////////////////////////////////////////////////////
///////////////////////      END OF FILE     ///////////////////////////////
////////////////////////////////////////////////////////////////////////////
