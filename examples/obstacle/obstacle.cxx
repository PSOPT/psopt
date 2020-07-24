//////////////////////////////////////////////////////////////////////////
////////////////        obstacle.cxx                 /////////////////////
//////////////////////////////////////////////////////////////////////////
////////////////           PSOPT  Example             ////////////////////
//////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////
//////// Title:         Obstacle avoidance problem        ////////////////
//////// Last modified: 05 January 2009                   ////////////////
//////// Reference:  PROPT User's Guide                   ////////////////
//////// (See PSOPT handbook for full reference)          ////////////////
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
    double V = 2.138;
    adouble theta = controls[ CINDEX(1) ];

    adouble dxdt = V*cos(theta);
    adouble dydt = V*sin(theta);

//     adouble dxdt, dydt;

//     get_state_derivative( &dxdt, 1, iphase, time, xad);

//     get_state_derivative( &dydt, 2, iphase, time, xad);

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

   adouble x    = states[ CINDEX(1) ];
   adouble y    = states[ CINDEX(2) ];


   adouble theta = controls[ CINDEX(1) ];

   double V = 2.138;

   adouble dxdt = V*cos(theta);
   adouble dydt = V*sin(theta);


   derivatives[ CINDEX(1) ] = dxdt;
   derivatives[ CINDEX(2) ] = dydt;


   path[ CINDEX(1) ] = pow(x-0.4,2.0) + pow(y-0.5,2.0);
   path[ CINDEX(2) ] = pow(x-0.8,2.0) + pow(y-1.5,2.0);

}

////////////////////////////////////////////////////////////////////////////
///////////////////  Define the events function ////////////////////////////
////////////////////////////////////////////////////////////////////////////

void events(adouble* e, adouble* initial_states, adouble* final_states,
            adouble* parameters,adouble& t0, adouble& tf, adouble* xad,
            int iphase, Workspace* workspace)

{
   adouble x0 = initial_states[ CINDEX(1) ];
   adouble y0 = initial_states[ CINDEX(2) ];
   adouble xf = final_states[   CINDEX(1) ];
   adouble yf = final_states[   CINDEX(2) ];

   e[ CINDEX(1) ] = x0;
   e[ CINDEX(2) ] = y0;
   e[ CINDEX(3) ] = xf;
   e[ CINDEX(4) ] = yf;

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

    problem.name        		= "Obstacle avoidance problem";
    problem.outfilename                 = "obstacle.txt";

////////////////////////////////////////////////////////////////////////////
////////////  Define problem level constants & do level 1 setup ////////////
////////////////////////////////////////////////////////////////////////////

    problem.nphases   			= 1;
    problem.nlinkages                   = 0;

    psopt_level1_setup(problem);

/////////////////////////////////////////////////////////////////////////////
/////////   Define phase related information & do level 2 setup  ////////////
/////////////////////////////////////////////////////////////////////////////

    problem.phases(1).nstates   		= 2;
    problem.phases(1).ncontrols 		= 1;
    problem.phases(1).nevents   		= 4;
    problem.phases(1).npath                     = 2;
    problem.phases(1).nodes                     = "[20]";

    psopt_level2_setup(problem, algorithm);


////////////////////////////////////////////////////////////////////////////
///////////////////  Enter problem bounds information //////////////////////
////////////////////////////////////////////////////////////////////////////

    double xL = 0.0;
    double yL = 0.0;
    double xU = 2.0;
    double yU = 2.0;

    double thetaL = -10.0;
    double thetaU = 10.0;

    double x0 = 0.0;
    double y0 = 0.0;
    double xf = 1.2;
    double yf = 1.6;


    problem.phases(1).bounds.lower.states(1) = xL;
    problem.phases(1).bounds.lower.states(2) = yL;

    problem.phases(1).bounds.upper.states(1) = xU;
    problem.phases(1).bounds.upper.states(2) = yU;


    problem.phases(1).bounds.lower.controls(1) = thetaL;
    problem.phases(1).bounds.upper.controls(1) = thetaU;

    problem.phases(1).bounds.lower.events(1) = x0;
    problem.phases(1).bounds.lower.events(2) = y0;
    problem.phases(1).bounds.lower.events(3) = xf;
    problem.phases(1).bounds.lower.events(4) = yf;


    problem.phases(1).bounds.upper.events(1) = x0;
    problem.phases(1).bounds.upper.events(2) = y0;
    problem.phases(1).bounds.upper.events(3) = xf;
    problem.phases(1).bounds.upper.events(4) = yf;

    problem.phases(1).bounds.lower.path(1) = 0.1;
    problem.phases(1).bounds.upper.path(1) = 100.0;

    problem.phases(1).bounds.lower.path(2) = 0.1;
    problem.phases(1).bounds.upper.path(2) = 100.0;


    problem.phases(1).bounds.lower.StartTime    = 0.0;
    problem.phases(1).bounds.upper.StartTime    = 0.0;

    problem.phases(1).bounds.lower.EndTime      = 1.0;
    problem.phases(1).bounds.upper.EndTime      = 1.0;



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

    int nnodes    			= 30;
    int ncontrols                       = problem.phases(1).ncontrols;
    int nstates                         = problem.phases(1).nstates;

    DMatrix u_guess    =  zeros(ncontrols,nnodes);
    DMatrix x_guess    =  zeros(nstates,nnodes);
    DMatrix time_guess =  linspace(0.0,1.0,nnodes);


    x_guess(1,colon()) = linspace(x0,xf,nnodes);
    x_guess(2,colon()) = linspace(y0,yf,nnodes);

    u_guess(1,colon()) = zeros(1,nnodes);

    problem.phases(1).guess.controls       = u_guess;
    problem.phases(1).guess.states         = x_guess;
    problem.phases(1).guess.time           = time_guess;


////////////////////////////////////////////////////////////////////////////
///////////////////  Enter algorithm options  //////////////////////////////
////////////////////////////////////////////////////////////////////////////
    algorithm.nlp_iter_max                = 1000;
    algorithm.nlp_tolerance               = 1.e-4;
    algorithm.nlp_method                  = "IPOPT";
    algorithm.scaling                     = "automatic";
    algorithm.derivatives                 = "automatic";
    algorithm.collocation_method          = "trapezoidal";
    algorithm.mesh_refinement             = "automatic";
    algorithm.ode_tolerance               = 1.0e-2;

////////////////////////////////////////////////////////////////////////////
///////////////////  Now call PSOPT to solve the problem   /////////////////
////////////////////////////////////////////////////////////////////////////

    psopt(solution, problem, algorithm);

////////////////////////////////////////////////////////////////////////////
///////////  Extract relevant variables from solution structure   //////////
////////////////////////////////////////////////////////////////////////////


    DMatrix states = solution.get_states_in_phase(1);
    DMatrix theta  = solution.get_controls_in_phase(1);
    DMatrix t      = solution.get_time_in_phase(1);
    DMatrix mu     = solution.get_dual_path_in_phase(1);
    DMatrix lambda = solution.get_dual_costates_in_phase(1);

    DMatrix x = states(1,colon());
    DMatrix y = states(2,colon());

////////////////////////////////////////////////////////////////////////////
///////////  Save solution data to files if desired ////////////////////////
////////////////////////////////////////////////////////////////////////////

    x.Save("obstacle_x.dat");
    y.Save("obstacle_y.dat");
    theta.Save("obstacle_theta.dat");
    t.Save("obstacle_t.dat");


////////////////////////////////////////////////////////////////////////////
///////////  Plot some results if desired (requires gnuplot) ///////////////
////////////////////////////////////////////////////////////////////////////

    DMatrix alpha = colon(0.0, pi/20, 2*pi);

    DMatrix xObs1 = sqrt(0.1)*cos(alpha) + 0.4;
    DMatrix yObs1 = sqrt(0.1)*sin(alpha) + 0.5;

    DMatrix xObs2 = sqrt(0.1)*cos(alpha) + 0.8;
    DMatrix yObs2 = sqrt(0.1)*sin(alpha) + 1.5;

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
