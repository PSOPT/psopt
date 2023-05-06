//////////////////////////////////////////////////////////////////////////
////////////////             wheat.cxx               /////////////////////
//////////////////////////////////////////////////////////////////////////
////////////////           PSOPT  Example            /////////////////////
//////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////
//////// Title:         Continuous Wheat Trading Model    ////////////////
//////// Last modified: 05 May 2023                       ////////////////
//////// Reference:     Ijiri and Thompson (1970)         ////////////////
//////// (See PSOPT handbook for full reference)          ////////////////
//////////////////////////////////////////////////////////////////////////
////////     Copyright (c) Victor M. Becerra, 2009        ////////////////
////////     Copyright (c) Eric Brown, 2023               ////////////////
//////////////////////////////////////////////////////////////////////////
//////// This is part of the PSOPT software library, which////////////////
//////// is distributed under the terms of the GNU Lesser ////////////////
//////// General Public License (LGPL)                    ////////////////
//////////////////////////////////////////////////////////////////////////

#include "psopt.h"

//////////////////////////////////////////////////////////////////////////
///////////////////  Define the end point (Mayer) cost function //////////
//////////////////////////////////////////////////////////////////////////

adouble pfunc(adouble &t)
{
    return 1.0 + 0.5 * sin(2.0 * 2.0*PSOPT::pi/12.0  * t + 3.0*PSOPT::pi/2.0);
}

adouble endpoint_cost(adouble* initial_states, adouble* final_states,
                      adouble* parameters,adouble& t0, adouble& tf,
                      adouble* xad, int iphase, Workspace* workspace)
{
   adouble c = final_states[0];
   adouble w = final_states[1];

   adouble P = pfunc(tf);

   return -(c + P*w);
}

//////////////////////////////////////////////////////////////////////////
///////////////////  Define the integrand (Lagrange) cost function  //////
//////////////////////////////////////////////////////////////////////////

adouble integrand_cost(adouble* states, adouble* controls, adouble* parameters,
                       adouble& time, adouble* xad, int iphase, Workspace* workspace)
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
   adouble cdot, wdot;

   adouble c = states[ 0 ];
   adouble w = states[ 1 ];

   adouble ub = controls[ 0 ];
   adouble us = controls[ 1 ];

   double R = 0.01;
   double S = 0.20;
   double fb = 0.01;
   double fs = 0.05;
   double A = 0.05;

   adouble p = pfunc(time);

   cdot = (R*c) - (S*w) - (p*ub) - (p*us) - (fb*ub) - (-fs*us) ;
   wdot =       - (A*w) +    ub  +    us                       ;

   derivatives[ 0 ] = cdot;
   derivatives[ 1 ] = wdot;
}

////////////////////////////////////////////////////////////////////////////
///////////////////  Define the events function ////////////////////////////
////////////////////////////////////////////////////////////////////////////

void events(adouble* e, adouble* initial_states, adouble* final_states,
            adouble* parameters,adouble& t0, adouble& tf, adouble* xad,
            int iphase, Workspace* workspace)

{
   adouble c0 = initial_states[ 0 ];
   adouble w0 = initial_states[ 1 ];

   e[ 0 ] = c0;
   e[ 1 ] = w0;
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

    problem.name        		        = "Wheat Problem";
    problem.outfilename                 = "wheat.txt";

////////////////////////////////////////////////////////////////////////////
////////////  Define problem level constants & do level 1 setup ////////////
////////////////////////////////////////////////////////////////////////////

    problem.nphases   		        	= 1;
    problem.nlinkages                   = 0;

    psopt_level1_setup(problem);

/////////////////////////////////////////////////////////////////////////////
/////////   Define phase related information & do level 2 setup  ////////////
/////////////////////////////////////////////////////////////////////////////

    problem.phases(1).nstates   		= 2;
    problem.phases(1).ncontrols 		= 2;
    problem.phases(1).nevents   		= 2;
    problem.phases(1).npath     		= 0;
    problem.phases(1).nodes       <<  32;

    psopt_level2_setup(problem, algorithm);

////////////////////////////////////////////////////////////////////////////
///////////////////  Declare MatrixXd objects to store results //////////////
////////////////////////////////////////////////////////////////////////////

    MatrixXd x, u, t;
    MatrixXd lambda, H;

////////////////////////////////////////////////////////////////////////////
///////////////////  Enter problem bounds information //////////////////////
////////////////////////////////////////////////////////////////////////////

    double cL = 0.0;
    double wL = 0.0;
    double cU = PSOPT::inf;
    double wU = PSOPT::inf;

    double uL = -1.0;
    double uH =  1.0;

    double c0 = 1.0;
    double w0 = 0.0;

    problem.phases(1).bounds.lower.states(0) = cL;
    problem.phases(1).bounds.lower.states(1) = wL;

    problem.phases(1).bounds.upper.states(0) = cU;
    problem.phases(1).bounds.upper.states(1) = wU;

    problem.phases(1).bounds.lower.controls(0) = 0.0;
    problem.phases(1).bounds.lower.controls(1) =  uL;
    problem.phases(1).bounds.upper.controls(0) =  uH;
    problem.phases(1).bounds.upper.controls(1) = 0.0;

    problem.phases(1).bounds.lower.events(0) = c0;
    problem.phases(1).bounds.lower.events(1) = w0;

    problem.phases(1).bounds.upper.events(0) = c0;
    problem.phases(1).bounds.upper.events(1) = w0;

    problem.phases(1).bounds.lower.StartTime    = 0.0;
    problem.phases(1).bounds.upper.StartTime    = 0.0;

    problem.phases(1).bounds.lower.EndTime      = 12.0;
    problem.phases(1).bounds.upper.EndTime      = 12.0;



////////////////////////////////////////////////////////////////////////////
///////////////////  Register problem functions  ///////////////////////////
////////////////////////////////////////////////////////////////////////////


    problem.integrand_cost 	= &integrand_cost;
    problem.endpoint_cost 	= &endpoint_cost;
    problem.dae             = &dae;
    problem.events 		    = &events;
    problem.linkages		= &linkages;

////////////////////////////////////////////////////////////////////////////
///////////////////  Define & register initial guess ///////////////////////
////////////////////////////////////////////////////////////////////////////

    int nnodes    			            = problem.phases(1).nodes(0);
    int ncontrols                       = problem.phases(1).ncontrols;
    int nstates                         = problem.phases(1).nstates;

    MatrixXd x_guess    =  zeros(nstates,nnodes);

    x_guess.row(0)  = c0*ones(1,nnodes);
    x_guess.row(1)  = w0*ones(1,nnodes);

    problem.phases(1).guess.controls       = zeros(ncontrols,nnodes);
    problem.phases(1).guess.states         = x_guess;
    problem.phases(1).guess.time           = linspace(0.0,12.0,nnodes);


////////////////////////////////////////////////////////////////////////////
///////////////////  Enter algorithm options  //////////////////////////////
////////////////////////////////////////////////////////////////////////////


    algorithm.nlp_iter_max                = 1000;
    algorithm.nlp_tolerance               = 1.e-4;
    algorithm.nlp_method                  = "IPOPT";
    algorithm.hessian                     = "exact";
    algorithm.scaling                     = "automatic";
    algorithm.derivatives                 = "automatic";
    algorithm.collocation_method          = "trapezoidal";
    algorithm.diff_matrix                 = "central-differences";
    algorithm.mesh_refinement             = "automatic";
    algorithm.switch_order                = 0;

////////////////////////////////////////////////////////////////////////////
///////////////////  Now call PSOPT to solve the problem   /////////////////
////////////////////////////////////////////////////////////////////////////

    psopt(solution, problem, algorithm);

////////////////////////////////////////////////////////////////////////////
///////////  Extract relevant variables from solution structure   //////////
////////////////////////////////////////////////////////////////////////////


    x      = solution.get_states_in_phase(1);
    u      = solution.get_controls_in_phase(1);
    t      = solution.get_time_in_phase(1);
    lambda = solution.get_dual_costates_in_phase(1);
    H      = solution.get_dual_hamiltonian_in_phase(1);


////////////////////////////////////////////////////////////////////////////
///////////  Save solution data to files if desired ////////////////////////
////////////////////////////////////////////////////////////////////////////

    Save(x, "x.dat");
    Save(u,"u.dat");
    Save(t,"t.dat");
    Save(lambda,"lambda.dat");
    Save(H,"H.dat");

////////////////////////////////////////////////////////////////////////////
///////////  Plot some results if desired (requires gnuplot) ///////////////
////////////////////////////////////////////////////////////////////////////

    plot(t,x,problem.name+": states", "time (s)", "states","c w");

    plot(t,u,problem.name+": controls","time (s)", "controls", "ub us");

    plot(t,lambda,problem.name+": lambda","time (s)", "lambda", "lambda_c lambda_w");

    plot(t,H,problem.name+": Hamiltonian", "time (s)", "Hamiltonian","H");

    plot(t,x,problem.name+": states", "time (s)", "states","c w",
                             "pdf", "commodity_states.pdf");

    plot(t,u,problem.name+": controls","time (s)", "controls", "ub us",
                             "pdf", "commodity_controls.pdf");

    plot(t,lambda,problem.name+": lambda", "time (s)", "lambda","lambda_c lambda_w",
         "pdf", "commodity_lambda.pdf");

    plot(t,H,problem.name+": Hamiltonian","time (s)", "Hamiltonian", "H",
         "pdf", "commodity_H.pdf");
}

////////////////////////////////////////////////////////////////////////////
///////////////////////      END OF FILE     ///////////////////////////////
////////////////////////////////////////////////////////////////////////////
