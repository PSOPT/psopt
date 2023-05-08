//////////////////////////////////////////////////////////////////////////
////////////////               dqc.cxx               /////////////////////
//////////////////////////////////////////////////////////////////////////
////////////////           PSOPT  Example            /////////////////////
//////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////
//////// Title:         Double-Integrator Quadratic Control (DQC)  ////
//////// Last modified: 05 May 2023                       ////////////////
//////// Reference:     Ross (2015)                       ////////////////
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

adouble endpoint_cost(adouble* initial_states, adouble* final_states,
                      adouble* parameters,adouble& t0, adouble& tf,
                      adouble* xad, int iphase, Workspace* workspace)
{
   return 0.0;
}

//////////////////////////////////////////////////////////////////////////
///////////////////  Define the integrand (Lagrange) cost function  //////
//////////////////////////////////////////////////////////////////////////

adouble integrand_cost(adouble* states, adouble* controls, adouble* parameters,
                     adouble& time, adouble* xad, int iphase, Workspace* workspace)
{
    adouble u = controls[0];

    adouble L = 0.5*(u*u);

    return L;
}


//////////////////////////////////////////////////////////////////////////
///////////////////  Define the DAE's ////////////////////////////////////
//////////////////////////////////////////////////////////////////////////

void dae(adouble* derivatives, adouble* path, adouble* states,
         adouble* controls, adouble* parameters, adouble& time,
         adouble* xad, int iphase, Workspace* workspace)
{
   adouble xdot, vdot;

   adouble x = states[ 0 ];
   adouble v = states[ 1 ];

   adouble u = controls[ 0 ];

   xdot = v;
   vdot = u;

   derivatives[ 0 ] = xdot;
   derivatives[ 1 ] = vdot;
}

////////////////////////////////////////////////////////////////////////////
///////////////////  Define the events function ////////////////////////////
////////////////////////////////////////////////////////////////////////////

void events(adouble* e, adouble* initial_states, adouble* final_states,
            adouble* parameters,adouble& t0, adouble& tf, adouble* xad,
            int iphase, Workspace* workspace)

{
   adouble x0 = initial_states[ 0 ];
   adouble v0 = initial_states[ 1 ];
   adouble xf = final_states[ 0 ];
   adouble vf = final_states[ 1 ];

   e[ 0 ] = x0;
   e[ 1 ] = v0;
   e[ 2 ] = xf;
   e[ 3 ] = vf;
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

    problem.name        		= "Double-Integrator Quadratic Control";
    problem.outfilename         = "dqc.txt";

////////////////////////////////////////////////////////////////////////////
////////////  Define problem level constants & do level 1 setup ////////////
////////////////////////////////////////////////////////////////////////////

    problem.nphases   			= 1;
    problem.nlinkages           = 0;

    psopt_level1_setup(problem);

/////////////////////////////////////////////////////////////////////////////
/////////   Define phase related information & do level 2 setup  ////////////
/////////////////////////////////////////////////////////////////////////////

    problem.phases(1).nstates   		= 2;
    problem.phases(1).ncontrols 		= 1;
    problem.phases(1).nevents   		= 4;
    problem.phases(1).npath     		= 0;
    problem.phases(1).nodes         << 50;

    psopt_level2_setup(problem, algorithm);

////////////////////////////////////////////////////////////////////////////
///////////////////  Declare MatrixXd objects to store results //////////////
////////////////////////////////////////////////////////////////////////////

    MatrixXd x, u, t;
    MatrixXd lambda, H;

////////////////////////////////////////////////////////////////////////////
///////////////////  Enter problem bounds information //////////////////////
////////////////////////////////////////////////////////////////////////////

    double xL = -PSOPT::inf;
    double vL = -PSOPT::inf;
    double xU =  PSOPT::inf;
    double vU =  PSOPT::inf;

    double uL = -PSOPT::inf;
    double uH =  PSOPT::inf;

    double x0 = 0.0;
    double v0 = 0.0;
    double xf = 1.0;
    double vf = 0.0;

    double t0 = 0.0;
    double tf = 1.0;

    problem.phases(1).bounds.lower.states(0) = xL;
    problem.phases(1).bounds.lower.states(1) = vL;

    problem.phases(1).bounds.upper.states(0) = xU;
    problem.phases(1).bounds.upper.states(1) = vU;

    problem.phases(1).bounds.lower.controls(0) = uL;
    problem.phases(1).bounds.upper.controls(0) = uH;

    problem.phases(1).bounds.lower.events(0) = x0;
    problem.phases(1).bounds.lower.events(1) = v0;
    problem.phases(1).bounds.lower.events(2) = xf;
    problem.phases(1).bounds.lower.events(3) = vf;


    problem.phases(1).bounds.upper.events(0) = x0;
    problem.phases(1).bounds.upper.events(1) = v0;
    problem.phases(1).bounds.upper.events(2) = xf;
    problem.phases(1).bounds.upper.events(3) = vf;

    problem.phases(1).bounds.lower.StartTime    = t0;
    problem.phases(1).bounds.upper.StartTime    = t0;

    problem.phases(1).bounds.lower.EndTime      = tf;
    problem.phases(1).bounds.upper.EndTime      = tf;



////////////////////////////////////////////////////////////////////////////
///////////////////  Register problem functions  ///////////////////////////
////////////////////////////////////////////////////////////////////////////


    problem.integrand_cost 	= &integrand_cost;
    problem.endpoint_cost 	= &endpoint_cost;
    problem.dae            	= &dae;
    problem.events 	 	    = &events;
    problem.linkages		= &linkages;

////////////////////////////////////////////////////////////////////////////
///////////////////  Define & register initial guess ///////////////////////
////////////////////////////////////////////////////////////////////////////

    int nnodes    			            = problem.phases(1).nodes(0);
    int ncontrols                       = problem.phases(1).ncontrols;
    int nstates                         = problem.phases(1).nstates;

    MatrixXd x_guess    =  zeros(nstates,nnodes);

    x_guess.row(0)  = x0*ones(1,nnodes);
    x_guess.row(1)  = v0*ones(1,nnodes);

    problem.phases(1).guess.controls       = zeros(ncontrols,nnodes);
    problem.phases(1).guess.states         = x_guess;
    problem.phases(1).guess.time           = linspace(t0,tf,nnodes);


////////////////////////////////////////////////////////////////////////////
///////////////////  Enter algorithm options  //////////////////////////////
////////////////////////////////////////////////////////////////////////////


    algorithm.nlp_iter_max                = 1000;
    algorithm.nlp_tolerance               = 1.e-4;
    algorithm.nlp_method                  = "IPOPT";
    algorithm.scaling                     = "automatic";
    algorithm.derivatives                 = "automatic";
//    algorithm.collocation_method          = "trapezoidal";
//    algorithm.mesh_refinement             = "automatic";
//    algorithm.switch_order                = 0;

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

    plot(t,x,problem.name+": states", "time (s)", "states","x v");

    plot(t,u,problem.name+": controls","time (s)", "controls", "u");

    plot(t,lambda,problem.name+": lambda","time (s)", "lambda", "lambda_x lambda_v");

    plot(t,H,problem.name+": Hamiltonian", "time (s)", "Hamiltonian","H");

    plot(t,x,problem.name+": states", "time (s)", "states","x v",
                             "pdf", "dqc_states.pdf");

    plot(t,u,problem.name+": controls","time (s)", "controls", "u",
                             "pdf", "dqc_controls.pdf");

    plot(t,lambda,problem.name+": lambda", "time (s)", "lambda","lambda_x lambda_v",
         "pdf", "dqc_lambda.pdf");

    plot(t,H,problem.name+": Hamiltonian","time (s)", "Hamiltonian", "H",
         "pdf", "dqc_H.pdf");
}

////////////////////////////////////////////////////////////////////////////
///////////////////////      END OF FILE     ///////////////////////////////
////////////////////////////////////////////////////////////////////////////
