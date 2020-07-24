//////////////////////////////////////////////////////////////////////////
//////////////////         ident1.cxx       ////////////////////////////
//////////////////////////////////////////////////////////////////////////
////////////////           PSOPT  Example             ////////////////////
//////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////
//////// Title:  DAE Index 3                            ////////////////
//////// Last modified:  07 June 2011                  ////////////////
//////// Reference:     Schittkowski (2002)        	  ////////////////
//////// (See PSOPT handbook for full reference)           ///////////////
//////////////////////////////////////////////////////////////////////////
////////     Copyright (c) Victor M. Becerra, 2011        ////////////////
//////////////////////////////////////////////////////////////////////////
//////// This is part of the PSOPT software library, which ///////////////
//////// is distributed under the terms of the GNU Lesser ////////////////
//////// General Public License (LGPL)                    ////////////////
//////////////////////////////////////////////////////////////////////////

#include "psopt.h"

//////////////////////////////////////////////////////////////////////////
///////////////////  Define the observation function //////////
//////////////////////////////////////////////////////////////////////////



void  observation_function( adouble* observations,
                            adouble* states, adouble* controls,
                            adouble* parameters, adouble& time, int k,
                            adouble* xad, int iphase, Workspace* workspace)
{
      observations[ CINDEX(1) ] = states[ CINDEX(1) ];
      observations[ CINDEX(2) ] = states[ CINDEX(2) ];

}



//////////////////////////////////////////////////////////////////////////
///////////////////  Define the DAE's ////////////////////////////////////
//////////////////////////////////////////////////////////////////////////

void dae(adouble* derivatives, adouble* path, adouble* states,
         adouble* controls, adouble* parameters, adouble& time,
         adouble* xad, int iphase, Workspace* workspace)
{

    // Variables
       adouble x1, x2, x3, x4, L, OMEGA, LAMBDA;
       adouble dx1, dx2, dx3, dx4;

    // Differential states
       x1 = states[CINDEX(1)];
       x2 = states[CINDEX(2)];
       x3 = states[CINDEX(3)];
       x4 = states[CINDEX(4)];

    // Algebraic variables
       LAMBDA = controls[CINDEX(1)];


    // Parameters
       L     = parameters[CINDEX(1)];
    // Differential equations

      dx1 = x3;

      dx2 = x4;

      dx3 = LAMBDA*x1;

      dx4 = LAMBDA*x2;

      derivatives[ CINDEX(1) ] = dx1;
      derivatives[ CINDEX(2) ] = dx2;
      derivatives[ CINDEX(3) ] = dx3;
      derivatives[ CINDEX(4) ] = dx4;


     // algebraic equation

      path[ CINDEX(1) ] = L*L - x1*x1 - x2*x2;

}

////////////////////////////////////////////////////////////////////////////
///////////////////  Define the events function ////////////////////////////
////////////////////////////////////////////////////////////////////////////

void events(adouble* e, adouble* initial_states, adouble* final_states,
            adouble* parameters,adouble& t0, adouble& tf, adouble* xad,
            int iphase, Workspace* workspace)
{
       // no events

       return;

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

    problem.name          		=       "DAE Index 3";
    problem.outfilename         =       "dae_i3.txt";

////////////////////////////////////////////////////////////////////////////
////////////  Define problem level constants & do level 1 setup ////////////
////////////////////////////////////////////////////////////////////////////

    problem.nphases   			        = 1;
    problem.nlinkages                   = 0;

    psopt_level1_setup(problem);


/////////////////////////////////////////////////////////////////////////////
/////////   Define phase related information & do level 2 setup /////////////
/////////////////////////////////////////////////////////////////////////////

    problem.phases(1).nstates   		= 4;
    problem.phases(1).ncontrols 		= 1;
    problem.phases(1).nevents   		= 0;
    problem.phases(1).npath     		= 1;
    problem.phases(1).nparameters       = 1;
    problem.phases(1).nodes    		    = "[10,20,30]";
    problem.phases(1).nobserved         = 2;
    problem.phases(1).nsamples          = 20;

    psopt_level2_setup(problem, algorithm);

////////////////////////////////////////////////////////////////////////////
////////////  Load data for parameter estimation                ////////////
////////////////////////////////////////////////////////////////////////////

   int iphase = 1;
   load_parameter_estimation_data(problem, iphase, "dae_i3.dat");

   problem.phases(1).observation_nodes.Print("observation nodes");
   problem.phases(1).observations.Print("observations");
   problem.phases(1).residual_weights.Print("weights");


////////////////////////////////////////////////////////////////////////////
///////////////////  Declare DMatrix objects to store results //////////////
////////////////////////////////////////////////////////////////////////////

    DMatrix x, u, p, t;

////////////////////////////////////////////////////////////////////////////
///////////////////  Enter problem bounds information //////////////////////
////////////////////////////////////////////////////////////////////////////


    problem.phases(1).bounds.lower.states(1) = -2.0;
    problem.phases(1).bounds.lower.states(2) = -2.0;
    problem.phases(1).bounds.lower.states(3) = -2.0;
    problem.phases(1).bounds.lower.states(4) = -2.0;

    problem.phases(1).bounds.upper.states(1) = 2.0;
    problem.phases(1).bounds.upper.states(2) = 2.0;
    problem.phases(1).bounds.upper.states(3) = 2.0;
    problem.phases(1).bounds.upper.states(4) = 2.0;

    problem.phases(1).bounds.lower.controls(1) = -10.0;
    problem.phases(1).bounds.upper.controls(1) =  10.0;

    problem.phases(1).bounds.lower.parameters(1)  = 0.0;
    problem.phases(1).bounds.upper.parameters(1)  = 5.0;


    problem.phases(1).bounds.lower.path(1)  = 0.0;
    problem.phases(1).bounds.upper.path(1)  = 0.0;

    problem.phases(1).bounds.lower.StartTime    = 0.5;
    problem.phases(1).bounds.upper.StartTime    = 0.5;

    problem.phases(1).bounds.lower.EndTime      = 10.0;
    problem.phases(1).bounds.upper.EndTime      = 10.0;

////////////////////////////////////////////////////////////////////////////
///////////////////  Register problem functions  ///////////////////////////
////////////////////////////////////////////////////////////////////////////

    problem.dae 		= &dae;
    problem.events 		= &events;
    problem.linkages		= &linkages;
    problem.observation_function = & observation_function;

////////////////////////////////////////////////////////////////////////////
///////////////////  Define & register initial guess ///////////////////////
////////////////////////////////////////////////////////////////////////////

    int nnodes =     (int) problem.phases(1).nsamples;
    DMatrix state_guess(4, nnodes);
    DMatrix control_guess(1,nnodes);
    DMatrix param_guess(1,1);

    state_guess(1,colon()) = problem.phases(1).observations(1,colon());
    state_guess(2,colon()) = problem.phases(1).observations(2,colon());
    state_guess(3,colon()) = ones(1,nnodes);
    state_guess(4,colon()) = ones(1,nnodes);

    control_guess(1,colon()) = zeros(1,nnodes);

    param_guess = 0.5;


    problem.phases(1).guess.states        = state_guess;
    problem.phases(1).guess.time          = problem.phases(1).observation_nodes;
    problem.phases(1).guess.parameters    = param_guess;
    problem.phases(1).guess.controls      = control_guess;



////////////////////////////////////////////////////////////////////////////
///////////////////  Enter algorithm options  //////////////////////////////
////////////////////////////////////////////////////////////////////////////

    algorithm.nlp_method                  = "IPOPT";
    algorithm.scaling                     = "automatic";
    algorithm.derivatives                 = "automatic";
    algorithm.collocation_method          = "Legendre";
    algorithm.jac_sparsity_ratio          = 0.50;

////////////////////////////////////////////////////////////////////////////
///////////////////  Now call PSOPT to solve the problem   //////////////////
////////////////////////////////////////////////////////////////////////////

    psopt(solution, problem, algorithm);

////////////////////////////////////////////////////////////////////////////
///////////  Extract relevant variables from solution structure   //////////
////////////////////////////////////////////////////////////////////////////

    x = solution.get_states_in_phase(1);
    u = solution.get_controls_in_phase(1);
    t = solution.get_time_in_phase(1);
    p = solution.get_parameters_in_phase(1);


////////////////////////////////////////////////////////////////////////////
///////////  Save solution data to files if desired ////////////////////////
////////////////////////////////////////////////////////////////////////////

    x.Save("x.dat");
    u.Save("u.dat");
    t.Save("t.dat");
    p.Print("Estimated parameter");


////////////////////////////////////////////////////////////////////////////
///////////  Plot some results if desired (requires gnuplot) ///////////////
////////////////////////////////////////////////////////////////////////////

     DMatrix tm;
     DMatrix ym;

     tm = problem.phases(1).observation_nodes;
     ym = problem.phases(1).observations;


     plot(t,x(1,colon()),tm,ym(1,colon()),problem.name, "time (s)", "state x1", "x1 yhat1");
     plot(t,x(2,colon()),tm,ym(2,colon()),problem.name, "time (s)", "state x2", "x2 yhat2");
     plot(t,u,problem.name, "time (s)", "algebraic state u", "u");

     plot(t,x(1,colon()),tm,ym(1,colon()),problem.name, "time (s)", "state x1", "x1 yhat1",
	  "pdf", "x1.pdf");
     plot(t,x(2,colon()),tm,ym(2,colon()),problem.name, "time (s)", "state x2", "x2 yhat2",
	  "pdf", "x2.pdf");
     plot(t,u,problem.name, "time (s)", "algebraic state lambda", "lambda", "pdf", "lambda.pdf");


}






////////////////////////////////////////////////////////////////////////////
///////////////////////      END OF FILE     ///////////////////////////////
////////////////////////////////////////////////////////////////////////////
