//////////////////////////////////////////////////////////////////////////
//////////////////         dae_i3.cxx       ////////////////////////////
//////////////////////////////////////////////////////////////////////////
////////////////           PSOPT  Example             ////////////////////
//////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////
//////// Title:  DAE Index 3                            ////////////////
//////// Last modified:  07 June 2011                  ////////////////
//////// Reference:     Schittkowski (2002)        	  ////////////////
//////// (See PSOPT handbook forf full reference)           ///////////////
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
      observations[ 0 ] = states[ 0 ];
      observations[ 1 ] = states[ 1 ];

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
       x1 = states[0];
       x2 = states[1];
       x3 = states[2];
       x4 = states[3];

    // Algebraic variables
       LAMBDA = controls[0];


    // Parameters
       L     = parameters[0];
    // Differential equations

      dx1 = x3;

      dx2 = x4;

      dx3 = LAMBDA*x1;

      dx4 = LAMBDA*x2;

      derivatives[ 0 ] = dx1;
      derivatives[ 1 ] = dx2;
      derivatives[ 2 ] = dx3;
      derivatives[ 3 ] = dx4;


     // algebraic equation

      path[ 0 ] = L*L - x1*x1 - x2*x2;

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

    problem.name          					=       "DAE Index 3";
    problem.outfilename         			=       "dae_i3.txt";

////////////////////////////////////////////////////////////////////////////
////////////  Define problem level constants & do level 1 setup ////////////
////////////////////////////////////////////////////////////////////////////

    problem.nphases   			        	= 1;
    problem.nlinkages                  = 0;

    psopt_level1_setup(problem);


/////////////////////////////////////////////////////////////////////////////
/////////   Define phase related information & do level 2 setup /////////////
/////////////////////////////////////////////////////////////////////////////

    problem.phases(1).nstates   			 = 4;
    problem.phases(1).ncontrols 			 = 1;
    problem.phases(1).nevents   			 = 0;
    problem.phases(1).npath     			 = 1;
    problem.phases(1).nparameters       = 1;
    problem.phases(1).nodes    		    << 30;
    problem.phases(1).nobserved         = 2;
    problem.phases(1).nsamples          = 20;

    psopt_level2_setup(problem, algorithm);

////////////////////////////////////////////////////////////////////////////
////////////  Load data for parameter estimation                ////////////
////////////////////////////////////////////////////////////////////////////

   int iphase = 1;
   load_parameter_estimation_data(problem, iphase, "../../../examples/dae_i3/dae_i3.dat");

   Print(problem.phases(1).observation_nodes, "observation nodes");
   Print(problem.phases(1).observations, "observations");
   Print(problem.phases(1).residual_weights, "weights");


////////////////////////////////////////////////////////////////////////////
///////////////////  Declare MatrixXd objects to store results //////////////
////////////////////////////////////////////////////////////////////////////

    MatrixXd x, u, p, t;

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

    problem.phases(1).bounds.lower.controls(0) = -10.0;
    problem.phases(1).bounds.upper.controls(0) =  10.0;

    problem.phases(1).bounds.lower.parameters(0)  = 0.0;
    problem.phases(1).bounds.upper.parameters(0)  = 5.0;


    problem.phases(1).bounds.lower.path(0)  = 0.0;
    problem.phases(1).bounds.upper.path(0)  = 0.0;

    problem.phases(1).bounds.lower.StartTime    = 0.5;
    problem.phases(1).bounds.upper.StartTime    = 0.5;

    problem.phases(1).bounds.lower.EndTime      = 10.0;
    problem.phases(1).bounds.upper.EndTime      = 10.0;

////////////////////////////////////////////////////////////////////////////
///////////////////  Register problem functions  ///////////////////////////
////////////////////////////////////////////////////////////////////////////

    problem.dae 						= &dae;
    problem.events 					= &events;
    problem.linkages					= &linkages;
    problem.observation_function = & observation_function;

////////////////////////////////////////////////////////////////////////////
///////////////////  Define & register initial guess ///////////////////////
////////////////////////////////////////////////////////////////////////////

    int nnodes =     (int) problem.phases(1).nsamples;
    
    MatrixXd state_guess(4, nnodes);
    MatrixXd control_guess(1,nnodes);
    MatrixXd param_guess(1,1);

    state_guess <<  	problem.phases(1).observations.row(0), 
    						problem.phases(1).observations.row(1), 
    						ones(1,nnodes),
    						ones(1,nnodes);

    control_guess = zeros(1,nnodes);

    param_guess << 0.5;


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

    Save(x,"x.dat");
    Save(u,"u.dat");
    Save(t,"t.dat");
    Print(p,"Estimated parameter");


////////////////////////////////////////////////////////////////////////////
///////////  Plot some results if desired (requires gnuplot) ///////////////
////////////////////////////////////////////////////////////////////////////

     MatrixXd tm;
     MatrixXd ym;

     tm = problem.phases(1).observation_nodes;
     ym = problem.phases(1).observations;


     plot(t,x.row(0) ,tm,ym.row(0) ,problem.name, "time (s)", "state x1", "x1 yhat1");
     plot(t,x.row(1) ,tm,ym.row(1) ,problem.name, "time (s)", "state x2", "x2 yhat2");
     plot(t,u,problem.name, "time (s)", "algebraic state u", "u");

     plot(t,x.row(0),tm,ym.row(0),problem.name, "time (s)", "state x1", "x1 yhat1",
	  "pdf", "x1.pdf");
     plot(t,x.row(1),tm,ym.row(1),problem.name, "time (s)", "state x2", "x2 yhat2",
	  "pdf", "x2.pdf");
     plot(t,u,problem.name, "time (s)", "algebraic state lambda", "lambda", "pdf", "lambda.pdf");


}






////////////////////////////////////////////////////////////////////////////
///////////////////////      END OF FILE     ///////////////////////////////
////////////////////////////////////////////////////////////////////////////
