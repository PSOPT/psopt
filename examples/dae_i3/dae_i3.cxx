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

using namespace std;

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

int main(int argc, char* argv[])
{
// Optional arguments select the integrated-residual variants of the comparison reported in
// the book, so that the table can be reproduced without editing the source:
//
//     ./dae_i3                       Legendre collocation on 30 nodes (the default)
//     ./dae_i3 <d> <delta> [in|out] [m]
//                                    integrated residual on 33 nodes, local order d,
//                                    residual box delta, holonomic constraint folded into
//                                    the residual ("in", the default) or left as a
//                                    pointwise path constraint ("out"), and m residual
//                                    points per element (default d+2)
//
// d must divide the number of intervals: 32 intervals admit d = 2, 4, 8, 16.
//
// The last argument is worth a word, because it is the pitfall this example was found to
// expose. The residual box binds only where the quadrature rule samples the residual. At
// d = 4 and delta = 1e-6 the reported maximum relative local error is 9.7e-2 with m = 4,
// 7.9e-4 with m = 5, and 1.5e-8 -- the box itself -- from m = 6 = d+2 upward: with too few
// points the element polynomial simply oscillates between them and the certificate is
// worthless. PSOPT now refuses m < d+2 for that reason.
    int    ir_d     = (argc > 1) ? atoi(argv[1]) : 0;
    double ir_delta = (argc > 2) ? atof(argv[2]) : 1.0e-6;
    bool   alg_in   = (argc > 3) ? (std::string(argv[3]) != "out") : true;
    int    ir_m     = (argc > 4) ? atoi(argv[4]) : 0;


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
    problem.phases(1).nodes    		    << ( (ir_d >= 2) ? 33 : 30 );
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

    // The pendulum length is bounded away from zero. With L allowed to reach zero the
    // problem has a degenerate stationary point -- x = 0 with L = 0 satisfies the
    // dynamics and the holonomic constraint exactly, fits none of the twenty
    // observations, and gives J = 20 -- and the integrated-residual solves fall into it
    // readily. Which of the two they return has turned out to depend on incidental
    // details of the NLP rather than on the transcription: it changed when a block of
    // inert constraint rows was written instead of left unwritten, it changes with the
    // residual box tolerance, and it changed again when the unread midpoint control
    // variables were removed from the decision vector. Excluding a pendulum of zero
    // length is a statement about the model, not about the answer, and it makes the
    // reported estimate a property of the method rather than of the barrier path.
    //
    // The value matters. At 0.1 the bound perturbs the collocation solve enough to stop
    // it converging -- it returns L = 1.000010 in restoration where it used to return
    // 1.0000000 cleanly -- and at 0.01 it is too close to zero to keep the
    // integrated-residual solves off the degenerate branch, which then simply stop at
    // the bound. At 0.05 the collocation solve is untouched and the residual solves
    // either return an admissible estimate or fail honestly.
    problem.phases(1).bounds.lower.parameters(0)  = 0.05;
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

    if (ir_d >= 2) {
        // The Nie-Kerrigan flexible-order local representation, over a Hermite-Simpson mesh:
        // the residual box is the accuracy specification and the local order d is what is
        // being varied.  With ir_include_path = "auto" the holonomic constraint joins the
        // residual, so the index-3 system is solved as posed rather than by index reduction.
        algorithm.collocation_method   = "Hermite-Simpson";
        algorithm.transcription_method = "integrated-residual";
        algorithm.ir_objective         = "cost";
        algorithm.ir_local_order       = ir_d;
        algorithm.ir_residual_bound    = ir_delta;
        algorithm.ir_residual_nodes    = (ir_m > 0) ? ir_m : (ir_d + 2);
        algorithm.ir_include_path      = alg_in ? "auto" : "none";
        algorithm.nlp_iter_max         = 3000;
    }


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
