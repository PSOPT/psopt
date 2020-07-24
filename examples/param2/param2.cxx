//////////////////////////////////////////////////////////////////////////
//////////////////         param2.cxx         ////////////////////////////
//////////////////////////////////////////////////////////////////////////
////////////////           PSOPT  Example             ////////////////////
//////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////
//////// Title:  Parameter estimation for ODE with 2 parameters //////////
//////// Last modified:  31 Jan 2014                      ////////////////
//////// Reference:     Li et al (2005)               	  ////////////////
//////// (See PSOPT handbook for full reference)           ///////////////
//////////////////////////////////////////////////////////////////////////
////////     Copyright (c) Victor M. Becerra, 2014        ////////////////
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
	int i;
	for (i=0; i<3; i++) {
            observations[ i ] = states[ i ];
	}
}



//////////////////////////////////////////////////////////////////////////
///////////////////  Define the DAE's ////////////////////////////////////
//////////////////////////////////////////////////////////////////////////

void dae(adouble* derivatives, adouble* path, adouble* states,
         adouble* controls, adouble* parameters, adouble& time,
         adouble* xad, int iphase, Workspace* workspace)
{

    // Variables
       adouble x1, x2, x3, p1, p2, t;

    // Differential states
       x1 = states[CINDEX(1)];
       x2 = states[CINDEX(2)];
       x3 = states[CINDEX(3)];

    // Parameters
       p1 = parameters[CINDEX(1)];
       p2 = parameters[CINDEX(2)];

       t= time;

	   ADMatrix M(3,3), f(3,1), dxdt(3,1), x(3,1), y(3,1);

	   x(1,1) = x1; x(2,1) = x2; x(3,1) = x3;

	   M(1,1) = p2-p1*cos(p2*t);  M(1,2) = 0.0; M(1,3) = p2 + p1*sin(p2*t);
	   M(2,1) = 0.0;              M(2,2) = p1;  M(2,3) = 0.0;
	   M(3,1) = -p2+p1*sin(p2*t); M(3,2) = 0.0; M(3,3) = p2 + p1*cos(p2*t);
	   
	   f(1,1) = exp(t)*(-1.0 + 19.0*( cos(t) - sin(t) ) );
	   f(2,1) = exp(t)*(-18.0);
	   f(3,1) = exp(t)*(1.0 - 19.0*( cos(t) + sin(t) ) );

       // Differential equations


           product_ad(M, x, &y);

	   sum_ad(y, f, &dxdt );

           derivatives[CINDEX(1)] = dxdt(1,1);

           derivatives[CINDEX(2)] = dxdt(2,1);

	   derivatives[CINDEX(3)] = dxdt(3,1);

}

////////////////////////////////////////////////////////////////////////////
///////////////////  Define the events function ////////////////////////////
////////////////////////////////////////////////////////////////////////////

void events(adouble* e, adouble* initial_states, adouble* final_states,
            adouble* parameters,adouble& t0, adouble& tf, adouble* xad,
            int iphase, Workspace* workspace)
{
   int i;
   double rhs = 1.0 + exp(pi);

   for (i=0; i<3; i++) {

      e[i]  = initial_states[i] + final_states[i] - rhs;

   }

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

    problem.name          	=       "Parameter estimation for ODE with two parameters";
    problem.outfilename         =       "param2.txt";

////////////////////////////////////////////////////////////////////////////
////////////  Define problem level constants & do level 1 setup ////////////
////////////////////////////////////////////////////////////////////////////

    problem.nphases   			        = 1;
    problem.nlinkages                   = 0;

    psopt_level1_setup(problem);


/////////////////////////////////////////////////////////////////////////////
/////////   Define phase related information & do level 2 setup /////////////
/////////////////////////////////////////////////////////////////////////////

    problem.phases(1).nstates   		= 3;
    problem.phases(1).ncontrols 		= 0;
    problem.phases(1).nevents   		= 3;
    problem.phases(1).npath     		= 0;
    problem.phases(1).nparameters               = 2;
    problem.phases(1).nodes    		        = 40;
    problem.phases(1).nobserved                 = 3;
    problem.phases(1).nsamples                  = 129;

    psopt_level2_setup(problem, algorithm);

////////////////////////////////////////////////////////////////////////////
////////////  Load data for parameter estimation                ////////////
////////////////////////////////////////////////////////////////////////////

   int iphase = 1;
   load_parameter_estimation_data(problem, iphase, "param2.dat");

   problem.phases(1).observation_nodes.Print("observation nodes");
   problem.phases(1).observations.Print("observations");
   problem.phases(1).residual_weights.Print("weights");


////////////////////////////////////////////////////////////////////////////
///////////////////  Declare DMatrix objects to store results //////////////
////////////////////////////////////////////////////////////////////////////

    DMatrix x, p, t;

////////////////////////////////////////////////////////////////////////////
///////////////////  Enter problem bounds information //////////////////////
////////////////////////////////////////////////////////////////////////////


    problem.phases(1).bounds.lower.states(1) = 0.0;
    problem.phases(1).bounds.lower.states(2) = 0.0;
    problem.phases(1).bounds.lower.states(3) = 0.0;


    problem.phases(1).bounds.upper.states(1) = 30.0;
    problem.phases(1).bounds.upper.states(2) = 30.0;
    problem.phases(1).bounds.upper.states(3) = 30.0;


    problem.phases(1).bounds.lower.events(1) = 0.0;
    problem.phases(1).bounds.lower.events(2) = 0.0;
    problem.phases(1).bounds.lower.events(3) = 0.0;


    problem.phases(1).bounds.upper.events(1) = 0.0;
    problem.phases(1).bounds.upper.events(2) = 0.0;
    problem.phases(1).bounds.upper.events(3) = 0.0;


    problem.phases(1).bounds.lower.parameters(1)  = 0.0;
    problem.phases(1).bounds.lower.parameters(2)  = 0.0;

    problem.phases(1).bounds.upper.parameters(1)  = 30.0;
    problem.phases(1).bounds.upper.parameters(2)  = 30.0;

    problem.phases(1).bounds.lower.StartTime    = 0.0;
    problem.phases(1).bounds.upper.StartTime    = 0.0;

    problem.phases(1).bounds.lower.EndTime      = pi;
    problem.phases(1).bounds.upper.EndTime      = pi;

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

    DMatrix state_guess(3, nnodes);
    DMatrix param_guess(2,1);


    state_guess(1,colon()) = linspace(1.0, exp(pi), nnodes );
    state_guess(2,colon()) = linspace(1.0, exp(pi), nnodes );
    state_guess(3,colon()) = linspace(1.0, exp(pi), nnodes );

    param_guess(1) = 19.0*1.5;
    param_guess(2) = 1.0*1.5;

    problem.phases(1).guess.states        = state_guess;
    problem.phases(1).guess.time          = linspace(0.0, pi, nnodes);
    problem.phases(1).guess.parameters    = param_guess;



////////////////////////////////////////////////////////////////////////////
///////////////////  Enter algorithm options  //////////////////////////////
////////////////////////////////////////////////////////////////////////////

    algorithm.nlp_method                  = "IPOPT";
    algorithm.scaling                     = "automatic";
    algorithm.derivatives                 = "automatic";
    algorithm.collocation_method          = "Legendre";
//    algorithm.mesh_refinement             = "automatic";
//    algorithm.defect_scaling              = "jacobian-based";
//    algorithm.nlp_iter_max                = 1000;
    algorithm.ode_tolerance               = 1.e-4;

////////////////////////////////////////////////////////////////////////////
///////////////////  Now call PSOPT to solve the problem   //////////////////
////////////////////////////////////////////////////////////////////////////

    psopt(solution, problem, algorithm);

////////////////////////////////////////////////////////////////////////////
///////////  Extract relevant variables from solution structure   //////////
////////////////////////////////////////////////////////////////////////////

    x = solution.get_states_in_phase(1);
    t = solution.get_time_in_phase(1);
    p = solution.get_parameters_in_phase(1);


////////////////////////////////////////////////////////////////////////////
///////////  Save solution data to files if desired ////////////////////////
////////////////////////////////////////////////////////////////////////////

    x.Save("x.dat");
    t.Save("t.dat");
    p.Print("Estimated parameters");


////////////////////////////////////////////////////////////////////////////
///////////  Plot some results if desired (requires gnuplot) ///////////////
////////////////////////////////////////////////////////////////////////////

     DMatrix tm;
     DMatrix ym;

     tm = problem.phases(1).observation_nodes;
     ym = problem.phases(1).observations;

     spplot(t,x(1,colon()),tm,ym(1,colon()),problem.name, "time (s)", "state x1", "x1 y1");
     spplot(t,x(2,colon()),tm,ym(2,colon()),problem.name, "time (s)", "state x2", "x2 y2");
     spplot(t,x(3,colon()),tm,ym(3,colon()),problem.name, "time (s)", "state x3", "x3 y3");


     spplot(t,x(1,colon()),tm,ym(1,colon()),problem.name, "time (s)", "state x1", "x1 y1", "pdf", "x1.pdf");
     spplot(t,x(2,colon()),tm,ym(2,colon()),problem.name, "time (s)", "state x2", "x2 y2", "pdf", "x2.pdf");
     spplot(t,x(3,colon()),tm,ym(3,colon()),problem.name, "time (s)", "state x3", "x3 y3", "pdf", "x3.pdf");


}






////////////////////////////////////////////////////////////////////////////
///////////////////////      END OF FILE     ///////////////////////////////
////////////////////////////////////////////////////////////////////////////
