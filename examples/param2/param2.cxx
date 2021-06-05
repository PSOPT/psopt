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

using namespace PSOPT;


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
       x1 = states[0];
       x2 = states[1];
       x3 = states[2];

    // Parameters
       p1 = parameters[0];
       p2 = parameters[1];

       t= time;

	   AutoDiffMatrix M(3,3), f(3,1), dxdt(3,1), x(3,1), y(3,1);

	   x(0) = x1; x(1) = x2; x(2) = x3;

	   M(0,0) = p2-p1*cos(p2*t);  M(0,1) = 0.0; M(0,2) = p2 + p1*sin(p2*t);
	   M(1,0) = 0.0;              M(1,1) = p1;  M(1,2) = 0.0;
	   M(2,0) = -p2+p1*sin(p2*t); M(2,1) = 0.0; M(2,2) = p2 + p1*cos(p2*t);
	   
	   f(0) = exp(t)*(-1.0 + 19.0*( cos(t) - sin(t) ) );
	   f(1) = exp(t)*(-18.0);
	   f(2) = exp(t)*(1.0 - 19.0*( cos(t) + sin(t) ) );

       // Differential equations


      product_ad(M, x, &y);

	   sum_ad(y, f, &dxdt );

      derivatives[0] = dxdt(0);

      derivatives[1] = dxdt(1);

	   derivatives[2] = dxdt(2);

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

    problem.name          			=       "Parameter estimation for ODE with two parameters";
    problem.outfilename         	=       "param2.txt";

////////////////////////////////////////////////////////////////////////////
////////////  Define problem level constants & do level 1 setup ////////////
////////////////////////////////////////////////////////////////////////////

    problem.nphases   			        = 1;
    problem.nlinkages                 = 0;

    psopt_level1_setup(problem);


/////////////////////////////////////////////////////////////////////////////
/////////   Define phase related information & do level 2 setup /////////////
/////////////////////////////////////////////////////////////////////////////

    problem.phases(1).nstates   			= 3;
    problem.phases(1).ncontrols 			= 0;
    problem.phases(1).nevents   			= 3;
    problem.phases(1).npath     			= 0;
    problem.phases(1).nparameters      = 2;
    problem.phases(1).nodes    		   << 40;
    problem.phases(1).nobserved        = 3;
    problem.phases(1).nsamples         = 129;

    psopt_level2_setup(problem, algorithm);

////////////////////////////////////////////////////////////////////////////
////////////  Load data for parameter estimation                ////////////
////////////////////////////////////////////////////////////////////////////

   int iphase = 1;
   
   load_parameter_estimation_data(problem, iphase, "../../../examples/param2/param2.dat");

   Print(problem.phases(1).observation_nodes,"observation nodes");
   Print(problem.phases(1).observations, "observations");
   Print(problem.phases(1).residual_weights, "weights");


////////////////////////////////////////////////////////////////////////////
///////////////////  Declare MatrixXd objects to store results //////////////
////////////////////////////////////////////////////////////////////////////

    MatrixXd x, p, t;

////////////////////////////////////////////////////////////////////////////
///////////////////  Enter problem bounds information //////////////////////
////////////////////////////////////////////////////////////////////////////


    problem.phases(1).bounds.lower.states(0) = 0.0;
    problem.phases(1).bounds.lower.states(1) = 0.0;
    problem.phases(1).bounds.lower.states(2) = 0.0;


    problem.phases(1).bounds.upper.states(0) = 30.0;
    problem.phases(1).bounds.upper.states(1) = 30.0;
    problem.phases(1).bounds.upper.states(2) = 30.0;


    problem.phases(1).bounds.lower.events(0) = 0.0;
    problem.phases(1).bounds.lower.events(1) = 0.0;
    problem.phases(1).bounds.lower.events(2) = 0.0;


    problem.phases(1).bounds.upper.events(0) = 0.0;
    problem.phases(1).bounds.upper.events(1) = 0.0;
    problem.phases(1).bounds.upper.events(2) = 0.0;


    problem.phases(1).bounds.lower.parameters(0)  = 0.0;
    problem.phases(1).bounds.lower.parameters(1)  = 0.0;

    problem.phases(1).bounds.upper.parameters(0)  = 30.0;
    problem.phases(1).bounds.upper.parameters(1)  = 30.0;

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

    MatrixXd state_guess(3, nnodes);
    MatrixXd param_guess(2,1);


    state_guess << linspace(1.0, exp(pi), nnodes ),
                   linspace(1.0, exp(pi), nnodes ),
                   linspace(1.0, exp(pi), nnodes );

    param_guess  << 19.0*1.5,
                    1.0*1.5;

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

    Save(x,"x.dat");
    Save(t,"t.dat");
    Print(p,"Estimated parameters");


////////////////////////////////////////////////////////////////////////////
///////////  Plot some results if desired (requires gnuplot) ///////////////
////////////////////////////////////////////////////////////////////////////

     MatrixXd tm;
     MatrixXd ym;

     tm = problem.phases(1).observation_nodes;
     ym = problem.phases(1).observations;

     spplot(t,x.row(0),tm,ym.row(0),problem.name, "time (s)", "state x1", "x1 y1");
     spplot(t,x.row(1),tm,ym.row(1),problem.name, "time (s)", "state x2", "x2 y2");
     spplot(t,x.row(2),tm,ym.row(2),problem.name, "time (s)", "state x3", "x3 y3");


     spplot(t,x.row(0),tm,ym.row(0),problem.name, "time (s)", "state x1", "x1 y1", "pdf", "x1.pdf");
     spplot(t,x.row(1),tm,ym.row(1),problem.name, "time (s)", "state x2", "x2 y2", "pdf", "x2.pdf");
     spplot(t,x.row(2),tm,ym.row(2),problem.name, "time (s)", "state x3", "x3 y3", "pdf", "x3.pdf");


}






////////////////////////////////////////////////////////////////////////////
///////////////////////      END OF FILE     ///////////////////////////////
////////////////////////////////////////////////////////////////////////////
