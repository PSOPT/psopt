//////////////////////////////////////////////////////////////////////////
////////////////        bryson_max_range.cxx         /////////////////////
//////////////////////////////////////////////////////////////////////////
////////////////           PSOPT  Example            /////////////////////
//////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////
//////// Title:         Bryson maximum range problem      ////////////////
//////// Last modified: 05 January 2009                   ////////////////
//////// Reference:     Bryson and Ho (1975)              ////////////////
//////// (See PSOPT handbook for full reference)          ////////////////
//////////////////////////////////////////////////////////////////////////
////////     Copyright (c) Victor M. Becerra, 2009        ////////////////
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
   adouble x = final_states[0];

   return (-x);
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
   adouble xdot, ydot, vdot;

   double g = 1.0;
   double a = 0.5*g;

   adouble x = states[ 0 ];
   adouble y = states[ 1 ];
   adouble v = states[ 2 ];

   adouble u1 = controls[ 0 ];
   adouble u2 = controls[ 1 ];

   xdot = v*u1;
   ydot = v*u2;
   vdot = a-g*u2;

   derivatives[ 0 ] = xdot;
   derivatives[ 1 ] = ydot;
   derivatives[ 2 ] = vdot;

   path[ 0 ] = (u1*u1) + (u2*u2);

}

////////////////////////////////////////////////////////////////////////////
///////////////////  Define the events function ////////////////////////////
////////////////////////////////////////////////////////////////////////////

void events(adouble* e, adouble* initial_states, adouble* final_states,
            adouble* parameters,adouble& t0, adouble& tf, adouble* xad,
            int iphase, Workspace* workspace)

{
   adouble x0 = initial_states[ 0 ];
   adouble y0 = initial_states[ 1 ];
   adouble v0 = initial_states[ 2 ];
   adouble xf = final_states[ 0 ];
   adouble yf = final_states[ 1 ];

   e[ 0 ] = x0;
   e[ 1 ] = y0;
   e[ 2 ] = v0;
   e[ 3 ] = yf;

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

    problem.name        		= "Bryson Maximum Range Problem";
    problem.outfilename                 = "brymr.txt";

////////////////////////////////////////////////////////////////////////////
////////////  Define problem level constants & do level 1 setup ////////////
////////////////////////////////////////////////////////////////////////////

    problem.nphases   			= 1;
    problem.nlinkages                   = 0;

    psopt_level1_setup(problem);

/////////////////////////////////////////////////////////////////////////////
/////////   Define phase related information & do level 2 setup  ////////////
/////////////////////////////////////////////////////////////////////////////

    problem.phases(1).nstates   		= 3;
    problem.phases(1).ncontrols 		= 2;
    problem.phases(1).nevents   		= 4;
    problem.phases(1).npath     		= 1;
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

    double xL = -10.0;
    double yL = -10.0;
    double vL = -10.0;
    double xU = 10.0;
    double yU = 10.0;
    double vU = 10.0;

    double u1L = -10.0;
    double u2L = -10.0;
    double u1U = 10.0;
    double u2U = 10.0;

    double x0 = 0.0;
    double y0 = 0.0;
    double v0 = 0.0;
    double yf = 0.1;


    problem.phases(1).bounds.lower.states(0) = xL;
    problem.phases(1).bounds.lower.states(1) = yL;
    problem.phases(1).bounds.lower.states(2) = vL;


    problem.phases(1).bounds.upper.states(0) = xU;
    problem.phases(1).bounds.upper.states(1) = yU;
    problem.phases(1).bounds.upper.states(2) = vU;


    problem.phases(1).bounds.lower.controls(0) = u1L;
    problem.phases(1).bounds.lower.controls(1) = u2L;
    problem.phases(1).bounds.upper.controls(0) = u1U;
    problem.phases(1).bounds.upper.controls(1) = u2U;

    problem.phases(1).bounds.lower.events(0) = x0;
    problem.phases(1).bounds.lower.events(1) = y0;
    problem.phases(1).bounds.lower.events(2) = v0;
    problem.phases(1).bounds.lower.events(3) = yf;


    problem.phases(1).bounds.upper.events(0) = x0;
    problem.phases(1).bounds.upper.events(1) = y0;
    problem.phases(1).bounds.upper.events(2) = v0;
    problem.phases(1).bounds.upper.events(3) = yf;

    problem.phases(1).bounds.upper.path(0) = 1.0;
    problem.phases(1).bounds.lower.path(0) = 1.0;



    problem.phases(1).bounds.lower.StartTime    = 0.0;
    problem.phases(1).bounds.upper.StartTime    = 0.0;

    problem.phases(1).bounds.lower.EndTime      = 2.0;
    problem.phases(1).bounds.upper.EndTime      = 2.0;



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

    int nnodes    			             = problem.phases(1).nodes(0);
    int ncontrols                       = problem.phases(1).ncontrols;
    int nstates                         = problem.phases(1).nstates;

    MatrixXd x_guess    =  zeros(nstates,nnodes);

    x_guess.row(0)  = x0*ones(1,nnodes);
    x_guess.row(1)  = y0*ones(1,nnodes);
    x_guess.row(2)  = v0*ones(1,nnodes);

    problem.phases(1).guess.controls       = zeros(ncontrols,nnodes);
    problem.phases(1).guess.states         = x_guess;
    problem.phases(1).guess.time           = linspace(0.0,2.0,nnodes);


////////////////////////////////////////////////////////////////////////////
///////////////////  Enter algorithm options  //////////////////////////////
////////////////////////////////////////////////////////////////////////////


    algorithm.nlp_iter_max                = 1000;
    algorithm.nlp_tolerance               = 1.e-4;
    algorithm.nlp_method                  = "IPOPT";
    algorithm.scaling                     = "automatic";
    algorithm.derivatives                 = "automatic";
//    algorithm.mesh_refinement             = "automatic";
    algorithm.collocation_method          = "trapezoidal";
//    algorithm.defect_scaling = "jacobian-based";
    algorithm.ode_tolerance               = 1.e-6;



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

    plot(t,x,problem.name+": states", "time (s)", "states","x y v");

    plot(t,u,problem.name+": controls","time (s)", "controls", "u_1 u_2");

    plot(t,x,problem.name+": states", "time (s)", "states","x y v",
                             "pdf", "brymr_states.pdf");

    plot(t,u,problem.name+": controls","time (s)", "controls", "u_1 u_2",
                             "pdf", "brymr_controls.pdf");
}

////////////////////////////////////////////////////////////////////////////
///////////////////////      END OF FILE     ///////////////////////////////
////////////////////////////////////////////////////////////////////////////
