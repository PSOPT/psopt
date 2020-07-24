//////////////////////////////////////////////////////////////////////////
////////////////        alpine_rider.cxx             /////////////////////
//////////////////////////////////////////////////////////////////////////
////////////////           PSOPT  Example            /////////////////////
//////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////
//////// Title:         Alp rider problem                 ////////////////
//////// Last modified: 09 February 2009                  ////////////////
//////// Reference:     Betts (2001)                      ////////////////
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
///////////////////  Define auxiliary function    ////////////////////////
//////////////////////////////////////////////////////////////////////////

adouble pk( adouble t, double a, double b)
{
   return exp(-b*(t-a)*(t-a));
}

//////////////////////////////////////////////////////////////////////////
///////////////////  Define the end point (Mayer) cost function //////////
//////////////////////////////////////////////////////////////////////////

adouble endpoint_cost(adouble* initial_states, adouble* final_states,
                      adouble* parameters,adouble& t0, adouble& tf,
                      adouble* xad, int iphase, Workspace* workspace)
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
    adouble x1 = states[CINDEX(1)];
    adouble x2 = states[CINDEX(2)];
    adouble x3 = states[CINDEX(3)];
    adouble x4 = states[CINDEX(4)];
    adouble u1 = controls[CINDEX(1)];
    adouble u2 = controls[CINDEX(2)];

    adouble L;

    L = 100.0*( x1*x1 + x2*x2 + x3*x3 + x4*x4 ) + 0.01*( u1*u1 + u2*u2 );

    return  L;
}


//////////////////////////////////////////////////////////////////////////
///////////////////  Define the DAE's ////////////////////////////////////
//////////////////////////////////////////////////////////////////////////

void dae(adouble* derivatives, adouble* path, adouble* states,
         adouble* controls, adouble* parameters, adouble& time,
         adouble* xad, int iphase, Workspace* workspace)
{
    adouble x1 = states[CINDEX(1)];
    adouble x2 = states[CINDEX(2)];
    adouble x3 = states[CINDEX(3)];
    adouble x4 = states[CINDEX(4)];
    adouble u1 = controls[CINDEX(1)];
    adouble u2 = controls[CINDEX(2)];
    adouble t  = time;

   derivatives[CINDEX(1)] = -10*x1 + u1 + u2;
   derivatives[CINDEX(2)] = -2*x2 + u1 + 2*u2;
   derivatives[CINDEX(3)] = -3*x3 + 5*x4 + u1 - u2;
   derivatives[CINDEX(4)] =  5*x3 - 3*x4 + u1 + 3*u2;

   path[0] = x1*x1 + x2*x2 + x3*x3 + x4*x4
             - 3*pk(t,3,12) - 3*pk(t,6,10) - 3*pk(t,10,6) - 8*pk(t,15,4)
             - 0.01;

}

////////////////////////////////////////////////////////////////////////////
///////////////////  Define the events function ////////////////////////////
////////////////////////////////////////////////////////////////////////////

void events(adouble* e, adouble* initial_states, adouble* final_states,
            adouble* parameters,adouble& t0, adouble& tf, adouble* xad,
            int iphase, Workspace* workspace)

{
   adouble x1i = initial_states[ CINDEX(1) ];
   adouble x2i = initial_states[ CINDEX(2) ];
   adouble x3i = initial_states[ CINDEX(3) ];
   adouble x4i = initial_states[ CINDEX(4) ];
   adouble x1f = final_states[ CINDEX(1) ];
   adouble x2f = final_states[ CINDEX(2) ];
   adouble x3f = final_states[ CINDEX(3) ];
   adouble x4f = final_states[ CINDEX(4) ];



   e[ CINDEX(1) ] = x1i;
   e[ CINDEX(2) ] = x2i;
   e[ CINDEX(3) ] = x3i;
   e[ CINDEX(4) ] = x4i;
   e[ CINDEX(5) ] = x1f;
   e[ CINDEX(6) ] = x2f;
   e[ CINDEX(7) ] = x3f;
   e[ CINDEX(8) ] = x4f;

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

    problem.name        		= "Alp rider problem";
    problem.outfilename                 = "alpine.txt";

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
    problem.phases(1).ncontrols 		= 2;
    problem.phases(1).nevents   		= 8;
    problem.phases(1).npath     		= 1;
    problem.phases(1).nodes                     = "[120]";



    psopt_level2_setup(problem, algorithm);

////////////////////////////////////////////////////////////////////////////
///////////////////  Enter problem bounds information //////////////////////
////////////////////////////////////////////////////////////////////////////


    problem.phases(1).bounds.lower.states = "[-4.0, -4.0, -4.0, -4.0]";

    problem.phases(1).bounds.upper.states = "[ 4.0,  4.0,  4.0,  4.0]";

    problem.phases(1).bounds.lower.controls = "[ -500.0, -500 ]";

    problem.phases(1).bounds.upper.controls = "[  500.0,  500 ]";

    problem.phases(1).bounds.lower.events =  "[2.0, 1.0, 2.0, 1.0, 2.0, 3.0, 1.0, -2.0]";

    problem.phases(1).bounds.upper.events = problem.phases(1).bounds.lower.events;

    problem.phases(1).bounds.upper.path         = 100.0;

    problem.phases(1).bounds.lower.path         = 0.0;

    problem.phases(1).bounds.lower.StartTime    = 0.0;

    problem.phases(1).bounds.upper.StartTime    = 0.0;

    problem.phases(1).bounds.lower.EndTime      = 20.0;

    problem.phases(1).bounds.upper.EndTime      = 20.0;


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

    int nnodes    			    = problem.phases(1).nodes(1);

    DMatrix x_guess    		            =  zeros(4,nnodes);

    x_guess(1,colon()) 			    = linspace(2,1,nnodes);
    x_guess(2,colon()) 			    = linspace(2,3,nnodes);
    x_guess(3,colon()) 			    = linspace(2,1,nnodes);
    x_guess(4,colon()) 			    = linspace(1,-2,nnodes);

    problem.phases(1).guess.controls        = zeros(2,nnodes);
    problem.phases(1).guess.states          = x_guess;
    problem.phases(1).guess.time            = linspace(0.0,20.0,nnodes+1);


////////////////////////////////////////////////////////////////////////////
///////////////////  Enter algorithm options  //////////////////////////////
////////////////////////////////////////////////////////////////////////////


    algorithm.nlp_iter_max                = 1000;
    algorithm.nlp_tolerance               = 1.e-6;
    algorithm.nlp_method                  = "IPOPT";
    algorithm.scaling                     = "automatic";
    algorithm.derivatives                 = "automatic";
    algorithm.jac_sparsity_ratio          = 0.20;
    algorithm.collocation_method          = "Legendre";
    algorithm.diff_matrix                 = "central-differences";
    algorithm.mesh_refinement             = "automatic";
    algorithm.mr_max_increment_factor     = 0.3;
    algorithm.defect_scaling              = "jacobian-based";






////////////////////////////////////////////////////////////////////////////
///////////////////  Now call PSOPT to solve the problem   //////////////////
////////////////////////////////////////////////////////////////////////////

    psopt(solution, problem, algorithm);

////////////////////////////////////////////////////////////////////////////
///////////  Extract relevant variables from solution structure   //////////
////////////////////////////////////////////////////////////////////////////


    DMatrix x = solution.get_states_in_phase(1);
    DMatrix u = solution.get_controls_in_phase(1);
    DMatrix t = solution.get_time_in_phase(1);


////////////////////////////////////////////////////////////////////////////
///////////  Save solution data to files if desired ////////////////////////
////////////////////////////////////////////////////////////////////////////

    x.Save("x.dat");
    u.Save("u.dat");
    t.Save("t.dat");


////////////////////////////////////////////////////////////////////////////
///////////  Plot some results if desired (requires gnuplot) ///////////////
////////////////////////////////////////////////////////////////////////////

    plot(t,x(1,colon()),problem.name+": state", "time (s)", "state","x1");

    plot(t,x(2,colon()),problem.name+": state", "time (s)", "state","x2");

    plot(t,x(3,colon()),problem.name+": state", "time (s)", "state","x3");

    plot(t,x(4,colon()),problem.name+": state", "time (s)", "state","x4");

    plot(t,u(1,colon()),problem.name+": control","time (s)", "control", "u1");

    plot(t,u(2,colon()),problem.name+": control","time (s)", "control", "u2");

    plot(t,x(1,colon()),problem.name+": state x1", "time (s)", "state","x1",
                                        "pdf", "alpine_state1.pdf");

    plot(t,x(2,colon()),problem.name+": state x2", "time (s)", "state","x2",
                                        "pdf", "alpine_state2.pdf");

    plot(t,x(3,colon()),problem.name+": state x3", "time (s)", "state","x2",
                                        "pdf", "alpine_state3.pdf");

    plot(t,x(4,colon()),problem.name+": state x4", "time (s)", "state","x2",
                                        "pdf", "alpine_state4.pdf");

    plot(t,u(1,colon()),problem.name+": control u1","time (s)", "control", "u1",
                                        "pdf", "alpine_control1.pdf");

    plot(t,u(2,colon()),problem.name+": control u1","time (s)", "control", "u2",
                                        "pdf", "alpine_control2.pdf");

}

////////////////////////////////////////////////////////////////////////////
///////////////////////      END OF FILE     ///////////////////////////////
////////////////////////////////////////////////////////////////////////////
