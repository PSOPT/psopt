//////////////////////////////////////////////////////////////////////////
////////////////              bioreactor.cxx         /////////////////////
//////////////////////////////////////////////////////////////////////////
////////////////           PSOPT Template            /////////////////////
//////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////
//////// Title: 	Lee-Ramirez bioreactor            ////////////////
//////// Last modified: 09 July 2009                      ////////////////
//////// Reference:     PROPT User Manual                 ////////////////
////////                                                  ////////////////
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
                      adouble* xad, int iphase, Workspace* workspace)
{
   adouble x1tf = final_states[CINDEX(1)];
   adouble x4tf = final_states[CINDEX(4)];

   return -(x1tf*x4tf);
}

//////////////////////////////////////////////////////////////////////////
///////////////////  Define the integrand (Lagrange) cost function  //////
//////////////////////////////////////////////////////////////////////////

adouble integrand_cost(adouble* states, adouble* controls, adouble* parameters,
                     adouble& time, adouble* xad, int iphase, Workspace* workspace)
{
     adouble u2 = controls[ CINDEX(2) ];

     adouble dot_u1, dot_u2;

     double Q = 0;

     get_control_derivative( &dot_u1, 1, iphase, time, xad, workspace);

     get_control_derivative( &dot_u2, 2, iphase, time, xad, workspace);

     double rho = 0.1/get_number_of_nodes(*workspace->problem, iphase);

     return   (Q*u2 + rho*(dot_u1*dot_u1 + dot_u2*dot_u2));

}


//////////////////////////////////////////////////////////////////////////
///////////////////  Define the DAE's ////////////////////////////////////
//////////////////////////////////////////////////////////////////////////


void dae(adouble* derivatives, adouble* path, adouble* states,
         adouble* controls, adouble* parameters, adouble& time,
         adouble* xad, int iphase, Workspace* workspace)
{

    adouble u1 = controls[ CINDEX(1) ];
    adouble u2 = controls[ CINDEX(2) ];

    adouble x1 = states[ CINDEX(1) ];
    adouble x2 = states[ CINDEX(2) ];
    adouble x3 = states[ CINDEX(3) ];
    adouble x4 = states[ CINDEX(4) ];
    adouble x5 = states[ CINDEX(5) ];
    adouble x6 = states[ CINDEX(6) ];
    adouble x7 = states[ CINDEX(7) ];


    adouble k1 = 0.09*x5/(0.034 + x5);

    adouble k2 = k1;

    adouble g1  = (x3/(14.35 + x3*(1.0+x3/111.5)))*(x6 + 0.22*x7/(0.22+x5));

    adouble Rfp = (0.233*x3/(14.35 + x3*(1.0+x3/111.5)))*((0.0005+x5)/(0.022+x5));

    derivatives[ CINDEX(1) ] =   u1 + u2;
    derivatives[ CINDEX(2) ] =  g1*x2 - (u1+u2)/x1*x2;
    derivatives[ CINDEX(3) ] =  100*u1/x1 - (u1+u2)/x1*x3 - g1/0.51*x2;
    derivatives[ CINDEX(4) ] =  Rfp*x2 - (u1+u2)/x1*x4;
    derivatives[ CINDEX(5) ] =  4*u2/x1 - (u1+u2)/x1*x5;
    derivatives[ CINDEX(6) ] =  -k1*x6;
    derivatives[ CINDEX(7) ] =  k2*(1-x7);


}

////////////////////////////////////////////////////////////////////////////
///////////////////  Define the events function ////////////////////////////
////////////////////////////////////////////////////////////////////////////

void events(adouble* e, adouble* initial_states, adouble* final_states,
            adouble* parameters,adouble& t0, adouble& tf, adouble* xad,
            int iphase, Workspace* workspace)

{

   int i;

   for(i=1;i<= 7; i++) {
          e[CINDEX(i)] = initial_states[CINDEX(i)];
   }



}



///////////////////////////////////////////////////////////////////////////
///////////////////  Define the phase linkages function ///////////////////
///////////////////////////////////////////////////////////////////////////

void linkages( adouble* linkages, adouble* xad, Workspace* workspace)
{
// Single phase problem

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

    problem.name        		= "Lee-Ramirez bioreactor";
    problem.outfilename                 = "bioreactor.txt";

////////////////////////////////////////////////////////////////////////////
////////////  Define problem level constants & do level 1 setup ////////////
////////////////////////////////////////////////////////////////////////////

    problem.nphases   			= 1;
    problem.nlinkages                   = 0;

    psopt_level1_setup(problem);

/////////////////////////////////////////////////////////////////////////////
/////////   Define phase related information & do level 2 setup  ////////////
/////////////////////////////////////////////////////////////////////////////

    problem.phases(1).nstates   		= 7;
    problem.phases(1).ncontrols 		= 2;
    problem.phases(1).nevents   		= 7;
    problem.phases(1).npath     		= 0;
    problem.phases(1).nodes                     = "[20 35 50]";


    psopt_level2_setup(problem, algorithm);



////////////////////////////////////////////////////////////////////////////
///////////////////  Enter problem bounds information //////////////////////
////////////////////////////////////////////////////////////////////////////

    int i;


    problem.phases(1).bounds.lower.states = -10.0*ones(7,1);
    problem.phases(1).bounds.upper.states = "[10.0 10.0 50.0 10.0 10.0 10.0 10.0]";


    problem.phases(1).bounds.lower.controls(1) = 0.0;
    problem.phases(1).bounds.lower.controls(2) = 0.0;



    problem.phases(1).bounds.upper.controls(1) = 1.0;
    problem.phases(1).bounds.upper.controls(2) = 1.0;




    problem.phases(1).bounds.lower.events = "[1.0 0.1 40.0 0.0 0.0 1.0 0.0]";
    problem.phases(1).bounds.upper.events = "[1.0 0.1 40.0 0.0 0.0 1.0 0.0]";


    problem.phases(1).bounds.lower.StartTime    = 0.0;
    problem.phases(1).bounds.upper.StartTime    = 0.0;

    problem.phases(1).bounds.lower.EndTime      = 10.0;
    problem.phases(1).bounds.upper.EndTime      = 10.0;


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

    int nnodes    			= problem.phases(1).nodes(1);
    int ncontrols                       = problem.phases(1).ncontrols;
    int nstates                         = problem.phases(1).nstates;

    DMatrix state_guess    =  zeros(nstates,nnodes);
    DMatrix control_guess  =  zeros(ncontrols,nnodes);
    DMatrix param_guess    =  zeros(0,0);
    DMatrix time_guess     =  linspace(0,10.0,nnodes);


    state_guess(1,colon()) = 1.0*ones(1,nnodes);
    state_guess(2,colon()) = 0.1*ones(1,nnodes);
    state_guess(3,colon()) = 40.0*ones(1,nnodes);
    state_guess(4,colon()) = 0.0*ones(1,nnodes);
    state_guess(5,colon()) = 0.0*ones(1,nnodes);
    state_guess(6,colon()) = 1.0*ones(1,nnodes);
    state_guess(7,colon()) = 0.0*ones(1,nnodes);



    control_guess(1,colon()) = time_guess/30.0;
    control_guess(2,colon()) = zeros(1,nnodes);


    auto_phase_guess(problem, control_guess, state_guess, param_guess, time_guess);


////////////////////////////////////////////////////////////////////////////
///////////////////  Enter algorithm options  //////////////////////////////
////////////////////////////////////////////////////////////////////////////


    algorithm.nlp_iter_max                = 1000;
    algorithm.nlp_tolerance               = 1.e-5;
    algorithm.nlp_method                  = "IPOPT";
    algorithm.scaling                     = "automatic";
    algorithm.derivatives                 = "automatic";
    algorithm.defect_scaling              = "jacobian-based";
    algorithm.diff_matrix                 = "central-differences";


////////////////////////////////////////////////////////////////////////////
///////////////////  Now call PSOPT to solve the problem   //////////////////
////////////////////////////////////////////////////////////////////////////

    psopt(solution, problem, algorithm);

////////////////////////////////////////////////////////////////////////////
///////////  Extract relevant variables from solution structure   //////////
////////////////////////////////////////////////////////////////////////////


    DMatrix x, u, t;

    x      = solution.get_states_in_phase(1);
    u      = solution.get_controls_in_phase(1);
    t      = solution.get_time_in_phase(1);


////////////////////////////////////////////////////////////////////////////
///////////  Save solution data to files if desired ////////////////////////
////////////////////////////////////////////////////////////////////////////

    x.Save("x.dat");
    u.Save("u.dat");
    t.Save("t.dat");


////////////////////////////////////////////////////////////////////////////
///////////  Plot some results if desired (requires gnuplot) ///////////////
////////////////////////////////////////////////////////////////////////////


    plot(t,u,problem.name+": control","time (s)", "control 1", "u1 u2");

    x(3,colon()) = x(3,colon())/10.0;

    plot(t,x,problem.name+": state","time (s)", "state", "x1 x2 x3/10 x4 x5 x6 x7");


    plot(t,u,problem.name+": control","time (s)", "control 1", "u1 u2",
	      "pdf", "bioreactor_controls.pdf");

    x(3,colon()) = x(3,colon())/10.0;

    plot(t,x,problem.name+": states","time (s)", "states", "x1 x2 x3/10 x4 x5 x6 x7",
	      "pdf", "bioreactor_states.pdf");




}

////////////////////////////////////////////////////////////////////////////
///////////////////////      END OF FILE     ///////////////////////////////
////////////////////////////////////////////////////////////////////////////
