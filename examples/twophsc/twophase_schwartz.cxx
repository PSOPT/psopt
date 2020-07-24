//////////////////////////////////////////////////////////////////////////
//////////////////      twophase_schwartz.cxx        /////////////////////
//////////////////////////////////////////////////////////////////////////
////////////////           PSOPT  Example             ////////////////////
//////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////
//////// Title:       Two phase Schwartz problem          ////////////////
//////// Last modified:         09 January 2009           ////////////////
//////// Reference:             PROPT Users Guide         ////////////////
//////// (See PSOPT handbook for full reference)           ///////////////
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
   adouble retval;
   adouble x1f, x2f;

   // Phase 1 cost


   if ( iphase ==  1)
   {
   	retval  = 0.0;

   }

   // Phase 2 cost



   if( iphase == 2 ) {

   	x1f = final_states[0];
   	x2f = final_states[1];

   	retval  = 5*(x1f*x1f + x2f*x2f);
   }

   return retval;

}


//////////////////////////////////////////////////////////////////////////
///////////////////  Define the integrand (Lagrange) cost function  //////
//////////////////////////////////////////////////////////////////////////

adouble integrand_cost(adouble* states, adouble* controls,
                       adouble* parameters, adouble& time, adouble* xad,
                       int iphase, Workspace* workspace)
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

   adouble x1 = states[0];
   adouble x2 = states[1];

   adouble u = controls[0];

   derivatives[0] = x2;
   derivatives[1] = u-0.1*(1+2*x1*x1)*x2;

    if(iphase==1) {

         path[0] =  1.0 - 9.0*((x1-1.0)*(x1-1.0))-((x2-0.4)*(x2-0.4))/(0.3*0.3);

    }

}

////////////////////////////////////////////////////////////////////////////
///////////////////  Define the events function ////////////////////////////
////////////////////////////////////////////////////////////////////////////

void events(adouble* e, adouble* initial_states, adouble* final_states,
            adouble* parameters,adouble& t0, adouble& tf, adouble* xad,
            int iphase, Workspace* workspace)
{
   adouble x1i = initial_states[0];
   adouble x2i = initial_states[1];


   if (iphase==1) {
   	e[0] = x1i;
   	e[1] = x2i;
   }


}


///////////////////////////////////////////////////////////////////////////
///////////////////  Define the phase linkages function ///////////////////
///////////////////////////////////////////////////////////////////////////

void linkages( adouble* linkages, adouble* xad, Workspace* workspace)
{
    int index = 0;
    auto_link( linkages, &index, xad, 1, 2, workspace );

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

    problem.name        		= "Two phase Schwartz problem";
    problem.outfilename                 = "twophsc.txt";

////////////////////////////////////////////////////////////////////////////
////////////  Define problem level constants & do level 1 setup ////////////
////////////////////////////////////////////////////////////////////////////

    problem.nphases   			= 2;
    problem.nlinkages                   = 3;

    psopt_level1_setup(problem);

/////////////////////////////////////////////////////////////////////////////
/////////   Define phase related information &  do level 2 setup ////////////
/////////////////////////////////////////////////////////////////////////////

    problem.phases(1).nstates   = 2;
    problem.phases(1).ncontrols = 1;
    problem.phases(1).nevents   = 2;
    problem.phases(1).npath     = 1;

    problem.phases(2).nstates   = 2;
    problem.phases(2).ncontrols = 1;
    problem.phases(2).nevents   = 0;
    problem.phases(2).npath     = 0;

    problem.phases(1).nodes     = 40;
    problem.phases(2).nodes     = 40;

    psopt_level2_setup(problem, algorithm);


////////////////////////////////////////////////////////////////////////////
///////////////////  Enter problem bounds information //////////////////////
////////////////////////////////////////////////////////////////////////////

    double x1L = -20.0;
    double x2L_phase1 = -0.8;
    double x2L_phase2 = -10.0;
    double uL  = -1.0;
    double x1U = 10.0;
    double x2U = 10.0;
    double uU  =  1.0;
    double hL  = -100.0;
    double hU  = 0.0;

    // Phase 0 bounds

    problem.phases(1).bounds.lower.states(1) = x1L;
    problem.phases(1).bounds.lower.states(2) = x2L_phase1;

    problem.phases(1).bounds.upper.states(1) = x1U;
    problem.phases(1).bounds.upper.states(2) = x2U;

    problem.phases(1).bounds.lower.controls(1) = uL;
    problem.phases(1).bounds.upper.controls(1) = uU;

    problem.phases(1).bounds.lower.events(1) = 1.0;
    problem.phases(1).bounds.lower.events(2) = 1.0;


    problem.phases(1).bounds.upper.events(1) = 1.0;
    problem.phases(1).bounds.upper.events(2) = 1.0;

    problem.phases(1).bounds.lower.path(1) = hL;
    problem.phases(1).bounds.upper.path(1) = hU;


    problem.phases(1).bounds.lower.StartTime    = 0.0;
    problem.phases(1).bounds.upper.StartTime    = 0.0;

    problem.phases(1).bounds.lower.EndTime      = 1.0;
    problem.phases(1).bounds.upper.EndTime      = 1.0;

    // Phase 1 bounds

    problem.phases(2).bounds.lower.states(1) = x1L;
    problem.phases(2).bounds.lower.states(2) = x2L_phase2;

    problem.phases(2).bounds.upper.states(1) = x1U;
    problem.phases(2).bounds.upper.states(2) = x2U;

    problem.phases(2).bounds.lower.controls(1) = -50.0;
    problem.phases(2).bounds.upper.controls(1) =  50.0;

    problem.phases(2).bounds.lower.StartTime    = 1.0;
    problem.phases(2).bounds.upper.StartTime    = 1.0;

    problem.phases(2).bounds.lower.EndTime      = 2.9;
    problem.phases(2).bounds.upper.EndTime      = 2.9;


////////////////////////////////////////////////////////////////////////////
///////////////////  Register problem functions  ///////////////////////////
////////////////////////////////////////////////////////////////////////////

    problem.integrand_cost 	= &integrand_cost;
    problem.endpoint_cost 	= &endpoint_cost;
    problem.dae		 	= &dae;
    problem.events 		= &events;
    problem.linkages		= &linkages;

////////////////////////////////////////////////////////////////////////////
///////////////////  Define & register initial guess ///////////////////////
////////////////////////////////////////////////////////////////////////////

    int iphase;


    DMatrix u0(1,40);
    DMatrix x0(2,40);


    DMatrix time_guess0    = linspace(0.0, 1.0 , 40);
    DMatrix time_guess1    = linspace(1.0, 2.9 , 40);


    iphase = 1;


    u0 = zeros(1,40);
    x0(1,colon()) = linspace(1.0,1.0, 40);
    x0(2,colon()) = linspace(1.0,1.0, 40);

    problem.phases(iphase).guess.controls = u0;
    problem.phases(iphase).guess.states   = x0;
    problem.phases(iphase).guess.time     = time_guess0;

    iphase = 2;

    u0 = zeros(1,40);
    x0(1,colon()) = linspace(1.0,1.0, 40);
    x0(2,colon()) = linspace(1.0,1.0, 40);

    problem.phases(iphase).guess.controls = u0;
    problem.phases(iphase).guess.states   = x0;
    problem.phases(iphase).guess.time     = time_guess1;

////////////////////////////////////////////////////////////////////////////
///////////////////  Enter algorithm options  //////////////////////////////
////////////////////////////////////////////////////////////////////////////

    algorithm.nlp_method                  = "IPOPT";
    algorithm.scaling                     = "automatic";
    algorithm.derivatives                 = "automatic";
    algorithm.nlp_iter_max                = 1000;
    algorithm.nlp_tolerance               = 1.e-6;



////////////////////////////////////////////////////////////////////////////
///////////////////  Now call PSOPT to solve the problem   /////////////////
////////////////////////////////////////////////////////////////////////////

    psopt(solution, problem, algorithm);

////////////////////////////////////////////////////////////////////////////
///////////  Extract relevant variables from solution structure   //////////
////////////////////////////////////////////////////////////////////////////

    DMatrix xphase1 = solution.get_states_in_phase(1);
    DMatrix uphase1 = solution.get_controls_in_phase(1);
    DMatrix tphase1 = solution.get_time_in_phase(1);

    DMatrix xphase2 = solution.get_states_in_phase(2);
    DMatrix uphase2 = solution.get_controls_in_phase(2);
    DMatrix tphase2 = solution.get_time_in_phase(2);

    DMatrix x = (xphase1 || xphase2);
    DMatrix u = (uphase1 || uphase2);
    DMatrix t = (tphase1 || tphase2);

////////////////////////////////////////////////////////////////////////////
///////////  Save solution data to files if desired ////////////////////////
////////////////////////////////////////////////////////////////////////////


    x.Save("x.dat");
    u.Save("u.dat");
    t.Save("t.dat");

////////////////////////////////////////////////////////////////////////////
///////////  Plot some results if desired (requires gnuplot) ///////////////
////////////////////////////////////////////////////////////////////////////

    plot(t,x,problem.name, "time (s)", "states", "x1 x2");

    plot(t,u,problem.name, "time (s)", "control", "u");

    plot(t,x,problem.name, "time (s)", "states", "x1 x2",
                           "pdf", "twophsc_states.pdf");

    plot(t,u,problem.name, "time (s)", "control", "u",
                            "pdf", "twophsc_control.pdf");


}

////////////////////////////////////////////////////////////////////////////
///////////////////////      END OF FILE     ///////////////////////////////
////////////////////////////////////////////////////////////////////////////

