//////////////////////////////////////////////////////////////////////////
//////////////////      twophase_robot.cxx        ////////////////////////
//////////////////////////////////////////////////////////////////////////
////////////////           PSOPT  Example             ////////////////////
//////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////
//////// Title:    Two phase path tracking robot          ////////////////
//////// Last modified:         09 January 2009           ////////////////
//////// Reference:             PROPT Users Guide         ////////////////
//////// (See PSOPT handbook for full reference)          ////////////////
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
   return 0.0;
}

//////////////////////////////////////////////////////////////////////////
///////////////////  Define the integrand (Lagrange) cost function  //////
//////////////////////////////////////////////////////////////////////////

adouble integrand_cost(adouble* states, adouble* controls,
                       adouble* parameters, adouble& time, adouble* xad,
                       int iphase, Workspace* workspace)
{
    adouble retval;

    adouble x1, x2, x3, x4;

    adouble x1ref, x2ref, x3ref, x4ref;

    double w1, w2, w3, w4;

    w1 = 100.0;
    w2 = 100.0;
    w3 = 500.0;
    w4 = 500.0;


   x1 = states[0];
   x2 = states[1];
   x3 = states[2];
   x4 = states[3];

   if (iphase==1) {
   	x1ref = time/2;

   	x2ref = 0.0;

   	x3ref = 0.5;

   	x4ref = 0.0;

   }

   if (iphase==2) {

   	x1ref =  0.5;

   	x2ref =  (time-1.0)/2.0;

   	x3ref =  0.0;

   	x4ref =  0.5;
   }

   retval = w1*pow(x1-x1ref,2)+w2*pow(x2-x2ref,2)+w3*pow(x3-x3ref,2)+w4*pow(x4-x4ref,2);

   return  retval;
}


//////////////////////////////////////////////////////////////////////////
///////////////////  Define the DAE's ////////////////////////////////////
//////////////////////////////////////////////////////////////////////////

void dae(adouble* derivatives, adouble* path, adouble* states,
         adouble* controls, adouble* parameters, adouble& time,
         adouble* xad, int iphase, Workspace* workspace)
{


   adouble  x3, x4, u1, u2;

   x3 = states[2];
   x4 = states[3];

   u1 = controls[0];
   u2 = controls[1];

   derivatives[0] = x3;
   derivatives[1] = x4;
   derivatives[2] = u1;
   derivatives[3] = u2;

}

////////////////////////////////////////////////////////////////////////////
///////////////////  Define the events function ////////////////////////////
////////////////////////////////////////////////////////////////////////////

void events(adouble* e, adouble* initial_states, adouble* final_states,
            adouble* parameters,adouble& t0, adouble& tf, adouble* xad,
            int iphase, Workspace* workspace)
{
   adouble x1i, x2i, x3i, x4i;
   adouble x1f, x2f, x3f, x4f;

   if (iphase == 1)
   {

      x1i = initial_states[0];
      x2i = initial_states[1];
      x3i = initial_states[2];
      x4i = initial_states[3];

      e[0] = x1i;
      e[1] = x2i;
      e[2] = x3i;
      e[3] = x4i;

   }

   else if (iphase == 2)
   {

      x1f = final_states[0];
      x2f = final_states[1];
      x3f = final_states[2];
      x4f = final_states[3];

      e[0] = x1f;
      e[1] = x2f;
      e[2] = x3f;
      e[3] = x4f;

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

    problem.name        		= "Two phase path tracking robot";
    problem.outfilename                 = "twophro.txt";

////////////////////////////////////////////////////////////////////////////
////////////  Define problem level constants & do level 1 setup ////////////
////////////////////////////////////////////////////////////////////////////

    problem.nphases   			      = 2;
    problem.nlinkages               = 5;

    psopt_level1_setup(problem);

/////////////////////////////////////////////////////////////////////////////
/////////   Define phase related information &  do level 2 setup  ///////////
/////////////////////////////////////////////////////////////////////////////

    problem.phases(1).nstates   		= 4;
    problem.phases(1).ncontrols 		= 2;
    problem.phases(1).nevents   		= 4;
    problem.phases(1).npath     		= 0;

    problem.phases(2).nstates   		= 4;
    problem.phases(2).ncontrols 		= 2;
    problem.phases(2).nevents   		= 4;
    problem.phases(2).npath     		= 0;

    problem.phases(1).nodes        	<< 30;
    problem.phases(2).nodes        	<< 30;

    psopt_level2_setup(problem, algorithm);


////////////////////////////////////////////////////////////////////////////
///////////////////  Enter problem bounds information //////////////////////
////////////////////////////////////////////////////////////////////////////

    double x1i = 0.0;
    double x2i = 0.0;
    double x3i = 0.5;
    double x4i = 0.0;

    double x1f = 0.5;
    double x2f = 0.5;
    double x3f = 0.0;
    double x4f = 0.5;

    // Phase 1 bounds

    problem.phases(1).bounds.lower.states(0) = -10.0;
    problem.phases(1).bounds.lower.states(1) = -10.0;
    problem.phases(1).bounds.lower.states(2) = -10.0;
    problem.phases(1).bounds.lower.states(3) = -10.0;

    problem.phases(1).bounds.upper.states(0) = 10.0;
    problem.phases(1).bounds.upper.states(1) = 10.0;
    problem.phases(1).bounds.upper.states(2) = 10.0;
    problem.phases(1).bounds.upper.states(3) = 10.0;


    problem.phases(1).bounds.lower.controls(0) = -10.0;
    problem.phases(1).bounds.upper.controls(0) =  10.0;
    problem.phases(1).bounds.lower.controls(1) = -10.0;
    problem.phases(1).bounds.upper.controls(1) =  10.0;

    problem.phases(1).bounds.lower.events(0) = x1i;
    problem.phases(1).bounds.lower.events(1) = x2i;
    problem.phases(1).bounds.lower.events(2) = x3i;
    problem.phases(1).bounds.lower.events(3) = x4i;


    problem.phases(1).bounds.upper.events(0) = x1i;
    problem.phases(1).bounds.upper.events(1) = x2i;
    problem.phases(1).bounds.upper.events(2) = x3i;
    problem.phases(1).bounds.upper.events(3) = x4i;


    problem.phases(1).bounds.lower.StartTime    = 0.0;
    problem.phases(1).bounds.upper.StartTime    = 0.0;

    problem.phases(1).bounds.lower.EndTime      = 1.0;
    problem.phases(1).bounds.upper.EndTime      = 1.0;


    // Phase 2 bounds


    problem.phases(2).bounds.lower.states(0) = -10.0;
    problem.phases(2).bounds.lower.states(1) = -10.0;
    problem.phases(2).bounds.lower.states(2) = -10.0;
    problem.phases(2).bounds.lower.states(3) = -10.0;

    problem.phases(2).bounds.upper.states(0) = 10.0;
    problem.phases(2).bounds.upper.states(1) = 10.0;
    problem.phases(2).bounds.upper.states(2) = 10.0;
    problem.phases(2).bounds.upper.states(3) = 10.0;

    problem.phases(2).bounds.lower.controls(0) = -10.0;
    problem.phases(2).bounds.upper.controls(0) =  10.0;
    problem.phases(2).bounds.lower.controls(1) = -10.0;
    problem.phases(2).bounds.upper.controls(1) =  10.0;

    problem.phases(2).bounds.lower.events(0) = x1f;
    problem.phases(2).bounds.lower.events(1) = x2f;
    problem.phases(2).bounds.lower.events(2) = x3f;
    problem.phases(2).bounds.lower.events(3) = x4f;


    problem.phases(2).bounds.upper.events(0) = x1f;
    problem.phases(2).bounds.upper.events(1) = x2f;
    problem.phases(2).bounds.upper.events(2) = x3f;
    problem.phases(2).bounds.upper.events(3) = x4f;

    problem.phases(2).bounds.lower.StartTime    = 1.0;
    problem.phases(2).bounds.upper.StartTime    = 1.0;

    problem.phases(2).bounds.lower.EndTime      = 2.0;
    problem.phases(2).bounds.upper.EndTime      = 2.0;


////////////////////////////////////////////////////////////////////////////
///////////////////  Register problem functions  ///////////////////////////
////////////////////////////////////////////////////////////////////////////


    problem.integrand_cost 		= &integrand_cost;
    problem.endpoint_cost 			= &endpoint_cost;
    problem.dae 						= &dae;
    problem.events 					= &events;
    problem.linkages					= &linkages;


////////////////////////////////////////////////////////////////////////////
///////////////////  Define & register initial guess ///////////////////////
////////////////////////////////////////////////////////////////////////////
    int iphase;

    MatrixXd u0(2,30);
    MatrixXd x0(4,30);

    MatrixXd time_guess0    = linspace(0.0, 1.0 , 30);
    MatrixXd time_guess1    = linspace(1.0, 2.0 , 30);

    iphase = 1;

    u0 = zeros(2,30);

    x0 << linspace(x1i,(x1i+x1f)/2, 30),
          linspace(x2i,(x2i+x2f)/2, 30),
          linspace(x3i,(x3i+x3f)/2, 30),
          linspace(x4i,(x4i+x4f)/2, 30);

    problem.phases(iphase).guess.controls = u0;
    problem.phases(iphase).guess.states   = x0;
    problem.phases(iphase).guess.time     = time_guess0;

    iphase = 2;

    u0 = zeros(2,30);

    x0  <<  linspace((x1i+x1f)/2, x1f, 30),
            linspace((x2i+x2f)/2, x2f, 30),
            linspace((x3i+x3f)/2, x3f, 30),
            linspace((x4i+x4f)/2, x4f, 30);

    problem.phases(iphase).guess.controls = u0;
    problem.phases(iphase).guess.states   = x0;
    problem.phases(iphase).guess.time     = time_guess1;

////////////////////////////////////////////////////////////////////////////
///////////////////  Enter algorithm options  //////////////////////////////
////////////////////////////////////////////////////////////////////////////


    algorithm.nlp_method                  	= "IPOPT";
    algorithm.scaling                     	= "automatic";
    algorithm.derivatives                 	= "automatic";
    algorithm.hessian                        = "exact";
    algorithm.nlp_iter_max                	= 1000;
    algorithm.nlp_tolerance               	= 1.e-6;



////////////////////////////////////////////////////////////////////////////
///////////////////  Now call PSOPT to solve the problem   /////////////////
////////////////////////////////////////////////////////////////////////////

    psopt(solution, problem, algorithm);

////////////////////////////////////////////////////////////////////////////
///////////  Extract relevant variables from solution structure   //////////
////////////////////////////////////////////////////////////////////////////

    MatrixXd xphase1 = solution.get_states_in_phase(1);
    MatrixXd uphase1 = solution.get_controls_in_phase(1);
    MatrixXd tphase1 = solution.get_time_in_phase(1);

    MatrixXd xphase2 = solution.get_states_in_phase(2);
    MatrixXd uphase2 = solution.get_controls_in_phase(2);
    MatrixXd tphase2 = solution.get_time_in_phase(2);

    MatrixXd x(4, length(tphase1)+length(tphase2));
    x  << xphase1 , xphase2;
    MatrixXd u(2, length(tphase1)+length(tphase2));
    u << uphase1 , uphase2;
    MatrixXd t(1, length(tphase1)+length(tphase2));  
    t << tphase1 , tphase2;

////////////////////////////////////////////////////////////////////////////
///////////  Save solution data to files if desired ////////////////////////
////////////////////////////////////////////////////////////////////////////

    Save(x,"x.dat");
    Save(u,"u.dat");
    Save(t,"t.dat");

////////////////////////////////////////////////////////////////////////////
///////////  Plot some results if desired (requires gnuplot) ///////////////
////////////////////////////////////////////////////////////////////////////

    plot(t,x,problem.name, "time (s)", "states");

    plot(t,u,problem.name, "time (s)", "control");

    plot(t,x,problem.name+": states", "time (s)", "states", "x1 x2 x3 x4",
                          "pdf", "twophro_states.pdf");

    plot(t,u,problem.name+": controls", "time (s)", "controls", "u1 u2",
                           "pdf", "twophro_controls.pdf" );

}

////////////////////////////////////////////////////////////////////////////
///////////////////////      END OF FILE     ///////////////////////////////
////////////////////////////////////////////////////////////////////////////

