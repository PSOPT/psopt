//////////////////////////////////////////////////////////////////////////
////////////////           steps.cxx                 /////////////////////
//////////////////////////////////////////////////////////////////////////
////////////////           PSOPT example             /////////////////////
//////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////
//////// Title: 	Steps problem                     ////////////////
//////// Last modified: 12 July 2009                      ////////////////
//////// Reference:      Gong, Farhoo, and Ross (2008)    ////////////////
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
   	return 0.0;
}

//////////////////////////////////////////////////////////////////////////
///////////////////  Define the integrand (Lagrange) cost function  //////
//////////////////////////////////////////////////////////////////////////

adouble integrand_cost(adouble* states, adouble* controls, adouble* parameters,
                     adouble& time, adouble* xad, int iphase, Workspace* workspace)
{
   adouble x = states[ 0 ];

   return (x);
}


//////////////////////////////////////////////////////////////////////////
///////////////////  Define the DAE's ////////////////////////////////////
//////////////////////////////////////////////////////////////////////////


void dae(adouble* derivatives, adouble* path, adouble* states,
         adouble* controls, adouble* parameters, adouble& time,
         adouble* xad, int iphase, Workspace* workspace)
{

    adouble u = controls[CINDEX(1)];

    derivatives[ 0 ]  = u;

}

////////////////////////////////////////////////////////////////////////////
///////////////////  Define the events function ////////////////////////////
////////////////////////////////////////////////////////////////////////////

void events(adouble* e, adouble* initial_states, adouble* final_states,
            adouble* parameters,adouble& t0, adouble& tf, adouble* xad,
            int iphase, Workspace* workspace)

{

    adouble x1_i 	= initial_states[ 0 ];

    adouble x1_f 	= final_states[   0 ];


   if ( iphase==1 ) {
		    e[ 0 ] =  x1_i;
   }
   else if ( iphase==3 ) {
		    e[ 0 ] =  x1_f;
   }

}



///////////////////////////////////////////////////////////////////////////
///////////////////  Define the phase linkages function ///////////////////
///////////////////////////////////////////////////////////////////////////

void linkages( adouble* linkages, adouble* xad, Workspace* workspace)
{

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
    MSdata msdata;

////////////////////////////////////////////////////////////////////////////
///////////////////  Register problem name  ////////////////////////////////
////////////////////////////////////////////////////////////////////////////

    problem.name        		          = "Steps problem";
    problem.outfilename                 = "steps.txt";

////////////////////////////////////////////////////////////////////////////
////////////  Define problem level constants & do level 1 setup ////////////
////////////////////////////////////////////////////////////////////////////

    msdata.nsegments       = 3;
    msdata.nstates         = 1;
    msdata.ncontrols       = 1;
    msdata.nparameters     = 0;
    msdata.npath           = 0;
    msdata.ninitial_events = 1;
    msdata.nfinal_events   = 1;
    msdata.nodes           << 20; // nodes per segment

    multi_segment_setup(problem, algorithm, msdata );

////////////////////////////////////////////////////////////////////////////
///////////////////  Enter problem bounds information //////////////////////
////////////////////////////////////////////////////////////////////////////


    problem.phases(1).bounds.lower.controls(0) = -1.0;
    problem.phases(1).bounds.upper.controls(0) =  1.0;
    problem.phases(1).bounds.lower.states(0) = 0.0;
    problem.phases(1).bounds.upper.states(0) = 5.0;
    problem.phases(1).bounds.lower.events(0) = 1.0;
    problem.phases(3).bounds.lower.events(0) = 1.0;


    problem.phases(1).bounds.upper.events=problem.phases(1).bounds.lower.events;
    problem.phases(3).bounds.upper.events=problem.phases(3).bounds.lower.events;



    problem.phases(1).bounds.lower.StartTime    = 0.0;
    problem.phases(1).bounds.upper.StartTime    = 0.0;

    problem.phases(3).bounds.lower.EndTime      = 3.0;
    problem.phases(3).bounds.upper.EndTime      = 3.0;

//    problem.bounds.lower.times = "[0.0, 1.0, 2.0, 3.0]";
//    problem.bounds.upper.times = "[0.0, 1.0, 2.0, 3.0]";

problem.bounds.lower.times.resize(1,4);
problem.bounds.upper.times.resize(1,4);

    problem.bounds.lower.times << 0.0, 1.0, 2.0, 3.0;
    problem.bounds.upper.times << 0.0, 1.0, 2.0, 3.0;

    auto_phase_bounds(problem);



////////////////////////////////////////////////////////////////////////////
///////////////////  Register problem functions  ///////////////////////////
////////////////////////////////////////////////////////////////////////////


    problem.integrand_cost 		= &integrand_cost;
    problem.endpoint_cost 			= &endpoint_cost;
    problem.dae             		= &dae;
    problem.events 					= &events;
    problem.linkages					= &linkages;

////////////////////////////////////////////////////////////////////////////
///////////////////  Define & register initial guess ///////////////////////
////////////////////////////////////////////////////////////////////////////

    int nnodes    							 = problem.phases(1).nodes(0);
    int ncontrols                       = problem.phases(1).ncontrols;
    int nstates                         = problem.phases(1).nstates;

    MatrixXd state_guess    				=  zeros(nstates,nnodes);
    MatrixXd control_guess  				=  zeros(ncontrols,nnodes);
    MatrixXd time_guess    				=  linspace(0.0,3.0,nnodes);
    MatrixXd param_guess;



    state_guess   = linspace(1.0, 1.0, nnodes);

    control_guess = zeros(1,nnodes);


    auto_phase_guess(problem, control_guess, state_guess, param_guess, time_guess);

////////////////////////////////////////////////////////////////////////////
///////////////////  Enter algorithm options  //////////////////////////////
////////////////////////////////////////////////////////////////////////////


    algorithm.nlp_iter_max                = 1000;
    algorithm.nlp_tolerance               = 1.e-6;
    algorithm.nlp_method                  = "IPOPT";
    algorithm.scaling                     = "automatic";
    algorithm.derivatives                 = "automatic";
    algorithm.hessian                     = "exact";
    algorithm.mesh_refinement             = "automatic";
    algorithm.ode_tolerance               = 1.e-5;


////////////////////////////////////////////////////////////////////////////
///////////////////  Now call PSOPT to solve the problem   //////////////////
////////////////////////////////////////////////////////////////////////////

    psopt(solution, problem, algorithm);

////////////////////////////////////////////////////////////////////////////
///////////  Extract relevant variables from solution structure   //////////
////////////////////////////////////////////////////////////////////////////


    MatrixXd x, u, t, x_ph1, u_ph1, t_ph1, x_ph2, u_ph2, t_ph2, x_ph3, u_ph3, t_ph3; 

    x_ph1      = solution.get_states_in_phase(1);
    u_ph1      = solution.get_controls_in_phase(1);
    t_ph1      = solution.get_time_in_phase(1);
    
    x_ph2      = solution.get_states_in_phase(2);
    u_ph2      = solution.get_controls_in_phase(2);
    t_ph2      = solution.get_time_in_phase(2);
    
    x_ph3      = solution.get_states_in_phase(3);
    u_ph3      = solution.get_controls_in_phase(3);
    t_ph3      = solution.get_time_in_phase(3);


    x.resize(1, length(t_ph1)+length(t_ph2)+length(t_ph3) );
    u.resize(1, length(t_ph1)+length(t_ph2)+length(t_ph3) ); 
    t.resize(1, length(t_ph1)+length(t_ph2)+length(t_ph3) );  
    
    x << x_ph2, x_ph2, x_ph3;
    u << u_ph1, u_ph2, u_ph3;
    t << t_ph1, t_ph2, t_ph3;


////////////////////////////////////////////////////////////////////////////
///////////  Save solution data to files if desired ////////////////////////
////////////////////////////////////////////////////////////////////////////

    Save(x,"x.dat");
    Save(u,"u.dat");
    Save(t,"t.dat");

////////////////////////////////////////////////////////////////////////////
///////////  Plot some results if desired (requires gnuplot) ///////////////
////////////////////////////////////////////////////////////////////////////

    plot(t,x,problem.name+": state","time (s)", "x", "x");

    plot(t,u,problem.name+": control","time (s)", "u", "u");

    plot(t,x,problem.name+": state","time (s)", "x", "x",
	                                "pdf", "steps_state.pdf");


    plot(t,u,problem.name+": control","time (s)", "u", "u",
	                                "pdf", "steps_control.pdf");


}

////////////////////////////////////////////////////////////////////////////
///////////////////////      END OF FILE     ///////////////////////////////
////////////////////////////////////////////////////////////////////////////
