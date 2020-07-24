//////////////////////////////////////////////////////////////////////////
////////////////             reorientation.cxx        ////////////////////
//////////////////////////////////////////////////////////////////////////
////////////////           PSOPT example            /////////////////////
//////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////
//////// Title: Reorientation of an asymmetric rigid body  ///////////////
//////// Last modified: 03 July 2010                      ////////////////
//////// Reference:     Fleming et al                     ////////////////
////////                                                  ////////////////
//////////////////////////////////////////////////////////////////////////
////////     Copyright (c) Victor M. Becerra, 2010        ////////////////
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
   return tf;
}

//////////////////////////////////////////////////////////////////////////
///////////////////  Define the integrand (Lagrange) cost function  //////
//////////////////////////////////////////////////////////////////////////

adouble integrand_cost(adouble* states, adouble* controls, adouble* parameters,
                     adouble& time, adouble* xad, int iphase, Workspace* workspace)
{
   return 0.0;
}


//////////////////////////////////////////////////////////////////////////
///////////////////  Define the DAE's ////////////////////////////////////
//////////////////////////////////////////////////////////////////////////


void dae(adouble* derivatives, adouble* path, adouble* states,
         adouble* controls, adouble* parameters, adouble& time,
         adouble* xad, int iphase, Workspace* workspace)
{

    adouble u1 		= controls[ CINDEX(1) ];
    adouble u2          = controls[ CINDEX(2) ];
    adouble u3          = controls[ CINDEX(3) ];
    adouble q4          = controls[ CINDEX(4) ];

    adouble q1 		= states[ CINDEX(1) ];
    adouble q2 		= states[ CINDEX(2) ];
    adouble q3 		= states[ CINDEX(3) ];
    adouble omega1 	= states[ CINDEX(4) ];
    adouble omega2      = states[ CINDEX(5) ];
    adouble omega3      = states[ CINDEX(6) ];

    double Ix = 5621.0;
    double Iy = 4547.0;
    double Iz = 2364.0;

    adouble dq1     = 0.5*( omega1*q4 - omega2*q3 + omega3*q2 );
    adouble dq2     = 0.5*( omega1*q3 + omega2*q4 - omega3*q1 );
    adouble dq3     = 0.5*(-omega1*q2 + omega2*q1 + omega3*q4 );
    adouble domega1 = u1/Ix - ((Iz-Iy)/Ix)*omega2*omega3;
    adouble domega2 = u2/Iy - ((Ix-Iz)/Iy)*omega1*omega3;
    adouble domega3 = u3/Iz - ((Iy-Ix)/Iz)*omega1*omega2;




    derivatives[ CINDEX(1) ] =   dq1;
    derivatives[ CINDEX(2) ] =   dq2;
    derivatives[ CINDEX(3) ] =   dq3;
    derivatives[ CINDEX(4) ] =   domega1;
    derivatives[ CINDEX(5) ] =   domega2;
    derivatives[ CINDEX(6) ] =   domega3;


    path[ CINDEX(1) ] = q1*q1 + q2*q2 + q3*q3 + q4*q4 - 1.0;

}

////////////////////////////////////////////////////////////////////////////
///////////////////  Define the events function ////////////////////////////
////////////////////////////////////////////////////////////////////////////

void events(adouble* e, adouble* initial_states, adouble* final_states,
            adouble* parameters,adouble& t0, adouble& tf, adouble* xad,
            int iphase, Workspace* workspace)

{

    adouble q1i 	= initial_states[ CINDEX(1) ];
    adouble q2i 	= initial_states[ CINDEX(2) ];
    adouble q3i 	= initial_states[ CINDEX(3) ];
    adouble omega1i 	= initial_states[ CINDEX(4) ];
    adouble omega2i     = initial_states[ CINDEX(5) ];
    adouble omega3i     = initial_states[ CINDEX(6) ];

    adouble initial_controls[4], final_controls[4], q4i, q4f;

    get_initial_controls( initial_controls, xad, iphase, workspace );

    get_final_controls(   final_controls  , xad, iphase, workspace );

    q4i = initial_controls[ CINDEX(4) ];

    q4f = final_controls[   CINDEX(4) ];

    adouble q1f 	= final_states[ CINDEX(1) ];
    adouble q2f 	= final_states[ CINDEX(2) ];
    adouble q3f 	= final_states[ CINDEX(3) ];
    adouble omega1f 	= final_states[ CINDEX(4) ];
    adouble omega2f     = final_states[ CINDEX(5) ];
    adouble omega3f     = final_states[ CINDEX(6) ];


    e[ CINDEX(1)  ] 	= q1i;
    e[ CINDEX(2)  ]     = q2i;
    e[ CINDEX(3)  ]     = q3i;
    e[ CINDEX(4)  ]   	= q4i;
    e[ CINDEX(5)  ]   	= omega1i;
    e[ CINDEX(6)  ]     = omega2i;
    e[ CINDEX(7)  ]  	= omega3i;

    e[ CINDEX(8)  ] 	= q1f;
    e[ CINDEX(9)  ]     = q2f;
    e[ CINDEX(10) ]     = q3f;
    e[ CINDEX(11) ]    	= q4f;
    e[ CINDEX(12) ]   	= omega1f;
    e[ CINDEX(13) ]     = omega2f;
    e[ CINDEX(14) ]  	= omega3f;


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

    problem.name        		= "Reorientation of a rigid body";
    problem.outfilename                 = "reorientation.txt";

////////////////////////////////////////////////////////////////////////////
////////////  Define problem level constants & do level 1 setup ////////////
////////////////////////////////////////////////////////////////////////////

    problem.nphases   			= 1;
    problem.nlinkages                   = 0;

    psopt_level1_setup(problem);

/////////////////////////////////////////////////////////////////////////////
/////////   Define phase related information & do level 2 setup  ////////////
/////////////////////////////////////////////////////////////////////////////

    problem.phases(1).nstates   		= 6;
    problem.phases(1).ncontrols 		= 4;
    problem.phases(1).nevents   		= 14;
    problem.phases(1).npath     		= 1;
    problem.phases(1).nodes                     = "[60]";


    psopt_level2_setup(problem, algorithm);



////////////////////////////////////////////////////////////////////////////
///////////////////  Enter problem bounds information //////////////////////
////////////////////////////////////////////////////////////////////////////


    problem.phases(1).bounds.lower.states = "[-1.0    -1.0    -1.0   -0.5    -0.5  -0.5]";
    problem.phases(1).bounds.upper.states = "[ 1.0     1.0     1.0    0.5     0.5   0.5]";


    problem.phases(1).bounds.lower.controls = "[-50.0  -50.0 -50.0  -1.0]";

    problem.phases(1).bounds.upper.controls = "[ 50.0   50.0  50.0   1.0]";


    problem.phases(1).bounds.lower.path(1) = 0.000;

    problem.phases(1).bounds.upper.path(1) =  0.000;

    DMatrix q0, qf, omega0, omegaf;

    double phi = 150.0*(pi/180.0);

    q0(1) = 0.0; q0(2) = 0.0; q0(3)=0.0; q0(4)=1.0;
    qf(1) = sin(phi/2.0); qf(2) = 0.0; qf(3) = 0.0; qf(4) = cos(phi/2.0);

    omega0 = zeros(1,3); omegaf = zeros(1,3);

    problem.phases(1).bounds.lower.events = q0 || omega0 || qf || omegaf;
    problem.phases(1).bounds.upper.events = problem.phases(1).bounds.lower.events;


    problem.phases(1).bounds.lower.StartTime    = 0.0;
    problem.phases(1).bounds.upper.StartTime    = 0.0;

    problem.phases(1).bounds.lower.EndTime      = 25.0;
    problem.phases(1).bounds.upper.EndTime      = 35.0;



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
    DMatrix control_guess  =  50.0*ones(ncontrols,nnodes);
    DMatrix time_guess     =  linspace(0.0,40,nnodes);


    state_guess(1, colon() ) = linspace( q0(1), qf(1), nnodes );
    state_guess(2, colon() ) = linspace( q0(2), qf(2), nnodes );
    state_guess(3, colon() ) = linspace( q0(3), qf(3), nnodes );
    control_guess(4, colon() )=linspace( q0(4), qf(4), nnodes );
    state_guess(4, colon() ) = linspace( omega0(1), omegaf(1), nnodes );
    state_guess(5, colon() ) = linspace( omega0(2), omegaf(2), nnodes );
    state_guess(6, colon() ) = linspace( omega0(3), omegaf(3), nnodes );

    problem.phases(1).guess.states   = state_guess;
    problem.phases(1).guess.controls = control_guess;
    problem.phases(1).guess.time     = time_guess;



////////////////////////////////////////////////////////////////////////////
///////////////////  Enter algorithm options  //////////////////////////////
////////////////////////////////////////////////////////////////////////////


    algorithm.nlp_iter_max                = 1000;
    algorithm.nlp_tolerance               = 1.e-6;
    algorithm.nlp_method                  = "IPOPT";
    algorithm.scaling                     = "automatic";
    algorithm.derivatives                 = "automatic";
    algorithm.collocation_method          = "trapezoidal";
    algorithm.mesh_refinement             = "automatic";
    algorithm.ode_tolerance               = 1.e-5;


////////////////////////////////////////////////////////////////////////////
///////////////////  Now call PSOPT to solve the problem   //////////////////
////////////////////////////////////////////////////////////////////////////

    psopt(solution, problem, algorithm);

////////////////////////////////////////////////////////////////////////////
///////////  Extract relevant variables from solution structure   //////////
////////////////////////////////////////////////////////////////////////////


    DMatrix states, controls, t, q1, q2, q3, q4, omega, u, q;

    states      = solution.get_states_in_phase(1);
    controls    = solution.get_controls_in_phase(1);
    t           = solution.get_time_in_phase(1);

    q1 = states(1,colon());
    q2 = states(2,colon());
    q3 = states(3,colon());
    q4 = controls(4, colon());

    q = q1 && q2 && q3 && q4;

    omega = states(colon(4,6), colon());

    u = controls(colon(1,3), colon());




////////////////////////////////////////////////////////////////////////////
///////////  Save solution data to files if desired ////////////////////////
////////////////////////////////////////////////////////////////////////////

    states.Save("states.dat");
    controls.Save("controls.dat");
    t.Save("t.dat");


////////////////////////////////////////////////////////////////////////////
///////////  Plot some results if desired (requires gnuplot) ///////////////
////////////////////////////////////////////////////////////////////////////



    multiplot(t,q,problem.name+": quarternion","time (s)", "q1 q2 q3 q4", "q1 q2 q3 q4", 2, 2);

    multiplot(t,u, problem.name+": controls", "time (s)", "u1 u2 u3", "u1 u2 u3");


    multiplot(t,q,problem.name+": quarternion","time (s)", "q1 q2 q3 q4", "q1 q2 q3 q4", 2, 2,
                "pdf", "reorientation_q.pdf");

    multiplot(t,u, problem.name+": controls", "time (s)", "u1 u2 u3", "u1 u2 u3", 3, 1,
               "pdf", "reorientation_u.pdf");


}

////////////////////////////////////////////////////////////////////////////
///////////////////////      END OF FILE     ///////////////////////////////
////////////////////////////////////////////////////////////////////////////
