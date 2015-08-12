//////////////////////////////////////////////////////////////////////////
////////////////        lambert.cxx               /////////////////////
//////////////////////////////////////////////////////////////////////////
////////////////           PSOPT  Example             /////////////////////
//////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
//////// Title: Lambert problem                           ////////////////
//////// Last modified: 27 December 2009                  ////////////////
//////// Reference:     Vallado (2001), page 470          ////////////////
//////// (See PSOPT handbook for full reference)           ////////////////
//////////////////////////////////////////////////////////////////////////
////////     Copyright (c) Victor M. Becerra, 2009        ////////////////
//////////////////////////////////////////////////////////////////////////
//////// This is part of the PSOPT software library, which ////////////////
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
    return  0.0;
}


//////////////////////////////////////////////////////////////////////////
///////////////////  Define the DAE's ////////////////////////////////////
//////////////////////////////////////////////////////////////////////////

void dae(adouble* derivatives, adouble* path, adouble* states,
         adouble* controls, adouble* parameters, adouble& time,
         adouble* xad, int iphase, Workspace* workspace)
{

   // Define constants:

   double G =   6.6720e-11; // Universal gravity constant [m^3/(kg s^2)]
   double Me =  5.9742e24;   // Mass of earth;[kg]

   double mu   = G*Me;       // [m^3/sec^2]

   ADMatrix r(3), v(3);

   // Extract individual variables

   r(1) = states[ CINDEX(1) ];
   r(2) = states[ CINDEX(2) ];
   r(3) = states[ CINDEX(3) ];

   v(1) = states[ CINDEX(4) ];
   v(2) = states[ CINDEX(5) ];
   v(3) = states[ CINDEX(6) ];


   ADMatrix rdd(3);

   adouble rr = sqrt( r(1)*r(1)+r(2)*r(2)+r(3)*r(3) );

   adouble r3 = pow(rr,3.0);

   rdd(1) = -mu*r(1)/r3;
   rdd(2) = -mu*r(2)/r3;
   rdd(3) = -mu*r(3)/r3;

   derivatives[ CINDEX(1) ] = v(1);
   derivatives[ CINDEX(2) ] = v(2);
   derivatives[ CINDEX(3) ] = v(3);
   derivatives[ CINDEX(4) ] = rdd(1);
   derivatives[ CINDEX(5) ] = rdd(2);
   derivatives[ CINDEX(6) ] = rdd(3);

}

////////////////////////////////////////////////////////////////////////////
///////////////////  Define the events function ////////////////////////////
////////////////////////////////////////////////////////////////////////////

void events(adouble* e, adouble* initial_states, adouble* final_states,
            adouble* parameters,adouble& t0, adouble& tf, adouble* xad,
            int iphase, Workspace* workspace)

{

   adouble ri1 = initial_states[ CINDEX(1) ];
   adouble ri2 = initial_states[ CINDEX(2) ];
   adouble ri3 = initial_states[ CINDEX(3) ];




   adouble rf1 = final_states[ CINDEX(1) ];
   adouble rf2 = final_states[ CINDEX(2) ];
   adouble rf3 = final_states[ CINDEX(3) ];

   e[ CINDEX(1) ] = ri1;
   e[ CINDEX(2) ] = ri2;
   e[ CINDEX(3) ] = ri3;

   e[ CINDEX(4) ] = rf1;
   e[ CINDEX(5) ] = rf2;
   e[ CINDEX(6) ] = rf3;

}



///////////////////////////////////////////////////////////////////////////
///////////////////  Define the phase linkages function ///////////////////
///////////////////////////////////////////////////////////////////////////

void linkages( adouble* linkages, adouble* xad, Workspace* workspace)
{
 //  auto_link_multiple(linkages, xad, N_PHASES);
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

    problem.name        		= "Lambert problem";
    problem.outfilename                 = "lambert.txt";

////////////////////////////////////////////////////////////////////////////
////////////  Define problem level constants & do level 1 setup ////////////
////////////////////////////////////////////////////////////////////////////

    problem.nphases   			= 1;
    problem.nlinkages                   = 0;

    psopt_level1_setup(problem);

/////////////////////////////////////////////////////////////////////////////
/////////   Define phase related information &  do level 2 setup ////////////
/////////////////////////////////////////////////////////////////////////////


    problem.phases(1).nstates   		= 6;
    problem.phases(1).ncontrols 		= 0;
    problem.phases(1).nevents   		= 6;
    problem.phases(1).nparameters               = 6;
    problem.phases(1).npath     		= 0;
    problem.phases(1).nodes              	= "[100]";

    psopt_level2_setup(problem, algorithm);


////////////////////////////////////////////////////////////////////////////
///////////////////  Enter problem bounds information //////////////////////
////////////////////////////////////////////////////////////////////////////

    double r1i = 15945.34e3; // m
    double r2i = 0.0;
    double r3i = 0.0;

    double r1f = 12214.83899e3; //m
    double r2f = 10249.46731e3; //m
    double r3f = 0.0;

    double TF = 76.0*60.0; // seconds

    problem.phases(1).bounds.lower.states(1) = -10*max(r1i,r1f);
    problem.phases(1).bounds.lower.states(2) = -10*max(r2i,r2f);
    problem.phases(1).bounds.lower.states(3) = -10*max(r3i,r3f);
    problem.phases(1).bounds.upper.states(1) = 10*max(r1i,r1f);
    problem.phases(1).bounds.upper.states(2) = 10*max(r2i,r2f);
    problem.phases(1).bounds.upper.states(3) = 10*max(r3i,r3f);

    problem.phases(1).bounds.lower.states(4) = -10*max(r1i,r1f)/TF;
    problem.phases(1).bounds.lower.states(5) = -10*max(r1i,r1f)/TF;;
    problem.phases(1).bounds.lower.states(6) = -10*max(r1i,r1f)/TF;;
    problem.phases(1).bounds.upper.states(4) = 10*max(r1i,r1f)/TF;
    problem.phases(1).bounds.upper.states(5) = 10*max(r2i,r2f)/TF;
    problem.phases(1).bounds.upper.states(6) = 10*max(r3i,r3f)/TF;


    problem.phases(1).bounds.lower.events(1) = r1i;
    problem.phases(1).bounds.upper.events(1) = r1i;

    problem.phases(1).bounds.lower.events(2) = r2i;
    problem.phases(1).bounds.upper.events(2) = r2i;

    problem.phases(1).bounds.lower.events(3) = r3i;
    problem.phases(1).bounds.upper.events(3) = r3i;

    problem.phases(1).bounds.lower.events(4) = r1f;
    problem.phases(1).bounds.upper.events(4) = r1f;

    problem.phases(1).bounds.lower.events(5) = r2f;
    problem.phases(1).bounds.upper.events(5) = r2f;

    problem.phases(1).bounds.lower.events(6) = r3f;
    problem.phases(1).bounds.upper.events(6) = r3f;


    problem.phases(1).bounds.lower.StartTime    = 0.0;
    problem.phases(1).bounds.upper.StartTime    = 0.0;

    problem.phases(1).bounds.lower.EndTime      = TF;
    problem.phases(1).bounds.upper.EndTime      = TF;


////////////////////////////////////////////////////////////////////////////
///////////////////  Enter problem names and units    //////////////////////
////////////////////////////////////////////////////////////////////////////


    problem.phases(1).name.states(1) = "x position";
    problem.phases(1).name.states(2) = "y position";
    problem.phases(1).name.states(3) = "z position";
    problem.phases(1).name.states(4) = "x velocity";
    problem.phases(1).name.states(5) = "y velocity";
    problem.phases(1).name.states(6) = "z velocity";

    problem.phases(1).units.states(1) = "m";
    problem.phases(1).units.states(2) = "m";
    problem.phases(1).units.states(3) = "m";
    problem.phases(1).units.states(4) = "m";
    problem.phases(1).units.states(5) = "m/s";
    problem.phases(1).units.states(6) = "m/s";




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

    int nnodes    			= 20;
    int ncontrols                       = problem.phases(1).ncontrols;
    int nstates                         = problem.phases(1).nstates;


    DMatrix x_guess    =  zeros(nstates,nnodes);
    DMatrix time_guess =  linspace(0.0,TF,nnodes);



    x_guess(1,colon()) = linspace(r1i,r1f, nnodes);
    x_guess(2,colon()) = linspace(r2i,r2f,nnodes);
    x_guess(3,colon()) = linspace(r3i,r3f,nnodes);
    x_guess(4,colon()) = linspace(r1i,r1f,nnodes)/TF;
    x_guess(5,colon()) = linspace(r2i,r2f,nnodes)/TF;
    x_guess(6,colon()) = linspace(r3i,r3f, nnodes)/TF;


    problem.phases(1).guess.states         = x_guess;
    problem.phases(1).guess.time           = time_guess;


////////////////////////////////////////////////////////////////////////////
///////////////////  Enter algorithm options  //////////////////////////////
////////////////////////////////////////////////////////////////////////////


    algorithm.nlp_iter_max                = 1000;
    algorithm.nlp_tolerance               = 1.e-6;
    algorithm.nlp_method                  = "IPOPT";
    algorithm.scaling                     = "automatic";
    algorithm.derivatives                 = "automatic";
    algorithm.defect_scaling              = "jacobian-based";
    algorithm.collocation_method          = "Hermite-Simpson";




////////////////////////////////////////////////////////////////////////////
///////////////////  Now call PSOPT to solve the problem   //////////////////
////////////////////////////////////////////////////////////////////////////

    psopt(solution, problem, algorithm);

////////////////////////////////////////////////////////////////////////////
///////////  Extract relevant variables from solution structure   //////////
////////////////////////////////////////////////////////////////////////////


    DMatrix x, u, t, xi, ui, ti;

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

    DMatrix r1 = x(1,colon());
    DMatrix r2 = x(2,colon());
    DMatrix r3 = x(3,colon());

    DMatrix v1 = x(4,colon());
    DMatrix v2 = x(5,colon());
    DMatrix v3 = x(6,colon());

    DMatrix vi(3), vf(3);

    vi(1) = v1(1);
    vi(2) = v2(1);
    vi(3) = v3(1);

    vf(1) = v1("end");
    vf(2) = v2("end");
    vf(3) = v3("end");

    vi.Print("Initial velocity vector [m/s]");

    vf.Print("Final velocity vector [m/s]");


    plot(t,r1,problem.name+": x", "time (s)", "x","x");

    plot(t,r2,problem.name+": y", "time (s)", "x","y");

    plot(t,r3,problem.name+": z", "time (s)", "z","z");

    plot(r1,r2,problem.name+": x-y", "x (m)", "y (m)","y");

    plot(r1,r2,problem.name+": x-y trajectory", "x (m)", "y (m)","y", "pdf", "lambert_xy.pdf");

    plot3(r1,r2,r3,problem.name+": trajectory","x (m)", "y (m)", "z (m)");



}

////////////////////////////////////////////////////////////////////////////
///////////////////////      END OF FILE     ///////////////////////////////
////////////////////////////////////////////////////////////////////////////
