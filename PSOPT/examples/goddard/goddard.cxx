//////////////////////////////////////////////////////////////////////////
//////////////////           goddard.cxx         /////////////////////////
//////////////////////////////////////////////////////////////////////////
////////////////           PSOPT  Example             ////////////////////
//////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////
//////// Title:         Goddard rocket maximum ascent     ////////////////
//////// Last modified: 05 January 2009                   ////////////////
//////// Reference:     Bryson (1999)                     ////////////////
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
    if (iphase==3) {
      adouble hf = final_states[1];
      return -hf;
    }
    else 
       return 0.0;
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
   adouble vdot, hdot, mdot;

   adouble v = states[0];
   adouble h = states[1];
   adouble m = states[2];

   adouble T = controls[0];

   double c        = 1580.9425;
   double sigma    = 5.4915e-5;
   double g        = 32.174;
   double h0       = 23800.0;

   adouble D;

   D = sigma*v*v*exp(-h/h0);

   adouble mg = (1.0+v/c)*D;
   adouble gg = mg/m;

   vdot = 1.0/m*(T-D)-g;
   hdot = v;
   mdot = -T/c;

   derivatives[0] = vdot;
   derivatives[1] = hdot;
   derivatives[2] = mdot;

   if (iphase ==2 ) {

 //     path[0] = T-D-mg- mg*(c*c*(1.0+v/c)/(h0*g)-1.0-2.0*c/v)/(1.0+4.0*(c/v)+2.0*pow(c/v,2.0));      
   
     path[0] = T-D-m*g- m*g*(c*c*(1.0+v/c)/(h0*g)-1.0-2.0*c/v)/(1.0+4.0*(c/v)+2.0*pow(c/v,2.0));
//     path[1] = m*g-(1.0+v/c)*sigma*v*v*exp(-h/h0);

   }
}

////////////////////////////////////////////////////////////////////////////
///////////////////  Define the events function ////////////////////////////
////////////////////////////////////////////////////////////////////////////

void events(adouble* e, adouble* initial_states, adouble* final_states, 
            adouble* parameters,adouble& t0, adouble& tf, adouble* xad, 
            int iphase, Workspace* workspace) 
{
   adouble vi = initial_states[0];
   adouble hi = initial_states[1];
   adouble mi = initial_states[2];
   adouble mf = final_states[2];

  
   if  (iphase == 1) {
     e[0] = vi;
     e[1] = hi;
     e[2] = mi;
   }

   if (iphase == 3) {
     e[0] = mf;
   }

}


///////////////////////////////////////////////////////////////////////////
///////////////////  Define the phase linkages function ///////////////////
///////////////////////////////////////////////////////////////////////////

void linkages( adouble* linkages, adouble* xad, Workspace* workspace)
{

    double c        = 1580.9425;
    double sigma    = 5.4915e-5;
    double g        = 32.174;
    double h0       = 23800.0;

    int index = 0;

   
    adouble xi[3], v, h, m;
    get_initial_states(xi, xad, 2, workspace);

    v  =  xi[0];
    h  =  xi[1];
    m  =  xi[2];

    auto_link(linkages, &index, xad, 1, 2, workspace );
    auto_link(linkages, &index, xad, 2, 3, workspace );

//    linkages[8] = m*g-(1.0+v/c)*sigma*v*v*exp(-h/h0);

    get_final_states(xi, xad, 2, workspace);

    v  =  xi[0];
    h  =  xi[1];
    m  =  xi[2];


    linkages[8] = m*g-(1.0+v/c)*sigma*v*v*exp(-h/h0);
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

    problem.name        		= "Goddard Rocket Maximum Ascent";
    problem.outfilename                 = "goddard.txt";

////////////////////////////////////////////////////////////////////////////
////////////  Declare problem level constants & setup phases  //////////////
////////////////////////////////////////////////////////////////////////////

    problem.nphases   			= 3;
    problem.nlinkages                   = 9;

    psopt_level1_setup(problem);


/////////////////////////////////////////////////////////////////////////////
/////////   Define phase related information &  setup PSOPT  /////////////////
/////////////////////////////////////////////////////////////////////////////

    int n1 = 20;
    int n2 = 10;
    int n3 = 30;


    problem.phases(1).nstates   		= 3;
    problem.phases(1).ncontrols 		= 1;
    problem.phases(1).nevents   		= 3;
    problem.phases(1).npath     		= 0;
    problem.phases(1).nodes                     = n1; 

    problem.phases(2).nstates   		= 3;
    problem.phases(2).ncontrols 		= 1;
    problem.phases(2).nevents   		= 0;
    problem.phases(2).npath     		= 1;
    problem.phases(2).nodes                     = n2; 

    problem.phases(3).nstates   		= 3;
    problem.phases(3).ncontrols 		= 1;
    problem.phases(3).nevents   		= 1;
    problem.phases(3).npath     		= 0;
    problem.phases(3).nodes                     = n3; 

    psopt_level2_setup(problem, algorithm);

////////////////////////////////////////////////////////////////////////////
///////////////////  Declare DMatrix objects to store results //////////////
////////////////////////////////////////////////////////////////////////////

    DMatrix x,u,t, x1, u1, t1, x2, u2, t2, x3, u3, t3;

////////////////////////////////////////////////////////////////////////////
///////////////////  Enter problem bounds information //////////////////////
////////////////////////////////////////////////////////////////////////////


    double v_L = 0;
    double h_L = 0.0;
    double m_L = 0.0;
    double v_U = 2000.0;
    double h_U = 30000.0;
    double m_U = 3.0;

    double T_L = 0.0;
    double T_U = 200.0;

    double v_i = 0.0;
    double h_i = 0.0;
    double m_i = 3.0;
    double m_f = 1.0;

    problem.phases(1).bounds.lower.states(1) = v_L;
    problem.phases(1).bounds.lower.states(2) = h_L;
    problem.phases(1).bounds.lower.states(3) = m_L;

    problem.phases(1).bounds.upper.states(1) = v_U;
    problem.phases(1).bounds.upper.states(2) = h_U;
    problem.phases(1).bounds.upper.states(3) = m_U;

    problem.phases(1).bounds.lower.controls(1) = T_L;
    problem.phases(1).bounds.upper.controls(1) = T_U;

    problem.phases(1).bounds.lower.events(1) = v_i;
    problem.phases(1).bounds.lower.events(2) = h_i;
    problem.phases(1).bounds.lower.events(3) = m_i;
   

    problem.phases(1).bounds.upper.events(1) = v_i;
    problem.phases(1).bounds.upper.events(2) = h_i;
    problem.phases(1).bounds.upper.events(3) = m_i;


    problem.phases(1).bounds.lower.StartTime    = 0.0;
    problem.phases(1).bounds.upper.StartTime    = 0.0;

    problem.phases(1).bounds.lower.EndTime      = 5.0;
    problem.phases(1).bounds.upper.EndTime      = 15.0;


    problem.phases(2).bounds.lower.states(1) = v_L;
    problem.phases(2).bounds.lower.states(2) = h_L;
    problem.phases(2).bounds.lower.states(3) = m_L;

    problem.phases(2).bounds.upper.states(1) = v_U;
    problem.phases(2).bounds.upper.states(2) = h_U;
    problem.phases(2).bounds.upper.states(3) = m_U;

    problem.phases(2).bounds.lower.controls(1) = T_L;
    problem.phases(2).bounds.upper.controls(1) = T_U;

    problem.phases(2).bounds.lower.path(1) = 0.0;
    problem.phases(2).bounds.upper.path(1) = 0.0;

 //   problem.phases(2).bounds.lower.path(2) = 0.0;
 //   problem.phases(2).bounds.upper.path(2) = 0.0;


    problem.phases(2).bounds.lower.StartTime    = 10.0;
    problem.phases(2).bounds.upper.StartTime    = 15.0;

    problem.phases(2).bounds.lower.EndTime      = 20.0;
    problem.phases(2).bounds.upper.EndTime      = 25.0;

    problem.phases(3).bounds.lower.states(1) = v_L;
    problem.phases(3).bounds.lower.states(2) = h_L;
    problem.phases(3).bounds.lower.states(3) = m_L;

    problem.phases(3).bounds.upper.states(1) = v_U;
    problem.phases(3).bounds.upper.states(2) = h_U;
    problem.phases(3).bounds.upper.states(3) = m_U;

    problem.phases(3).bounds.lower.controls(1) = T_L;
    problem.phases(3).bounds.upper.controls(1) = T_U;

   problem.phases(3).bounds.lower.events(1) = m_f;
   problem.phases(3).bounds.upper.events(1) = m_f;


    problem.phases(3).bounds.lower.StartTime    = 20.0;
    problem.phases(3).bounds.upper.StartTime    = 25.0;

    problem.phases(3).bounds.lower.EndTime      = 40.0;
    problem.phases(3).bounds.upper.EndTime      = 50.0;

////////////////////////////////////////////////////////////////////////////
///////////////////  Register problem functions  ///////////////////////////
////////////////////////////////////////////////////////////////////////////


    problem.integrand_cost 	= &integrand_cost;
    problem.endpoint_cost 	= &endpoint_cost;
    problem.dae 		= &dae;
    problem.events 		= &events;
    problem.linkages		= &linkages;

////////////////////////////////////////////////////////////////////////////
///////////////////  Define & register initial guess ///////////////////////
////////////////////////////////////////////////////////////////////////////


    DMatrix x0(3,n1);

    x0(1,colon()) = linspace(v_i,v_i, n1);
    x0(2,colon()) = linspace(h_i,h_i, n1);
    x0(3,colon()) = linspace(m_i,m_i, n1);

    problem.phases(1).guess.controls       = T_U*ones(1,n1);
    problem.phases(1).guess.states         = x0;
    problem.phases(1).guess.time           = linspace(0.0, 15.0, n1); 

    x0.Resize(3,n2);

    x0(1,colon()) = linspace(v_i,v_i, n2);
    x0(2,colon()) = linspace(h_i,h_i, n2);
    x0(3,colon()) = linspace(m_i,m_i, n2);


    problem.phases(2).guess.controls       = zeros(1,n2);
    problem.phases(2).guess.states         = x0;
    problem.phases(2).guess.time           = linspace(15.0, 20.0, 20); 

    x0.Resize(3,n3);

    x0(1,colon()) = linspace(v_i,v_i, n3);
    x0(2,colon()) = linspace(h_i,h_i, n3);
    x0(3,colon()) = linspace(m_i,m_i, n3);

    problem.phases(3).guess.controls       = zeros(1,n3);
    problem.phases(3).guess.states         = x0;
    problem.phases(3).guess.time           = linspace(20.0, 40.0, n3); 


////////////////////////////////////////////////////////////////////////////
///////////////////  Enter algorithm options  //////////////////////////////
////////////////////////////////////////////////////////////////////////////
    

    algorithm.nlp_method                  = "IPOPT";
    algorithm.scaling                     = "automatic";
    algorithm.derivatives                 = "automatic";
    algorithm.nlp_iter_max                = 1000;
    algorithm.nlp_tolerance               = 1.e-6;
    algorithm.collocation_method          = "trapezoidal";
    algorithm.mesh_refinement             = "automatic";
    algorithm.mr_max_iterations           = 5;


////////////////////////////////////////////////////////////////////////////
///////////////////  Now call PSOPT to solve the problem   //////////////////
////////////////////////////////////////////////////////////////////////////

    psopt(solution, problem, algorithm);

////////////////////////////////////////////////////////////////////////////
///////////  Extract relevant variables from solution structure   //////////
////////////////////////////////////////////////////////////////////////////

    x1 = solution.get_states_in_phase(1);
    u1 = solution.get_controls_in_phase(1);
    t1 = solution.get_time_in_phase(1);

    x2 = solution.get_states_in_phase(2);
    u2 = solution.get_controls_in_phase(2);
    t2 = solution.get_time_in_phase(2);

    x3 = solution.get_states_in_phase(3);
    u3 = solution.get_controls_in_phase(3);
    t3 = solution.get_time_in_phase(3);

    x = x1 || x2 || x3;

    u = u1 || u2 || u3;
 
    t = t1 || t2 || t3;


////////////////////////////////////////////////////////////////////////////
///////////  Save solution data to files if desired ////////////////////////
////////////////////////////////////////////////////////////////////////////

    x.Save("x.dat");
    u.Save("u.dat");
    t.Save("t.dat");

////////////////////////////////////////////////////////////////////////////
///////////  Plot some results if desired (requires gnuplot) ///////////////
////////////////////////////////////////////////////////////////////////////

    plot(t,x,problem.name, "time (s)", "states", "v h m");

    plot(t,u,problem.name, "time (s)", "control", "T");

    plot(t,x,problem.name, "time (s)", "states", "v h m",
                           "pdf", "goddard_states.pdf");

    plot(t,u,problem.name, "time (s)", "control", "T",
                           "pdf", "goddard_control.pdf");

}

////////////////////////////////////////////////////////////////////////////
///////////////////////      END OF FILE     ///////////////////////////////
////////////////////////////////////////////////////////////////////////////

