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
   return (0);
}

//////////////////////////////////////////////////////////////////////////
///////////////////  Define the integrand (Lagrange) cost function  //////
//////////////////////////////////////////////////////////////////////////

adouble integrand_cost(adouble* states, adouble* controls, adouble* parameters,
                     adouble& time, adouble* xad, int iphase, Workspace* workspace)
{
    adouble u = controls[0];
    return  0.5*u*u;
}


//////////////////////////////////////////////////////////////////////////
///////////////////  Define the DAE's ////////////////////////////////////
//////////////////////////////////////////////////////////////////////////

void dae(adouble* derivatives, adouble* path, adouble* states,
         adouble* controls, adouble* parameters, adouble& time,
         adouble* xad, int iphase, Workspace* workspace)
{
   adouble xdot, vdot;

   double g = 1.0;
   double a = 0.5*g;

   adouble x = states[ CINDEX(1) ];
   adouble v = states[ CINDEX(2) ];

   adouble u = controls[ CINDEX(1) ];


   xdot = v;
   vdot = u;

   derivatives[ CINDEX(1) ] = xdot;
   derivatives[ CINDEX(2) ] = vdot;

}

////////////////////////////////////////////////////////////////////////////
///////////////////  Define the events function ////////////////////////////
////////////////////////////////////////////////////////////////////////////

void events(adouble* e, adouble* initial_states, adouble* final_states,
            adouble* parameters,adouble& t0, adouble& tf, adouble* xad,
            int iphase, Workspace* workspace)

{
   adouble x0 = initial_states[ CINDEX(1) ];
   adouble v0 = initial_states[ CINDEX(2) ];
   adouble xf = final_states[ CINDEX(1) ];
   adouble vf = final_states[ CINDEX(2) ];

   e[ CINDEX(1) ] = x0;
   e[ CINDEX(2) ] = v0;
   e[ CINDEX(3) ] = xf;
   e[ CINDEX(4) ] = vf;

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

    problem.name        		= "Breakwell Problem";
    problem.outfilename                 = "breakwell.txt";

////////////////////////////////////////////////////////////////////////////
////////////  Define problem level constants & do level 1 setup ////////////
////////////////////////////////////////////////////////////////////////////

    problem.nphases   			= 1;
    problem.nlinkages                   = 0;

    psopt_level1_setup(problem);

/////////////////////////////////////////////////////////////////////////////
/////////   Define phase related information & do level 2 setup  ////////////
/////////////////////////////////////////////////////////////////////////////

    problem.phases(1).nstates   		= 2;
    problem.phases(1).ncontrols 		= 1;
    problem.phases(1).nevents   		= 4;
    problem.phases(1).npath     		= 0;
    problem.phases(1).nodes             = "[200]";

    psopt_level2_setup(problem, algorithm);

////////////////////////////////////////////////////////////////////////////
///////////////////  Declare DMatrix objects to store results //////////////
////////////////////////////////////////////////////////////////////////////

    DMatrix x, u, t;
    DMatrix lambda, H;

////////////////////////////////////////////////////////////////////////////
///////////////////  Enter problem bounds information //////////////////////
////////////////////////////////////////////////////////////////////////////

    double xL = -2.0;
    double vL = -2.0;
    double xU = 0.1;
    double vU = 2.0;

    double uL = -10.0;
    double uU = 10.0;


    double x0 = 0.0;
    double v0 = 1.0;
    double xf = 0.0;
    double vf = -1.0;


    problem.phases(1).bounds.lower.states(1) = xL;
    problem.phases(1).bounds.lower.states(2) = vL;


    problem.phases(1).bounds.upper.states(1) = xU;
    problem.phases(1).bounds.upper.states(2) = vU;


    problem.phases(1).bounds.lower.controls(1) = uL;
    problem.phases(1).bounds.upper.controls(1) = uU;

    problem.phases(1).bounds.lower.events(1) = x0;
    problem.phases(1).bounds.lower.events(2) = v0;
    problem.phases(1).bounds.lower.events(3) = xf;
    problem.phases(1).bounds.lower.events(4) = vf;


    problem.phases(1).bounds.upper.events(1) = x0;
    problem.phases(1).bounds.upper.events(2) = v0;
    problem.phases(1).bounds.upper.events(3) = xf;
    problem.phases(1).bounds.upper.events(4) = vf;


    problem.phases(1).bounds.lower.StartTime    = 0.0;
    problem.phases(1).bounds.upper.StartTime    = 0.0;

    problem.phases(1).bounds.lower.EndTime      = 1.0;
    problem.phases(1).bounds.upper.EndTime      = 1.0;



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

    DMatrix x_guess    =  zeros(nstates,nnodes);

    x_guess(1,colon()) = x0*ones(1,nnodes);
    x_guess(2,colon()) = v0*ones(1,nnodes);


    problem.phases(1).guess.controls       = zeros(ncontrols,nnodes);
    problem.phases(1).guess.states         = x_guess;
    problem.phases(1).guess.time           = linspace(0.0,1.0,nnodes);


////////////////////////////////////////////////////////////////////////////
///////////////////  Enter algorithm options  //////////////////////////////
////////////////////////////////////////////////////////////////////////////


    algorithm.nlp_iter_max                = 1000;
    algorithm.nlp_tolerance               = 1.e-6;
    algorithm.nlp_method                  = "IPOPT";
    algorithm.scaling                     = "automatic";
    algorithm.derivatives                 = "automatic";
//    algorithm.hessian                     = "exact";
//    algorithm.mesh_refinement             = "automatic";
    algorithm.collocation_method = "Hermite-Simpson";
//    algorithm.diff_matrix                  = "central-differences";
//    algorithm.defect_scaling = "jacobian-based";
    algorithm.nlp_tolerance               = 1.e-6;



////////////////////////////////////////////////////////////////////////////
///////////////////  Now call PSOPT to solve the problem   /////////////////
////////////////////////////////////////////////////////////////////////////

    psopt(solution, problem, algorithm);

////////////////////////////////////////////////////////////////////////////
///////////  Extract relevant variables from solution structure   //////////
////////////////////////////////////////////////////////////////////////////


    DMatrix muE;

    x      = solution.get_states_in_phase(1);
    u      = solution.get_controls_in_phase(1);
    t      = solution.get_time_in_phase(1);
    lambda = solution.get_dual_costates_in_phase(1);
    H      = solution.get_dual_hamiltonian_in_phase(1);
    muE    = solution.get_dual_events_in_phase(1);


////////////////////////////////////////////////////////////////////////////
///////////  Compute the analytical solution, which is valid when the //////
///////////  state constraint x<=l, with 0<=l<=1/6. In this case l=0.1 /////
///////////  The analytical solution is given (with some errors) in    /////
///////////  the book by Bryson and Ho (1975), pages 121-122.          /////
////////////////////////////////////////////////////////////////////////////


    double l = 0.1;

    double t1 = 3*l;
    double t2 = 1.0-3*l;

    int nn = length(t);

    DMatrix ua(1,nn), xa(2,nn), pa(2,nn);

    for(int i=1;i <=nn;i++) {

        if (t(i)<t1) {
           ua(1,i) = -2.0/(3.0*l)*(1.0-t(i)/(3.0*l));
           xa(1,i) = l*( 1.0 - pow( (1.0-t(i)/(3.0*l)), 3.0) ) ;
           xa(2,i) = pow( 1.0 - t(i)/(3.0*l), 2.0);
           pa(1,i) = 2.0/(9.0*l*l);
           pa(2,i) = 2.0/(3.0*l)*(1.0-t(i)/(3*l));

        }
        else if (t(i)>=t1 && t(i)<t2) {

           ua(1,i) = 0.0;
           xa(1,i) = l;
           xa(2,i) = 0.0;
           pa(1,i) = 0.0;
           pa(2,i) = 0.0;

       }
       else if (t(i)>=t2) {
           ua(1,i) = -2.0/(3*l)*(1.0-(1.0-t(i))/(3.0*l));
           xa(1,i) = l*(1.0 - pow( (1.0-(1.0-t(i))/(3.0*l)), 3.0));
           xa(2,i) = -pow(1.0-(1.0-t(i))/(3.0*l), 2.0) ;
           pa(1,i) = -2.0/(9.0*l*l);
           pa(2,i) = 2.0/(3.0*l)*(1.0-(1.0-t(i))/(3*l));

       }


    }


////////////////////////////////////////////////////////////////////////////
///////////  Save solution data to files if desired ////////////////////////
////////////////////////////////////////////////////////////////////////////

    x.Save("x.dat");
    u.Save("u.dat");
    t.Save("t.dat");
    lambda.Save("p.dat");
    H.Save("H.dat");

////////////////////////////////////////////////////////////////////////////
///////////  Plot some results if desired (requires gnuplot) ///////////////
////////////////////////////////////////////////////////////////////////////

    plot(t,x,t,xa,problem.name+": states", "time (s)", "states","x v xa va");

    plot(t,u,t,ua,problem.name+": controls","time (s)", "control", "u ua");

    plot(t,lambda,t,pa, problem.name+": costates","time (s)", "costates", "l_x l_v la_x la_v");

    plot(t,H,problem.name+": Hamiltonian","time (s)", "H", "H");


    plot(t,x,t,xa,problem.name+": states", "time (s)", "states","x v xa va", "pdf", "breakwell_states.pdf");

    plot(t,u,t,ua,problem.name+": control","time (s)", "control", "u ua", "pdf", "breakwell_control.pdf");

    plot(t,lambda,t,pa, problem.name+": costates","time (s)", "costates", "l_x l_v la_x la_v", "pdf", "breakwell_costates.pdf");


}

////////////////////////////////////////////////////////////////////////////
///////////////////////      END OF FILE     ///////////////////////////////
////////////////////////////////////////////////////////////////////////////
