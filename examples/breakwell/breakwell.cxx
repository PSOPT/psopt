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

   adouble x = states[ 0 ];
   adouble v = states[ 1 ];

   adouble u = controls[ 0 ];


   xdot = v;
   vdot = u;

   derivatives[ 0 ] = xdot;
   derivatives[ 1 ] = vdot;

}

////////////////////////////////////////////////////////////////////////////
///////////////////  Define the events function ////////////////////////////
////////////////////////////////////////////////////////////////////////////

void events(adouble* e, adouble* initial_states, adouble* final_states,
            adouble* parameters,adouble& t0, adouble& tf, adouble* xad,
            int iphase, Workspace* workspace)

{
   adouble x0 = initial_states[ 0 ];
   adouble v0 = initial_states[ 1 ];
   adouble xf = final_states[ 0 ];
   adouble vf = final_states[ 1 ];

   e[ 0 ] = x0;
   e[ 1 ] = v0;
   e[ 2 ] = xf;
   e[ 3 ] = vf;

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

    problem.name        		       = "Breakwell Problem";
    problem.outfilename              = "breakwell.txt";

////////////////////////////////////////////////////////////////////////////
////////////  Define problem level constants & do level 1 setup ////////////
////////////////////////////////////////////////////////////////////////////

    problem.nphases   			      = 1;
    problem.nlinkages               = 0;

    psopt_level1_setup(problem);

/////////////////////////////////////////////////////////////////////////////
/////////   Define phase related information & do level 2 setup  ////////////
/////////////////////////////////////////////////////////////////////////////

    problem.phases(1).nstates   		= 2;
    problem.phases(1).ncontrols 		= 1;
    problem.phases(1).nevents   		= 4;
    problem.phases(1).npath     		= 0;
    problem.phases(1).nodes(0)      = 200;

    psopt_level2_setup(problem, algorithm);

////////////////////////////////////////////////////////////////////////////
///////////////////  Declare DMatrix objects to store results //////////////
////////////////////////////////////////////////////////////////////////////

    MatrixXd x, u, t;
    MatrixXd lambda, H;

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

    problem.phases(1).bounds.lower.states   << xL, vL;

    problem.phases(1).bounds.upper.states   << xU, vU;

    problem.phases(1).bounds.lower.controls <<  uL;
    
    problem.phases(1).bounds.upper.controls <<  uU;

    problem.phases(1).bounds.lower.events   <<  x0, v0, xf, vf;

    problem.phases(1).bounds.upper.events   <<  x0, v0, xf, vf;

    problem.phases(1).bounds.lower.StartTime    = 0.0;
    
    problem.phases(1).bounds.upper.StartTime    = 0.0;

    problem.phases(1).bounds.lower.EndTime      = 1.0;

    problem.phases(1).bounds.upper.EndTime      = 1.0;



////////////////////////////////////////////////////////////////////////////
///////////////////  Register problem functions  ///////////////////////////
////////////////////////////////////////////////////////////////////////////


    problem.integrand_cost 	                 = &integrand_cost;
    problem.endpoint_cost 	                    = &endpoint_cost;
    problem.dae             	                 = &dae;
    problem.events 		                       = &events;
    problem.linkages		                       = &linkages;

////////////////////////////////////////////////////////////////////////////
///////////////////  Define & register initial guess ///////////////////////
////////////////////////////////////////////////////////////////////////////

    int nnodes    			             = problem.phases(1).nodes(0);
    int ncontrols                       = problem.phases(1).ncontrols;
    int nstates                         = problem.phases(1).nstates;

    MatrixXd x_guess                    =  zeros(nstates,nnodes);

    x_guess.row(0)                      = x0*ones(1,nnodes);
    x_guess.row(1)                      = v0*ones(1,nnodes);


    problem.phases(1).guess.controls    = zeros(ncontrols,nnodes);
    problem.phases(1).guess.states      = x_guess;
    problem.phases(1).guess.time        = linspace(0.0,1.0,nnodes);
    
    Save(problem.phases(1).guess.time, "guess_time.dat");
    Save(problem.phases(1).guess.states, "guess_states.dat");


////////////////////////////////////////////////////////////////////////////
///////////////////  Enter algorithm options  //////////////////////////////
////////////////////////////////////////////////////////////////////////////


    algorithm.nlp_iter_max                = 1000;
    algorithm.nlp_tolerance               = 1.e-6;
    algorithm.nlp_method                  = "IPOPT";
    algorithm.scaling                     = "automatic";
    algorithm.derivatives                 = "automatic";
//    algorithm.hessian                   = "exact";
//    algorithm.mesh_refinement           = "automatic";
    algorithm.collocation_method          = "Legendre";
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

    MatrixXd ua(1,nn), xa(2,nn), pa(2,nn);

    for(int i=0;i <nn;i++) {

        if (t(i)<t1) {
           ua(0,i) = -2.0/(3.0*l)*(1.0-t(i)/(3.0*l));
           xa(0,i) = l*( 1.0 - pow( (1.0-t(i)/(3.0*l)), 3.0) ) ;
           xa(1,i) = pow( 1.0 - t(i)/(3.0*l), 2.0);
           pa(0,i) = 2.0/(9.0*l*l);
           pa(1,i) = 2.0/(3.0*l)*(1.0-t(i)/(3*l));

        }
        else if (t(i)>=t1 && t(i)<t2) {

           ua(0,i) = 0.0;
           xa(0,i) = l;
           xa(1,i) = 0.0;
           pa(0,i) = 0.0;
           pa(1,i) = 0.0;

       }
       else if (t(i)>=t2) {
           ua(0,i) = -2.0/(3*l)*(1.0-(1.0-t(i))/(3.0*l));
           xa(0,i) = l*(1.0 - pow( (1.0-(1.0-t(i))/(3.0*l)), 3.0));
           xa(1,i) = -pow(1.0-(1.0-t(i))/(3.0*l), 2.0) ;
           pa(0,i) = -2.0/(9.0*l*l);
           pa(1,i) = 2.0/(3.0*l)*(1.0-(1.0-t(i))/(3*l));

       }


    }


////////////////////////////////////////////////////////////////////////////
///////////  Save solution data to files if desired ////////////////////////
////////////////////////////////////////////////////////////////////////////

    Save(x,"x.dat");
    Save(u,"u.dat");
    Save(t,"t.dat");
    Save(lambda,"p.dat");
    Save(H,"H.dat");

////////////////////////////////////////////////////////////////////////////
///////////  Plot some results if desired (requires gnuplot) ///////////////
////////////////////////////////////////////////////////////////////////////

    plot(t,xa,problem.name+": analytical states", "time (s)", "states","xa va");

    plot(t,ua,problem.name+": analytical controls","time (s)", "control", "ua");

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
