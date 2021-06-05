//////////////////////////////////////////////////////////////////////////
////////////////              heat.cxx               /////////////////////
//////////////////////////////////////////////////////////////////////////
////////////////           PSOPT Template            /////////////////////
//////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////
//////// Title: 	Heat difussion process            ////////////////
//////// Last modified: 09 July 2009                      ////////////////
//////// Reference:     Betts (2001)             	  ////////////////
////////                                                  ////////////////
//////////////////////////////////////////////////////////////////////////
////////     Copyright (c) Victor M. Becerra, 2009        ////////////////
//////////////////////////////////////////////////////////////////////////
//////// This is part of the PSOPT software library, which ////////////////
//////// is distributed under the terms of the GNU Lesser ////////////////
//////// General Public License (LGPL)                    ////////////////
//////////////////////////////////////////////////////////////////////////

#include "psopt.h"

using namespace PSOPT;



//////////////////////////////////////////////////////////////////////////
///////////////////  Define the end point (Mayer) cost function //////////
//////////////////////////////////////////////////////////////////////////


#define N_DISCRETIZATION  10


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

     int N = N_DISCRETIZATION;
 
     double gamma 	= 1.0e-3; 

     double rho 		= -1.0;

     adouble yd = 2.0-exp(rho*time);
 
     adouble yN = states[ N-1 ];

     adouble v1 = controls[ 0 ];

     return   0.5*(  pow(yN-yd,2.0) + gamma*pow(v1,2.0) );
} 


//////////////////////////////////////////////////////////////////////////
///////////////////  Define the DAE's ////////////////////////////////////
//////////////////////////////////////////////////////////////////////////


void dae(adouble* derivatives, adouble* path, adouble* states, 
         adouble* controls, adouble* parameters, adouble& time, 
         adouble* xad, int iphase, Workspace* workspace)
{

   adouble u;

   double a1 		= 4.0;
   double a2 		= 1.0;
   double a3 		= 4.0;
   double a4 		= -1.0;
   double rho 		= -1.0;
   double T 		= 0.5;
   double g 		= 1.0;

   int N 		= N_DISCRETIZATION;

   int i;

   double delta  	= 1.0/(N-1);

  

   adouble* y 		= states;

   adouble v1 		= controls[ 0 ];
   adouble v2 		= controls[ 1 ];
   adouble v3 		= controls[ 2 ];

   adouble y1 		= y[ 0 ];
   adouble y2 		= y[ 1 ];
   adouble yN 		= y[ N-1 ];
   adouble yNm1 	= y[ N-2 ];

   adouble x1 		=   0;
   adouble xN 		=   1;

   adouble q1 =  (rho*(a1+2*a2) + pi*pi*(a3+2*a4))*exp(rho*time)*cos(pi*x1) 
                     - a4*pi*pi*exp(2*rho*time) + (2*a4*pi*pi+rho*a2)*exp(2*rho*time)*pow(cos(pi*x1),2.0);
   adouble qN =  (rho*(a1+2*a2) + pi*pi*(a3+2*a4))*exp(rho*time)*cos(pi*xN) 
                     - a4*pi*pi*exp(2*rho*time) + (2*a4*pi*pi+rho*a2)*exp(2*rho*time)*pow(cos(pi*xN),2.0);

   derivatives[ 0 ] = 1.0/(a1+a2*y1)*(q1 + 
                              (1.0/pow(delta,2.0))*(a3+a4*y1)*(y2-2*y1+v2)+a4*pow( (y2-v2)/(2*delta), 2.0) );
  
   for(i=1;i<(N-1);i++) {
       adouble yi   = y[i];
       adouble yim1 = y[i-1];
       adouble yip1 = y[i+1];
       adouble xi =   ( (double) (i) ) /( (double) (N-1) );
       adouble qi = (rho*(a1+2*a2) + pi*pi*(a3+2*a4))*exp(rho*time)*cos(pi*xi) 
                     - a4*pi*pi*exp(2*rho*time) + (2*a4*pi*pi+rho*a2)*exp(2*rho*time)*pow(cos(pi*xi),2.0);

       derivatives[i]=1.0/(a1+a2*yi)*(qi + (1.0/pow(delta,2.0))*
                              (a3+a4*yi)*(yip1-2*yi+yim1)+a4*pow( (yip1-yim1)/(2*delta), 2.0) );
   }

   derivatives[ N-1 ] = 1.0/(a1+a2*yN)*(qN + (1.0/pow(delta,2.0))*
                              (a3+a4*yN)*(v3-2*yN+yNm1)+a4*pow( (v3-yNm1)/(2*delta), 2.0) );

   path[ 0 ] =  g*(y1-v1) - (1.0/(2*delta))*(a3+a4*y1)*(y2-v2); 

   path[ 1 ] = (1.0/(2*delta))*(a3+a4*yN)*(v3-yNm1);

}

////////////////////////////////////////////////////////////////////////////
///////////////////  Define the events function ////////////////////////////
////////////////////////////////////////////////////////////////////////////

void events(adouble* e, adouble* initial_states, adouble* final_states, 
            adouble* parameters,adouble& t0, adouble& tf, adouble* xad, 
            int iphase, Workspace* workspace) 

{

   int i;

   int N = N_DISCRETIZATION;

   for(i=0;i< N; i++) {
          e[ i ] = initial_states[ i ];
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

   int N = N_DISCRETIZATION;

////////////////////////////////////////////////////////////////////////////
///////////////////  Declare key structures ////////////////////////////////
////////////////////////////////////////////////////////////////////////////

    Alg  algorithm;
    Sol  solution;
    Prob problem;

////////////////////////////////////////////////////////////////////////////
///////////////////  Register problem name  ////////////////////////////////
////////////////////////////////////////////////////////////////////////////

    problem.name        		= "Heat diffusion process";
    problem.outfilename                 = "heat.txt";

////////////////////////////////////////////////////////////////////////////
////////////  Define problem level constants & do level 1 setup ////////////
////////////////////////////////////////////////////////////////////////////

    problem.nphases   					= 1;
    problem.nlinkages               = 0;

    psopt_level1_setup(problem);

/////////////////////////////////////////////////////////////////////////////
/////////   Define phase related information & do level 2 setup  ////////////
/////////////////////////////////////////////////////////////////////////////

    problem.phases(1).nstates   		= N;
    problem.phases(1).ncontrols 		= 3;
    problem.phases(1).nevents   		= N;
    problem.phases(1).npath     		= 2;
    problem.phases(1).nodes         = (RowVectorXi(3)<< 20, 50, 60 ).finished();            //  "[20 50 60]"; 


    psopt_level2_setup(problem, algorithm);



////////////////////////////////////////////////////////////////////////////
///////////////////  Enter problem bounds information //////////////////////
////////////////////////////////////////////////////////////////////////////

    int i, j;


    problem.phases(1).bounds.lower.states = zeros(N,1);
    problem.phases(1).bounds.upper.states = 3.0*ones(N,1);


    problem.phases(1).bounds.lower.controls <<  0.0, 0.0, 0.0;

    problem.phases(1).bounds.upper.controls << 0.1, 3.0, 3.0;


    for(i = 0; i< N; i++ ) {
        double xi =    ( (double) (i) ) /( (double) (N-1) );
    	double yI = 2.0+cos(pi*xi);

    	problem.phases(1).bounds.lower.events(i) = yI;
    	problem.phases(1).bounds.upper.events(i) = yI;

    }

    problem.phases(1).bounds.lower.path << 0.0, 0.0;

    problem.phases(1).bounds.upper.path << 0.0, 0.0;

    problem.phases(1).bounds.lower.StartTime    = 0.0;
    problem.phases(1).bounds.upper.StartTime    = 0.0;

    problem.phases(1).bounds.lower.EndTime      = 0.5;
    problem.phases(1).bounds.upper.EndTime      = 0.5;


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

    int nnodes    						    = problem.phases(1).nodes(0);
    int ncontrols                       = problem.phases(1).ncontrols;
    int nstates                         = problem.phases(1).nstates;

    MatrixXd state_guess    =  zeros(nstates,nnodes);
    MatrixXd control_guess  =  zeros(ncontrols,nnodes);
    MatrixXd param_guess    =  zeros(0,0);
    MatrixXd time_guess     =  linspace(0,0.5,nnodes);

    for(i = 0; i< N; i++ ) {
        double xi =    ( (double) (i) ) /( (double) (N-1) );
    	double yI = 2.0+cos(pi*xi);
    	state_guess.row(i) = yI*ones(1,nnodes);
    }

    control_guess <<  linspace(0.1,0.1, nnodes),
                      state_guess.row(0), 
                      state_guess.row(N-1);


    auto_phase_guess(problem, control_guess, state_guess, param_guess, time_guess);


////////////////////////////////////////////////////////////////////////////
///////////////////  Enter algorithm options  //////////////////////////////
////////////////////////////////////////////////////////////////////////////

    algorithm.nlp_iter_max                = 1000;
    algorithm.nlp_tolerance               = 1.e-6;
    algorithm.nlp_method                  = "IPOPT";
    algorithm.scaling                     = "automatic";
    algorithm.derivatives                 = "automatic";
    algorithm.collocation_method          = "trapezoidal";
//    algorithm.defect_scaling              = "jacobian-based";    

////////////////////////////////////////////////////////////////////////////
///////////////////  Now call PSOPT to solve the problem   //////////////////
////////////////////////////////////////////////////////////////////////////

    psopt(solution, problem, algorithm);

////////////////////////////////////////////////////////////////////////////
///////////  Extract relevant variables from solution structure   //////////
////////////////////////////////////////////////////////////////////////////


    MatrixXd y, u, t;

    y      = solution.get_states_in_phase(1);
    u      = solution.get_controls_in_phase(1);
    t      = solution.get_time_in_phase(1);


////////////////////////////////////////////////////////////////////////////
///////////  Save solution data to files if desired ////////////////////////
////////////////////////////////////////////////////////////////////////////

    Save(y,"y.dat");
    Save(u,"u.dat");
    Save(t,"t.dat");    

    MatrixXd x(N,1);

    for(i = 0; i< N; i++ ) {
        x(i) =    ( (double) (i) ) /( (double) (N-1) );
    }

////////////////////////////////////////////////////////////////////////////
///////////  Plot some results if desired (requires gnuplot) ///////////////
////////////////////////////////////////////////////////////////////////////


    plot(t,u.row(0),problem.name+": control","time (s)", "control", "u");

    plot(t,u.row(0),problem.name+": control","time (s)", "control", "v1", 
			"pdf", "heat_control.pdf");

    plot(t,y.row(0) ,problem.name+": state","time (s)", "state", "x1");

    plot(t,y.row(0) ,problem.name+": state","time (s)", "state", "x1",
			"pdf", "heat_state1.pdf");

    plot(t,y.row(N-1) ,problem.name+": state","time (s)", "state", "xN");

    plot(t,y.row(N-1) ,problem.name+": state","time (s)", "state", "xN",
			"pdf", "heat_stateN.pdf");

    surf(x, t,  y, problem.name, "x", "t", "h"); 

    surf(x, t,  y, problem.name, "z", "t", "x", "pdf", "heat_surf.pdf"); 

}

////////////////////////////////////////////////////////////////////////////
///////////////////////      END OF FILE     ///////////////////////////////
////////////////////////////////////////////////////////////////////////////
