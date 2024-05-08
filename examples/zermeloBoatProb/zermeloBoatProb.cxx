//////////////////////////////////////////////////////////////////////////
////////////////        zermeloBoatProb.cxx         /////////////////////
//////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////
//////// Title:         zermelo Boat Nav Prob problem      ////////////////
//////// Last modified: 24 March 2024                   ////////////////
//////// Reference:     Seywald, Kevin, and Hans Seywald. "Desensitized optimal control." AIAA SciTech 2019 forum. 2019.////////////////
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
   adouble x = final_states[0];

   return (-x);
}

//////////////////////////////////////////////////////////////////////////
///////////////////  Define the integrand (Lagrange) cost function  //////
//////////////////////////////////////////////////////////////////////////

adouble integrand_cost(adouble* states, adouble* controls, adouble* parameters,
                     adouble& time, adouble* xad, int iphase, Workspace* workspace)
{
    adouble x1  = states[0];
    adouble x2  = states[1];
    adouble x3  = states[2];
    adouble s11 = states[3];
    adouble s12 = states[4];
    adouble s13 = states[5];
    adouble s21 = states[6];
    adouble s22 = states[7];
    adouble s23 = states[8];
    adouble u1  = controls[0];
    adouble k1  = controls[1];
    adouble k2  = controls[1];

    adouble L;

    L = 50.0*( s13*s13 + s23*s23 ) + 10.0*( k1*k1 + k2*k2 );

    return  L;
}


//////////////////////////////////////////////////////////////////////////
///////////////////  Define the DAE's ////////////////////////////////////
//////////////////////////////////////////////////////////////////////////

void dae(adouble* derivatives, adouble* path, adouble* states,
         adouble* controls, adouble* parameters, adouble& time,
         adouble* xad, int iphase, Workspace* workspace)
{
   adouble x1dot, x2dot, x3dot, s11dot, s12dot, s13dot, s21dot, s22dot, s23dot;
   
    adouble x1  = states[0];
    adouble x2  = states[1];
    adouble x3  = states[2];
    adouble s11 = states[3];
    adouble s12 = states[4];
    adouble s13 = states[5];
    adouble s21 = states[6];
    adouble s22 = states[7];
    adouble s23 = states[8];
    adouble u1  = controls[0];
    adouble k1  = controls[1];
    adouble k2  = controls[1];

   x1dot = cos(u1)+x3*x2;
   x2dot = sin(u1);
   x3dot = 0;
   s11dot = -s11*(-k1*sin(u1))-s12*(k1*cos(u1));
   s12dot = -s11*(-k2*sin(u1)+x3)-s12*(k2*cos(u1));
   s13dot = -s11*(x2);
   s21dot = -s21*(-k1*sin(u1))-s22*(k1*cos(u1));
   s22dot = -s21*(-k2*sin(u1)+x3)-s22*(k2*cos(u1));
   s23dot = -s21*(x2);



   derivatives[ 0 ] = x1dot;
   derivatives[ 1 ] = x2dot;
   derivatives[ 2 ] = x3dot;
   derivatives[ 3 ] = s11dot;
   derivatives[ 4 ] = s12dot;
   derivatives[ 5 ] = s13dot;
   derivatives[ 6 ] = s21dot;
   derivatives[ 7 ] = s22dot;
   derivatives[ 8 ] = s23dot;


}

////////////////////////////////////////////////////////////////////////////
///////////////////  Define the events function ////////////////////////////
////////////////////////////////////////////////////////////////////////////

void events(adouble* e, adouble* initial_states, adouble* final_states,
            adouble* parameters,adouble& t0, adouble& tf, adouble* xad,
            int iphase, Workspace* workspace)

{
   adouble x10  = initial_states[ 0 ];
   adouble x20  = initial_states[ 1 ];
   adouble x30  = initial_states[ 2 ];
   adouble x2f  = final_states[ 1 ];
   adouble s11f = final_states[ 3 ];
   adouble s12f = final_states[ 4 ];
   adouble s13f = final_states[ 5 ];
   adouble s21f = final_states[ 6 ];
   adouble s22f = final_states[ 7 ];
   adouble s23f = final_states[ 8 ];

   e[ 0 ] = x10;
   e[ 1 ] = x20;
   e[ 2 ] = x30;
   e[ 3 ] = x2f;
   e[ 4 ] = s11f;
   e[ 5 ] = s12f;
   e[ 6 ] = s13f;
   e[ 7 ] = s21f;
   e[ 8 ] = s22f;
   e[ 9 ] = s23f;

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

    problem.name        		= "Zermelo Boat Navigation Problem";
    problem.outfilename         = "zermeloBoatNavProb.txt";

////////////////////////////////////////////////////////////////////////////
////////////  Define problem level constants & do level 1 setup ////////////
////////////////////////////////////////////////////////////////////////////

    problem.nphases   			= 1;
    problem.nlinkages           = 0;

    psopt_level1_setup(problem);

/////////////////////////////////////////////////////////////////////////////
/////////   Define phase related information & do level 2 setup  ////////////
/////////////////////////////////////////////////////////////////////////////

    problem.phases(1).nstates   		= 9;
    problem.phases(1).ncontrols 		= 3;
    problem.phases(1).nevents   		= 10;
    problem.phases(1).npath     		= 0;
    problem.phases(1).nodes            << 10;

    psopt_level2_setup(problem, algorithm);

////////////////////////////////////////////////////////////////////////////
///////////////////  Declare MatrixXd objects to store results //////////////
////////////////////////////////////////////////////////////////////////////

    MatrixXd x, u, t;
    MatrixXd lambda, H;

////////////////////////////////////////////////////////////////////////////
///////////////////  Enter problem bounds information //////////////////////
////////////////////////////////////////////////////////////////////////////

    double x1L  =   0; 
    double x2L  = -10;
    double x3L  =  10;
    double s11L = -10;
    double s12L = -10;
    double s13L = -10;
    double s21L = -10;
    double s22L = -10;
    double s23L = -10;

    double x1U  =  10; 
    double x2U  =  10;
    double x3U  =  10;
    double s11U =  10;
    double s12U =  10;
    double s13U =  10;
    double s21U =  10;
    double s22U =  10;
    double s23U =  10;

    double u1L = -10.0;
    double k1L = -10.0;
    double k2L = -10.0;
    double u1U =  10.0;
    double k1U =  10.0;
    double k2U =  10.0;

    double x10  = 0;
    double x20  = 0;
    double x30  = 10;
    double x2f  = 0;
    double s11f = 1;
    double s12f = 0;
    double s13f = 0;
    double s21f = 0;
    double s22f = 1;
    double s23f = 0;


    problem.phases(1).bounds.lower.states(0) = x1L ; 
    problem.phases(1).bounds.lower.states(1) = x2L ;
    problem.phases(1).bounds.lower.states(2) = x3L ;
    problem.phases(1).bounds.lower.states(3) = s11L;
    problem.phases(1).bounds.lower.states(4) = s12L;
    problem.phases(1).bounds.lower.states(5) = s13L;
    problem.phases(1).bounds.lower.states(6) = s21L;
    problem.phases(1).bounds.lower.states(7) = s22L;
    problem.phases(1).bounds.lower.states(8) = s23L;


    problem.phases(1).bounds.upper.states(0) = x1U; 
    problem.phases(1).bounds.upper.states(1) = x2U;
    problem.phases(1).bounds.upper.states(2) = x3U;
    problem.phases(1).bounds.upper.states(3) = s11U;
    problem.phases(1).bounds.upper.states(4) = s12U;
    problem.phases(1).bounds.upper.states(5) = s13U;
    problem.phases(1).bounds.upper.states(6) = s21U;
    problem.phases(1).bounds.upper.states(7) = s22U;
    problem.phases(1).bounds.upper.states(8) = s23U;


    problem.phases(1).bounds.lower.controls(0) = u1L;
    problem.phases(1).bounds.lower.controls(1) = k1L;
    problem.phases(1).bounds.lower.controls(2) = k2L;
    problem.phases(1).bounds.upper.controls(0) = u1U;
    problem.phases(1).bounds.upper.controls(1) = k1U;
    problem.phases(1).bounds.upper.controls(2) = k2U;

    problem.phases(1).bounds.lower.events(0) = x10 ;
    problem.phases(1).bounds.lower.events(1) = x20 ;
    problem.phases(1).bounds.lower.events(2) = x30 ;
    problem.phases(1).bounds.lower.events(3) = x2f ;
    problem.phases(1).bounds.lower.events(4) = s11f;
    problem.phases(1).bounds.lower.events(5) = s12f;
    problem.phases(1).bounds.lower.events(6) = s13f;
    problem.phases(1).bounds.lower.events(7) = s21f;
    problem.phases(1).bounds.lower.events(8) = s22f;
    problem.phases(1).bounds.lower.events(9) = s23f;


    problem.phases(1).bounds.upper.events(0) = x10 ;
    problem.phases(1).bounds.upper.events(1) = x20 ;
    problem.phases(1).bounds.upper.events(2) = x30 ;
    problem.phases(1).bounds.upper.events(3) = x2f ;
    problem.phases(1).bounds.upper.events(4) = s11f;
    problem.phases(1).bounds.upper.events(5) = s12f;
    problem.phases(1).bounds.upper.events(6) = s13f;
    problem.phases(1).bounds.upper.events(7) = s21f;
    problem.phases(1).bounds.upper.events(8) = s22f;
    problem.phases(1).bounds.upper.events(9) = s23f;



    problem.phases(1).bounds.lower.StartTime    = 0.0;
    problem.phases(1).bounds.upper.StartTime    = 0.0;

    problem.phases(1).bounds.lower.EndTime      = 1.0;
    problem.phases(1).bounds.upper.EndTime      = 1.0;



////////////////////////////////////////////////////////////////////////////
///////////////////  Register problem functions  ///////////////////////////
////////////////////////////////////////////////////////////////////////////


    problem.integrand_cost 	= &integrand_cost;
    problem.endpoint_cost 	= &endpoint_cost;
    problem.dae             = &dae;
    problem.events 		    = &events;
    problem.linkages		= &linkages;

////////////////////////////////////////////////////////////////////////////
///////////////////  Define & register initial guess ///////////////////////
////////////////////////////////////////////////////////////////////////////

    int nnodes    			            = problem.phases(1).nodes(0);
    int ncontrols                       = problem.phases(1).ncontrols;
    int nstates                         = problem.phases(1).nstates;

    MatrixXd x_guess    =  zeros(nstates,nnodes);

    x_guess.row(0)  = x10*ones(1,nnodes);
    x_guess.row(1)  = x20*ones(1,nnodes);
    x_guess.row(2)  = x30*ones(1,nnodes);
    x_guess.row(3)  = s11f*ones(1,nnodes);
    x_guess.row(4)  = s12f*ones(1,nnodes);
    x_guess.row(5)  = s13f*ones(1,nnodes);
    x_guess.row(6)  = s21f*ones(1,nnodes);
    x_guess.row(7)  = s22f*ones(1,nnodes);
    x_guess.row(8)  = s23f*ones(1,nnodes);

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
    algorithm.jac_sparsity_ratio          = 0.20;
    algorithm.collocation_method          = "Legendre";
    algorithm.diff_matrix                 = "central-differences";
    algorithm.mesh_refinement             = "automatic";
    algorithm.mr_max_increment_factor     = 0.3;
    algorithm.mr_max_iterations           = 3;
    algorithm.defect_scaling              = "jacobian-based";



////////////////////////////////////////////////////////////////////////////
///////////////////  Now call PSOPT to solve the problem   /////////////////
////////////////////////////////////////////////////////////////////////////

    psopt(solution, problem, algorithm);

////////////////////////////////////////////////////////////////////////////
///////////  Extract relevant variables from solution structure   //////////
////////////////////////////////////////////////////////////////////////////


    x      = solution.get_states_in_phase(1);
    u      = solution.get_controls_in_phase(1);
    t      = solution.get_time_in_phase(1);
    lambda = solution.get_dual_costates_in_phase(1);
    H      = solution.get_dual_hamiltonian_in_phase(1);


////////////////////////////////////////////////////////////////////////////
///////////  Save solution data to files if desired ////////////////////////
////////////////////////////////////////////////////////////////////////////

    Save(x, "x.dat");
    Save(u,"u.dat");
    Save(t,"t.dat");
    Save(lambda,"lambda.dat");
    Save(H,"H.dat");

////////////////////////////////////////////////////////////////////////////
///////////  Plot some results if desired (requires gnuplot) ///////////////
////////////////////////////////////////////////////////////////////////////

    plot(t,x.row(0),problem.name+": state", "time (s)", "state","x1");

    plot(t,x.row(1),problem.name+": state", "time (s)", "state","x2");

    plot(t,x.row(5),problem.name+": sensitivity", "time (s)", "sensitivity","s13");

    plot(t,x.row(8),problem.name+": sensitivity", "time (s)", "sensitivity","s23");

    plot(t,u.row(0),problem.name+": control","time (s)", "control", "u1");

    plot(t,u.row(1),problem.name+": gain","time (s)", "gain", "k1");

    plot(t,u.row(2),problem.name+": gain","time (s)", "gain", "k2");

    plot(t,x.row(0),problem.name+": state x1", "time (s)", "state","x1",
                                        "pdf", "zermeloBoatProbResults_x1.pdf");

    plot(t,x.row(1),problem.name+": state x2", "time (s)", "state","x2",
                                        "pdf", "zermeloBoatProbResults_x2.pdf");

    plot(t,x.row(5),problem.name+": sensitivity", "time (s)", "sensitivity","s13",
                                        "pdf", "zermeloBoatProbResults_13.pdf");

    plot(t,x.row(8),problem.name+": sensitivity", "time (s)", "sensitivity","s23",
                                        "pdf", "zermeloBoatProbResults_s23.pdf");

    plot(t,u.row(0),problem.name+": control u1","time (s)", "control", "u1",
                                        "pdf", "zermeloBoatProbResults_u1.pdf");

    plot(t,u.row(1),problem.name+": gain","time (s)", "gain", "k1",
                                        "pdf", "zermeloBoatProbResults_k1.pdf");

    plot(t,u.row(2),problem.name+": gain","time (s)", "gain", "k2",
                                        "pdf", "zermeloBoatProbResults_k2.pdf");
}

////////////////////////////////////////////////////////////////////////////
///////////////////////      END OF FILE     ///////////////////////////////
////////////////////////////////////////////////////////////////////////////
