//////////////////////////////////////////////////////////////////////////
//////////////////       shuttle_reentry1.cxx       //////////////////////
//////////////////////////////////////////////////////////////////////////
////////////////           PSOPT  Example             ////////////////////
//////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////
//////// Title:         Shuttle reentry problem           ////////////////
//////// Last modified: 05 January 2009                   ////////////////
//////// Reference:     Betts (2001)                      ////////////////
//////// (See PSOPT handbook for full reference)          ////////////////
//////////////////////////////////////////////////////////////////////////
////////     Copyright (c) Victor M. Becerra, 2009        ////////////////
//////////////////////////////////////////////////////////////////////////
//////// This is part of the PSOPT software library, which ///////////////
//////// is distributed under the terms of the GNU Lesser  ///////////////
//////// General Public License (LGPL)                     ///////////////
//////////////////////////////////////////////////////////////////////////

#include "psopt.h"

//////////////////////////////////////////////////////////////////////////
//////////////////    Define some local macros      //////////////////////
//////////////////////////////////////////////////////////////////////////


#define H_INDX   	1
#define PHI_INDX 	2
#define THETA_INDX 	3
#define V_INDX		4
#define GAMMA_INDX	5
#define PSI_INDX	6

#define ALPHA_INDX   	1
#define BETA_INDX	2

#define DEG2RAD(x)  (3.141592653589793*(x)/180.0)


//////////////////////////////////////////////////////////////////////////
///////////////////  Define the end point (Mayer) cost function //////////
//////////////////////////////////////////////////////////////////////////

adouble endpoint_cost(adouble* initial_states, adouble* final_states,
                      adouble* parameters,adouble& t0, adouble& tf,
                      adouble* xad, int iphase, Workspace* workspace)
{
    adouble theta = final_states[THETA_INDX-1];

    return (-theta*180/3.141592653589793);
}

//////////////////////////////////////////////////////////////////////////
///////////////////  Define the integrand (Lagrange) cost function  //////
//////////////////////////////////////////////////////////////////////////

adouble integrand_cost(adouble* states, adouble* controls,
                       adouble* parameters,
                       adouble& time, adouble* xad,
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

	adouble alpha, beta;


	adouble alt    = states[H_INDX-1];
	adouble lon    = states[PHI_INDX-1];
	adouble lat    = states[THETA_INDX-1];
	adouble vel    = states[V_INDX-1];
	adouble gamma  = states[GAMMA_INDX-1];
	adouble azi    = states[PSI_INDX-1];


	alpha = controls[ALPHA_INDX-1];
	beta  = controls[BETA_INDX-1];


	double pi      = 3.141592653589793;
	double cr2d    = 180.0/pi;
	double weight  = 203000.0;
	double cm2w    = 32.174;
	double cea     = 20902900.0;
	double mu      = 0.14076539e17;
	double rho0    = 0.002378;
	double href    = 23800.0;
	double cl0     = -0.20704;
	double cl1     = 0.029244;
	double cd0     = 0.07854;
	double cd1     = -6.1592e-3;
	double cd2     = 6.21408e-4;
	double sref    = 2690.0;


	double mass    = weight/cm2w;
	adouble sgamma  = sin(gamma);
	adouble cgamma  = cos(gamma);
	adouble sazi    = sin(azi);
	adouble cazi    = cos(azi);
	adouble sbeta   = sin(beta);
	adouble cbeta   = cos(beta);
	adouble slat    = sin(lat);
	adouble clat    = cos(lat);

	adouble alphad  = cr2d*alpha;
	adouble radius  = cea+alt;
	adouble grav    = mu/pow(radius,2);

	adouble rhodns  = rho0*exp(-alt/href);
	adouble dynp    = 0.5*(rhodns*pow(vel,2));
	adouble subl    = cl0+cl1*alphad;
	adouble subd    = cd0+((cd1+cd2*alphad)*alphad);
	adouble drag    = (dynp*subd)*sref;
	adouble lift    = (dynp*subl)*sref;
	adouble vrelg   = (vel/radius)-(grav/vel);

	adouble d_alt_dt   = vel*sgamma;
	adouble d_lon_dt   = ((vel*cgamma)*sazi)/(radius*clat);
	adouble d_lat_dt   = ((vel*cgamma)*cazi)/radius;
	adouble d_vel_dt   = (-(drag/mass)-(grav*sgamma));
	adouble d_gamma_dt = ((lift*cbeta)/(mass*vel))+(cgamma*vrelg);
	adouble d_azi_dt   = ((lift*sbeta)/((mass*vel)*cgamma))
                              + ((vel*cgamma)*(sazi*slat)/(radius*clat));

	derivatives[H_INDX-1] 	= d_alt_dt      ;
	derivatives[PHI_INDX-1] 	= d_lon_dt      ;
	derivatives[THETA_INDX-1] 	= d_lat_dt      ;
	derivatives[V_INDX-1] 	= d_vel_dt      ;
	derivatives[GAMMA_INDX-1] 	= d_gamma_dt    ;
	derivatives[PSI_INDX-1] 	= d_azi_dt      ;

//    adouble alpha_hat = 180.0*alpha/pi;

//    adouble qa  = c0  + c1*alpha_hat + c2*pow(alpha_hat,2) + c3*pow(alpha_hat,3);

//    adouble qr = (17700.0*sqrt(rhodns))*pow(0.0001*v,3.07);

//    adouble q = qa*qr;

//    path[0] = q;


}

////////////////////////////////////////////////////////////////////////////
///////////////////  Define the events function ////////////////////////////
////////////////////////////////////////////////////////////////////////////

void events(adouble* e, adouble* initial_states, adouble* final_states,
            adouble* parameters,adouble& t0, adouble& tf, adouble* xad,
            int iphase, Workspace* workspace)
{
   adouble h0 		= initial_states[H_INDX-1];
   adouble phi0 	= initial_states[PHI_INDX-1];
   adouble theta0 	= initial_states[THETA_INDX-1];
   adouble v0 		= initial_states[V_INDX-1];
   adouble gamma0 	= initial_states[GAMMA_INDX-1];
   adouble psi0         = initial_states[PSI_INDX-1];

   adouble hf           = final_states[H_INDX-1];
   adouble vf           = final_states[V_INDX-1];
   adouble gammaf       = final_states[GAMMA_INDX-1];

   e[H_INDX-1] 		= h0;
   e[PHI_INDX-1] 	= phi0;
   e[THETA_INDX-1]	= theta0;
   e[V_INDX-1] 		= v0;
   e[GAMMA_INDX-1]	= gamma0;
   e[PSI_INDX-1] 	= psi0;

   e[6] 		= hf;
   e[7] 		= vf;
   e[8] 		= gammaf;
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

    problem.name        		= "Shuttle re-entry problem";
    problem.outfilename                 = "shuttle.txt";

////////////////////////////////////////////////////////////////////////////
////////////  Define problem level constants & do level 1 setup ////////////
////////////////////////////////////////////////////////////////////////////

    problem.nphases   			= 1;
    problem.nlinkages                   = 0;

    psopt_level1_setup(problem);

/////////////////////////////////////////////////////////////////////////////
/////////   Define phase related information & do level 2 setup /////////////
/////////////////////////////////////////////////////////////////////////////

    problem.phases(1).nstates   		= 6;
    problem.phases(1).ncontrols 		= 2;
    problem.phases(1).nevents   		= 9;
    problem.phases(1).npath    		        = 0;
    problem.phases(1).nodes                     = "[60 80]";
    problem.phases(1).zero_cost_integrand       = true;

    psopt_level2_setup(problem, algorithm);


////////////////////////////////////////////////////////////////////////////
///////////////////  Declare DMatrix objects to store results //////////////
////////////////////////////////////////////////////////////////////////////

    DMatrix x, u, t;

////////////////////////////////////////////////////////////////////////////
///////////////////  Enter problem bounds information //////////////////////
////////////////////////////////////////////////////////////////////////////

    double hL     = 0.0;
    double hU     = 300000.0;
    double phiL   = DEG2RAD(-45.0);
    double phiU   = DEG2RAD( 45.0);
    double vL     = 1000.0  ;
    double vU     = 40000.0     ;
    double thetaL = DEG2RAD(-89.0) ;
    double thetaU = DEG2RAD( 89.0) ;
    double gammaL = DEG2RAD(-89.0) ;
    double gammaU = DEG2RAD(89.0)  ;
    double psiL   = DEG2RAD(-180.0);
    double psiU   = DEG2RAD(180.0) ;

    double alphaL = DEG2RAD(-89.0) ;
    double alphaU = DEG2RAD( 89.0) ;
    double betaL  = DEG2RAD(-90.0) ;
    double betaU  = DEG2RAD( 1.0)  ;

    double qU     = 70.0;
    double qL     = -INF;


    double h0      =  260000.0     ;
    double v0      =  25600.0       ;
    double phi0    =  DEG2RAD(-0.5*75.3153) ;
    double gamma0  =  DEG2RAD(-1.0);
    double theta0  =  0.0         ;
    double psi0    =  DEG2RAD(90.0) ;

    double hf      = 80000.0     ;
    double vf      = 2500.0      ;
    double gammaf  = DEG2RAD(-5.0);


    double pi = 3.141592653589793;


    problem.phases(1).bounds.lower.states(H_INDX) 	= hL;
    problem.phases(1).bounds.lower.states(PHI_INDX) 	= phiL;
    problem.phases(1).bounds.lower.states(THETA_INDX) 	= thetaL;
    problem.phases(1).bounds.lower.states(V_INDX) 	= vL;
    problem.phases(1).bounds.lower.states(GAMMA_INDX) 	= gammaL;
    problem.phases(1).bounds.lower.states(PSI_INDX) 	= psiL;


    problem.phases(1).bounds.upper.states(H_INDX) = hU;
    problem.phases(1).bounds.upper.states(PHI_INDX) = phiU;
    problem.phases(1).bounds.upper.states(THETA_INDX) = thetaU;
    problem.phases(1).bounds.upper.states(V_INDX) = vU;
    problem.phases(1).bounds.upper.states(GAMMA_INDX) = gammaU;
    problem.phases(1).bounds.upper.states(PSI_INDX) = psiU;

    problem.phases(1).bounds.lower.controls(ALPHA_INDX) = alphaL;
    problem.phases(1).bounds.upper.controls(ALPHA_INDX) = alphaU;
    problem.phases(1).bounds.lower.controls(BETA_INDX) = betaL;
    problem.phases(1).bounds.upper.controls(BETA_INDX) = betaU;


    problem.phases(1).bounds.lower.events(H_INDX) 	= h0;
    problem.phases(1).bounds.lower.events(PHI_INDX) 	= phi0;
    problem.phases(1).bounds.lower.events(THETA_INDX) 	= theta0;
    problem.phases(1).bounds.lower.events(V_INDX) 	= v0;
    problem.phases(1).bounds.lower.events(GAMMA_INDX) 	= gamma0;
    problem.phases(1).bounds.lower.events(PSI_INDX) 	= psi0;
    problem.phases(1).bounds.lower.events(7) 		= hf;
    problem.phases(1).bounds.lower.events(8) 		= vf;
    problem.phases(1).bounds.lower.events(9) 		= gammaf;


    problem.phases(1).bounds.upper.events(H_INDX) 	= h0;
    problem.phases(1).bounds.upper.events(PHI_INDX) 	= phi0;
    problem.phases(1).bounds.upper.events(THETA_INDX) 	= theta0;
    problem.phases(1).bounds.upper.events(V_INDX) 	= v0;
    problem.phases(1).bounds.upper.events(GAMMA_INDX) 	= gamma0;
    problem.phases(1).bounds.upper.events(PSI_INDX) 	= psi0;
    problem.phases(1).bounds.upper.events(7) 		= hf;
    problem.phases(1).bounds.upper.events(8) 		= vf;
    problem.phases(1).bounds.upper.events(9) 		= gammaf;


    problem.phases(1).bounds.lower.path(1)              = qL;
    problem.phases(1).bounds.upper.path(1) 		= qU;


    problem.phases(1).bounds.lower.StartTime    = 0.0;
    problem.phases(1).bounds.upper.StartTime    = 0.0;

    problem.phases(1).bounds.lower.EndTime      = 100.0;
    problem.phases(1).bounds.upper.EndTime      = 4000.0;


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

    int nodes = (int) problem.phases(1).nodes(1);

    DMatrix u_guess =  ( zeros(1,nodes+1) && (DEG2RAD(1.0)*ones(1,nodes+1) ) );
    DMatrix x_guess = zeros(6,nodes+1);
    DMatrix time_guess;

    x_guess(H_INDX,colon()) 	= linspace(h0, hf, nodes+1);
    x_guess(PHI_INDX,colon()) 	= -0.5*DEG2RAD(90.0)*ones(1,nodes+1);
    x_guess(THETA_INDX,colon()) = DEG2RAD(-89.0)*ones(1,nodes+1);
    x_guess(V_INDX,colon()) 	= linspace(v0,vf, nodes+1);
    x_guess(GAMMA_INDX,colon()) = linspace(gamma0, gammaf, nodes+1);
    x_guess(PSI_INDX,colon()) 	= linspace(pi/2, -pi/2, nodes+1);

    time_guess = linspace(0.0, 1000.0, nodes+1);

    problem.phases(1).guess.controls = u_guess;
    problem.phases(1).guess.states   = x_guess;
    problem.phases(1).guess.time     = time_guess;



////////////////////////////////////////////////////////////////////////////
///////////////////  Enter algorithm options  //////////////////////////////
////////////////////////////////////////////////////////////////////////////



    algorithm.nlp_method                  = "IPOPT";
    algorithm.nlp_iter_max                = 1000;
    algorithm.nlp_tolerance               = 5e-6;
    algorithm.scaling                     = "automatic";
    algorithm.derivatives 	          = "automatic";
    algorithm.collocation_method          = "trapezoidal";
    algorithm.mesh_refinement             = "automatic";

////////////////////////////////////////////////////////////////////////////
///////////////////  Now call PSOPT to solve the problem   /////////////////
////////////////////////////////////////////////////////////////////////////

    psopt(solution, problem, algorithm);

////////////////////////////////////////////////////////////////////////////
///////////  Extract relevant variables from solution structure   //////////
////////////////////////////////////////////////////////////////////////////

    x = solution.get_states_in_phase(1);
    u = solution.get_controls_in_phase(1);
    t = solution.get_time_in_phase(1);

////////////////////////////////////////////////////////////////////////////
///////////  Save solution data to files if desired ////////////////////////
////////////////////////////////////////////////////////////////////////////

    x.Save("x.dat");
    u.Save("u.dat");
    t.Save("t.dat");

////////////////////////////////////////////////////////////////////////////
///////////  Plot some results if desired (requires gnuplot) ///////////////
////////////////////////////////////////////////////////////////////////////

    DMatrix h     = x(1,colon());
    DMatrix phi   = x(2,colon());
    DMatrix theta = x(3,colon());
    DMatrix v     = x(4,colon());
    DMatrix gamma = x(5,colon());
    DMatrix psi   = x(6,colon());
    DMatrix alpha = u(1,colon());
    DMatrix beta  = u(2,colon());


    plot(t,x(1,colon()),problem.name, "time (s)", "x1","altitude");

    plot(t,u(1,colon()),problem.name,"time (s)", "alpha");

    plot(t,u(2,colon()),problem.name,"time (s)", "beta");

    plot(t,h, problem.name+": altitude", "time (s)", "h (ft)",
            "altitude","pdf", "shutt_alt.pdf" );

    plot(t,phi, problem.name+": longitude", "time (s)", "phi (rad)",
            "longitude","pdf", "shutt_lon.pdf" );

    plot(t,theta, problem.name+": latitude", "time (s)", "theta (rad)",
            "latitude","pdf", "shutt_lat.pdf" );

    plot(t,v, problem.name+": velocity", "time (s)", "v (ft/s)",
            "velocity","pdf", "shutt_vel.pdf" );

    plot(t,gamma, problem.name+": flight path angle", "time (s)", "gamma (rad)",
            "fpa","pdf", "shutt_fpa.pdf" );

    plot(t,psi, problem.name+": azimuth", "time (s)", "psi (rad)",
            "azi","pdf", "shutt_azi.pdf" );

    plot(t,alpha, problem.name+": alpha", "time (s)", "alpha (rad)",
            "alpha","pdf", "shutt_alpha.pdf" );

    plot(t,beta, problem.name+": beta", "time (s)", "beta (rad)",
            "beta","pdf", "shutt_beta.pdf" );


}

////////////////////////////////////////////////////////////////////////////
///////////////////////      END OF FILE     ///////////////////////////////
////////////////////////////////////////////////////////////////////////////

