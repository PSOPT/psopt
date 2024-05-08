//////////////////////////////////////////////////////////////////////////
////////////////        aerocaptureDOC.cxx         /////////////////////
//////////////////////////////////////////////////////////////////////////
////////////////           PSOPT Application         /////////////////////
//////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////
//////// Title:         DOC Aerocapture Problem      ////////////////
//////// Last modified: 27 March 2024                   ////////////////
//////////////////////////////////////////////////////////////////////////
////////     Copyright (c) Pardhasai Chadalavada        ////////////////
//////////////////////////////////////////////////////////////////////////
//////// This is part of the PSOPT software library, which////////////////
//////// is distributed under the terms of the GNU Lesser ////////////////
//////// General Public License (LGPL)                    ////////////////
//////////////////////////////////////////////////////////////////////////

#include "psopt.h"
#include <math.h>

struct Constants {
  double mu;
  double cd;
  double cl;
  double m;
  double J2;
  double r0;
  double S_ref;
  double rho0;
  double Omega;
  double h_0;
  double Scale_Height;
  double Re;
  double g0;
  double target_ap;
  double target_pp;
};

typedef struct Constants Constants_;

//////////////////////////////////////////////////////////////////////////
///////////////////  Define the end point (Mayer) cost function //////////
//////////////////////////////////////////////////////////////////////////

adouble endpoint_cost(adouble* initial_states, adouble* final_states,
                      adouble* parameters,adouble& t0, adouble& tf,
                      adouble* xad, int iphase, Workspace* workspace)
{

   Constants_& CONSTANTS = *( (Constants_ *) workspace->problem->user_data );
    adouble fval;
      adouble r_exit = final_states[0];
         adouble v_exit_x = final_states[3]*sin(final_states[4]);
         adouble v_exit_y = final_states[3]*cos(final_states[4])*cos(final_states[5]);
         adouble v_exit_z = final_states[3]*cos(final_states[4])*sin(final_states[5]);
         adouble v_exit = sqrt(v_exit_x*v_exit_x + v_exit_y*v_exit_y+v_exit_z*v_exit_z);
         adouble gamma_exit_d_1 = final_states[3]*cos(final_states[4]);
         adouble gamma_exit_d_2 = 2*CONSTANTS.Omega*final_states[0]*final_states[3]*cos(final_states[4])*cos(final_states[5])*cos(final_states[2]);
         adouble gamma_exit_d_3 = CONSTANTS.Omega*final_states[0]*cos(final_states[2]);
         adouble gamma_exit = atan(final_states[3]*final_states[4]/(sqrt(gamma_exit_d_1*gamma_exit_d_1+gamma_exit_d_2+gamma_exit_d_3*gamma_exit_d_3)));
    

       
   if (iphase < 2)
        
        fval = 0.0;

   if (iphase == 2)
        
        fval = sqrt((2*CONSTANTS.mu))*(sqrt(1/CONSTANTS.target_ap-1/(CONSTANTS.target_ap+ CONSTANTS.target_pp))-sqrt(1/CONSTANTS.target_ap-(2*CONSTANTS.mu-r_exit*pow(v_exit,2))/(2*CONSTANTS.mu*r_exit)));

   return (fval);
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

   Constants_& CONSTANTS = *( (Constants_ *) workspace->problem->user_data );

   adouble r     = states[ 0 ];
   adouble theta = states[ 1 ];
   adouble phi   = states[ 2 ];
   adouble v     = states[ 3 ];
   adouble gamma = states[ 4 ];
   adouble azim  = states[ 5 ];
   adouble p0    = states[ 6 ];

   adouble bank  = controls[ 0 ];
   
    //Get Atmospheric Density:   
    //Back out Current Altitude:
    adouble h = r - CONSTANTS.Re;

    //Current Scale Height:
    double B = 1/(CONSTANTS.Scale_Height);
    adouble p = p0 * exp(-B * (h - CONSTANTS.h_0)); //Atmospheric Denisty 
    
    //Get Aerodynamic forces:
    adouble Lift = .5 * p * pow(v,2) * CONSTANTS.cl * CONSTANTS.S_ref/CONSTANTS.m;
    adouble Drag  = .5 * p * pow(v,2) * CONSTANTS.cd * CONSTANTS.S_ref/CONSTANTS.m;
    

    //Gravitational ptential terms
    adouble gr = (CONSTANTS.mu/pow(r,2))*(1+CONSTANTS.J2*pow((CONSTANTS.r0/r),2)*(1.5-4.5*pow(sin(phi),2)));
    adouble g_phi = (CONSTANTS.mu/pow(r,2))*(CONSTANTS.J2*pow((CONSTANTS.r0/r),2)*(3*sin(phi)*cos(phi)));
    
    adouble rdot     = v * sin(gamma);
    adouble thetadot = (v * cos(gamma) * sin(azim))/(r * cos(phi));
    adouble phidot   = (v*cos(gamma) * cos(azim))/r;
    adouble vdot     = -Drag-gr*sin(gamma)-g_phi*(cos(gamma*cos(azim)))+pow(CONSTANTS.Omega,2)*r*cos(phi)*(sin(gamma)*cos(phi)-cos(gamma)*sin(phi)*cos(azim));
    adouble gammadot = 1/(v)*(Lift*cos(bank)+(pow((v),2)/r-gr)* cos(gamma)+g_phi*(sin(gamma)*cos(azim))+ 2*CONSTANTS.Omega*(v)*cos(phi)*sin(azim)+pow(CONSTANTS.Omega,2)*r*cos(phi)*(cos(gamma)*cos(phi)+sin(gamma)*sin(phi)*cos(azim)));
    adouble azimdot  = 1/v*(Lift*sin(bank)/cos(gamma)+ (pow(v,2)/r)*cos(gamma)*sin(azim)*tan(phi) +g_phi*(sin(azim)/cos(gamma))-2*CONSTANTS.Omega*v*(tan(gamma)*cos(azim)*cos(phi) - sin(phi)) + (pow(CONSTANTS.Omega,2)*r/cos(gamma))*sin(azim)*cos(phi)*sin(phi));
    adouble p0dot    = 0;
 
    derivatives[ 0 ] = rdot;
    derivatives[ 1 ] = thetadot;
    derivatives[ 2 ] = phidot;
    derivatives[ 3 ] = vdot;
    derivatives[ 4 ] = gammadot;
    derivatives[ 5 ] = azimdot ;
    derivatives[ 6 ] = p0dot;

}

////////////////////////////////////////////////////////////////////////////
///////////////////  Define the events function ////////////////////////////
////////////////////////////////////////////////////////////////////////////

void events(adouble* e, adouble* initial_states, adouble* final_states,
            adouble* parameters,adouble& t0, adouble& tf, adouble* xad,
            int iphase, Workspace* workspace)

{

   Constants_& CONSTANTS = *( (Constants_ *) workspace->problem->user_data );
  
   int j;

   if(iphase==1) {
       // These events are related to the initial state conditions in phase 1
       for(j=0;j<7;j++) e[j] = initial_states[j];
   }
   if (iphase==2) {
          adouble r_exit = final_states[0];
         adouble v_exit_x = final_states[3]*sin(final_states[4]);
         adouble v_exit_y = final_states[3]*cos(final_states[4])*cos(final_states[5]);
         adouble v_exit_z = final_states[3]*cos(final_states[4])*sin(final_states[5]);
         adouble v_exit = sqrt(v_exit_x*v_exit_x + v_exit_y*v_exit_y+v_exit_z*v_exit_z);
         adouble gamma_exit_d_1 = final_states[3]*cos(final_states[4]);
         adouble gamma_exit_d_2 = 2*CONSTANTS.Omega*final_states[0]*final_states[3]*cos(final_states[4])*cos(final_states[5])*cos(final_states[2]);
         adouble gamma_exit_d_3 = CONSTANTS.Omega*final_states[0]*cos(final_states[2]);
         adouble gamma_exit = atan(final_states[3]*final_states[4]/(sqrt(gamma_exit_d_1*gamma_exit_d_1+gamma_exit_d_2+gamma_exit_d_3*gamma_exit_d_3)));
    
        
        adouble a = CONSTANTS.mu/((2*CONSTANTS.mu/r_exit)-pow(v_exit,2));
        adouble r_a =CONSTANTS.mu/((2*CONSTANTS.mu/r_exit)-pow(v_exit,2))*(1+sqrt(1-(pow(v_exit,2)*pow(r_exit,2)*pow(cos(gamma_exit),2))/(CONSTANTS.mu*CONSTANTS.mu/((2*CONSTANTS.mu/r_exit)-pow(v_exit,2))))); //Apogee Raduis
         e[0] = final_states[0];
         e[1] = CONSTANTS.mu/((2*CONSTANTS.mu/r_exit)-pow(v_exit,2))*(1+sqrt(1-(pow(v_exit,2)*pow(r_exit,2)*pow(cos(gamma_exit),2))/(CONSTANTS.mu*CONSTANTS.mu/((2*CONSTANTS.mu/r_exit)-pow(v_exit,2)))));
    }

}



///////////////////////////////////////////////////////////////////////////
///////////////////  Define the phase linkages function ///////////////////
///////////////////////////////////////////////////////////////////////////

void linkages( adouble* linkages, adouble* xad, Workspace* workspace)
{
    int index=0;
    auto_link(linkages, &index, xad, 1, 2, workspace );
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

    problem.name        		= "DOC Aerocapture Problem";
    problem.outfilename         = "aerocaptureDOC.txt";

////////////////////////////////////////////////////////////////////////////
///////////////////  Declare an instance of Constants structure /////////////
////////////////////////////////////////////////////////////////////////////

    
    Constants_ CONSTANTS;    
    
    problem.user_data = (void*) &CONSTANTS;


////////////////////////////////////////////////////////////////////////////
////////////  Define problem level constants & do level 1 setup ////////////
////////////////////////////////////////////////////////////////////////////

    problem.nphases   			= 2;
    problem.nlinkages           = 8; // States and time continuity constraints

    psopt_level1_setup(problem);

/////////////////////////////////////////////////////////////////////////////
/////////   Define phase related information & do level 2 setup  ////////////
/////////////////////////////////////////////////////////////////////////////

    problem.phases(1).nstates   		= 7;
    problem.phases(1).ncontrols 		= 2;
    problem.phases(1).nevents   		= 7;

    problem.phases(2).nstates   		= 7;
    problem.phases(2).ncontrols 		= 2;
    problem.phases(2).nevents   		= 2;

    problem.phases(1).nodes      << 400; 
    problem.phases(2).nodes      << 400; 
    

    psopt_level2_setup(problem, algorithm);

////////////////////////////////////////////////////////////////////////////
///////////////////  Declare MatrixXd objects to store results //////////////
////////////////////////////////////////////////////////////////////////////

    MatrixXd x, u, t, H;

////////////////////////////////////////////////////////////////////////////
///////////////////  Initialize CONSTANTS and //////////////////////////////
///////////////////  declare local variables  //////////////////////////////
////////////////////////////////////////////////////////////////////////////


    CONSTANTS.mu = 3.986012e14;       // Gravitational parameter (m^3/s^2)
    CONSTANTS.cd = 1.37;
    CONSTANTS.cl = 0.37;
    CONSTANTS.m = 8983;
    CONSTANTS.J2 = 0.00108263;
    CONSTANTS.r0 = 6378 *1000;
    CONSTANTS.S_ref = 19.86;
    CONSTANTS.rho0 = 3;
    CONSTANTS.Omega = 7.29211585e-5; // Earth rotation rate (rad/s);
    CONSTANTS.h_0 = 0;
    CONSTANTS.Scale_Height = 8.5*1000;
    CONSTANTS.Re = 6378 *1000;
    CONSTANTS.g0 = 9.81;
    CONSTANTS.target_ap = (200 + 6378) * 1000;
    CONSTANTS.target_pp = (200 + 6378) * 1000;

////////////////////////////////////////////////////////////////////////////
///////////////////  Enter problem bounds information //////////////////////
////////////////////////////////////////////////////////////////////////////
    
    int iphase;
    
    // Phase 1 bounds

    iphase =  1;
    double r_min     = 6378 *1000;
    double theta_min = -2*M_PI;
    double phi_min   = -2*M_PI;
    double v_min     = 5000;
    double gamma_min = -2*M_PI;
    double azim_min  = -2*M_PI;
    double p0_min    = CONSTANTS.rho0;

    double r_max     = 1.5 * 121900 + 6378 *1000; 
    double theta_max = 2*M_PI;
    double phi_max   = 2*M_PI;
    double v_max     = 1.5 * 11055;
    double gamma_max = 2*M_PI;
    double azim_max  = 2*M_PI;
    double p0_max    = CONSTANTS.rho0;

    double bank_min  = 15*M_PI/180;
    double bank_max  = 90*M_PI/180;
    double k1_min  = 0;
    double k1_max  = 0;

    problem.phases(iphase).bounds.lower.states   << r_min, theta_min, phi_min, v_min, gamma_min, azim_min, p0_min;
    problem.phases(iphase).bounds.upper.states   << r_max, theta_max, phi_max, v_max, gamma_max, azim_max, p0_max;

    problem.phases(iphase).bounds.lower.controls << bank_min, k1_min;
    problem.phases(iphase).bounds.upper.controls << bank_max, k1_max;

    // The following bounds fix the initial state conditions in phase 0.
    
    double r_0     = 121900 + 6378 *1000;
    double theta_0 = 4.2368;
    double phi_0   = -0.8145;
    double v_0     = 11055;
    double gamma_0 = -0.097;// -0.1031;
    double azim_0  = 0;
    double p0_0    = CONSTANTS.rho0;
    problem.phases(iphase).bounds.lower.events << r_0, theta_0, phi_0, v_0, gamma_0, azim_0, p0_0;  
    problem.phases(iphase).bounds.upper.events << r_0, theta_0, phi_0, v_0, gamma_0, azim_0, p0_0;
    
    problem.phases(iphase).bounds.lower.StartTime    = 0.0;
    problem.phases(iphase).bounds.upper.StartTime    = 0.0;

    problem.phases(iphase).bounds.lower.EndTime    = 85.0;
    problem.phases(iphase).bounds.upper.EndTime    = 125.0; 

    // problem.bounds.lower.linkage(0)= -1e-3;
    // problem.bounds.upper.linkage(0)= 1e-3;
    
    // Phase 2 bounds

    iphase =  2;

    double r_f = 121900 + 6378 *1000 ;
    double ra_f_min = CONSTANTS.target_ap-5;
    double ra_f_max = CONSTANTS.target_ap+5;

    bank_min  = 90*M_PI/180;
    bank_max  = 165*M_PI/180; 


    problem.phases(iphase).bounds.lower.states   << r_min, theta_min, phi_min, v_min, gamma_min, azim_min, p0_min;
    problem.phases(iphase).bounds.upper.states   << r_max, theta_max, phi_max, v_max, gamma_max, azim_max, p0_max;


    problem.phases(iphase).bounds.lower.controls << bank_min, k1_min;;
    problem.phases(iphase).bounds.upper.controls << bank_max, k1_max;;


    problem.phases(iphase).bounds.lower.events   << r_f, ra_f_min;
    problem.phases(iphase).bounds.upper.events   << r_f, ra_f_max;   

    problem.phases(iphase).bounds.lower.StartTime    = 85.0;
    problem.phases(iphase).bounds.upper.StartTime    = 125.0;
    
    problem.phases(iphase).bounds.lower.EndTime    = 200.0;
    problem.phases(iphase).bounds.upper.EndTime    = 400.0; 
    

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
    iphase = 1; 

    int nnodes    			            = problem.phases(iphase).nodes(0);
    int ncontrols                       = problem.phases(iphase).ncontrols;
    int nstates                         = problem.phases(iphase).nstates;

    MatrixXd x_guess    =  zeros(nstates,nnodes);

    x_guess.row(0)  = r_0*ones(1,nnodes);
    x_guess.row(1)  = theta_0*ones(1,nnodes);
    x_guess.row(2)  = phi_0*ones(1,nnodes);
    x_guess.row(3)  = v_0*ones(1,nnodes);
    x_guess.row(4)  = gamma_0*ones(1,nnodes);
    x_guess.row(5)  = azim_0*ones(1,nnodes);
    x_guess.row(6)  = p0_0*ones(1,nnodes);

    problem.phases(iphase).guess.controls       = 20*M_PI/180*ones(ncontrols,nnodes);
    problem.phases(iphase).guess.states         = x_guess;
    problem.phases(iphase).guess.time           = linspace(0.0,110.0,nnodes);
    

    iphase = 2; 

    nnodes    			            = problem.phases(iphase).nodes(0);
    ncontrols                       = problem.phases(iphase).ncontrols;
    nstates                         = problem.phases(iphase).nstates;

    x_guess    =  zeros(nstates,nnodes);

    x_guess.row(0)  = r_0*ones(1,nnodes);
    x_guess.row(1)  = theta_0*ones(1,nnodes);
    x_guess.row(2)  = phi_0*ones(1,nnodes);
    x_guess.row(3)  = v_0*ones(1,nnodes);
    x_guess.row(4)  = gamma_0*ones(1,nnodes);
    x_guess.row(5)  = azim_0*ones(1,nnodes);
    x_guess.row(6)  = p0_0*ones(1,nnodes);

    problem.phases(iphase).guess.controls       = 90*M_PI/180*ones(ncontrols,nnodes);
    problem.phases(iphase).guess.states         = x_guess;
    problem.phases(iphase).guess.time           = linspace(110,300,nnodes);


////////////////////////////////////////////////////////////////////////////
///////////////////  Enter algorithm options  //////////////////////////////
////////////////////////////////////////////////////////////////////////////


    algorithm.nlp_iter_max                = 3000;
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

// ////////////////////////////////////////////////////////////////////////////
// ///////////  Extract relevant variables from solution structure   //////////
// ////////////////////////////////////////////////////////////////////////////


    MatrixXd x_ph1, x_ph2, u_ph1, u_ph2;
    MatrixXd t_ph1, t_ph2;

    x_ph1 = solution.get_states_in_phase(1);
    x_ph2 = solution.get_states_in_phase(2);
    
    u_ph1 = solution.get_controls_in_phase(1);
    u_ph2 = solution.get_controls_in_phase(2);
    
    t_ph1 = solution.get_time_in_phase(1);
    t_ph2 = solution.get_time_in_phase(2);

    x.resize(7, x_ph1.cols()+ x_ph2.cols());
    u.resize(2, u_ph1.cols() + u_ph2.cols());
    t.resize(1, t_ph1.cols()+ t_ph2.cols());

    x << x_ph1, x_ph2; 
    u << u_ph1, u_ph2;
    t << t_ph1, t_ph2;

// ////////////////////////////////////////////////////////////////////////////
// ///////////  Save solution data to files if desired ////////////////////////
// ////////////////////////////////////////////////////////////////////////////

    Save(x, "x.dat");
    Save(u,"u.dat");
    Save(t,"t.dat");

// ////////////////////////////////////////////////////////////////////////////
// ///////////  Plot some results if desired (requires gnuplot) ///////////////
// ////////////////////////////////////////////////////////////////////////////
     MatrixXd r, v;
    
    r = x.block(0,0,1,x.cols())/1000; 
    r.array() -= CONSTANTS.Re/1000;
    v = x.block(3,0,1,x.cols())/1000; 

    plot(t,r,problem.name, "time (s)", "Raduis (km)");

    plot(t,v,problem.name, "time (s)", "Velovity (km/sec)");

    plot(v,r,problem.name, "Velovity (km/sec)","Raduis (km)");

    plot(t,u*180/M_PI,problem.name,"time (s)", "Bank (deg)");
    
    plot(t,u*180/M_PI,problem.name,"time (s)", "u (dimensionless)", "Bank",
                               "pdf", "bank_control.pdf");

}
////////////////////////////////////////////////////////////////////////////
///////////////////////      END OF FILE     ///////////////////////////////
////////////////////////////////////////////////////////////////////////////
