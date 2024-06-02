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
  double bank_min;
  double bank_max;
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
         adouble v_exit_y = final_states[3]*cos(final_states[4])*cos(M_PI/2-final_states[5])+final_states[0]*CONSTANTS.Omega*cos(final_states[2]);
         adouble v_exit_z = final_states[3]*cos(final_states[4])*sin(M_PI/2-final_states[5]);
         adouble v_exit = sqrt(v_exit_x*v_exit_x + v_exit_y*v_exit_y+v_exit_z*v_exit_z);
         adouble gamma_exit_d_1 = final_states[3]*cos(final_states[4]);
         adouble gamma_exit_d_2 = 2*CONSTANTS.Omega*final_states[0]*final_states[3]*cos(final_states[4])*cos(M_PI/2-final_states[5])*cos(final_states[2]);
         adouble gamma_exit_d_3 = CONSTANTS.Omega*final_states[0]*cos(final_states[2]);
         adouble gamma_exit = atan(final_states[3]*sin(final_states[4])/(sqrt(gamma_exit_d_1*gamma_exit_d_1+gamma_exit_d_2+gamma_exit_d_3*gamma_exit_d_3)));
    

       
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
    double q = 1e-3;
    double r = 1e2;
    adouble runCost =q*( states[13]*states[13]+states[20]*states[20]) + r *(controls[1]*controls[1] + controls[2]*controls[2] +controls[3]*controls[3]+controls[4]*controls[4]+controls[5]*controls[5]+controls[6]*controls[6]+controls[7]*controls[7]);
    return  runCost;
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
    adouble k1  = controls[ 1 ];
    adouble k2  = controls[ 2 ];
    adouble k3  = controls[ 3 ];
    adouble k4  = controls[ 4 ];
    adouble k5  = controls[ 5 ];
    adouble k6  = controls[ 6 ];
    adouble k7  = controls[ 7 ];
   
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
    adouble g = (CONSTANTS.mu/r*r);
    
    adouble rdot     = v * sin(gamma);
    adouble thetadot = (v * cos(gamma) * sin(azim))/(r * cos(phi));
    adouble phidot   = (v*cos(gamma) * cos(azim))/r;
    adouble vdot     = -Drag-gr*sin(gamma)-g_phi*(cos(gamma*cos(azim)))+pow(CONSTANTS.Omega,2)*r*cos(phi)*(sin(gamma)*cos(phi)-cos(gamma)*sin(phi)*cos(azim));
    adouble gammadot = 1/(v)*(Lift*cos(bank)+(pow((v),2)/r-gr)* cos(gamma)+g_phi*(sin(gamma)*cos(azim))+ 2*CONSTANTS.Omega*(v)*cos(phi)*sin(azim)+pow(CONSTANTS.Omega,2)*r*cos(phi)*(cos(gamma)*cos(phi)+sin(gamma)*sin(phi)*cos(azim)));
    adouble azimdot  = 1/v*(Lift*sin(bank)/cos(gamma)+ (pow(v,2)/r)*cos(gamma)*sin(azim)*tan(phi) +g_phi*(sin(azim)/cos(gamma))-2*CONSTANTS.Omega*v*(tan(gamma)*cos(azim)*cos(phi) - sin(phi)) + (pow(CONSTANTS.Omega,2)*r/cos(gamma))*sin(azim)*cos(phi)*sin(phi));
    adouble p0dot    = 0;
 
   // Sensititvity matrix elements diff eqs
    adouble a41    = states[ 7 ];
    adouble a42    = states[ 8 ];
    adouble a43    = states[ 9 ];
    adouble a44    = states[ 10 ];
    adouble a45    = states[ 11 ];
    adouble a46    = states[ 12 ];
    adouble a47    = states[ 13 ];
    adouble a51    = states[ 14 ];
    adouble a52    = states[ 15 ];
    adouble a53    = states[ 16 ];
    adouble a54    = states[ 17 ];
    adouble a55    = states[ 18 ];
    adouble a56    = states[ 19 ];
    adouble a57    = states[ 20 ];
    
    

    double Omega    = CONSTANTS.Omega;
    double m        = CONSTANTS.m;
    double J2       = CONSTANTS.J2;
    double mu       = CONSTANTS.mu;
    double r0       = CONSTANTS.r0;
    double Cl       = CONSTANTS.cl;
    double Cd       = CONSTANTS.cd;
    double S_ref    = CONSTANTS.S_ref;
    double Re       = CONSTANTS.Re;
    double h_0      = CONSTANTS.h_0;
    double bank_min = CONSTANTS.bank_min;
    double bank_max = CONSTANTS.bank_max;
    

    adouble a41_dot= a46*((1.0/(r*r*r*r*r)*(m*(r*r*r)*(v*v)*pow(cos(gamma),2.0)*sin(azim)*sin(phi)*2.0+(Omega*Omega)*m*(r*r*r*r*r)*sin(azim)*sin(phi)*(pow(sin(phi),2.0)-1.0)*2.0-J2*(Re*Re)*m*mu*sin(azim)*sin(phi)*(pow(sin(phi),2.0)-1.0)*2.4E+1+B*Cl*S_ref*p0*(r*r*r*r*r)*(v*v)*exp(B*(Re+h_0-r))*cos(phi)*sin(bank)))/(m*v*cos(gamma)*cos(phi)*2.0)+(Cl*S_ref*k1*p0*v*exp(B*(Re+h_0-r))*cos(bank)*(bank-bank_min)*(bank-bank_max)*1.0/pow(bank_min-bank_max,2.0)*2.0)/(m*cos(gamma)))-a44*((Omega*Omega)*cos(phi)*(cos(phi)*sin(gamma)-cos(azim)*cos(gamma)*sin(phi))-mu*1.0/(r*r*r)*sin(gamma)*(J2*(Re*Re)*1.0/(r*r)*(pow(sin(phi),2.0)*(9.0/2.0)-3.0/2.0)-1.0)*2.0-J2*(Re*Re)*mu*1.0/(r*r*r*r*r)*sin(gamma)*(pow(sin(phi),2.0)*(9.0/2.0)-3.0/2.0)*2.0+(B*Cd*S_ref*p0*(v*v)*exp(B*(Re+h_0-r)))/(m*2.0)+J2*(Re*Re)*mu*1.0/(r*r*r*r*r)*cos(phi)*sin(phi)*cos(gamma*cos(azim))*1.2E+1)+a45*((cos(gamma)*(1.0/(r*r)*(v*v)+mu*1.0/(r*r*r)*(J2*(Re*Re)*1.0/(r*r)*(pow(sin(phi),2.0)*(9.0/2.0)-3.0/2.0)-1.0)*2.0+J2*(Re*Re)*mu*1.0/(r*r*r*r*r)*(pow(sin(phi),2.0)*(9.0/2.0)-3.0/2.0)*2.0)-(Omega*Omega)*cos(phi)*(cos(gamma)*cos(phi)+cos(azim)*sin(gamma)*sin(phi))+(B*Cl*S_ref*p0*(v*v)*exp(B*(Re+h_0-r))*cos(bank))/(m*2.0)+J2*(Re*Re)*mu*1.0/(r*r*r*r*r)*cos(azim)*cos(phi)*sin(gamma)*sin(phi)*1.2E+1)/v-(Cl*S_ref*k1*p0*v*exp(B*(Re+h_0-r))*sin(bank)*(bank-bank_min)*(bank-bank_max)*1.0/pow(bank_min-bank_max,2.0)*2.0)/m)+a43*1.0/(r*r)*v*cos(azim)*cos(gamma)+(a42*1.0/(r*r)*v*cos(gamma)*sin(azim))/cos(phi);
    adouble a42_dot= (Cl*S_ref*k2*p0*v*exp(B*(Re+h_0-r))*(a46*cos(bank)-a45*cos(gamma)*sin(bank))*(bank-bank_min)*(bank-bank_max)*1.0/pow(bank_min-bank_max,2.0)*2.0)/(m*cos(gamma));
    adouble a43_dot= -a46*((1.0/(r*r*r*r)*1.0/pow(cos(phi),2.0)*((r*r*r)*(v*v)*pow(cos(gamma),2.0)*sin(azim)-(Omega*Omega)*(r*r*r*r*r)*pow(cos(phi),2.0)*sin(azim)+(Omega*Omega)*(r*r*r*r*r)*pow(cos(phi),4.0)*sin(azim)*2.0-J2*(Re*Re)*mu*pow(cos(phi),2.0)*sin(azim)*3.0+J2*(Re*Re)*mu*pow(cos(phi),4.0)*sin(azim)*6.0+Omega*(r*r*r*r)*v*cos(gamma)*pow(cos(phi),3.0)*2.0+Omega*(r*r*r*r)*v*cos(azim)*pow(cos(phi),2.0)*sin(gamma)*sin(phi)*2.0))/(v*cos(gamma))-(Cl*S_ref*k3*p0*v*exp(B*(Re+h_0-r))*cos(bank)*(bank-bank_min)*(bank-bank_max)*1.0/pow(bank_min-bank_max,2.0)*2.0)/(m*cos(gamma)))+a45*((1.0/(r*r*r*r)*((Omega*Omega)*(r*r*r*r*r)*cos(azim)*sin(gamma)+(Omega*Omega)*(r*r*r*r*r)*cos(gamma)*cos(phi)*sin(phi)*2.0+J2*(Re*Re)*mu*cos(azim)*sin(gamma)*3.0+Omega*(r*r*r*r)*v*sin(azim)*sin(phi)*2.0-(Omega*Omega)*(r*r*r*r*r)*cos(azim)*pow(cos(phi),2.0)*sin(gamma)*2.0-J2*(Re*Re)*mu*cos(gamma)*cos(phi)*sin(phi)*9.0-J2*(Re*Re)*mu*cos(azim)*pow(cos(phi),2.0)*sin(gamma)*6.0))/v-(Cl*S_ref*k3*p0*v*exp(B*(Re+h_0-r))*sin(bank)*(bank-bank_min)*(bank-bank_max)*1.0/pow(bank_min-bank_max,2.0)*2.0)/m)-a44*1.0/(r*r*r*r)*((Omega*Omega)*(r*r*r*r*r)*cos(azim)*cos(gamma)+J2*(Re*Re)*mu*cos(gamma*cos(azim))*3.0-(Omega*Omega)*(r*r*r*r*r)*cos(phi)*sin(gamma)*sin(phi)*2.0-(Omega*Omega)*(r*r*r*r*r)*cos(azim)*cos(gamma)*pow(cos(phi),2.0)*2.0-J2*(Re*Re)*mu*pow(cos(phi),2.0)*cos(gamma*cos(azim))*6.0+J2*(Re*Re)*mu*cos(phi)*sin(gamma)*sin(phi)*9.0)-(a42*v*cos(gamma)*1.0/pow(cos(phi),2.0)*sin(azim)*sin(phi))/r;
    adouble a44_dot= -a45*(-1.0/(v*v)*(cos(gamma)*((v*v)/r+mu*1.0/(r*r)*(J2*(Re*Re)*1.0/(r*r)*(pow(sin(phi),2.0)*(9.0/2.0)-3.0/2.0)-1.0))+(Omega*Omega)*r*cos(phi)*(cos(gamma)*cos(phi)+cos(azim)*sin(gamma)*sin(phi))+Omega*v*cos(phi)*sin(azim)*2.0+(Cl*S_ref*p0*(v*v)*exp(B*(Re+h_0-r))*cos(bank))/(m*2.0)+J2*(Re*Re)*mu*1.0/(r*r*r*r)*cos(azim)*cos(phi)*sin(gamma)*sin(phi)*3.0)+(Omega*cos(phi)*sin(azim)*2.0+(v*cos(gamma)*2.0)/r+(Cl*S_ref*p0*v*exp(B*(Re+h_0-r))*cos(bank))/m)/v+(Cl*S_ref*k4*p0*v*exp(B*(Re+h_0-r))*sin(bank)*(bank-bank_min)*(bank-bank_max)*1.0/pow(bank_min-bank_max,2.0)*2.0)/m)-a41*sin(gamma)+a46*(-(Omega*(sin(phi)-cos(azim)*cos(phi)*tan(gamma))*2.0+(v*cos(gamma)*sin(azim)*tan(phi)*2.0)/r+(Cl*S_ref*p0*v*exp(B*(Re+h_0-r))*sin(bank))/(m*cos(gamma)))/v+1.0/(v*v)*(Omega*v*(sin(phi)-cos(azim)*cos(phi)*tan(gamma))*2.0+((v*v)*cos(gamma)*sin(azim)*tan(phi))/r+((Omega*Omega)*r*cos(phi)*sin(azim)*sin(phi))/cos(gamma)+(J2*(Re*Re)*mu*1.0/(r*r*r*r)*cos(phi)*sin(azim)*sin(phi)*3.0)/cos(gamma)+(Cl*S_ref*p0*(v*v)*exp(B*(Re+h_0-r))*sin(bank))/(m*cos(gamma)*2.0))+(Cl*S_ref*k4*p0*v*exp(B*(Re+h_0-r))*cos(bank)*(bank-bank_min)*(bank-bank_max)*1.0/pow(bank_min-bank_max,2.0)*2.0)/(m*cos(gamma)))-(a43*cos(azim)*cos(gamma))/r-(a42*cos(gamma)*sin(azim))/(r*cos(phi))+(Cd*S_ref*a44*p0*v*exp(B*(Re+h_0-r)))/m;
    adouble a45_dot= -a46*((1.0/(r*r*r*r)*1.0/pow(cos(gamma),2.0)*(Omega*m*(r*r*r*r)*v*cos(azim)*pow(cos(phi),2.0)*-4.0+(Omega*Omega)*m*(r*r*r*r*r)*pow(cos(phi),2.0)*sin(azim)*sin(gamma)*sin(phi)*2.0-m*(r*r*r)*(v*v)*pow(cos(gamma),2.0)*sin(azim)*sin(gamma)*sin(phi)*2.0+J2*(Re*Re)*m*mu*pow(cos(phi),2.0)*sin(azim)*sin(gamma)*sin(phi)*6.0+Cl*S_ref*p0*(r*r*r*r)*(v*v)*exp(B*(Re+h_0-r))*cos(phi)*sin(bank)*sin(gamma)))/(m*v*cos(phi)*2.0)-(Cl*S_ref*k5*p0*v*exp(B*(Re+h_0-r))*cos(bank)*(bank-bank_min)*(bank-bank_max)*1.0/pow(bank_min-bank_max,2.0)*2.0)/(m*cos(gamma)))+a45*((sin(gamma)*((v*v)/r+mu*1.0/(r*r)*(J2*(Re*Re)*1.0/(r*r)*(pow(sin(phi),2.0)*(9.0/2.0)-3.0/2.0)-1.0))+(Omega*Omega)*r*cos(phi)*(cos(phi)*sin(gamma)-cos(azim)*cos(gamma)*sin(phi))-J2*(Re*Re)*mu*1.0/(r*r*r*r)*cos(azim)*cos(gamma)*cos(phi)*sin(phi)*3.0)/v-(Cl*S_ref*k5*p0*v*exp(B*(Re+h_0-r))*sin(bank)*(bank-bank_min)*(bank-bank_max)*1.0/pow(bank_min-bank_max,2.0)*2.0)/m)-a44*((Omega*Omega)*r*cos(phi)*(cos(gamma)*cos(phi)+cos(azim)*sin(gamma)*sin(phi))+mu*1.0/(r*r)*cos(gamma)*(J2*(Re*Re)*1.0/(r*r)*(pow(sin(phi),2.0)*(9.0/2.0)-3.0/2.0)-1.0)+J2*(Re*Re)*mu*1.0/(r*r*r*r)*cos(azim)*cos(phi)*sin(phi)*sin(gamma*cos(azim))*3.0)-a41*v*cos(gamma)+(a43*v*cos(azim)*sin(gamma))/r+(a42*v*sin(azim)*sin(gamma))/(r*cos(phi));
    adouble a46_dot= a45*((1.0/(r*r*r*r)*cos(phi)*(Omega*(r*r*r*r)*v*cos(azim)*-2.0+(Omega*Omega)*(r*r*r*r*r)*sin(azim)*sin(gamma)*sin(phi)+J2*(Re*Re)*mu*sin(azim)*sin(gamma)*sin(phi)*3.0))/v-(Cl*S_ref*k6*p0*v*exp(B*(Re+h_0-r))*sin(bank)*(bank-bank_min)*(bank-bank_max)*1.0/pow(bank_min-bank_max,2.0)*2.0)/m)-a46*((1.0/(r*r*r*r)*((Omega*Omega)*(r*r*r*r*r)*cos(azim)*pow(cos(phi),2.0)*sin(phi)+(r*r*r)*(v*v)*cos(azim)*pow(cos(gamma),2.0)*sin(phi)+J2*(Re*Re)*mu*cos(azim)*pow(cos(phi),2.0)*sin(phi)*3.0+Omega*(r*r*r*r)*v*pow(cos(phi),2.0)*sin(azim)*sin(gamma)*2.0))/(v*cos(gamma)*cos(phi))-(Cl*S_ref*k6*p0*v*exp(B*(Re+h_0-r))*cos(bank)*(bank-bank_min)*(bank-bank_max)*1.0/pow(bank_min-bank_max,2.0)*2.0)/(m*cos(gamma)))+(a43*v*cos(gamma)*sin(azim))/r-a44*1.0/(r*r*r*r)*cos(phi)*sin(azim)*sin(phi)*((Omega*Omega)*(r*r*r*r*r)*cos(gamma)-J2*(Re*Re)*gamma*mu*sin(gamma*cos(azim))*3.0)-(a42*v*cos(azim)*cos(gamma))/(r*cos(phi));
    adouble a47_dot= -a46*((Cl*S_ref*v*exp(B*(Re+h_0-r))*sin(bank))/(m*cos(gamma)*2.0)-(Cl*S_ref*k7*p0*v*exp(B*(Re+h_0-r))*cos(bank)*(bank-bank_min)*(bank-bank_max)*1.0/pow(bank_min-bank_max,2.0)*2.0)/(m*cos(gamma)))-a45*((Cl*S_ref*v*exp(B*(Re+h_0-r))*cos(bank))/(m*2.0)+(Cl*S_ref*k7*p0*v*exp(B*(Re+h_0-r))*sin(bank)*(bank-bank_min)*(bank-bank_max)*1.0/pow(bank_min-bank_max,2.0)*2.0)/m)+(Cd*S_ref*a44*(v*v)*exp(B*(Re+h_0-r)))/(m*2.0);
    adouble a51_dot= a56*((1.0/(r*r*r*r*r)*(m*(r*r*r)*(v*v)*pow(cos(gamma),2.0)*sin(azim)*sin(phi)*2.0+(Omega*Omega)*m*(r*r*r*r*r)*sin(azim)*sin(phi)*(pow(sin(phi),2.0)-1.0)*2.0-J2*(Re*Re)*m*mu*sin(azim)*sin(phi)*(pow(sin(phi),2.0)-1.0)*2.4E+1+B*Cl*S_ref*p0*(r*r*r*r*r)*(v*v)*exp(B*(Re+h_0-r))*cos(phi)*sin(bank)))/(m*v*cos(gamma)*cos(phi)*2.0)+(Cl*S_ref*k1*p0*v*exp(B*(Re+h_0-r))*cos(bank)*(bank-bank_min)*(bank-bank_max)*1.0/pow(bank_min-bank_max,2.0)*2.0)/(m*cos(gamma)))-a54*((Omega*Omega)*cos(phi)*(cos(phi)*sin(gamma)-cos(azim)*cos(gamma)*sin(phi))-mu*1.0/(r*r*r)*sin(gamma)*(J2*(Re*Re)*1.0/(r*r)*(pow(sin(phi),2.0)*(9.0/2.0)-3.0/2.0)-1.0)*2.0-J2*(Re*Re)*mu*1.0/(r*r*r*r*r)*sin(gamma)*(pow(sin(phi),2.0)*(9.0/2.0)-3.0/2.0)*2.0+(B*Cd*S_ref*p0*(v*v)*exp(B*(Re+h_0-r)))/(m*2.0)+J2*(Re*Re)*mu*1.0/(r*r*r*r*r)*cos(phi)*sin(phi)*cos(gamma*cos(azim))*1.2E+1)+a55*((cos(gamma)*(1.0/(r*r)*(v*v)+mu*1.0/(r*r*r)*(J2*(Re*Re)*1.0/(r*r)*(pow(sin(phi),2.0)*(9.0/2.0)-3.0/2.0)-1.0)*2.0+J2*(Re*Re)*mu*1.0/(r*r*r*r*r)*(pow(sin(phi),2.0)*(9.0/2.0)-3.0/2.0)*2.0)-(Omega*Omega)*cos(phi)*(cos(gamma)*cos(phi)+cos(azim)*sin(gamma)*sin(phi))+(B*Cl*S_ref*p0*(v*v)*exp(B*(Re+h_0-r))*cos(bank))/(m*2.0)+J2*(Re*Re)*mu*1.0/(r*r*r*r*r)*cos(azim)*cos(phi)*sin(gamma)*sin(phi)*1.2E+1)/v-(Cl*S_ref*k1*p0*v*exp(B*(Re+h_0-r))*sin(bank)*(bank-bank_min)*(bank-bank_max)*1.0/pow(bank_min-bank_max,2.0)*2.0)/m)+a53*1.0/(r*r)*v*cos(azim)*cos(gamma)+(a52*1.0/(r*r)*v*cos(gamma)*sin(azim))/cos(phi);
    adouble a52_dot= (Cl*S_ref*k2*p0*v*exp(B*(Re+h_0-r))*(a56*cos(bank)-a55*cos(gamma)*sin(bank))*(bank-bank_min)*(bank-bank_max)*1.0/pow(bank_min-bank_max,2.0)*2.0)/(m*cos(gamma));
    adouble a53_dot= -a56*((1.0/(r*r*r*r)*1.0/pow(cos(phi),2.0)*((r*r*r)*(v*v)*pow(cos(gamma),2.0)*sin(azim)-(Omega*Omega)*(r*r*r*r*r)*pow(cos(phi),2.0)*sin(azim)+(Omega*Omega)*(r*r*r*r*r)*pow(cos(phi),4.0)*sin(azim)*2.0-J2*(Re*Re)*mu*pow(cos(phi),2.0)*sin(azim)*3.0+J2*(Re*Re)*mu*pow(cos(phi),4.0)*sin(azim)*6.0+Omega*(r*r*r*r)*v*cos(gamma)*pow(cos(phi),3.0)*2.0+Omega*(r*r*r*r)*v*cos(azim)*pow(cos(phi),2.0)*sin(gamma)*sin(phi)*2.0))/(v*cos(gamma))-(Cl*S_ref*k3*p0*v*exp(B*(Re+h_0-r))*cos(bank)*(bank-bank_min)*(bank-bank_max)*1.0/pow(bank_min-bank_max,2.0)*2.0)/(m*cos(gamma)))+a55*((1.0/(r*r*r*r)*((Omega*Omega)*(r*r*r*r*r)*cos(azim)*sin(gamma)+(Omega*Omega)*(r*r*r*r*r)*cos(gamma)*cos(phi)*sin(phi)*2.0+J2*(Re*Re)*mu*cos(azim)*sin(gamma)*3.0+Omega*(r*r*r*r)*v*sin(azim)*sin(phi)*2.0-(Omega*Omega)*(r*r*r*r*r)*cos(azim)*pow(cos(phi),2.0)*sin(gamma)*2.0-J2*(Re*Re)*mu*cos(gamma)*cos(phi)*sin(phi)*9.0-J2*(Re*Re)*mu*cos(azim)*pow(cos(phi),2.0)*sin(gamma)*6.0))/v-(Cl*S_ref*k3*p0*v*exp(B*(Re+h_0-r))*sin(bank)*(bank-bank_min)*(bank-bank_max)*1.0/pow(bank_min-bank_max,2.0)*2.0)/m)-a54*1.0/(r*r*r*r)*((Omega*Omega)*(r*r*r*r*r)*cos(azim)*cos(gamma)+J2*(Re*Re)*mu*cos(gamma*cos(azim))*3.0-(Omega*Omega)*(r*r*r*r*r)*cos(phi)*sin(gamma)*sin(phi)*2.0-(Omega*Omega)*(r*r*r*r*r)*cos(azim)*cos(gamma)*pow(cos(phi),2.0)*2.0-J2*(Re*Re)*mu*pow(cos(phi),2.0)*cos(gamma*cos(azim))*6.0+J2*(Re*Re)*mu*cos(phi)*sin(gamma)*sin(phi)*9.0)-(a52*v*cos(gamma)*1.0/pow(cos(phi),2.0)*sin(azim)*sin(phi))/r;
    adouble a54_dot= -a55*(-1.0/(v*v)*(cos(gamma)*((v*v)/r+mu*1.0/(r*r)*(J2*(Re*Re)*1.0/(r*r)*(pow(sin(phi),2.0)*(9.0/2.0)-3.0/2.0)-1.0))+(Omega*Omega)*r*cos(phi)*(cos(gamma)*cos(phi)+cos(azim)*sin(gamma)*sin(phi))+Omega*v*cos(phi)*sin(azim)*2.0+(Cl*S_ref*p0*(v*v)*exp(B*(Re+h_0-r))*cos(bank))/(m*2.0)+J2*(Re*Re)*mu*1.0/(r*r*r*r)*cos(azim)*cos(phi)*sin(gamma)*sin(phi)*3.0)+(Omega*cos(phi)*sin(azim)*2.0+(v*cos(gamma)*2.0)/r+(Cl*S_ref*p0*v*exp(B*(Re+h_0-r))*cos(bank))/m)/v+(Cl*S_ref*k4*p0*v*exp(B*(Re+h_0-r))*sin(bank)*(bank-bank_min)*(bank-bank_max)*1.0/pow(bank_min-bank_max,2.0)*2.0)/m)-a51*sin(gamma)+a56*(-(Omega*(sin(phi)-cos(azim)*cos(phi)*tan(gamma))*2.0+(v*cos(gamma)*sin(azim)*tan(phi)*2.0)/r+(Cl*S_ref*p0*v*exp(B*(Re+h_0-r))*sin(bank))/(m*cos(gamma)))/v+1.0/(v*v)*(Omega*v*(sin(phi)-cos(azim)*cos(phi)*tan(gamma))*2.0+((v*v)*cos(gamma)*sin(azim)*tan(phi))/r+((Omega*Omega)*r*cos(phi)*sin(azim)*sin(phi))/cos(gamma)+(J2*(Re*Re)*mu*1.0/(r*r*r*r)*cos(phi)*sin(azim)*sin(phi)*3.0)/cos(gamma)+(Cl*S_ref*p0*(v*v)*exp(B*(Re+h_0-r))*sin(bank))/(m*cos(gamma)*2.0))+(Cl*S_ref*k4*p0*v*exp(B*(Re+h_0-r))*cos(bank)*(bank-bank_min)*(bank-bank_max)*1.0/pow(bank_min-bank_max,2.0)*2.0)/(m*cos(gamma)))-(a53*cos(azim)*cos(gamma))/r-(a52*cos(gamma)*sin(azim))/(r*cos(phi))+(Cd*S_ref*a54*p0*v*exp(B*(Re+h_0-r)))/m;
    adouble a55_dot= -a56*((1.0/(r*r*r*r)*1.0/pow(cos(gamma),2.0)*(Omega*m*(r*r*r*r)*v*cos(azim)*pow(cos(phi),2.0)*-4.0+(Omega*Omega)*m*(r*r*r*r*r)*pow(cos(phi),2.0)*sin(azim)*sin(gamma)*sin(phi)*2.0-m*(r*r*r)*(v*v)*pow(cos(gamma),2.0)*sin(azim)*sin(gamma)*sin(phi)*2.0+J2*(Re*Re)*m*mu*pow(cos(phi),2.0)*sin(azim)*sin(gamma)*sin(phi)*6.0+Cl*S_ref*p0*(r*r*r*r)*(v*v)*exp(B*(Re+h_0-r))*cos(phi)*sin(bank)*sin(gamma)))/(m*v*cos(phi)*2.0)-(Cl*S_ref*k5*p0*v*exp(B*(Re+h_0-r))*cos(bank)*(bank-bank_min)*(bank-bank_max)*1.0/pow(bank_min-bank_max,2.0)*2.0)/(m*cos(gamma)))+a55*((sin(gamma)*((v*v)/r+mu*1.0/(r*r)*(J2*(Re*Re)*1.0/(r*r)*(pow(sin(phi),2.0)*(9.0/2.0)-3.0/2.0)-1.0))+(Omega*Omega)*r*cos(phi)*(cos(phi)*sin(gamma)-cos(azim)*cos(gamma)*sin(phi))-J2*(Re*Re)*mu*1.0/(r*r*r*r)*cos(azim)*cos(gamma)*cos(phi)*sin(phi)*3.0)/v-(Cl*S_ref*k5*p0*v*exp(B*(Re+h_0-r))*sin(bank)*(bank-bank_min)*(bank-bank_max)*1.0/pow(bank_min-bank_max,2.0)*2.0)/m)-a54*((Omega*Omega)*r*cos(phi)*(cos(gamma)*cos(phi)+cos(azim)*sin(gamma)*sin(phi))+mu*1.0/(r*r)*cos(gamma)*(J2*(Re*Re)*1.0/(r*r)*(pow(sin(phi),2.0)*(9.0/2.0)-3.0/2.0)-1.0)+J2*(Re*Re)*mu*1.0/(r*r*r*r)*cos(azim)*cos(phi)*sin(phi)*sin(gamma*cos(azim))*3.0)-a51*v*cos(gamma)+(a53*v*cos(azim)*sin(gamma))/r+(a52*v*sin(azim)*sin(gamma))/(r*cos(phi));
    adouble a56_dot= a55*((1.0/(r*r*r*r)*cos(phi)*(Omega*(r*r*r*r)*v*cos(azim)*-2.0+(Omega*Omega)*(r*r*r*r*r)*sin(azim)*sin(gamma)*sin(phi)+J2*(Re*Re)*mu*sin(azim)*sin(gamma)*sin(phi)*3.0))/v-(Cl*S_ref*k6*p0*v*exp(B*(Re+h_0-r))*sin(bank)*(bank-bank_min)*(bank-bank_max)*1.0/pow(bank_min-bank_max,2.0)*2.0)/m)-a56*((1.0/(r*r*r*r)*((Omega*Omega)*(r*r*r*r*r)*cos(azim)*pow(cos(phi),2.0)*sin(phi)+(r*r*r)*(v*v)*cos(azim)*pow(cos(gamma),2.0)*sin(phi)+J2*(Re*Re)*mu*cos(azim)*pow(cos(phi),2.0)*sin(phi)*3.0+Omega*(r*r*r*r)*v*pow(cos(phi),2.0)*sin(azim)*sin(gamma)*2.0))/(v*cos(gamma)*cos(phi))-(Cl*S_ref*k6*p0*v*exp(B*(Re+h_0-r))*cos(bank)*(bank-bank_min)*(bank-bank_max)*1.0/pow(bank_min-bank_max,2.0)*2.0)/(m*cos(gamma)))+(a53*v*cos(gamma)*sin(azim))/r-a54*1.0/(r*r*r*r)*cos(phi)*sin(azim)*sin(phi)*((Omega*Omega)*(r*r*r*r*r)*cos(gamma)-J2*(Re*Re)*gamma*mu*sin(gamma*cos(azim))*3.0)-(a52*v*cos(azim)*cos(gamma))/(r*cos(phi));
    adouble a57_dot= -a56*((Cl*S_ref*v*exp(B*(Re+h_0-r))*sin(bank))/(m*cos(gamma)*2.0)-(Cl*S_ref*k7*p0*v*exp(B*(Re+h_0-r))*cos(bank)*(bank-bank_min)*(bank-bank_max)*1.0/pow(bank_min-bank_max,2.0)*2.0)/(m*cos(gamma)))-a55*((Cl*S_ref*v*exp(B*(Re+h_0-r))*cos(bank))/(m*2.0)+(Cl*S_ref*k7*p0*v*exp(B*(Re+h_0-r))*sin(bank)*(bank-bank_min)*(bank-bank_max)*1.0/pow(bank_min-bank_max,2.0)*2.0)/m)+(Cd*S_ref*a54*(v*v)*exp(B*(Re+h_0-r)))/(m*2.0);


    derivatives[ 0 ] = rdot;
    derivatives[ 1 ] = thetadot;
    derivatives[ 2 ] = phidot;
    derivatives[ 3 ] = vdot;
    derivatives[ 4 ] = gammadot;
    derivatives[ 5 ] = azimdot ;
    derivatives[ 6 ] = p0dot;

    derivatives[ 7 ]  = a41_dot;  
    derivatives[ 8 ]  = a42_dot;    
    derivatives[ 9 ]  = a43_dot;    
    derivatives[ 10 ] = a44_dot;
    derivatives[ 11 ] = a45_dot;
    derivatives[ 12 ] = a46_dot;
    derivatives[ 13 ] = a47_dot;

    derivatives[ 14 ] = a51_dot;
    derivatives[ 15 ] = a52_dot;
    derivatives[ 16 ] = a53_dot;
    derivatives[ 17 ] = a54_dot;
    derivatives[ 18 ] = a55_dot;
    derivatives[ 19 ] = a56_dot;
    derivatives[ 20 ] = a57_dot;

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
         adouble v_exit_y = final_states[3]*cos(final_states[4])*cos(M_PI/2-final_states[5])+final_states[0]*CONSTANTS.Omega*cos(final_states[2]);
         adouble v_exit_z = final_states[3]*cos(final_states[4])*sin(M_PI/2-final_states[5]);
         adouble v_exit = sqrt(v_exit_x*v_exit_x + v_exit_y*v_exit_y+v_exit_z*v_exit_z);
         adouble gamma_exit_d_1 = final_states[3]*cos(final_states[4]);
         adouble gamma_exit_d_2 = 2*CONSTANTS.Omega*final_states[0]*final_states[3]*cos(final_states[4])*cos(M_PI/2-final_states[5])*cos(final_states[2]);
         adouble gamma_exit_d_3 = CONSTANTS.Omega*final_states[0]*cos(final_states[2]);
         adouble gamma_exit = atan(final_states[3]*sin(final_states[4])/(sqrt(gamma_exit_d_1*gamma_exit_d_1+gamma_exit_d_2+gamma_exit_d_3*gamma_exit_d_3)));
    
        
         e[0] = final_states[0];
         e[1] = CONSTANTS.mu/((2*CONSTANTS.mu/r_exit)-pow(v_exit,2))*(1+sqrt(1-(pow(v_exit,2)*pow(r_exit,2)*pow(cos(gamma_exit),2))/(CONSTANTS.mu*CONSTANTS.mu/((2*CONSTANTS.mu/r_exit)-pow(v_exit,2)))));
         e[2] = final_states[7]; //Sensitivity matrix final states
         e[3] = final_states[8];
         e[4] = final_states[9];
         e[5] = final_states[10];
         e[6] = final_states[11];
         e[7] = final_states[12];
         e[8] = final_states[13];

          e[9] = final_states[14]; 
         e[10] = final_states[15];
         e[11] = final_states[16];
         e[12] = final_states[17];
         e[13] = final_states[18];
         e[14] = final_states[19];
         e[15] = final_states[20];

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
    problem.nlinkages           = 22; // States and time continuity constraints

    psopt_level1_setup(problem);

/////////////////////////////////////////////////////////////////////////////
/////////   Define phase related information & do level 2 setup  ////////////
/////////////////////////////////////////////////////////////////////////////

    problem.phases(1).nstates   		= 21;
    problem.phases(1).ncontrols 		= 8;
    problem.phases(1).nevents   		= 7;

    problem.phases(2).nstates   		= 21;
    problem.phases(2).ncontrols 		= 8;
    problem.phases(2).nevents   		= 16;

    problem.phases(1).nodes      << 50; 
    problem.phases(2).nodes      << 50; 
    

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
    CONSTANTS.rho0 = 1.2092;
    CONSTANTS.Omega = 7.29211585e-5; // Earth rotation rate (rad/s);
    CONSTANTS.h_0 = 0;
    CONSTANTS.Scale_Height = 8.5*1000;
    CONSTANTS.Re = 6378 *1000;
    CONSTANTS.g0 = 9.81;
    CONSTANTS.target_ap = (200 + 6378) * 1000;
    CONSTANTS.target_pp = (200 + 6378) * 1000;
    CONSTANTS.bank_max  = 165*M_PI/180;
    CONSTANTS.bank_min  = 15*M_PI/180;


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
    double bank_max  = 165*M_PI/180;
    double k1_min   =-1e-5;
    double k1_max   = 1e-5; 
    double k2_min   =-1e-5;
    double k2_max   = 1e-5; 
    double k3_min   =-1e-5;
    double k3_max   = 1e-5; 
    double k4_min   =-1e-5;
    double k4_max   = 1e-5; 
    double k5_min   =-1e-5;
    double k5_max   = 1e-5; 
    double k6_min   =-1e-5;
    double k6_max   = 1e-5; 
    double k7_min   =-1e-5;
    double k7_max   = 1e-5; 


    double a41_min  =10*-1.01528000000000;
    double a42_min  =10*-0.431156000000000;
    double a43_min  =10*-101.275000000000;
    double a44_min  =10*-13.0960000000000;
    double a45_min  =10*-920965;
    double a46_min  =10*-5049.03000000000;
    double a47_min  =10*-10161.9000000000;
    double a51_min  =10*-4.40308000000000e-05;
    double a52_min  =10*-3.73145000000000e-05;
    double a53_min  =10*-0.00251662000000000;
    double a54_min  =10*-0.000372154000000000;
    double a55_min  =10*-32.9471000000000;
    double a56_min  =10*-0.318531000000000;
    double a57_min  =10*-0.517590000000000;
    
    double a41_max  =10*1.01528000000000;
    double a42_max  =10*0.431156000000000;
    double a43_max  =10*101.275000000000;
    double a44_max  =10*13.0960000000000;
    double a45_max  =10*920965;
    double a46_max  =10*5049.03000000000;
    double a47_max  =10*10161.9000000000;
    double a51_max  =10*4.40308000000000e-05;
    double a52_max  =10*3.73145000000000e-05;
    double a53_max  =10*0.00251662000000000;
    double a54_max  =10*0.000372154000000000;
    double a55_max  =10*32.9471000000000;
    double a56_max  =10*0.318531000000000;
    double a57_max  =10*0.517590000000000;
    
    problem.phases(iphase).bounds.lower.states   << r_min, theta_min, phi_min, v_min, gamma_min, azim_min, p0_min, a41_min, a42_min, a43_min, a44_min, a45_min, a46_min, a47_min, a51_min, a52_min, a53_min, a54_min, a55_min, a56_min, a57_min;
    problem.phases(iphase).bounds.upper.states   << r_max, theta_max, phi_max, v_max, gamma_max, azim_max, p0_max, a41_max, a42_max, a43_max, a44_max, a45_max, a46_max, a47_max, a51_max, a52_max, a53_max, a54_max, a55_max, a56_max, a57_max;

    problem.phases(iphase).bounds.lower.controls << bank_min, k1_min, k2_min, k3_min, k4_min, k5_min, k6_min, k7_min;
    problem.phases(iphase).bounds.upper.controls << bank_max, k1_max, k2_max, k3_max, k4_max, k5_max, k6_max, k7_max;

    // The following bounds fix the initial state conditions in phase 0.
    
    double r_0     = 121900 + 6378 *1000;
    double theta_0 = 4.2368;
    double phi_0   = -0.8145;
    double v_0     = 11055;
    double gamma_0 = -0.1031;
    double azim_0  = 0;
    double p0_0    = CONSTANTS.rho0;
    problem.phases(iphase).bounds.lower.events << r_0, theta_0, phi_0, v_0, gamma_0, azim_0, p0_0;  
    problem.phases(iphase).bounds.upper.events << r_0, theta_0, phi_0, v_0, gamma_0, azim_0, p0_0;
    
    problem.phases(iphase).bounds.lower.StartTime    = 0.0;
    problem.phases(iphase).bounds.upper.StartTime    = 0.0;

    problem.phases(iphase).bounds.lower.EndTime    = 25;
    problem.phases(iphase).bounds.upper.EndTime    = 175; 

    
    // Phase 2 bounds

    iphase =  2;

    double r_f = 121900 + 6378 *1000 ;
    double ra_f_min = CONSTANTS.target_ap-500;
    double ra_f_max = CONSTANTS.target_ap+500;

    bank_min  = CONSTANTS.bank_min;
    bank_max  = CONSTANTS.bank_max;


    double a41_0    = 0;
    double a42_0    = 0;
    double a43_0    = 0;
    double a44_0    = 1;
    double a45_0    = 0;
    double a46_0    = 0;
    double a47_0    = 0;
    double a51_0    = 0;
    double a52_0    = 0;
    double a53_0    = 0;
    double a54_0    = 0;
    double a55_0    = 1;
    double a56_0    = 0;
    double a57_0    = 0;
    


    problem.phases(iphase).bounds.lower.states   << r_min, theta_min, phi_min, v_min, gamma_min, azim_min, p0_min, a41_min, a42_min, a43_min, a44_min, a45_min, a46_min, a47_min, a51_min, a52_min, a53_min, a54_min, a55_min, a56_min, a57_min;
    problem.phases(iphase).bounds.upper.states   << r_max, theta_max, phi_max, v_max, gamma_max, azim_max, p0_max, a41_max, a42_max, a43_max, a44_max, a45_max, a46_max, a47_max, a51_max, a52_max, a53_max, a54_max, a55_max, a56_max, a57_max;

    problem.phases(iphase).bounds.lower.controls << bank_min, k1_min, k2_min, k3_min, k4_min, k5_min, k6_min, k7_min;
    problem.phases(iphase).bounds.upper.controls << bank_max, k1_max, k2_max, k3_max, k4_max, k5_max, k6_max, k7_max;


    problem.phases(iphase).bounds.lower.events   << r_f, ra_f_min, a41_0, a42_0, a43_0, a44_0, a45_0, a46_0, a47_0, a51_0, a52_0, a53_0, a54_0, a55_0, a56_0, a57_0;
    problem.phases(iphase).bounds.upper.events   << r_f, ra_f_max, a41_0, a42_0, a43_0, a44_0, a45_0, a46_0, a47_0, a51_0, a52_0, a53_0, a54_0, a55_0, a56_0, a57_0;   
    problem.phases(iphase).bounds.lower.StartTime    = 25;
    problem.phases(iphase).bounds.upper.StartTime    = 175;
    
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
    
    x_guess.row(7)  = 2*ones(1,nnodes);
    x_guess.row(8)  = 2*ones(1,nnodes);
    x_guess.row(9)  = 2*ones(1,nnodes);
    x_guess.row(10)  =2*ones(1,nnodes);
    x_guess.row(11)  =2*ones(1,nnodes);
    x_guess.row(12)  =2*ones(1,nnodes);
    x_guess.row(13)  =2*ones(1,nnodes);

    x_guess.row(14)  = 2* ones(1,nnodes);
    x_guess.row(15)  = 2* ones(1,nnodes);
    x_guess.row(16)  = 2* ones(1,nnodes);
    x_guess.row(17)  = 2*ones(1,nnodes);
    x_guess.row(18)  = 2*ones(1,nnodes);
    x_guess.row(19)  = 2*ones(1,nnodes);
    x_guess.row(20)  = 2*ones(1,nnodes);

    MatrixXd u_guess    =  zeros(ncontrols,nnodes);
    u_guess.row(0)  = 20*M_PI/180*ones(1,nnodes);
    u_guess.row(1)  =1e-6*ones(1,nnodes);
    u_guess.row(2)  =1e-6*ones(1,nnodes);
    u_guess.row(3)  =1e-6*ones(1,nnodes);
    u_guess.row(4)  =1e-6*ones(1,nnodes);
    u_guess.row(5)  =1e-6*ones(1,nnodes);
    u_guess.row(6)  =1e-6*ones(1,nnodes);
    u_guess.row(7)  =1e-6*ones(1,nnodes);
    problem.phases(iphase).guess.controls       = u_guess;
    problem.phases(iphase).guess.states         = x_guess;
    problem.phases(iphase).guess.time           = linspace(0.0,110.0,nnodes);
    

    iphase = 2; 

    nnodes    			            = problem.phases(iphase).nodes(0);
    ncontrols                       = problem.phases(iphase).ncontrols;
    nstates                         = problem.phases(iphase).nstates;

    u_guess    =  zeros(ncontrols,nnodes);
    u_guess.row(0)  = 165*M_PI/180*ones(1,nnodes);
    u_guess.row(1)  =1e-6*ones(1,nnodes);
    u_guess.row(2)  =1e-6*ones(1,nnodes);
    u_guess.row(3)  =1e-6*ones(1,nnodes);
    u_guess.row(4)  =1e-6*ones(1,nnodes);
    u_guess.row(5)  =1e-6*ones(1,nnodes);
    u_guess.row(6)  =1e-6*ones(1,nnodes);
    u_guess.row(7)  =1e-6*ones(1,nnodes);
    problem.phases(iphase).guess.controls       = u_guess;
    problem.phases(iphase).guess.states         = x_guess;
    problem.phases(iphase).guess.time           = linspace(110,300,nnodes);



////////////////////////////////////////////////////////////////////////////
///////////////////  Enter algorithm options  //////////////////////////////
////////////////////////////////////////////////////////////////////////////


    algorithm.nlp_iter_max                = 5000;
    algorithm.nlp_tolerance               = 1e-1;
    algorithm.nlp_method                  = "IPOPT";
    algorithm.scaling                     = "automatic";
    algorithm.derivatives                 = "automatic";
    algorithm.jac_sparsity_ratio          = 0.20;
    algorithm.collocation_method          = "Legendre";
    algorithm.diff_matrix                 = "central-differences";
    algorithm.mesh_refinement             = "automatic";
    algorithm.mr_max_increment_factor     = 0.3;
    algorithm.mr_max_iterations           = 1;
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

    x.resize(nstates, x_ph1.cols()+ x_ph2.cols());
    u.resize(ncontrols, u_ph1.cols() + u_ph2.cols());
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
     MatrixXd r, v, b, k;
    
    r = x.block(0,0,1,x.cols())/1000; 
    r.array() -= CONSTANTS.Re/1000;
    v = x.block(3,0,1,x.cols())/1000; 

    b = u.block(0,0,1,u.cols()); 
    k = u.block(1,0,7,u.cols());

    plot(t,r,problem.name, "time (s)", "Raduis (km)");

    plot(t,v,problem.name, "time (s)", "Velovity (km/sec)");

    plot(v,r,problem.name, "Velovity (km/sec)","Raduis (km)");

    plot(t,b*180/M_PI,problem.name,"time (s)", "Bank (deg)");
    
    plot(t,b*180/M_PI,problem.name,"time (s)", "u (dimensionless)", "Bank",
                               "pdf", "bank_control.pdf");

    plot(t,k,problem.name,"time (s)", "Gains");
    
    plot(t,k,problem.name,"time (s)", "u (dimensionless)", "Gains",
                               "pdf", "gains_control.pdf");

}
////////////////////////////////////////////////////////////////////////////
///////////////////////      END OF FILE     ///////////////////////////////
////////////////////////////////////////////////////////////////////////////
