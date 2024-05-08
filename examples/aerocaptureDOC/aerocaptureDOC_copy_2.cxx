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
#include <iostream>
#include <boost/numeric/odeint.hpp>
#include <fstream>

using namespace std;
using namespace boost::numeric::odeint;

// Define your state type
using state_type = std::vector<double>;

// Function to convert a vector of vectors into an Eigen matrix
int vectorToinit(const std::vector<state_type>& vec, const std::vector<state_type>& vec2, Eigen::MatrixXd *x_guess_1, Eigen::MatrixXd *x_guess_2 ) {
    

    size_t rows = vec.size();
    size_t cols = vec[0].size(); // Assuming all state vectors have the same size
    size_t rows2 = vec2.size();
    size_t cols2 = vec2[0].size(); // Assuming all state vectors have the same size

    Eigen::MatrixXd mat(4, 4);
    Eigen::MatrixXd mat2(4, 4);
    Eigen::MatrixXd matf(4, 4);
    Eigen::MatrixXd lam(4, 4);
    Eigen::MatrixXd lam2(4, 4);
    
    for (size_t i = 0; i < 4; ++i) {
        for (size_t j = 0; j < 4; ++j) {
            matf(i, j) = vec2[rows2-1][7+4*i+j];
        }
    }

    for (size_t l = 0; l < rows; ++l) {
        for (size_t i = 0; i < 4; ++i) {
            for (size_t j = 0; j < 4; ++j) {
                mat(i, j) = vec[l][7+4*i+j];
            }
        }
         lam = matf * mat.inverse();
        

        (*x_guess_1)(0,l)  = vec[l][0];
        (*x_guess_1)(1,l)  = vec[l][1];
        (*x_guess_1)(2,l)  = vec[l][2];
        (*x_guess_1)(3,l)  = vec[l][3];
        (*x_guess_1)(4,l)  = vec[l][4];
        (*x_guess_1)(5,l)  = vec[l][5];
        (*x_guess_1)(6,l)  = vec[l][6];
        (*x_guess_1)(7,l)   = lam(0,0);
        (*x_guess_1)(8,l)   = lam(0,1);
        (*x_guess_1)(9,l)   = lam(0,2);
        (*x_guess_1)(10,l)  = lam(0,3);
        (*x_guess_1)(11,l)  = lam(1,0);
        (*x_guess_1)(12,l)  = lam(1,1);
        (*x_guess_1)(13,l)  = lam(1,2);
        (*x_guess_1)(14,l)  = lam(1,3);
        (*x_guess_1)(15,l)  = lam(2,0);
        (*x_guess_1)(16,l)  = lam(2,1);
        (*x_guess_1)(17,l)  = lam(2,2);
        (*x_guess_1)(18,l)  = lam(2,3);
    }

    for (size_t l = 0; l < rows2; ++l) {
        for (size_t i = 0; i < 4; ++i) {
            for (size_t j = 0; j < 4; ++j) {
                mat2(i, j) = vec2[l][7+4*i+j];
            }
        }
    
        lam2 = matf * mat2.inverse();

        // if (l==rows2-1){
        //     std::ofstream outFile3("m_f.dat");
        //     if (outFile3.is_open()) {
        //         outFile3 <<  matf;
        //         outFile3 <<  endl;
        //         outFile3 <<  lam2;
        //         outFile3 <<  endl;
        //         outFile3.close();
        //     }
        // }

        (*x_guess_2)(0,l)  = vec2[l][0];
        (*x_guess_2)(1,l)  = vec2[l][1];
        (*x_guess_2)(2,l)  = vec2[l][2];
        (*x_guess_2)(3,l)  = vec2[l][3];
        (*x_guess_2)(4,l)  = vec2[l][4];
        (*x_guess_2)(5,l)  = vec2[l][5];
        (*x_guess_2)(6,l)  = vec2[l][6];
        (*x_guess_2)(7,l)   = lam2(0,0);
        (*x_guess_2)(8,l)   = lam2(0,1);
        (*x_guess_2)(9,l)   = lam2(0,2);
        (*x_guess_2)(10,l)  = lam2(0,3);
        (*x_guess_2)(11,l)  = lam2(1,0);
        (*x_guess_2)(12,l)  = lam2(1,1);
        (*x_guess_2)(13,l)  = lam2(1,2);
        (*x_guess_2)(14,l)  = lam2(1,3);
        (*x_guess_2)(15,l)  = lam2(2,0);
        (*x_guess_2)(16,l)  = lam2(2,1);
        (*x_guess_2)(17,l)  = lam2(2,2);
        (*x_guess_2)(18,l)  = lam2(2,3);
    }

    return 0;
}

struct push_back_state_and_time
{
    typedef vector<double> state_type;
    std::vector< state_type >& m_states;
    std::vector< double >& m_times;

    push_back_state_and_time( std::vector< state_type > &states , std::vector< double > &times )
    : m_states( states ) , m_times( times ) { }

    void operator()( const state_type &x , double t )
    {
        m_states.push_back( x );
        m_times.push_back( t );

        // cout << "Time: " << t << ", Solution: ";
        // cout << endl;
    }
};


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

typedef struct Constants Constants_;

void rhs(const vector<double>& y, vector<double>& dydt, const double t, Constants_& CONSTANTS, double bank) {

   double r     = y[ 0 ];
   double theta = y[ 1 ];
   double phi   = y[ 2 ];
   double v     = y[ 3 ];
   double gamma = y[ 4 ];
   double azim  = y[ 5 ];
   double p0    = y[ 6 ];

    //Get Atmospheric Density:   
    //Back out Current Altitude:
    double h = r - CONSTANTS.Re;

    //Current Scale Height:
    double B = 1/(CONSTANTS.Scale_Height);
    double p = p0 * exp(-B * (h - CONSTANTS.h_0)); //Atmospheric Denisty 
    
    //Get Aerodynamic forces:
    double Lift = .5 * p * pow(v,2) * CONSTANTS.cl * CONSTANTS.S_ref/CONSTANTS.m;
    double Drag  = .5 * p * pow(v,2) * CONSTANTS.cd * CONSTANTS.S_ref/CONSTANTS.m;
    

    //Gravitational ptential terms
    double gr = (CONSTANTS.mu/pow(r,2))*(1+CONSTANTS.J2*pow((CONSTANTS.r0/r),2)*(1.5-4.5*pow(sin(phi),2)));
    double g_phi = (CONSTANTS.mu/pow(r,2))*(CONSTANTS.J2*pow((CONSTANTS.r0/r),2)*(3*sin(phi)*cos(phi)));
    double g = (CONSTANTS.mu/r*r);
    
    double rdot     = v * sin(gamma);
    double thetadot = (v * cos(gamma) * sin(azim))/(r * cos(phi));
    double phidot   = (v*cos(gamma) * cos(azim))/r;
    double vdot     = -Drag-gr*sin(gamma)-g_phi*(cos(gamma*cos(azim)))+pow(CONSTANTS.Omega,2)*r*cos(phi)*(sin(gamma)*cos(phi)-cos(gamma)*sin(phi)*cos(azim));
    double gammadot = 1/(v)*(Lift*cos(bank)+(pow((v),2)/r-gr)* cos(gamma)+g_phi*(sin(gamma)*cos(azim))+ 2*CONSTANTS.Omega*(v)*cos(phi)*sin(azim)+pow(CONSTANTS.Omega,2)*r*cos(phi)*(cos(gamma)*cos(phi)+sin(gamma)*sin(phi)*cos(azim)));
    double azimdot  = 1/v*(Lift*sin(bank)/cos(gamma)+ (pow(v,2)/r)*cos(gamma)*sin(azim)*tan(phi) +g_phi*(sin(azim)/cos(gamma))-2*CONSTANTS.Omega*v*(tan(gamma)*cos(azim)*cos(phi) - sin(phi)) + (pow(CONSTANTS.Omega,2)*r/cos(gamma))*sin(azim)*cos(phi)*sin(phi));
    double p0dot    = 0;

    // Sensititvity matrix elements diff eqs
    double a11    = y[ 7 ];
    double a12    = y[ 8 ];
    double a13    = y[ 9 ];
    double a14    = y[ 10 ];

    double a21    = y[ 11 ];
    double a22    = y[ 12 ];
    double a23    = y[ 13 ];
    double a24    = y[ 14 ];

    double a31    = y[ 15 ];
    double a32    = y[ 16 ];
    double a33    = y[ 17 ];
    double a34    = y[ 18 ];
    

    double a41    = y[ 19 ];
    double a42    = y[ 20 ];
    double a43    = y[ 21 ];
    double a44    = y[ 22 ];
    

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
    double k = 1e-2;
    double bank_min = CONSTANTS.bank_min;
    double bank_max = CONSTANTS.bank_max;
    double eta      = ((bank_max-bank)*(bank-bank_min))/((bank_max-bank_min)*(bank_max-bank_min));
    


    double s11 = 0;
    double s12 = sin(gamma);
    double s13 = v*cos(gamma);
    double s14 = 0;
    // Done

    double s41 = (2*mu*sin(gamma))/(r * r * r) + (B*Cd*S_ref*p0*v * v*exp(B*(Re + h_0 - r)))/(2*m);
    double s42 = -(Cd*S_ref*p0*v*exp(B*(Re + h_0 - r)))/m;
    double s43 = -(mu*cos(gamma))/(r*r);
    double s44 = -(Cd * S_ref * v * v * exp(B * (Re + h_0 - r))) / (2 * m);

    double s51 = - (cos(gamma)*(v*v/(r*r) - (2*mu)/(r*r*r)) + (B*Cl*S_ref*p0*v*v*exp(B*(Re + h_0 - r))*cos(bank))/(2*m))/v - eta*k*(Cl*S_ref*p0*v*exp(B*(Re + h_0 - r))*sin(bank))/(2*m);
    double s52 =(2*m*mu*cos(gamma) + 2*m*r*v*v*cos(gamma) + Cl*S_ref*p0*r*r*v*v*exp(B*(Re + h_0 - r))*cos(bank) -  eta*k*Cl*S_ref*p0*r*r*v*v*v*exp(B*(Re + h_0 - r))*sin(bank))/(2*m*r*r*v*v);
    double s53 = - (sin(gamma)*(v*v/r - mu/(r*r)))/v -  eta*k* (Cl*S_ref*p0*v*exp(B*(Re + h_0 - r))*sin(bank))/(2*m);
    double s54 = (Cl*S_ref*v*exp(B*(Re + h_0 - r))*(cos(bank)))/(2*m);

    
    double s71 = 0;
    double s72 = 0;
    double s73 = 0;
    double s74 = 0;
    
    int size = 4;
    MatrixXd a_matrix = zeros(size, size);
    MatrixXd s_matrix = zeros(size, size);
    MatrixXd aDot_matrix1 = zeros(size, size);
    
    // Initialize matrices with some values (you can use your own initialization)
    a_matrix << a11, a12, a13, a14, 
                a21, a22, a23, a24, 
                a31, a32, a33, a34, 
                a41, a42, a43, a44; 
    
    s_matrix << s11, s12, s13, s14, 
                s41, s42, s43, s44, 
                s51, s52, s53, s54, 
                s71, s72, s73, s74; 

    // Perform matrix multiplication
    aDot_matrix1 = a_matrix * s_matrix;
    
    dydt[ 0 ] = rdot;
    dydt[ 1 ] = thetadot;
    dydt[ 2 ] = phidot;
    dydt[ 3 ] = vdot;
    dydt[ 4 ] = gammadot;
    dydt[ 5 ] = azimdot ;
    dydt[ 6 ] = p0dot;

    dydt[ 7 ]  = aDot_matrix1(0,0);  
    dydt[ 8 ]  = aDot_matrix1(0,1);    
    dydt[ 9 ]  = aDot_matrix1(0,2);    
    dydt[ 10 ] = aDot_matrix1(0,3);

    dydt[ 11 ] = aDot_matrix1(1,0);
    dydt[ 12 ] = aDot_matrix1(1,1);
    dydt[ 13 ] = aDot_matrix1(1,2);
    dydt[ 14 ] = aDot_matrix1(1,3);


 
    dydt[ 15 ] = aDot_matrix1(2,0);
    dydt[ 16 ] = aDot_matrix1(2,1);
    dydt[ 17 ] = aDot_matrix1(2,2);
    dydt[ 18 ] = aDot_matrix1(2,3);



    dydt[ 19 ] = aDot_matrix1(3,0);
    dydt[ 20 ] = aDot_matrix1(3,1);
    dydt[ 21 ] = aDot_matrix1(3,2);
    dydt[ 22 ] = aDot_matrix1(3,3);

}

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
    adouble v_exit = final_states[3];
    adouble gamma_exit = final_states[4];

    adouble s_r_gam0 = initial_states[9];
    adouble s_v_gam0 = initial_states[13];
    adouble s_gam_gam0 = initial_states[17];
    double l = 1e-2;

       
   if (iphase < 2)
        
        fval = l*(s_r_gam0*s_r_gam0 + s_v_gam0*s_v_gam0) + 1e2*s_gam_gam0*s_gam_gam0;

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
    adouble k1  = controls[ 1 ];
    adouble k2  = controls[ 2 ];
    adouble k3  = controls[ 3 ];
   
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
 
    // Sensititvity matrix elements diff eqs
    adouble a11    = states[ 7 ];
    adouble a14    = states[ 8 ];
    adouble a15    = states[ 9 ];
    adouble a17    = states[ 10 ];
    adouble a41    = states[ 11 ];
    adouble a44    = states[ 12 ];
    adouble a45    = states[ 13 ];
    adouble a47    = states[ 14 ];
    adouble a51    = states[ 15 ];
    adouble a54    = states[ 16 ];
    adouble a55    = states[ 17 ];
    adouble a57    = states[ 18 ];


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
    adouble eta      = ((bank_max-bank)*(bank-bank_min))/((bank_max-bank_min)*(bank_max-bank_min));
    


    adouble s11 = 0;
    adouble s14 = sin(gamma);
    adouble s15 = v*cos(gamma);
    adouble s17 = 0;
    // Done


    adouble s41 = (2*mu*sin(gamma))/(r * r * r) + (B*Cd*S_ref*p0*v * v*exp(B*(Re + h_0 - r)))/(2*m);
    adouble s44 = -(Cd*S_ref*p0*v*exp(B*(Re + h_0 - r)))/m;
    adouble s45 = -(mu*cos(gamma))/(r*r);
    adouble s47 = -(Cd * S_ref * v * v * exp(B * (Re + h_0 - r))) / (2 * m);

    adouble s51 = - (cos(gamma)*(v*v/(r*r) - (2*mu)/(r*r*r)) + (B*Cl*S_ref*p0*v*v*exp(B*(Re + h_0 - r))*cos(bank))/(2*m))/v- eta*k1*(Cl*S_ref*p0*v*exp(B*(Re + h_0 - r))*sin(bank))/(2*m);
    adouble s54 =(2*m*mu*cos(gamma) + 2*m*r*v*v*cos(gamma) + Cl*S_ref*p0*r*r*v*v*exp(B*(Re + h_0 - r))*cos(bank) - eta*k2*Cl*S_ref*p0*r*r*v*v*v*exp(B*(Re + h_0 - r))*sin(bank))/(2*m*r*r*v*v);
    adouble s55 = - (sin(gamma)*(v*v/r - mu/(r*r)))/v- eta*k3*(Cl*S_ref*p0*v*exp(B*(Re + h_0 - r))*sin(bank))/(2*m);
    adouble s57 = (Cl*S_ref*v*exp(B*(Re + h_0 - r))*(cos(bank)))/(2*m);

    adouble s71 = 0;
    adouble s74 = 0;
    adouble s75 = 0;
    adouble s77 = 0;



    adouble a11_dot    = -1*(a11*s11+a14*s41+a15*s51+a17*s71);
    adouble a14_dot    = -1*(a11*s14+a14*s44+a15*s54+a17*s74);
    adouble a15_dot    = -1*(a11*s15+a14*s45+a15*s55+a17*s75);
    adouble a17_dot    = -1*(a11*s17+a14*s47+a15*s57+a17*s77);

    adouble a41_dot    = -1*(a41*s11+a44*s41+a45*s51+a47*s71);
    adouble a44_dot    = -1*(a41*s14+a44*s44+a45*s54+a47*s74);
    adouble a45_dot    = -1*(a41*s15+a44*s45+a45*s55+a47*s75);
    adouble a47_dot    = -1*(a41*s17+a44*s47+a45*s57+a47*s77);
    
    adouble a51_dot    = -1*(a51*s11+a54*s41+a55*s51+a57*s71);
    adouble a54_dot    = -1*(a51*s14+a54*s44+a55*s54+a57*s74);
    adouble a55_dot    = -1*(a51*s15+a54*s45+a55*s55+a57*s75);
    adouble a57_dot    = -1*(a51*s17+a54*s47+a55*s57+a57*s77);

    derivatives[ 0 ] = rdot;
    derivatives[ 1 ] = thetadot;
    derivatives[ 2 ] = phidot;
    derivatives[ 3 ] = vdot;
    derivatives[ 4 ] = gammadot;
    derivatives[ 5 ] = azimdot ;
    derivatives[ 6 ] = p0dot;

    derivatives[ 7 ]  = a11_dot;  
    derivatives[ 8 ]  = a14_dot;    
    derivatives[ 9 ]  = a15_dot;    
    derivatives[ 10 ] = a17_dot;
    derivatives[ 11 ] = a41_dot;
    derivatives[ 12 ] = a44_dot;
    derivatives[ 13 ] = a45_dot;
    derivatives[ 14 ] = a47_dot;
    derivatives[ 15 ] = a51_dot;
    derivatives[ 16 ] = a54_dot;
    derivatives[ 17 ] = a55_dot;
    derivatives[ 18 ] = a57_dot;

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
        adouble v_exit = final_states[3];
        adouble gamma_exit = final_states[4];
        
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
    problem.nlinkages           = 20; // States and time continuity constraints

    psopt_level1_setup(problem);

/////////////////////////////////////////////////////////////////////////////
/////////   Define phase related information & do level 2 setup  ////////////
/////////////////////////////////////////////////////////////////////////////

    problem.phases(1).nstates   		= 19;
    problem.phases(1).ncontrols 		= 4;
    problem.phases(1).nevents   		= 7;

    problem.phases(2).nstates   		= 19;
    problem.phases(2).ncontrols 		= 4;
    problem.phases(2).nevents   		= 14;

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
    double p0_min    = -2*CONSTANTS.rho0;

    double r_max     = 1.5 * 121900 + 6378 *1000; 
    double theta_max = 2*M_PI;
    double phi_max   = 2*M_PI;
    double v_max     = 1.5 * 11055;
    double gamma_max = 2*M_PI;
    double azim_max  = 2*M_PI;
    double p0_max    = 2*CONSTANTS.rho0;

    double bank_min  = 15*M_PI/180;
    double bank_max  = 90*M_PI/180;
    double k1_min  = -1;
    double k1_max  =  1;
    double k2_min  = -1;
    double k2_max  =  1;
    double k3_min  = -1;
    double k3_max  =  1;

    double a11_min  =  -1e7;
    double a14_min  =  -1e7;
    double a15_min  =  -1e7;
    double a17_min  =  -1e7;
    double a41_min  =  -1e7;
    double a44_min  =  -1e7;
    double a45_min  =  -1e7;
    double a47_min  =  -1e7;
    double a51_min  =  -1e7;
    double a54_min  =  -1e7;
    double a55_min  =  -1e7;
    double a57_min  =  -1e7;

    double a11_max  =  1e7;
    double a14_max  =  1e7;
    double a15_max  =  1e7;
    double a17_max  =  1e7;
    double a41_max  =  1e7;
    double a44_max  =  1e7;
    double a45_max  =  1e7;
    double a47_max  =  1e7;
    double a51_max  =  1e7;
    double a54_max  =  1e7;
    double a55_max  =  1e7;
    double a57_max  =  1e7;

    problem.phases(iphase).bounds.lower.states   << r_min, theta_min, phi_min, v_min, gamma_min, azim_min, p0_min, a11_min, a14_min, a15_min, a17_min, a41_min, a44_min, a45_min, a47_min, a51_min, a54_min, a55_min, a57_min;
    problem.phases(iphase).bounds.upper.states   << r_max, theta_max, phi_max, v_max, gamma_max, azim_max, p0_max, a11_max, a14_max, a15_max, a17_max, a41_max, a44_max, a45_max, a47_max, a51_max, a54_max, a55_max, a57_max;

    problem.phases(iphase).bounds.lower.controls << bank_min, k1_min, k2_min, k3_min;
    problem.phases(iphase).bounds.upper.controls << bank_max, k1_max, k2_max, k3_max;

    // The following bounds fix the initial state conditions in phase 0.
    
    double r_0     = 121900 + 6378 *1000;
    double theta_0 = 4.2368;
    double phi_0   = -0.8145;
    double v_0     = 11055;
    double gamma_0 = -0.1031;
    double azim_0  = 0;
    double p0_0    = 1.2092;
    problem.phases(iphase).bounds.lower.events << r_0, theta_0, phi_0, v_0, gamma_0, azim_0, p0_0;  
    problem.phases(iphase).bounds.upper.events << r_0, theta_0, phi_0, v_0, gamma_0, azim_0, p0_0;
    
    problem.phases(iphase).bounds.lower.StartTime    = 0.0;
    problem.phases(iphase).bounds.upper.StartTime    = 0.0;

    problem.phases(iphase).bounds.lower.EndTime    = 85.0;
    problem.phases(iphase).bounds.upper.EndTime    = 150.0; 

    
    // Phase 2 bounds

    iphase =  2;

    double r_f = 121900 + 6378 *1000 ;
    double ra_f_min = CONSTANTS.target_ap-5;
    double ra_f_max = CONSTANTS.target_ap+5;

    bank_min  = CONSTANTS.bank_min;
    bank_max  = CONSTANTS.bank_max;

    double a11_0 = 1;
    double a14_0 = 0;
    double a15_0 = 0;
    double a17_0 = 0;
    double a41_0 = 0;
    double a44_0 = 1;
    double a45_0 = 0;
    double a47_0 = 0;
    double a51_0 = 0;
    double a54_0 = 0;
    double a55_0 = 1;
    double a57_0 = 0;


    problem.phases(iphase).bounds.lower.states   << r_min, theta_min, phi_min, v_min, gamma_min, azim_min, p0_min, a11_min, a14_min, a15_min, a17_min, a41_min, a44_min, a45_min, a47_min, a51_min, a54_min, a55_min, a57_min;
    problem.phases(iphase).bounds.upper.states   << r_max, theta_max, phi_max, v_max, gamma_max, azim_max, p0_max, a11_max, a14_max, a15_max, a17_max, a41_max, a44_max, a45_max, a47_max, a51_max, a54_max, a55_max, a57_max;


    problem.phases(iphase).bounds.lower.controls << bank_min, k1_min, k2_min, k3_min;
    problem.phases(iphase).bounds.upper.controls << bank_max, k1_max, k2_max, k3_max;


    problem.phases(iphase).bounds.lower.events   << r_f, ra_f_min, a11_0, a14_0, a15_0, a17_0, a41_0, a44_0, a45_0, a47_0, a51_0, a54_0, a55_0, a57_0;
    problem.phases(iphase).bounds.upper.events   << r_f, ra_f_max, a11_0, a14_0, a15_0, a17_0, a41_0, a44_0, a45_0, a47_0, a51_0, a54_0, a55_0, a57_0;   

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

    // MatrixXd x_guess    =  zeros(nstates,nnodes);

    // x_guess.row(0)  = r_0*ones(1,nnodes);
    // x_guess.row(1)  = theta_0*ones(1,nnodes);
    // x_guess.row(2)  = phi_0*ones(1,nnodes);
    // x_guess.row(3)  = v_0*ones(1,nnodes);
    // x_guess.row(4)  = gamma_0*ones(1,nnodes);
    // x_guess.row(5)  = azim_0*ones(1,nnodes);
    // x_guess.row(6)  = p0_0*ones(1,nnodes);
    // x_guess.row(7)   = a11_0*ones(1,nnodes);
    // x_guess.row(8)   = a14_0*ones(1,nnodes);
    // x_guess.row(9)   = a15_0*ones(1,nnodes);
    // x_guess.row(10)  = a17_0*ones(1,nnodes);
    // x_guess.row(11)  = a41_0*ones(1,nnodes);
    // x_guess.row(12)  = a44_0*ones(1,nnodes);
    // x_guess.row(13)  = a45_0*ones(1,nnodes);
    // x_guess.row(14)  = a47_0*ones(1,nnodes);
    // x_guess.row(15)  = a51_0*ones(1,nnodes);
    // x_guess.row(16)  = a54_0*ones(1,nnodes);
    // x_guess.row(17)  = a55_0*ones(1,nnodes);
    // x_guess.row(18)  = a57_0*ones(1,nnodes);
    
    MatrixXd u_guess    =  zeros(ncontrols,nnodes);
    u_guess.row(0)  = 15*M_PI/180*ones(1,nnodes);
    u_guess.row(1)  = 1e-2*ones(1,nnodes);
    u_guess.row(2)  = 1e-2*ones(1,nnodes);
    u_guess.row(3)  = 1e-2*ones(1,nnodes);
    problem.phases(iphase).guess.controls       = u_guess;
    // problem.phases(iphase).guess.states         = x_guess;
    problem.phases(iphase).guess.time           = linspace(0.0,110.0,nnodes);
    

    iphase = 2; 

    nnodes    			            = problem.phases(iphase).nodes(0);
    ncontrols                       = problem.phases(iphase).ncontrols;
    nstates                         = problem.phases(iphase).nstates;

    // x_guess    =  zeros(nstates,nnodes);

    // x_guess.row(0)  = r_0*ones(1,nnodes);
    // x_guess.row(1)  = theta_0*ones(1,nnodes);
    // x_guess.row(2)  = phi_0*ones(1,nnodes);
    // x_guess.row(3)  = v_0*ones(1,nnodes);
    // x_guess.row(4)  = gamma_0*ones(1,nnodes);
    // x_guess.row(5)  = azim_0*ones(1,nnodes);
    // x_guess.row(6)  = p0_0*ones(1,nnodes);
    // x_guess.row(7)   = a11_0*ones(1,nnodes);
    // x_guess.row(8)   = a14_0*ones(1,nnodes);
    // x_guess.row(9)   = a15_0*ones(1,nnodes);
    // x_guess.row(10)  = a17_0*ones(1,nnodes);
    // x_guess.row(11)  = a41_0*ones(1,nnodes);
    // x_guess.row(12)  = a44_0*ones(1,nnodes);
    // x_guess.row(13)  = a45_0*ones(1,nnodes);
    // x_guess.row(14)  = a47_0*ones(1,nnodes);
    // x_guess.row(15)  = a51_0*ones(1,nnodes);
    // x_guess.row(16)  = a54_0*ones(1,nnodes);
    // x_guess.row(17)  = a55_0*ones(1,nnodes);
    // x_guess.row(18)  = a57_0*ones(1,nnodes);
    
    u_guess    =  zeros(ncontrols,nnodes);
    u_guess.row(0)  = 90*M_PI/180*ones(1,nnodes);
    u_guess.row(1)  = 1e-2*ones(1,nnodes);
    u_guess.row(2)  = 1e-2*ones(1,nnodes);
    u_guess.row(3)  = 1e-2*ones(1,nnodes);
    problem.phases(iphase).guess.controls       = u_guess;
    // problem.phases(iphase).guess.states         = x_guess;
    problem.phases(iphase).guess.time           = linspace(110,300,nnodes);


    // Define the state type for Boost ODEINT
    typedef vector<double> state_type;
    state_type y(23);

    y[0] = r_0    ;
    y[1] = theta_0;
    y[2] = phi_0  ;
    y[3] = v_0    ;
    y[4] = gamma_0;
    y[5] = azim_0 ;
    y[6] = p0_0   ;
    double smallInit = 1;
    y[7] = smallInit; y[8] = 0.0; y[9] = 0.0; y[10] = 0.0;
    y[11] = 0.0; y[12] =smallInit; y[13] = 0.0; y[14] = 0.0; 
    y[15] = 0.0; y[16] = 0.0; y[17] = smallInit; y[18] = 0.0;
    y[19] = 0.0; y[20] = 0.0; y[21] = 0.0; y[22] = smallInit; 

    double t_start = 0.0;
    double t_end = 110.0;
    double dt = (t_end-t_start)/(nnodes-1);

    // Define the ODE stepper (e.g., Runge-Kutta 4th order)
    typedef runge_kutta4<state_type> stepper_type;

    double bank = 15*M_PI/180;

    auto wrapped_rhs = [&](const vector<double>& y, vector<double>& dydt, const double t) {
        rhs(y, dydt, t, CONSTANTS, bank);
    };

    vector<state_type> x_vec_1;
    vector<state_type> x_vec_2;
    vector<double> times_1;
    vector<double> times_2;

    // Create the push_back_state_and_time observer
    push_back_state_and_time obs(x_vec_1, times_1);


    // Integrate the ODE
    integrate_const(stepper_type(), wrapped_rhs ,
        y , t_start , t_end , dt ,
        obs );
    

    t_start = 110.0;
    t_end = 300.0;
    dt = (t_end-t_start)/(nnodes-1);
    
    size_t rows = x_vec_1.size();
    y = x_vec_1.back();

    // Create the push_back_state_and_time observer
    push_back_state_and_time obs2(x_vec_2, times_2);
    
    bank = 90*M_PI/180;
    // Integrate the ODE
    integrate_const(stepper_type(), wrapped_rhs ,
        y , t_start , t_end , dt ,
        obs2 );

    
    MatrixXd x_guess_1    =  zeros(nstates,nnodes);
    MatrixXd x_guess_2    =  zeros(nstates,nnodes);
    // Convert the vector of vectors to an initial guess
    vectorToinit(x_vec_1, x_vec_2, &x_guess_1, &x_guess_2);
    
    // Save the vector to a dat file
    std::ofstream outFile("x_int_sol.dat");
    if (outFile.is_open()) {
        for (const auto& row : x_vec_1) {
            for (const auto& value : row) {
                outFile << value << " ";
            }
            outFile << std::endl;           
        }
        for (const auto& row : x_vec_2) {
            for (const auto& value : row) {
                outFile << value << " ";
            }
            outFile << std::endl; // Move to the next line after each row
        }
        outFile.close();
        std::cout << "x_int_sol saved to output.dat" << std::endl;
    } else {
        std::cerr << "Error: Unable to open file for writing!" << std::endl;
        return 1;
    }
   


    // Save the matrix to a dat file
    std::ofstream outFile2("x_guess_phaase1.dat");
    if (outFile2.is_open()) {
        outFile2 <<  x_guess_1;
        outFile2 << endl;
        outFile2 << x_guess_2;
        outFile2.close();
        std::cout << "x_guess saved to output.dat" << std::endl;
    } else {
        std::cerr << "Error: Unable to open file for writing!" << std::endl;
        return 1;
    }

    problem.phases(iphase).guess.states         = x_guess_2;

    iphase = 1; 
    problem.phases(iphase).guess.states         = x_guess_1;




////////////////////////////////////////////////////////////////////////////
///////////////////  Enter algorithm options  //////////////////////////////
////////////////////////////////////////////////////////////////////////////


    algorithm.nlp_iter_max                = 5000;
    algorithm.nlp_tolerance               = 1e-6;
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
    k = u.block(1,0,3,u.cols());

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
