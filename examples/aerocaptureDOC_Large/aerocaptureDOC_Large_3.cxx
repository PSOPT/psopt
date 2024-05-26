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

// Function to check if all elements in a column are zero
bool is_column_zero(const MatrixXd &matrix, int col) {
    for (int i = 0; i < matrix.rows(); ++i) {
        if (matrix(i, col) != 0) {
            return false; // Column is not all zeros
        }
    }
    return true; // Column is all zeros
}

// Function to replace a column with the previous column
void replace_column_with_previous(MatrixXd &matrix, int col) {
    if (col > 0) {
        matrix.col(col) = matrix.col(col - 1);
    }
}

// Function to convert a vector of vectors into an Eigen matrix
int vectorToinit(const std::vector<state_type>& vec, const std::vector<state_type>& vec2, Eigen::MatrixXd *x_guess_1, Eigen::MatrixXd *x_guess_2 ) {
    

    size_t rows = vec.size();
    size_t cols = vec[0].size(); // Assuming all state vectors have the same size
    size_t rows2 = vec2.size();
    size_t cols2 = vec2[0].size(); // Assuming all state vectors have the same size

    Eigen::MatrixXd mat(7, 7);
    Eigen::MatrixXd mat2(7, 7);
    Eigen::MatrixXd matf(7, 7);
    Eigen::MatrixXd lam(7, 7);
    Eigen::MatrixXd lam2(7, 7);
    
    for (size_t i = 0; i < 7; ++i) {
        for (size_t j = 0; j < 7; ++j) {
            matf(i, j) = vec2[rows2-1][7+7*i+j];
        }
    }

    for (size_t l = 0; l < rows; ++l) {
        for (size_t i = 0; i < 7; ++i) {
            for (size_t j = 0; j < 7; ++j) {
                mat(i, j) = vec[l][7+7*i+j];
            }
        }
         lam = matf * mat.inverse();
        
        if (l==rows2-1){
            std::ofstream outFile3("m_f.dat");
            if (outFile3.is_open()) {
                outFile3 <<  matf;
                outFile3 <<  endl;
                outFile3 <<  lam2;
                outFile3 <<  endl;
                outFile3.close();
            }
        }

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
        (*x_guess_1)(11,l)  = lam(0,4);
        (*x_guess_1)(12,l)  = lam(0,5);
        (*x_guess_1)(13,l)  = lam(0,6);
        (*x_guess_1)(14,l)  = lam(1,0);
        (*x_guess_1)(15,l)  = lam(1,1);
        (*x_guess_1)(16,l)  = lam(1,2);
        (*x_guess_1)(17,l)  = lam(1,3);
        (*x_guess_1)(18,l)  = lam(1,4);
        (*x_guess_1)(19,l)  = lam(1,5);
        (*x_guess_1)(20,l)  = lam(1,6);
        (*x_guess_1)(21,l)  = lam(2,0);
        (*x_guess_1)(22,l)  = lam(2,1);
        (*x_guess_1)(23,l)  = lam(2,2);
        (*x_guess_1)(24,l)  = lam(2,3);
        (*x_guess_1)(25,l)  = lam(2,4);
        (*x_guess_1)(26,l)  = lam(2,5);
        (*x_guess_1)(27,l)  = lam(2,6);
        (*x_guess_1)(28,l)  = lam(3,0);
        (*x_guess_1)(29,l)  = lam(3,1);
        (*x_guess_1)(30,l)  = lam(3,2);
        (*x_guess_1)(31,l)  = lam(3,3);
        (*x_guess_1)(32,l)  = lam(3,4);
        (*x_guess_1)(33,l)  = lam(3,5);
        (*x_guess_1)(34,l)  = lam(3,6);
        (*x_guess_1)(35,l)  = lam(4,0);
        (*x_guess_1)(36,l)  = lam(4,1);
        (*x_guess_1)(37,l)  = lam(4,2);
        (*x_guess_1)(38,l)  = lam(4,3);
        (*x_guess_1)(39,l)  = lam(4,4);
        (*x_guess_1)(40,l)  = lam(4,5);
        (*x_guess_1)(41,l)  = lam(4,6);
        (*x_guess_1)(42,l)  = lam(5,0);
        (*x_guess_1)(43,l)  = lam(5,1);
        (*x_guess_1)(44,l)  = lam(5,2);
        (*x_guess_1)(45,l)  = lam(5,3);
        (*x_guess_1)(46,l)  = lam(5,4);
        (*x_guess_1)(47,l)  = lam(5,5);
        (*x_guess_1)(48,l)  = lam(5,6);
        (*x_guess_1)(49,l)  = lam(6,0);
        (*x_guess_1)(50,l)  = lam(6,1);
        (*x_guess_1)(51,l)  = lam(6,2);
        (*x_guess_1)(52,l)  = lam(6,3);
        (*x_guess_1)(53,l)  = lam(6,4);
        (*x_guess_1)(54,l)  = lam(6,5);
        (*x_guess_1)(55,l)  = lam(6,6);
    }

    for (size_t l = 0; l < rows2; ++l) {
        for (size_t i = 0; i < 7; ++i) {
            for (size_t j = 0; j < 7; ++j) {
                mat2(i, j) = vec2[l][7+7*i+j];
            }
        }
    
        lam2 = matf * mat2.inverse();

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
        (*x_guess_2)(11,l)  = lam2(0,4);
        (*x_guess_2)(12,l)  = lam2(0,5);
        (*x_guess_2)(13,l)  = lam2(0,6);
        (*x_guess_2)(14,l)  = lam2(1,0);
        (*x_guess_2)(15,l)  = lam2(1,1);
        (*x_guess_2)(16,l)  = lam2(1,2);
        (*x_guess_2)(17,l)  = lam2(1,3);
        (*x_guess_2)(18,l)  = lam2(1,4);
        (*x_guess_2)(19,l)  = lam2(1,5);
        (*x_guess_2)(20,l)  = lam2(1,6);
        (*x_guess_2)(21,l)  = lam2(2,0);
        (*x_guess_2)(22,l)  = lam2(2,1);
        (*x_guess_2)(23,l)  = lam2(2,2);
        (*x_guess_2)(24,l)  = lam2(2,3);
        (*x_guess_2)(25,l)  = lam2(2,4);
        (*x_guess_2)(26,l)  = lam2(2,5);
        (*x_guess_2)(27,l)  = lam2(2,6);
        (*x_guess_2)(28,l)  = lam2(3,0);
        (*x_guess_2)(29,l)  = lam2(3,1);
        (*x_guess_2)(30,l)  = lam2(3,2);
        (*x_guess_2)(31,l)  = lam2(3,3);
        (*x_guess_2)(32,l)  = lam2(3,4);
        (*x_guess_2)(33,l)  = lam2(3,5);
        (*x_guess_2)(34,l)  = lam2(3,6);
        (*x_guess_2)(35,l)  = lam2(4,0);
        (*x_guess_2)(36,l)  = lam2(4,1);
        (*x_guess_2)(37,l)  = lam2(4,2);
        (*x_guess_2)(38,l)  = lam2(4,3);
        (*x_guess_2)(39,l)  = lam2(4,4);
        (*x_guess_2)(40,l)  = lam2(4,5);
        (*x_guess_2)(41,l)  = lam2(4,6);
        (*x_guess_2)(42,l)  = lam2(5,0);
        (*x_guess_2)(43,l)  = lam2(5,1);
        (*x_guess_2)(44,l)  = lam2(5,2);
        (*x_guess_2)(45,l)  = lam2(5,3);
        (*x_guess_2)(46,l)  = lam2(5,4);
        (*x_guess_2)(47,l)  = lam2(5,5);
        (*x_guess_2)(48,l)  = lam2(5,6);
        (*x_guess_2)(49,l)  = lam2(6,0);
        (*x_guess_2)(50,l)  = lam2(6,1);
        (*x_guess_2)(51,l)  = lam2(6,2);
        (*x_guess_2)(52,l)  = lam2(6,3);
        (*x_guess_2)(53,l)  = lam2(6,4);
        (*x_guess_2)(54,l)  = lam2(6,5);
        (*x_guess_2)(55,l)  = lam2(6,6);
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
    double a15    = y[ 11 ];
    double a16    = y[ 12 ];
    double a17    = y[ 13 ];

    double a21    = y[ 14 ];
    double a22    = y[ 15 ];
    double a23    = y[ 16 ];
    double a24    = y[ 17 ];
    double a25    = y[ 18 ];
    double a26    = y[ 19 ];
    double a27    = y[ 20 ];

    double a31    = y[ 21 ];
    double a32    = y[ 22 ];
    double a33    = y[ 23 ];
    double a34    = y[ 24 ];
    double a35    = y[ 25 ];
    double a36    = y[ 26 ];
    double a37    = y[ 27 ];

    double a41    = y[ 28 ];
    double a42    = y[ 29 ];
    double a43    = y[ 30 ];
    double a44    = y[ 31 ];
    double a45    = y[ 32 ];
    double a46    = y[ 33 ];
    double a47    = y[ 34 ];

    double a51    = y[ 35 ];
    double a52    = y[ 36 ];
    double a53    = y[ 37 ];
    double a54    = y[ 38 ];
    double a55    = y[ 39 ];
    double a56    = y[ 40 ];
    double a57    = y[ 41 ];

    double a61    = y[ 42 ];
    double a62    = y[ 43 ];
    double a63    = y[ 44 ];
    double a64    = y[ 45 ];
    double a65    = y[ 46 ];
    double a66    = y[ 47 ];
    double a67    = y[ 48 ];

    double a71    = y[ 49 ];
    double a72    = y[ 50 ];
    double a73    = y[ 51 ];
    double a74    = y[ 52 ];
    double a75    = y[ 53 ];
    double a76    = y[ 54 ];
    double a77    = y[ 55 ];
    

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
    double eta      = 4*((bank_max-bank)*(bank-bank_min))/((bank_max-bank_min)*(bank_max-bank_min));
    


double s11 = 0.0;
double s12 = 0.0;
double s13 = 0.0;
double s14 = sin(gamma);
double s15 = v*cos(gamma);
double s16 = 0.0;
double s17 = 0.0;
double s21 = -(1.0/(r*r)*v*cos(gamma)*sin(azim))/cos(phi);
double s22 = 0.0;
double s23 = (v*cos(gamma)*1.0/pow(cos(phi),2.0)*sin(azim)*sin(phi))/r;
double s24 = (cos(gamma)*sin(azim))/(r*cos(phi));
double s25 = -(v*sin(azim)*sin(gamma))/(r*cos(phi));
double s26 = (v*cos(azim)*cos(gamma))/(r*cos(phi));
double s27 = 0.0;
double s31 = -1.0/(r*r)*v*cos(azim)*cos(gamma);
double s32 = 0.0;
double s33 = 0.0;
double s34 = (cos(azim)*cos(gamma))/r;
double s35 = -(v*cos(azim)*sin(gamma))/r;
double s36 = -(v*cos(gamma)*sin(azim))/r;
double s37 = 0.0;
double s41 = (Omega*Omega)*cos(phi)*(cos(phi)*sin(gamma)-cos(azim)*cos(gamma)*sin(phi))-mu*1.0/(r*r*r)*sin(gamma)*(J2*(Re*Re)*1.0/(r*r)*(pow(sin(phi),2.0)*(9.0/2.0)-3.0/2.0)-1.0)*2.0-J2*(Re*Re)*mu*1.0/(r*r*r*r*r)*sin(gamma)*(pow(sin(phi),2.0)*(9.0/2.0)-3.0/2.0)*2.0+(B*Cd*S_ref*p0*(v*v)*exp(B*(Re+h_0-r)))/(m*2.0)+J2*(Re*Re)*mu*1.0/(r*r*r*r*r)*cos(phi)*sin(phi)*cos(gamma*cos(azim))*1.2E+1;
double s42 = 0.0;
double s43 = 1.0/(r*r*r*r)*((Omega*Omega)*(r*r*r*r*r)*cos(azim)*cos(gamma)+J2*(Re*Re)*mu*cos(gamma*cos(azim))*3.0-(Omega*Omega)*(r*r*r*r*r)*cos(phi)*sin(gamma)*sin(phi)*2.0-(Omega*Omega)*(r*r*r*r*r)*cos(azim)*cos(gamma)*pow(cos(phi),2.0)*2.0-J2*(Re*Re)*mu*pow(cos(phi),2.0)*cos(gamma*cos(azim))*6.0+J2*(Re*Re)*mu*cos(phi)*sin(gamma)*sin(phi)*9.0);
double s44 = -(Cd*S_ref*p0*v*exp(B*(Re+h_0-r)))/m;
double s45 = (Omega*Omega)*r*cos(phi)*(cos(gamma)*cos(phi)+cos(azim)*sin(gamma)*sin(phi))+mu*1.0/(r*r)*cos(gamma)*(J2*(Re*Re)*1.0/(r*r)*(pow(sin(phi),2.0)*(9.0/2.0)-3.0/2.0)-1.0)+J2*(Re*Re)*mu*1.0/(r*r*r*r)*cos(azim)*cos(phi)*sin(phi)*sin(gamma*cos(azim))*3.0;
double s46 = 1.0/(r*r*r*r)*cos(phi)*sin(azim)*sin(phi)*((Omega*Omega)*(r*r*r*r*r)*cos(gamma)-J2*(Re*Re)*gamma*mu*sin(gamma*cos(azim))*3.0);
double s47 = (Cd*S_ref*(v*v)*exp(B*(Re+h_0-r))*(-1.0/2.0))/m;
double s51 = -(cos(gamma)*(1.0/(r*r)*(v*v)+mu*1.0/(r*r*r)*(J2*(Re*Re)*1.0/(r*r)*(pow(sin(phi),2.0)*(9.0/2.0)-3.0/2.0)-1.0)*2.0+J2*(Re*Re)*mu*1.0/(r*r*r*r*r)*(pow(sin(phi),2.0)*(9.0/2.0)-3.0/2.0)*2.0)-(Omega*Omega)*cos(phi)*(cos(gamma)*cos(phi)+cos(azim)*sin(gamma)*sin(phi))+(B*Cl*S_ref*p0*(v*v)*exp(B*(Re+h_0-r))*cos(bank))/(m*2.0)+J2*(Re*Re)*mu*1.0/(r*r*r*r*r)*cos(azim)*cos(phi)*sin(gamma)*sin(phi)*1.2E+1)/v;
double s52 = 0.0;
double s53 = -(1.0/(r*r*r*r)*((Omega*Omega)*(r*r*r*r*r)*cos(azim)*sin(gamma)+(Omega*Omega)*(r*r*r*r*r)*cos(gamma)*cos(phi)*sin(phi)*2.0+J2*(Re*Re)*mu*cos(azim)*sin(gamma)*3.0+Omega*(r*r*r*r)*v*sin(azim)*sin(phi)*2.0-(Omega*Omega)*(r*r*r*r*r)*cos(azim)*pow(cos(phi),2.0)*sin(gamma)*2.0-J2*(Re*Re)*mu*cos(gamma)*cos(phi)*sin(phi)*9.0-J2*(Re*Re)*mu*cos(azim)*pow(cos(phi),2.0)*sin(gamma)*6.0))/v;
double s54 = (1.0/(r*r*r*r)*1.0/(v*v)*(m*(r*r*r)*(v*v)*cos(gamma)*2.0+m*mu*(r*r)*cos(gamma)*2.0-J2*(Re*Re)*m*mu*cos(gamma)*6.0-(Omega*Omega)*m*(r*r*r*r*r)*cos(gamma)*pow(cos(phi),2.0)*2.0+J2*(Re*Re)*m*mu*cos(gamma)*pow(cos(phi),2.0)*9.0-(Omega*Omega)*m*(r*r*r*r*r)*cos(azim)*cos(phi)*sin(gamma)*sin(phi)*2.0+Cl*S_ref*p0*(r*r*r*r)*(v*v)*exp(B*(Re+h_0-r))*cos(bank)-J2*(Re*Re)*m*mu*cos(azim)*cos(phi)*sin(gamma)*sin(phi)*6.0))/(m*2.0);
double s55 = -(sin(gamma)*((v*v)/r+mu*1.0/(r*r)*(J2*(Re*Re)*1.0/(r*r)*(pow(sin(phi),2.0)*(9.0/2.0)-3.0/2.0)-1.0))+(Omega*Omega)*r*cos(phi)*(cos(phi)*sin(gamma)-cos(azim)*cos(gamma)*sin(phi))-J2*(Re*Re)*mu*1.0/(r*r*r*r)*cos(azim)*cos(gamma)*cos(phi)*sin(phi)*3.0)/v;
double s56 = -(1.0/(r*r*r*r)*cos(phi)*(Omega*(r*r*r*r)*v*cos(azim)*-2.0+(Omega*Omega)*(r*r*r*r*r)*sin(azim)*sin(gamma)*sin(phi)+J2*(Re*Re)*mu*sin(azim)*sin(gamma)*sin(phi)*3.0))/v;
double s57 = (Cl*S_ref*v*exp(B*(Re+h_0-r))*cos(bank))/(m*2.0);
double s61 = (1.0/(r*r*r*r*r)*(m*(r*r*r)*(v*v)*pow(cos(gamma),2.0)*sin(azim)*sin(phi)*2.0+(Omega*Omega)*m*(r*r*r*r*r)*sin(azim)*sin(phi)*(pow(sin(phi),2.0)-1.0)*2.0-J2*(Re*Re)*m*mu*sin(azim)*sin(phi)*(pow(sin(phi),2.0)-1.0)*2.4E+1+B*Cl*S_ref*p0*(r*r*r*r*r)*(v*v)*exp(B*(Re+h_0-r))*cos(phi)*sin(bank))*(-1.0/2.0))/(m*v*cos(gamma)*cos(phi));
double s62 = 0.0;
double s63 = (1.0/(r*r*r*r)*1.0/pow(cos(phi),2.0)*((r*r*r)*(v*v)*pow(cos(gamma),2.0)*sin(azim)-(Omega*Omega)*(r*r*r*r*r)*pow(cos(phi),2.0)*sin(azim)+(Omega*Omega)*(r*r*r*r*r)*pow(cos(phi),4.0)*sin(azim)*2.0-J2*(Re*Re)*mu*pow(cos(phi),2.0)*sin(azim)*3.0+J2*(Re*Re)*mu*pow(cos(phi),4.0)*sin(azim)*6.0+Omega*(r*r*r*r)*v*cos(gamma)*pow(cos(phi),3.0)*2.0+Omega*(r*r*r*r)*v*cos(azim)*pow(cos(phi),2.0)*sin(gamma)*sin(phi)*2.0))/(v*cos(gamma));
double s64 = (1.0/(r*r*r*r)*1.0/(v*v)*(m*(r*r*r)*(v*v)*pow(cos(gamma),2.0)*sin(azim)*sin(phi)*2.0+(Omega*Omega)*m*(r*r*r*r*r)*sin(azim)*sin(phi)*(pow(sin(phi),2.0)-1.0)*2.0+J2*(Re*Re)*m*mu*sin(azim)*sin(phi)*(pow(sin(phi),2.0)-1.0)*6.0+Cl*S_ref*p0*(r*r*r*r)*(v*v)*exp(B*(Re+h_0-r))*cos(phi)*sin(bank)))/(m*cos(gamma)*cos(phi)*2.0);
double s65 = (1.0/(r*r*r*r)*1.0/pow(cos(gamma),2.0)*(Omega*m*(r*r*r*r)*v*cos(azim)*pow(cos(phi),2.0)*-4.0+(Omega*Omega)*m*(r*r*r*r*r)*pow(cos(phi),2.0)*sin(azim)*sin(gamma)*sin(phi)*2.0-m*(r*r*r)*(v*v)*pow(cos(gamma),2.0)*sin(azim)*sin(gamma)*sin(phi)*2.0+J2*(Re*Re)*m*mu*pow(cos(phi),2.0)*sin(azim)*sin(gamma)*sin(phi)*6.0+Cl*S_ref*p0*(r*r*r*r)*(v*v)*exp(B*(Re+h_0-r))*cos(phi)*sin(bank)*sin(gamma)))/(m*v*cos(phi)*2.0);
double s66 = (1.0/(r*r*r*r)*((Omega*Omega)*(r*r*r*r*r)*cos(azim)*pow(cos(phi),2.0)*sin(phi)+(r*r*r)*(v*v)*cos(azim)*pow(cos(gamma),2.0)*sin(phi)+J2*(Re*Re)*mu*cos(azim)*pow(cos(phi),2.0)*sin(phi)*3.0+Omega*(r*r*r*r)*v*pow(cos(phi),2.0)*sin(azim)*sin(gamma)*2.0))/(v*cos(gamma)*cos(phi));
double s67 = (Cl*S_ref*v*exp(B*(Re+h_0-r))*sin(bank))/(m*cos(gamma)*2.0);
double s71 = 0.0;
double s72 = 0.0;
double s73 = 0.0;
double s74 = 0.0;
double s75 = 0.0;
double s76 = 0.0;
double s77 = 0.0;

    
    
    int size = 7;
    MatrixXd a_matrix = zeros(size, size);
    MatrixXd s_matrix = zeros(size, size);
    MatrixXd aDot_matrix1 = zeros(size, size);
    
    // Initialize matrices with some values (you can use your own initialization)
    a_matrix << a11, a12, a13, a14, a15, a16, a17,
                a21, a22, a23, a24, a25, a26, a27,
                a31, a32, a33, a34, a35, a36, a37,
                a41, a42, a43, a44, a45, a46, a47,
                a51, a52, a53, a54, a55, a56, a57,
                a61, a62, a63, a64, a65, a66, a67,
                a71, a72, a73, a74, a75, a76, a77;
    
    s_matrix << s11, s12, s13, s14, s15, s16, s17,
                s21, s22, s23, s24, s25, s26, s27,
                s31, s32, s33, s34, s35, s36, s37,
                s41, s42, s43, s44, s45, s36, s47,
                s51, s52, s53, s54, s55, s56, s57,
                s61, s62, s63, s64, s65, s66, s67,
                s71, s72, s73, s74, s75, s76, s77;

    // Perform matrix multiplication
    aDot_matrix1 = s_matrix * a_matrix;
    
    dydt[ 0 ]  = rdot;
    dydt[ 1 ]  = thetadot;
    dydt[ 2 ]  = phidot;
    dydt[ 3 ]  = vdot;
    dydt[ 4 ]  = gammadot;
    dydt[ 5 ]  = azimdot ;
    dydt[ 6 ]  = p0dot;
    dydt[ 7 ]  = aDot_matrix1(0,0);  
    dydt[ 8 ]  = aDot_matrix1(0,1);    
    dydt[ 9 ]  = aDot_matrix1(0,2);    
    dydt[ 10 ] = aDot_matrix1(0,3);
    dydt[ 11 ] = aDot_matrix1(0,4);
    dydt[ 12 ] = aDot_matrix1(0,5);
    dydt[ 13 ] = aDot_matrix1(0,6);
    dydt[ 14 ] = aDot_matrix1(1,0);
    dydt[ 15 ] = aDot_matrix1(1,1);
    dydt[ 16 ] = aDot_matrix1(1,2);
    dydt[ 17 ] = aDot_matrix1(1,3);
    dydt[ 18 ] = aDot_matrix1(1,4);
    dydt[ 19 ] = aDot_matrix1(1,5);
    dydt[ 20 ] = aDot_matrix1(1,6);
    dydt[ 21 ] = aDot_matrix1(2,0);
    dydt[ 22 ] = aDot_matrix1(2,1);
    dydt[ 23 ] = aDot_matrix1(2,2);
    dydt[ 24 ] = aDot_matrix1(2,3);
    dydt[ 25 ] = aDot_matrix1(2,4);
    dydt[ 26 ] = aDot_matrix1(2,5);
    dydt[ 27 ] = aDot_matrix1(2,6);
    dydt[ 28 ] = aDot_matrix1(3,0);
    dydt[ 29 ] = aDot_matrix1(3,1);
    dydt[ 30 ] = aDot_matrix1(3,2);
    dydt[ 31 ] = aDot_matrix1(3,3);
    dydt[ 32 ] = aDot_matrix1(3,4);
    dydt[ 33 ] = aDot_matrix1(3,5);
    dydt[ 34 ] = aDot_matrix1(3,6);
    dydt[ 35 ] = aDot_matrix1(4,0);
    dydt[ 36 ] = aDot_matrix1(4,1);
    dydt[ 37 ] = aDot_matrix1(4,2);
    dydt[ 38 ] = aDot_matrix1(4,3);
    dydt[ 39 ] = aDot_matrix1(4,4);
    dydt[ 40 ] = aDot_matrix1(4,5);
    dydt[ 41 ] = aDot_matrix1(4,6);
    dydt[ 42 ] = aDot_matrix1(5,0);
    dydt[ 43 ] = aDot_matrix1(5,1);
    dydt[ 44 ] = aDot_matrix1(5,2);
    dydt[ 45 ] = aDot_matrix1(5,3);
    dydt[ 46 ] = aDot_matrix1(5,4);
    dydt[ 46 ] = aDot_matrix1(5,5);
    dydt[ 48 ] = aDot_matrix1(5,6);
    dydt[ 49 ] = aDot_matrix1(6,0);
    dydt[ 50 ] = aDot_matrix1(6,1);
    dydt[ 51 ] = aDot_matrix1(6,2);
    dydt[ 52 ] = aDot_matrix1(6,3);
    dydt[ 53 ] = aDot_matrix1(6,4);
    dydt[ 54 ] = aDot_matrix1(6,5);
    dydt[ 55 ] = aDot_matrix1(6,6);

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
    double q = 1e-2;
    adouble runCost =q*( states[34]*states[34]+states[41]*states[41]);// + r *(controls[1]*controls[1] + controls[2]*controls[2] +controls[3]*controls[3]+controls[4]*controls[4]+controls[5]*controls[5]+controls[6]*controls[6]+controls[7]*controls[7]);
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
    adouble a11    = states[ 7 ];
    adouble a12    = states[ 8 ];
    adouble a13    = states[ 9 ];
    adouble a14    = states[ 10 ];
    adouble a15    = states[ 11 ];
    adouble a16    = states[ 12 ];
    adouble a17    = states[ 13 ];
    adouble a21    = states[ 14 ];
    adouble a22    = states[ 15 ];
    adouble a23    = states[ 16 ];
    adouble a24    = states[ 17 ];
    adouble a25    = states[ 18 ];
    adouble a26    = states[ 19 ];
    adouble a27    = states[ 20 ];
    adouble a31    = states[ 21 ];
    adouble a32    = states[ 22 ];
    adouble a33    = states[ 23 ];
    adouble a34    = states[ 24 ];
    adouble a35    = states[ 25 ];
    adouble a36    = states[ 26 ];
    adouble a37    = states[ 27 ];
    adouble a41    = states[ 28 ];
    adouble a42    = states[ 29 ];
    adouble a43    = states[ 30 ];
    adouble a44    = states[ 31 ];
    adouble a45    = states[ 32 ];
    adouble a46    = states[ 33 ];
    adouble a47    = states[ 34 ];
    adouble a51    = states[ 35 ];
    adouble a52    = states[ 36 ];
    adouble a53    = states[ 37 ];
    adouble a54    = states[ 38 ];
    adouble a55    = states[ 39 ];
    adouble a56    = states[ 40 ];
    adouble a57    = states[ 41 ];
    adouble a61    = states[ 42 ];
    adouble a62    = states[ 43 ];
    adouble a63    = states[ 44 ];
    adouble a64    = states[ 45 ];
    adouble a65    = states[ 46 ];
    adouble a66    = states[ 47 ];
    adouble a67    = states[ 48 ];
    adouble a71    = states[ 49 ];
    adouble a72    = states[ 50 ];
    adouble a73    = states[ 51 ];
    adouble a74    = states[ 52 ];
    adouble a75    = states[ 53 ];
    adouble a76    = states[ 54 ];
    adouble a77    = states[ 55 ];
    

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
    
adouble s11 = 0.0;
adouble s12 = 0.0;
adouble s13 = 0.0;
adouble s14 = sin(gamma);
adouble s15 = v*cos(gamma);
adouble s16 = 0.0;
adouble s17 = 0.0;
adouble s21 = -(1.0/(r*r)*v*cos(gamma)*sin(azim))/cos(phi);
adouble s22 = 0.0;
adouble s23 = (v*cos(gamma)*1.0/pow(cos(phi),2.0)*sin(azim)*sin(phi))/r;
adouble s24 = (cos(gamma)*sin(azim))/(r*cos(phi));
adouble s25 = -(v*sin(azim)*sin(gamma))/(r*cos(phi));
adouble s26 = (v*cos(azim)*cos(gamma))/(r*cos(phi));
adouble s27 = 0.0;
adouble s31 = -1.0/(r*r)*v*cos(azim)*cos(gamma);
adouble s32 = 0.0;
adouble s33 = 0.0;
adouble s34 = (cos(azim)*cos(gamma))/r;
adouble s35 = -(v*cos(azim)*sin(gamma))/r;
adouble s36 = -(v*cos(gamma)*sin(azim))/r;
adouble s37 = 0.0;
adouble s41 = (Omega*Omega)*cos(phi)*(cos(phi)*sin(gamma)-cos(azim)*cos(gamma)*sin(phi))-mu*1.0/(r*r*r)*sin(gamma)*(J2*(Re*Re)*1.0/(r*r)*(pow(sin(phi),2.0)*(9.0/2.0)-3.0/2.0)-1.0)*2.0-J2*(Re*Re)*mu*1.0/(r*r*r*r*r)*sin(gamma)*(pow(sin(phi),2.0)*(9.0/2.0)-3.0/2.0)*2.0+(B*Cd*S_ref*p0*(v*v)*exp(B*(Re+h_0-r)))/(m*2.0)+J2*(Re*Re)*mu*1.0/(r*r*r*r*r)*cos(phi)*sin(phi)*cos(gamma*cos(azim))*1.2E+1;
adouble s42 = 0.0;
adouble s43 = 1.0/(r*r*r*r)*((Omega*Omega)*(r*r*r*r*r)*cos(azim)*cos(gamma)+J2*(Re*Re)*mu*cos(gamma*cos(azim))*3.0-(Omega*Omega)*(r*r*r*r*r)*cos(phi)*sin(gamma)*sin(phi)*2.0-(Omega*Omega)*(r*r*r*r*r)*cos(azim)*cos(gamma)*pow(cos(phi),2.0)*2.0-J2*(Re*Re)*mu*pow(cos(phi),2.0)*cos(gamma*cos(azim))*6.0+J2*(Re*Re)*mu*cos(phi)*sin(gamma)*sin(phi)*9.0);
adouble s44 = -(Cd*S_ref*p0*v*exp(B*(Re+h_0-r)))/m;
adouble s45 = (Omega*Omega)*r*cos(phi)*(cos(gamma)*cos(phi)+cos(azim)*sin(gamma)*sin(phi))+mu*1.0/(r*r)*cos(gamma)*(J2*(Re*Re)*1.0/(r*r)*(pow(sin(phi),2.0)*(9.0/2.0)-3.0/2.0)-1.0)+J2*(Re*Re)*mu*1.0/(r*r*r*r)*cos(azim)*cos(phi)*sin(phi)*sin(gamma*cos(azim))*3.0;
adouble s46 = 1.0/(r*r*r*r)*cos(phi)*sin(azim)*sin(phi)*((Omega*Omega)*(r*r*r*r*r)*cos(gamma)-J2*(Re*Re)*gamma*mu*sin(gamma*cos(azim))*3.0);
adouble s47 = (Cd*S_ref*(v*v)*exp(B*(Re+h_0-r))*(-1.0/2.0))/m;
adouble s51 = -(cos(gamma)*(1.0/(r*r)*(v*v)+mu*1.0/(r*r*r)*(J2*(Re*Re)*1.0/(r*r)*(pow(sin(phi),2.0)*(9.0/2.0)-3.0/2.0)-1.0)*2.0+J2*(Re*Re)*mu*1.0/(r*r*r*r*r)*(pow(sin(phi),2.0)*(9.0/2.0)-3.0/2.0)*2.0)-(Omega*Omega)*cos(phi)*(cos(gamma)*cos(phi)+cos(azim)*sin(gamma)*sin(phi))+(B*Cl*S_ref*p0*(v*v)*exp(B*(Re+h_0-r))*cos(bank))/(m*2.0)+J2*(Re*Re)*mu*1.0/(r*r*r*r*r)*cos(azim)*cos(phi)*sin(gamma)*sin(phi)*1.2E+1)/v;
adouble s52 = 0.0;
adouble s53 = -(1.0/(r*r*r*r)*((Omega*Omega)*(r*r*r*r*r)*cos(azim)*sin(gamma)+(Omega*Omega)*(r*r*r*r*r)*cos(gamma)*cos(phi)*sin(phi)*2.0+J2*(Re*Re)*mu*cos(azim)*sin(gamma)*3.0+Omega*(r*r*r*r)*v*sin(azim)*sin(phi)*2.0-(Omega*Omega)*(r*r*r*r*r)*cos(azim)*pow(cos(phi),2.0)*sin(gamma)*2.0-J2*(Re*Re)*mu*cos(gamma)*cos(phi)*sin(phi)*9.0-J2*(Re*Re)*mu*cos(azim)*pow(cos(phi),2.0)*sin(gamma)*6.0))/v;
adouble s54 = (1.0/(r*r*r*r)*1.0/(v*v)*(m*(r*r*r)*(v*v)*cos(gamma)*2.0+m*mu*(r*r)*cos(gamma)*2.0-J2*(Re*Re)*m*mu*cos(gamma)*6.0-(Omega*Omega)*m*(r*r*r*r*r)*cos(gamma)*pow(cos(phi),2.0)*2.0+J2*(Re*Re)*m*mu*cos(gamma)*pow(cos(phi),2.0)*9.0-(Omega*Omega)*m*(r*r*r*r*r)*cos(azim)*cos(phi)*sin(gamma)*sin(phi)*2.0+Cl*S_ref*p0*(r*r*r*r)*(v*v)*exp(B*(Re+h_0-r))*cos(bank)-J2*(Re*Re)*m*mu*cos(azim)*cos(phi)*sin(gamma)*sin(phi)*6.0))/(m*2.0);
adouble s55 = -(sin(gamma)*((v*v)/r+mu*1.0/(r*r)*(J2*(Re*Re)*1.0/(r*r)*(pow(sin(phi),2.0)*(9.0/2.0)-3.0/2.0)-1.0))+(Omega*Omega)*r*cos(phi)*(cos(phi)*sin(gamma)-cos(azim)*cos(gamma)*sin(phi))-J2*(Re*Re)*mu*1.0/(r*r*r*r)*cos(azim)*cos(gamma)*cos(phi)*sin(phi)*3.0)/v;
adouble s56 = -(1.0/(r*r*r*r)*cos(phi)*(Omega*(r*r*r*r)*v*cos(azim)*-2.0+(Omega*Omega)*(r*r*r*r*r)*sin(azim)*sin(gamma)*sin(phi)+J2*(Re*Re)*mu*sin(azim)*sin(gamma)*sin(phi)*3.0))/v;
adouble s57 = (Cl*S_ref*v*exp(B*(Re+h_0-r))*cos(bank))/(m*2.0);
adouble s61 = (1.0/(r*r*r*r*r)*(m*(r*r*r)*(v*v)*pow(cos(gamma),2.0)*sin(azim)*sin(phi)*2.0+(Omega*Omega)*m*(r*r*r*r*r)*sin(azim)*sin(phi)*(pow(sin(phi),2.0)-1.0)*2.0-J2*(Re*Re)*m*mu*sin(azim)*sin(phi)*(pow(sin(phi),2.0)-1.0)*2.4E+1+B*Cl*S_ref*p0*(r*r*r*r*r)*(v*v)*exp(B*(Re+h_0-r))*cos(phi)*sin(bank))*(-1.0/2.0))/(m*v*cos(gamma)*cos(phi));
adouble s62 = 0.0;
adouble s63 = (1.0/(r*r*r*r)*1.0/pow(cos(phi),2.0)*((r*r*r)*(v*v)*pow(cos(gamma),2.0)*sin(azim)-(Omega*Omega)*(r*r*r*r*r)*pow(cos(phi),2.0)*sin(azim)+(Omega*Omega)*(r*r*r*r*r)*pow(cos(phi),4.0)*sin(azim)*2.0-J2*(Re*Re)*mu*pow(cos(phi),2.0)*sin(azim)*3.0+J2*(Re*Re)*mu*pow(cos(phi),4.0)*sin(azim)*6.0+Omega*(r*r*r*r)*v*cos(gamma)*pow(cos(phi),3.0)*2.0+Omega*(r*r*r*r)*v*cos(azim)*pow(cos(phi),2.0)*sin(gamma)*sin(phi)*2.0))/(v*cos(gamma));
adouble s64 = (1.0/(r*r*r*r)*1.0/(v*v)*(m*(r*r*r)*(v*v)*pow(cos(gamma),2.0)*sin(azim)*sin(phi)*2.0+(Omega*Omega)*m*(r*r*r*r*r)*sin(azim)*sin(phi)*(pow(sin(phi),2.0)-1.0)*2.0+J2*(Re*Re)*m*mu*sin(azim)*sin(phi)*(pow(sin(phi),2.0)-1.0)*6.0+Cl*S_ref*p0*(r*r*r*r)*(v*v)*exp(B*(Re+h_0-r))*cos(phi)*sin(bank)))/(m*cos(gamma)*cos(phi)*2.0);
adouble s65 = (1.0/(r*r*r*r)*1.0/pow(cos(gamma),2.0)*(Omega*m*(r*r*r*r)*v*cos(azim)*pow(cos(phi),2.0)*-4.0+(Omega*Omega)*m*(r*r*r*r*r)*pow(cos(phi),2.0)*sin(azim)*sin(gamma)*sin(phi)*2.0-m*(r*r*r)*(v*v)*pow(cos(gamma),2.0)*sin(azim)*sin(gamma)*sin(phi)*2.0+J2*(Re*Re)*m*mu*pow(cos(phi),2.0)*sin(azim)*sin(gamma)*sin(phi)*6.0+Cl*S_ref*p0*(r*r*r*r)*(v*v)*exp(B*(Re+h_0-r))*cos(phi)*sin(bank)*sin(gamma)))/(m*v*cos(phi)*2.0);
adouble s66 = (1.0/(r*r*r*r)*((Omega*Omega)*(r*r*r*r*r)*cos(azim)*pow(cos(phi),2.0)*sin(phi)+(r*r*r)*(v*v)*cos(azim)*pow(cos(gamma),2.0)*sin(phi)+J2*(Re*Re)*mu*cos(azim)*pow(cos(phi),2.0)*sin(phi)*3.0+Omega*(r*r*r*r)*v*pow(cos(phi),2.0)*sin(azim)*sin(gamma)*2.0))/(v*cos(gamma)*cos(phi));
adouble s67 = (Cl*S_ref*v*exp(B*(Re+h_0-r))*sin(bank))/(m*cos(gamma)*2.0);
adouble s71 = 0.0;
adouble s72 = 0.0;
adouble s73 = 0.0;
adouble s74 = 0.0;
adouble s75 = 0.0;
adouble s76 = 0.0;
adouble s77 = 0.0;


    



    adouble a11_dot    = -1*(a11*s11+a12*s21+a13*s31+a14*s41+a15*s51+a16*s61+a17*s71);   
    adouble a12_dot    = -1*(a11*s12+a12*s22+a13*s32+a14*s42+a15*s52+a16*s62+a17*s72);
    adouble a13_dot    = -1*(a11*s13+a12*s23+a13*s33+a14*s43+a15*s53+a16*s63+a17*s73);
    adouble a14_dot    = -1*(a11*s14+a12*s24+a13*s34+a14*s44+a15*s54+a16*s64+a17*s74);
    adouble a15_dot    = -1*(a11*s15+a12*s25+a13*s35+a14*s45+a15*s55+a16*s65+a17*s75);
    adouble a16_dot    = -1*(a11*s16+a12*s26+a13*s36+a14*s46+a15*s56+a16*s66+a17*s76);
    adouble a17_dot    = -1*(a11*s17+a12*s27+a13*s37+a14*s47+a15*s57+a16*s67+a17*s77);

    adouble a21_dot    = -1*(a21*s11+a22*s21+a23*s31+a24*s41+a25*s51+a26*s61+a27*s71);   
    adouble a22_dot    = -1*(a21*s12+a22*s22+a23*s32+a24*s42+a25*s52+a26*s62+a27*s72);
    adouble a23_dot    = -1*(a21*s13+a22*s23+a23*s33+a24*s43+a25*s53+a26*s63+a27*s73);
    adouble a24_dot    = -1*(a21*s14+a22*s24+a23*s34+a24*s44+a25*s54+a26*s64+a27*s74);
    adouble a25_dot    = -1*(a21*s15+a22*s25+a23*s35+a24*s45+a25*s55+a26*s65+a27*s75);
    adouble a26_dot    = -1*(a21*s16+a22*s26+a23*s36+a24*s46+a25*s56+a26*s66+a27*s76);
    adouble a27_dot    = -1*(a21*s17+a22*s27+a23*s37+a24*s47+a25*s57+a26*s67+a27*s77);

    adouble a31_dot    = -1*(a31*s11+a32*s21+a33*s31+a34*s41+a35*s51+a36*s61+a37*s71);   
    adouble a32_dot    = -1*(a31*s12+a32*s22+a33*s32+a34*s42+a35*s52+a36*s62+a37*s72);
    adouble a33_dot    = -1*(a31*s13+a32*s23+a33*s33+a34*s43+a35*s53+a36*s63+a37*s73);
    adouble a34_dot    = -1*(a31*s14+a32*s24+a33*s34+a34*s44+a35*s54+a36*s64+a37*s74);
    adouble a35_dot    = -1*(a31*s15+a32*s25+a33*s35+a34*s45+a35*s55+a36*s65+a37*s75);
    adouble a36_dot    = -1*(a31*s16+a32*s26+a33*s36+a34*s46+a35*s56+a36*s66+a37*s76);
    adouble a37_dot    = -1*(a31*s17+a32*s27+a33*s37+a34*s47+a35*s57+a36*s67+a37*s77);

    adouble a41_dot    = -1*(a41*s11+a42*s21+a43*s31+a44*s41+a45*s51+a46*s61+a47*s71);   
    adouble a42_dot    = -1*(a41*s12+a42*s22+a43*s32+a44*s42+a45*s52+a46*s62+a47*s72);
    adouble a43_dot    = -1*(a41*s13+a42*s23+a43*s33+a44*s43+a45*s53+a46*s63+a47*s73);
    adouble a44_dot    = -1*(a41*s14+a42*s24+a43*s34+a44*s44+a45*s54+a46*s64+a47*s74);
    adouble a45_dot    = -1*(a41*s15+a42*s25+a43*s35+a44*s45+a45*s55+a46*s65+a47*s75);
    adouble a46_dot    = -1*(a41*s16+a42*s26+a43*s36+a44*s46+a45*s56+a46*s66+a47*s76);
    adouble a47_dot    = -1*(a41*s17+a42*s27+a43*s37+a44*s47+a45*s57+a46*s67+a47*s77);

    adouble a51_dot    = -1*(a51*s11+a52*s21+a53*s31+a54*s41+a55*s51+a56*s61+a57*s71);   
    adouble a52_dot    = -1*(a51*s12+a52*s22+a53*s32+a54*s42+a55*s52+a56*s62+a57*s72);
    adouble a53_dot    = -1*(a51*s13+a52*s23+a53*s33+a54*s43+a55*s53+a56*s63+a57*s73);
    adouble a54_dot    = -1*(a51*s14+a52*s24+a53*s34+a54*s44+a55*s54+a56*s64+a57*s74);
    adouble a55_dot    = -1*(a51*s15+a52*s25+a53*s35+a54*s45+a55*s55+a56*s65+a57*s75);
    adouble a56_dot    = -1*(a51*s16+a52*s26+a53*s36+a54*s46+a55*s56+a56*s66+a57*s76);
    adouble a57_dot    = -1*(a51*s17+a52*s27+a53*s37+a54*s47+a55*s57+a56*s67+a57*s77);

    adouble a61_dot    = -1*(a61*s11+a62*s21+a63*s31+a64*s41+a65*s51+a66*s61+a67*s71);   
    adouble a62_dot    = -1*(a61*s12+a62*s22+a63*s32+a64*s42+a65*s52+a66*s62+a67*s72);
    adouble a63_dot    = -1*(a61*s13+a62*s23+a63*s33+a64*s43+a65*s53+a66*s63+a67*s73);
    adouble a64_dot    = -1*(a61*s14+a62*s24+a63*s34+a64*s44+a65*s54+a66*s64+a67*s74);
    adouble a65_dot    = -1*(a61*s15+a62*s25+a63*s35+a64*s45+a65*s55+a66*s65+a67*s75);
    adouble a66_dot    = -1*(a61*s16+a62*s26+a63*s36+a64*s46+a65*s56+a66*s66+a67*s76);
    adouble a67_dot    = -1*(a61*s17+a62*s27+a63*s37+a64*s47+a65*s57+a66*s67+a67*s77);

    adouble a71_dot    = -1*(a71*s11+a72*s21+a73*s31+a74*s41+a75*s51+a76*s61+a77*s71);   
    adouble a72_dot    = -1*(a71*s12+a72*s22+a73*s32+a74*s42+a75*s52+a76*s62+a77*s72);
    adouble a73_dot    = -1*(a71*s13+a72*s23+a73*s33+a74*s43+a75*s53+a76*s63+a77*s73);
    adouble a74_dot    = -1*(a71*s14+a72*s24+a73*s34+a74*s44+a75*s54+a76*s64+a77*s74);
    adouble a75_dot    = -1*(a71*s15+a72*s25+a73*s35+a74*s45+a75*s55+a76*s65+a77*s75);
    adouble a76_dot    = -1*(a71*s16+a72*s26+a73*s36+a74*s46+a75*s56+a76*s66+a77*s76);
    adouble a77_dot    = -1*(a71*s17+a72*s27+a73*s37+a74*s47+a75*s57+a76*s67+a77*s77);

    derivatives[ 0 ] = rdot;
    derivatives[ 1 ] = thetadot;
    derivatives[ 2 ] = phidot;
    derivatives[ 3 ] = vdot;
    derivatives[ 4 ] = gammadot;
    derivatives[ 5 ] = azimdot ;
    derivatives[ 6 ] = p0dot;

    derivatives[ 7 ]  = a11_dot;  
    derivatives[ 8 ]  = a12_dot;    
    derivatives[ 9 ]  = a13_dot;    
    derivatives[ 10 ] = a14_dot;
    derivatives[ 11 ] = a15_dot;
    derivatives[ 12 ] = a16_dot;
    derivatives[ 13 ] = a17_dot;

    derivatives[ 14 ] = a21_dot;
    derivatives[ 15 ] = a22_dot;
    derivatives[ 16 ] = a23_dot;
    derivatives[ 17 ] = a24_dot;
    derivatives[ 18 ] = a25_dot;
    derivatives[ 19 ] = a26_dot;
    derivatives[ 20 ] = a27_dot;

    derivatives[ 21 ] = a31_dot;
    derivatives[ 22 ] = a32_dot;
    derivatives[ 23 ] = a33_dot;
    derivatives[ 24 ] = a34_dot;
    derivatives[ 25 ] = a35_dot;
    derivatives[ 26 ] = a36_dot;
    derivatives[ 27 ] = a37_dot;

    derivatives[ 28 ] = a41_dot;
    derivatives[ 29 ] = a42_dot;
    derivatives[ 30 ] = a43_dot;
    derivatives[ 31 ] = a44_dot;
    derivatives[ 32 ] = a45_dot;
    derivatives[ 33 ] = a46_dot;
    derivatives[ 34 ] = a47_dot;

    derivatives[ 35 ] = a51_dot;
    derivatives[ 36 ] = a52_dot;
    derivatives[ 37 ] = a53_dot;
    derivatives[ 38 ] = a54_dot;
    derivatives[ 39 ] = a55_dot;
    derivatives[ 40 ] = a56_dot;
    derivatives[ 41 ] = a57_dot;

    derivatives[ 42 ] = a61_dot;
    derivatives[ 43 ] = a62_dot;
    derivatives[ 44 ] = a63_dot;
    derivatives[ 45 ] = a64_dot;
    derivatives[ 46 ] = a65_dot;
    derivatives[ 47 ] = a66_dot;
    derivatives[ 48 ] = a67_dot;

    derivatives[ 49 ] = a71_dot;
    derivatives[ 50 ] = a72_dot;
    derivatives[ 51 ] = a73_dot;
    derivatives[ 52 ] = a74_dot;
    derivatives[ 53 ] = a75_dot;
    derivatives[ 54 ] = a76_dot;
    derivatives[ 55 ] = a77_dot;

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

         e[16] = final_states[21]; 
         e[17] = final_states[22];
         e[18] = final_states[23];
         e[19] = final_states[24];
         e[20] = final_states[25];
         e[21] = final_states[26];
         e[22] = final_states[27];

         e[23] = final_states[28]; 
         e[24] = final_states[29];
         e[25] = final_states[30];
         e[26] = final_states[31];
         e[27] = final_states[32];
         e[28] = final_states[33];
         e[29] = final_states[34];

         e[30] = final_states[35]; 
         e[31] = final_states[36];
         e[32] = final_states[37];
         e[33] = final_states[38];
         e[34] = final_states[39];
         e[35] = final_states[40];
         e[36] = final_states[41];

         e[37] = final_states[42]; 
         e[38] = final_states[43];
         e[39] = final_states[44];
         e[40] = final_states[45];
         e[41] = final_states[46];
         e[42] = final_states[47];
         e[43] = final_states[48];

         e[44] = final_states[49]; 
         e[45] = final_states[50];
         e[46] = final_states[51];
         e[47] = final_states[52];
         e[48] = final_states[53];
         e[49] = final_states[54];
         e[50] = final_states[55];

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
    problem.nlinkages           = 57; // States and time continuity constraints

    psopt_level1_setup(problem);

/////////////////////////////////////////////////////////////////////////////
/////////   Define phase related information & do level 2 setup  ////////////
/////////////////////////////////////////////////////////////////////////////

    problem.phases(1).nstates   		= 56;
    problem.phases(1).ncontrols 		= 2;
    problem.phases(1).nevents   		= 7;

    problem.phases(2).nstates   		= 56;
    problem.phases(2).ncontrols 		= 2;
    problem.phases(2).nevents   		= 51;

    problem.phases(1).nodes      << 80; 
    problem.phases(2).nodes      << 80; 
    

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
    double k1_min  =0;
    double k1_max  =0;

    double a11_min  =-32.6431000000000;
    double a12_min  =-0;
    double a13_min  =-936.482000000000;
    double a14_min  =-1072.18000000000;
    double a15_min  =-40603300;
    double a16_min  =-539011;
    double a17_min  =-157301;
    double a21_min  =-3.77131000000000e-06;
    double a22_min  =-10;
    double a23_min  =-0.0530223000000000;
    double a24_min  =-5.24883000000000e-05;
    double a25_min  =-1.74820000000000;
    double a26_min  =-5.02219000000000;
    double a27_min  =-0.0266189000000000;
    double a31_min  =-3.89337000000000e-05;
    double a32_min  =-0;
    double a33_min  =-10;
    double a34_min  =-0.000646268000000000;
    double a35_min  =-37.3379000000000;
    double a36_min  =-0.173007000000000;
    double a37_min  =-0.276693000000000;
    double a41_min  =-1.11368000000000;
    double a42_min  =-0;
    double a43_min  =-108.815000000000;
    double a44_min  =-17.5451000000000;
    double a45_min  =-1050090;
    double a46_min  =-6608.41000000000;
    double a47_min  =-7820.33000000000;
    double a51_min  =-3.81511000000000e-05;
    double a52_min  =-0;
    double a53_min  =-0.00226876000000000;
    double a54_min  =-0.00109013000000000;
    double a55_min  =-50.9192000000000;
    double a56_min  =-0.627949000000000;
    double a57_min  =-0.265072000000000;
    double a61_min  =-1.56492000000000e-05;
    double a62_min  =-0;
    double a63_min  =-0.365060000000000;
    double a64_min  =-0.000235869000000000;
    double a65_min  =-8.01958000000000;
    double a66_min  =-10.0574000000000;
    double a67_min  =-0.109882000000000;
    double a71_min  =-0;
    double a72_min  =-0;
    double a73_min  =-0;
    double a74_min  =-0;
    double a75_min  =-0;
    double a76_min  =-0;
    double a77_min  =-10;

    double a11_max  =32.6431000000000;
    double a12_max  =0;
    double a13_max  =936.482000000000;
    double a14_max  =1072.18000000000;
    double a15_max  =40603300;
    double a16_max  =539011;
    double a17_max  =157301;
    double a21_max  =3.77131000000000e-06;
    double a22_max  =10;
    double a23_max  =0.0530223000000000;
    double a24_max  =5.24883000000000e-05;
    double a25_max  =1.74820000000000;
    double a26_max  =5.02219000000000;
    double a27_max  =0.0266189000000000;
    double a31_max  =3.89337000000000e-05;
    double a32_max  =0;
    double a33_max  =10;
    double a34_max  =0.000646268000000000;
    double a35_max  =37.3379000000000;
    double a36_max  =0.173007000000000;
    double a37_max  =0.276693000000000;
    double a41_max  =1.11368000000000;
    double a42_max  =0;
    double a43_max  =108.815000000000;
    double a44_max  =17.5451000000000;
    double a45_max  =1050090;
    double a46_max  =6608.41000000000;
    double a47_max  =7820.33000000000;
    double a51_max  =3.81511000000000e-05;
    double a52_max  =0;
    double a53_max  =0.00226876000000000;
    double a54_max  =0.00109013000000000;
    double a55_max  =50.9192000000000;
    double a56_max  =0.627949000000000;
    double a57_max  =0.265072000000000;
    double a61_max  =1.56492000000000e-05;
    double a62_max  =0;
    double a63_max  =0.365060000000000;
    double a64_max  =0.000235869000000000;
    double a65_max  =8.01958000000000;
    double a66_max  =10.0574000000000;
    double a67_max  =0.109882000000000;
    double a71_max  =0;
    double a72_max  =0;
    double a73_max  =0;
    double a74_max  =0;
    double a75_max  =0;
    double a76_max  =0;
    double a77_max  =10;

    problem.phases(iphase).bounds.lower.states   << r_min, theta_min, phi_min, v_min, gamma_min, azim_min, p0_min, a11_min, a12_min, a13_min, a14_min, a15_min, a16_min, a17_min, a21_min, a22_min, a23_min, a24_min, a25_min, a26_min, a27_min, a31_min, a32_min, a33_min, a34_min, a35_min, a36_min, a37_min, a41_min, a42_min, a43_min, a44_min, a45_min, a46_min, a47_min, a51_min, a52_min, a53_min, a54_min, a55_min, a56_min, a57_min, a61_min, a62_min, a63_min, a64_min, a65_min, a66_min, a67_min, a71_min, a72_min, a73_min, a74_min, a75_min, a76_min, a77_min;
    problem.phases(iphase).bounds.upper.states   << r_max, theta_max, phi_max, v_max, gamma_max, azim_max, p0_max, a11_max, a12_max, a13_max, a14_max, a15_max, a16_max, a17_max, a21_max, a22_max, a23_max, a24_max, a25_max, a26_max, a27_max, a31_max, a32_max, a33_max, a34_max, a35_max, a36_max, a37_max, a41_max, a42_max, a43_max, a44_max, a45_max, a46_max, a47_max, a51_max, a52_max, a53_max, a54_max, a55_max, a56_max, a57_max, a61_max, a62_max, a63_max, a64_max, a65_max, a66_max, a67_max, a71_max, a72_max, a73_max, a74_max, a75_max, a76_max, a77_max;

    problem.phases(iphase).bounds.lower.controls << bank_min, k1_min;
    problem.phases(iphase).bounds.upper.controls << bank_max, k1_max;

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

    double a11_0    = 1;
    double a12_0    = 0;
    double a13_0    = 0;
    double a14_0    = 0;
    double a15_0    = 0;
    double a16_0    = 0;
    double a17_0    = 0;
    double a21_0    = 0;
    double a22_0    = 1;
    double a23_0    = 0;
    double a24_0    = 0;
    double a25_0    = 0;
    double a26_0    = 0;
    double a27_0    = 0;
    double a31_0    = 0;
    double a32_0    = 0;
    double a33_0    = 1;
    double a34_0    = 0;
    double a35_0    = 0;
    double a36_0    = 0;
    double a37_0    = 0;
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
    double a61_0    = 0;
    double a62_0    = 0;
    double a63_0    = 0;
    double a64_0    = 0;
    double a65_0    = 0;
    double a66_0    = 1;
    double a67_0    = 0;
    double a71_0    = 0;
    double a72_0    = 0;
    double a73_0    = 0;
    double a74_0    = 0;
    double a75_0    = 0;
    double a76_0    = 0;
    double a77_0    = 1;


    problem.phases(iphase).bounds.lower.states   << r_min, theta_min, phi_min, v_min, gamma_min, azim_min, p0_min, a11_min, a12_min, a13_min, a14_min, a15_min, a16_min, a17_min, a21_min, a22_min, a23_min, a24_min, a25_min, a26_min, a27_min, a31_min, a32_min, a33_min, a34_min, a35_min, a36_min, a37_min, a41_min, a42_min, a43_min, a44_min, a45_min, a46_min, a47_min, a51_min, a52_min, a53_min, a54_min, a55_min, a56_min, a57_min, a61_min, a62_min, a63_min, a64_min, a65_min, a66_min, a67_min, a71_min, a72_min, a73_min, a74_min, a75_min, a76_min, a77_min;
    problem.phases(iphase).bounds.upper.states   << r_max, theta_max, phi_max, v_max, gamma_max, azim_max, p0_max, a11_max, a12_max, a13_max, a14_max, a15_max, a16_max, a17_max, a21_max, a22_max, a23_max, a24_max, a25_max, a26_max, a27_max, a31_max, a32_max, a33_max, a34_max, a35_max, a36_max, a37_max, a41_max, a42_max, a43_max, a44_max, a45_max, a46_max, a47_max, a51_max, a52_max, a53_max, a54_max, a55_max, a56_max, a57_max, a61_max, a62_max, a63_max, a64_max, a65_max, a66_max, a67_max, a71_max, a72_max, a73_max, a74_max, a75_max, a76_max, a77_max;

    problem.phases(iphase).bounds.lower.controls << bank_min, k1_min;
    problem.phases(iphase).bounds.upper.controls << bank_max, k1_max;


    problem.phases(iphase).bounds.lower.events   << r_f, ra_f_min, a11_0, a12_0, a13_0, a14_0, a15_0, a16_0, a17_0, a21_0, a22_0, a23_0, a24_0, a25_0, a26_0, a27_0, a31_0, a32_0, a33_0, a34_0, a35_0, a36_0, a37_0, a41_0, a42_0, a43_0, a44_0, a45_0, a46_0, a47_0, a51_0, a52_0, a53_0, a54_0, a55_0, a56_0, a57_0, a61_0, a62_0, a63_0, a64_0, a65_0, a66_0, a67_0, a71_0, a72_0, a73_0, a74_0, a75_0, a76_0, a77_0;
    problem.phases(iphase).bounds.upper.events   << r_f, ra_f_max, a11_0, a12_0, a13_0, a14_0, a15_0, a16_0, a17_0, a21_0, a22_0, a23_0, a24_0, a25_0, a26_0, a27_0, a31_0, a32_0, a33_0, a34_0, a35_0, a36_0, a37_0, a41_0, a42_0, a43_0, a44_0, a45_0, a46_0, a47_0, a51_0, a52_0, a53_0, a54_0, a55_0, a56_0, a57_0, a61_0, a62_0, a63_0, a64_0, a65_0, a66_0, a67_0, a71_0, a72_0, a73_0, a74_0, a75_0, a76_0, a77_0;   
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
    
    
    MatrixXd u_guess    =  zeros(ncontrols,nnodes);
    u_guess.row(0)  = 20*M_PI/180*ones(1,nnodes);
    u_guess.row(1)  =0*ones(1,nnodes);

    problem.phases(iphase).guess.controls       = u_guess;
    problem.phases(iphase).guess.time           = linspace(0.0,110.0,nnodes);
    

    iphase = 2; 

    nnodes    			            = problem.phases(iphase).nodes(0);
    ncontrols                       = problem.phases(iphase).ncontrols;
    nstates                         = problem.phases(iphase).nstates;

    u_guess    =  zeros(ncontrols,nnodes);
    u_guess.row(0)  = 20*M_PI/180*ones(1,nnodes);
    u_guess.row(1)  =0*ones(1,nnodes);

    problem.phases(iphase).guess.controls       = u_guess;
    problem.phases(iphase).guess.time           = linspace(110,300,nnodes);


    // Define the state type for Boost ODEINT
    typedef vector<double> state_type;
    state_type y(56);

    y[0] = r_0    ;
    y[1] = theta_0;
    y[2] = phi_0  ;
    y[3] = v_0    ;
    y[4] = gamma_0;
    y[5] = azim_0 ;
    y[6] = p0_0   ;
    double smallInit = 1;
    y[7] = smallInit; y[8] = 0.0; y[9] = 0.0; y[10] = 0.0; y[11] = 0.0; y[12] =0.0; y[13] = 0.0; 
     y[14] = 0.0; y[15] = smallInit; y[16] = 0.0; y[17] = 0.0; y[18] = 0.0; y[19] =0.0; y[20] = 0.0;
     y[21] = 0.0; y[22] = 0.0; y[23] = smallInit; y[24] = 0.0; y[25] = 0.0; y[26] =0.0; y[27] = 0.0;
     y[28] = 0.0; y[29] = 0.0; y[30] = 0.0; y[31] =smallInit; y[32] = 0.0; y[33] =0.0; y[34] = 0.0;
     y[35] = 0.0; y[36] = 0.0; y[37] = 0.0; y[38] = 0.0; y[39] = smallInit; y[40] =0.0; y[41] = 0.0;
     y[42] = 0.0; y[43] = 0.0; y[44] = 0.0; y[45] = 0.0; y[46] = 0.0; y[47] =smallInit; y[48] = 0.0;
     y[49] = 0.0; y[50] = 0.0; y[51] = 0.0; y[52] = 0.0; y[53] = 0.0; y[54] =0.0; y[55] = smallInit;
    

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

    // Process each column starting from the second one
    for (int j = 1; j < x_guess_1.cols(); ++j) {
        if (is_column_zero(x_guess_1, j)) {
            replace_column_with_previous(x_guess_1, j);
        }
    }

    // Process each column starting from the second one
    for (int j = 1; j < x_guess_2.cols(); ++j) {
        if (is_column_zero(x_guess_2, j)) {
            replace_column_with_previous(x_guess_2, j);
        }
    }


    
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
    

    plot(t,r,problem.name, "time (s)", "Raduis (km)");

    plot(t,v,problem.name, "time (s)", "Velovity (km/sec)");

    plot(v,r,problem.name, "Velovity (km/sec)","Raduis (km)");

    plot(t,b*180/M_PI,problem.name,"time (s)", "Bank (deg)");
    
    plot(t,b*180/M_PI,problem.name,"time (s)", "u (dimensionless)", "Bank",
                               "pdf", "bank_control.pdf");

}
////////////////////////////////////////////////////////////////////////////
///////////////////////      END OF FILE     ///////////////////////////////
////////////////////////////////////////////////////////////////////////////
