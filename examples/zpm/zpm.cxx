//////////////////////////////////////////////////////////////////////////
////////////////              zpm.cxx                /////////////////////
//////////////////////////////////////////////////////////////////////////
////////////////           PSOPT Example             /////////////////////
//////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////
//////// Title: 	Zero propellant maneouvre problem  ///////////////
//////// Last modified: 09 November 2009                  ////////////////
//////// Reference:     S.A. Bhatt (2007), Masters Thesis,////////////////
////////                Rice University, Houston TX       ////////////////
//////////////////////////////////////////////////////////////////////////
////////     Copyright (c) Victor M. Becerra, 2009        ////////////////
//////////////////////////////////////////////////////////////////////////
//////// This is part of the PSOPT software library, which ///////////////
//////// is distributed under the terms of the GNU Lesser ////////////////
//////// General Public License (LGPL)                    ////////////////
//////////////////////////////////////////////////////////////////////////

#include "psopt.h"

using namespace PSOPT;


// Set CASE below to 1: 6000 s maneouvre, or to 2: 7200 maneouvre
#define CASE 2



struct Constants {
	MatrixXd J;
	double n;
	double Kp;
	double Kd;
	double hmax;
};

typedef struct Constants Constants_;

static Constants_ CONSTANTS;


void Tfun( adouble T[][3], adouble *q )
{

adouble q1 = q[0];
adouble q2 = q[1];
adouble q3 = q[2];
adouble q4 = q[3];

T[0][0] = -q2 ; T[0][1] = -q3; T[0][2] = -q4;

T[1][0] =  q1 ; T[1][1] = -q4; T[1][2] =  q3;

T[2][0] =  q4 ; T[2][1] =  q1; T[2][2] = -q2;

T[3][0] = -q3;  T[3][1] =  q2; T[3][2] =  q1;

}


void compute_omega0(adouble* omega0, adouble* q)
{
// This function computes the angular speed in the rotating LVLH reference frame
	int i;
	double n = CONSTANTS.n;
	adouble C2[3];

   adouble q1 = q[0];
   adouble q2 = q[1];
   adouble q3 = q[2];
   adouble q4 = q[3];

	C2[ 0 ] = 2*(q2*q3 + q1*q4);
	C2[ 1 ] = 1.0-2.0*(q2*q2+q4*q4);
	C2[ 2 ] = 2*(q3*q4-q1*q2);

	for (i=0;i<3;i++)   omega0[i] = -n*C2[i];

}

void compute_control_torque(adouble* u, adouble* q, adouble* qc, adouble* omega )
{
// This function computes the control torque
//

double Kp = CONSTANTS.Kp; // Proportional gain
double Kd = CONSTANTS.Kd; // Derivative gain
double n  = CONSTANTS.n;  // Orbital rotation rate [rad/s]


int i, j;

MatrixXd& J = CONSTANTS.J;


adouble T[4][3];

Tfun( T, q );

adouble Tc[4][3];

Tfun( Tc, qc );

adouble epsilon_tilde[3];

for(i=0;i<3;i++) {
	epsilon_tilde[i] = 0.0;
	for(j=0;j<4;j++) {
		epsilon_tilde[i] += 2*Tc[j][i]*q[j];
	}
}

adouble omega_c[3];

compute_omega0( omega_c, qc );

adouble omega_tilde[3];

for(i=0;i<3;i++) {
	omega_tilde[i] = omega[i]-omega_c[i];
}

adouble uaux[3];

for(i=0;i<3;i++) {
	uaux[i]= Kp*epsilon_tilde[i]+Kd*omega_tilde[i];
}

product_ad( J, uaux, 3, u );




}

void quarternion2Euler( MatrixXd& phi, MatrixXd& theta, MatrixXd& psi, MatrixXd& q)
{
// This function finds the Euler angles given the quarternion vector
//
	long N = q.cols();
	MatrixXd q0; q0 = q.row(0);
 	MatrixXd q1; q1 = q.row(1);
 	MatrixXd q2; q2 = q.row(2);
	MatrixXd q3; q3 = q.row(3);

	phi.resize(1,N);
	theta.resize(1,N);
	psi.resize(1,N);

   	for(int i=0;i<N;i++) {
   	      phi(i)=atan2( 2*(q0(i)*q1(i) + q2(i)*q3(i)), 1.0-2.0*(q1(i)*q1(i)+q2(i)*q2(i)) );
            theta(i)=asin( 2*(q0(i)*q2(i)-q3(i)*q1(i)) );
            psi(i) = atan2( 2*(q0(i)*q3(i)+q1(i)*q2(i)), 1.0-2.0*(q2(i)*q2(i)+q3(i)*q3(i)) );
	   }
}



void compute_aerodynamic_torque(adouble* tau_aero, adouble& time )
{
// This function approximates the aerodynamic torque by using the model and
// parameters given in the following reference:
// A. Chun Lee (2003) "Robust Momemtum Manager Controller for Space Station Applications".
// Master of Arts Thesis, Rice University.
//

	double alpha1[3] = {1.0, 4.0, 1.0};
   double alpha2[3] = {1.0, 2.0, 1.0};
   double alpha3[3] = {0.5, 0.5, 0.5};
	adouble t = time;

	double phi1 = 0.0;
	double phi2 = 0.0;

	double n = CONSTANTS.n;

 	for(int i=0;i<3;i++) {
                // Aerodynamic torque in [lb-ft]
		tau_aero[i] = alpha1[i] + alpha2[i]*sin( n*t + phi1 ) + alpha3[i]*sin( 2*n*t + phi2);
	}
}


//////////////////////////////////////////////////////////////////////////
///////////////////  Define the end point (Mayer) cost function //////////
//////////////////////////////////////////////////////////////////////////


adouble endpoint_cost(adouble* initial_states, adouble* final_states,
                      adouble* parameters,adouble& t0, adouble& tf,
                      adouble* xad, int iphase, Workspace* workspace)
{

	double end_point_weight = 0.1;

   adouble gamma     = parameters[ 0 ];

	return (end_point_weight*gamma);

}

//////////////////////////////////////////////////////////////////////////
///////////////////  Define the integrand (Lagrange) cost function  //////
//////////////////////////////////////////////////////////////////////////

adouble integrand_cost(adouble* states, adouble* controls, adouble* parameters,
                     adouble& time, adouble* xad, int iphase, Workspace* workspace)
{
		double running_cost_weight = 1.0;

		adouble q[4]; // quarternion vector
		adouble u[3]; // control torque

		q[0] = states[ 0 ];
		q[1] = states[ 1 ];
		q[2] = states[ 2 ];
		q[3] = states[ 3 ];

		adouble omega[3]; // angular rate vector

		omega[0] = states[ 4 ];
		omega[1] = states[ 5 ];
		omega[2] = states[ 6 ];


		adouble qc[4]; // control vector

		qc[0] = controls[ 0 ];
		qc[1] = controls[ 1 ];
		qc[2] = controls[ 2 ];
		qc[3] = controls[ 3 ];

		compute_control_torque(u,q,qc,omega);

		return running_cost_weight*dot(u,u,3);

}

//////////////////////////////////////////////////////////////////////////
///////////////////  Define the DAE's ////////////////////////////////////
//////////////////////////////////////////////////////////////////////////
void dae(adouble* derivatives, adouble* path, adouble* states,
         adouble* controls, adouble* parameters, adouble& time,
         adouble* xad, int iphase, Workspace* workspace)
{

int i,j;


double n = CONSTANTS.n; // Orbital rotation rate [rad/s]

adouble q[4]; // quarternion vector

q[0] = states[ 0 ];
q[1] = states[ 1 ];
q[2] = states[ 2 ];
q[3] = states[ 3 ];

adouble omega[3]; // angular rate vector

omega[0] = states[ 4 ];
omega[1] = states[ 5 ];
omega[2] = states[ 6 ];

adouble h[3]; // momentum vector

h[0] = states[ 7 ];
h[1] = states[ 8 ];
h[2] = states[ 9 ];

adouble qc[4]; // control vector

qc[0] = controls[ 0 ];
qc[1] = controls[ 1 ];
qc[2] = controls[ 2 ];
qc[3] = controls[ 3 ];

adouble C2[3], C3[3];

adouble u[3];

adouble gamma;


gamma = parameters[ 0 ];


// Inertia matrix in slug-ft^2

MatrixXd& J = CONSTANTS.J;

MatrixXd Jinv; Jinv = J.inverse();

adouble q1 = q[0];
adouble q2 = q[1];
adouble q3 = q[2];
adouble q4 = q[3];

C3[ 0 ] = 2*(q2*q4 - q1*q3);
C3[ 1 ] = 2*(q3*q4 + q1*q2);
C3[ 2 ] = 1.0-2.0*(q2*q2 + q3*q3);



adouble T[4][3];

Tfun( T, q );

adouble qdot[4];

adouble omega0[3];

compute_omega0( omega0, q);

// Quarternion attitude kinematics

for(j=0;j<4;j++) {

	qdot[j]=0;

	for(i=0;i<3;i++) {

  		qdot[j] += 0.5*T[j][i]*(omega[i]-omega0[i]);

	}

}


adouble Jomega[3];

product_ad( J, omega, 3, Jomega );

adouble omegaCrossJomega[3];

cross(omega,Jomega, omegaCrossJomega);

adouble F[3];

// Compute the torque disturbances:

adouble tau_grav[3], tau_aero[3];

adouble v1[3];

for(i=0;i<3;i++) {
  v1[i] = 3*pow(n,2)*C3[i];
}

adouble JC3[3];

product_ad( J, C3, 3, JC3 );

//gravity gradient torque

cross( v1, JC3, tau_grav );

//Aerodynamic torque
compute_aerodynamic_torque(tau_aero, time );

for(i=0;i<3;i++) {
  // Uncomment this section to ignore the aerodynamic disturbance torque
  tau_aero[i] = 0.0;
}
adouble tau_d[3];

for (i=0;i<3;i++) {

   tau_d[i] =  tau_grav[i] + tau_aero[i];

}

compute_control_torque(u, q, qc, omega );


for (i=0;i<3;i++) {
   F[i]     =  tau_d[i] - omegaCrossJomega[i] - u[i];
}

adouble omega_dot[3];

// Rotational dynamics

product_ad( Jinv, F, 3, omega_dot );

adouble OmegaCrossH[3];

cross( omega, h , OmegaCrossH );

adouble hdot[3];


//Momemtum derivative

for(i=0; i<3; i++) {
   hdot[i] = u[i] - OmegaCrossH[i];
}



derivatives[0] =       qdot[ 0 ];
derivatives[1] =       qdot[ 1 ];
derivatives[2] =       qdot[ 2 ];
derivatives[3] =       qdot[ 3 ];
derivatives[4] =       omega_dot[ 0 ];
derivatives[5] =       omega_dot[ 1 ];
derivatives[6] =       omega_dot[ 2 ];
derivatives[7] =       hdot[ 0 ];
derivatives[8] =       hdot[ 1 ];
derivatives[9]=        hdot[ 2 ];


path[ 0 ] =   dot( q,  q,  4);
path[ 1 ] =   dot( qc, qc, 4);
path[ 2 ] =   dot( h, h, 3 ) - gamma;   //  <= 0
path[ 3 ] =   dot( hdot, hdot,3); // <= hdotmax^2,

}







////////////////////////////////////////////////////////////////////////////
///////////////////  Define the events function ////////////////////////////
////////////////////////////////////////////////////////////////////////////

void events(adouble* e, adouble* initial_states, adouble* final_states,
            adouble* parameters,adouble& t0, adouble& tf, adouble* xad,
            int iphase, Workspace* workspace)

{

adouble q1_i        = initial_states[0];
adouble q2_i        = initial_states[1];
adouble q3_i        = initial_states[2];
adouble q4_i        = initial_states[3];
adouble omega1_i    = initial_states[4];
adouble omega2_i    = initial_states[5];
adouble omega3_i    = initial_states[6];
adouble h1_i        = initial_states[7];
adouble h2_i        = initial_states[8];
adouble h3_i        = initial_states[9];

adouble q1_f        = final_states[0];
adouble q2_f        = final_states[1];
adouble q3_f        = final_states[2];
adouble q4_f        = final_states[3];
adouble omega1_f    = final_states[4];
adouble omega2_f    = final_states[5];
adouble omega3_f    = final_states[6];
adouble h1_f        = final_states[7];
adouble h2_f        = final_states[8];
adouble h3_f        = final_states[9];

// Initial conditions

e[ 0 ] = q1_i;
e[ 1 ] = q2_i;
e[ 2 ] = q3_i;
e[ 3 ] = q4_i;
e[ 4 ] = omega1_i;
e[ 5 ] = omega2_i;
e[ 6 ] = omega3_i;
e[ 7 ] = h1_i;
e[ 8 ] = h2_i;
e[ 9 ] = h3_i;

// Final conditions

 e[ 10 ] = q1_f;
 e[ 11 ] = q2_f;
 e[ 12 ] = q3_f;
 e[ 13 ] = q4_f;
 e[ 14 ] = omega1_f;
 e[ 15 ] = omega2_f;
 e[ 16 ] = omega3_f;
 e[ 17 ] = h1_f;
 e[ 18 ] = h2_f;
 e[ 19 ] = h3_f;


}



///////////////////////////////////////////////////////////////////////////
///////////////////  Define the phase linkages function ///////////////////
///////////////////////////////////////////////////////////////////////////

void linkages( adouble* linkages, adouble* xad, Workspace* workspace)
{
    // Single phase
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



   CONSTANTS.Kp = 0.000128; // Proportional gain
   CONSTANTS.Kd = 0.015846; // Derivative gain

   double hmax; // maximum momentum magnitude  in [ft-lbf-sec]


   if (CASE==1) {
	CONSTANTS.n = 1.1461E-3; // Orbital rotation rate [rad/s]
	hmax = 4*3600.0; // 4 CMG's
   }
   else if (CASE==2) {
	CONSTANTS.n = 1.1475E-3;
	hmax = 3*3600.0; // 3 CMG's
   }

CONSTANTS.hmax = hmax;

   MatrixXd& J = CONSTANTS.J;

   J.resize(3,3);

   // Inertia matrix in slug-ft^2

   if (CASE==1) {
	J(0,0) = 17834580.0 ; J(0,1)= 2787992.0; J(0,2)= 2873636.0;
	J(1,0) = 2787992.0  ; J(1,1)= 2773815.0; J(1,2)= -863810.0;
	J(2,0) = 28736361.0 ; J(2,1)= -863810.0; J(2,2)=  38030467.0;
   }

   else if (CASE==2) {
	J(0,0) = 18836544.0 ; J(0,1)= 3666370.0;  J(0,2)= 2965301.0;
	J(1,0) = 3666370.0  ; J(1,1)= 27984088.0; J(1,2)= -1129004.0;
	J(2,0) = 2965301.0  ; J(2,1)= -1129004.0; J(2,2)=  39442649.0;
   }



////////////////////////////////////////////////////////////////////////////
///////////////////  Register problem name  ////////////////////////////////
////////////////////////////////////////////////////////////////////////////

    problem.name        		= "Zero Propellant Maneouvre of the ISS";
    problem.outfilename                 = "zpm.txt";

////////////////////////////////////////////////////////////////////////////
////////////  Define problem level constants & do level 1 setup ////////////
////////////////////////////////////////////////////////////////////////////

    problem.nphases   							= 1;
    problem.nlinkages                   	= 0;

    psopt_level1_setup(problem);
    
/////////////////////////////////////////////////////////////////////////////
/////////   Define phase related information & do level 2 setup  ////////////
/////////////////////////////////////////////////////////////////////////////


    problem.phases(1).nstates   				= 10;
    problem.phases(1).ncontrols 				= 4;
    problem.phases(1).nevents   	        	= 20;
    problem.phases(1).npath     				= 4;
    problem.phases(1).nodes               = (RowVectorXi(5) << 20, 30, 40, 50, 60).finished();   // <<  20, 30, 40, 50, 60;


    problem.phases(1).nparameters         = 1;


    psopt_level2_setup(problem, algorithm);

////////////////////////////////////////////////////////////////////////////
///////////////////  Enter problem bounds information //////////////////////
////////////////////////////////////////////////////////////////////////////

   // Control bounds

    problem.phases(1).bounds.lower.controls(0) =  -1.0;
    problem.phases(1).bounds.lower.controls(1) =  -1.0;
    problem.phases(1).bounds.lower.controls(2) =  -1.0;
    problem.phases(1).bounds.lower.controls(3) =  -1.0;


    problem.phases(1).bounds.upper.controls(0) =   1.0;
    problem.phases(1).bounds.upper.controls(1) =   1.0;
    problem.phases(1).bounds.upper.controls(2) =   1.0;
    problem.phases(1).bounds.upper.controls(3) =   1.0;

    // state bounds

    problem.phases(1).bounds.lower.states(0) =  -1.0;
    problem.phases(1).bounds.lower.states(1) =  -0.2;
    problem.phases(1).bounds.lower.states(2) =  -0.2;
    problem.phases(1).bounds.lower.states(3) =  -1.0;
    problem.phases(1).bounds.lower.states(4) =  -1.E-2;
    problem.phases(1).bounds.lower.states(5) =  -1.E-2;
    problem.phases(1).bounds.lower.states(6) =  -1.E-2;
    problem.phases(1).bounds.lower.states(7) =  -8000.0;
    problem.phases(1).bounds.lower.states(8) =  -8000.0;
    problem.phases(1).bounds.lower.states(9) =  -8000.0;

    problem.phases(1).bounds.upper.states(0) =  1.0;
    problem.phases(1).bounds.upper.states(1) =  0.2;
    problem.phases(1).bounds.upper.states(2) =  0.2;
    problem.phases(1).bounds.upper.states(3) =  1.0;
    problem.phases(1).bounds.upper.states(4) =  1.E-2;
    problem.phases(1).bounds.upper.states(5) =  1.E-2;
    problem.phases(1).bounds.upper.states(6) =  1.E-2;
    problem.phases(1).bounds.upper.states(7) =  8000.0;
    problem.phases(1).bounds.upper.states(8) =  8000.0;
    problem.phases(1).bounds.upper.states(9) =  8000.0;

   // Parameter bound


    problem.phases(1).bounds.lower.parameters(0) =  0.0;
    problem.phases(1).bounds.upper.parameters(0) =  hmax*hmax;


    // Event bounds

    // Initial / Final  condition values:
    MatrixXd q_i(4,1),  q_f(4,1), omega_i(3,1), omega_f(3,1), h_i(3,1), h_f(3,1);

    adouble q_ad[4], omega_ad[3];

// Initial conditions
    if (CASE==1) {
    	q_i(0) =  0.98966;
    	q_i(1) =  0.02690;
    	q_i(2) = -0.08246;
    	q_i(3) =  0.11425;

        q_ad[ 0 ]=q_i(0);
        q_ad[ 1 ]=q_i(1);
        q_ad[ 2 ]=q_i(2);
        q_ad[ 3 ]=q_i(3);
        compute_omega0( omega_ad, q_ad);
        omega_i(0) = omega_ad[ 0 ].value();
        omega_i(1) = omega_ad[ 1 ].value();
        omega_i(2) = omega_ad[ 2 ].value();

 //   	omega_i(1) = -2.5410E-4;
 //   	omega_i(2) = -1.1145E-3;
 //   	omega_i(3) =  8.2609E-5;

    	h_i(0) =  -496.0;
    	h_i(1) =  -175.0;
    	h_i(2) = -3892.0;
    }

    else if (CASE==2) {
	   q_i(0) =  0.98996;
    	q_i(1) =  0.02650;
    	q_i(2) = -0.07891;
    	q_i(3) =  0.11422;
        q_ad[ 0 ]=q_i(0);
        q_ad[ 1 ]=q_i(1);
        q_ad[ 2 ]=q_i(2);
        q_ad[ 3 ]=q_i(3);
        compute_omega0( omega_ad, q_ad);
        omega_i(0) = omega_ad[ 0 ].value();
        omega_i(1) = omega_ad[ 1 ].value();
        omega_i(2) = omega_ad[ 2 ].value();

    	h_i(0) =  1000.0;
    	h_i(1) =  -500.0;
    	h_i(2) = -4200.0;
    }

    // Final conditions
    q_f(0) =  0.70531;
    q_f(1) = -0.06201;
    q_f(2) = -0.03518;
    q_f(3) = -0.70531;

    q_ad[ 0 ]=q_f(0);
    q_ad[ 1 ]=q_f(1);
    q_ad[ 2 ]=q_f(2);
    q_ad[ 3 ]=q_f(3);
    compute_omega0( omega_ad, q_ad);
    omega_f(0) = omega_ad[ 0 ].value();
    omega_f(1) = omega_ad[ 1 ].value();
    omega_f(2) = omega_ad[ 2 ].value();

//    omega_f(1) = 1.1353E-3;
//    omega_f(2) = 3.0062E-6;
//    omega_f(3) = -1.5713E-4;

    h_f(0) = -9.0;
    h_f(1) = -3557.0;
    h_f(2) =  -135.0;

    double DQ = 0.0001;
    double DWF = 0.0;
    double DHF = 0.0;
    
    problem.phases(1).bounds.lower.events(0)  = q_i(0)-DQ;
    problem.phases(1).bounds.lower.events(1)  = q_i(1)-DQ;
    problem.phases(1).bounds.lower.events(2)  = q_i(2)-DQ;
    problem.phases(1).bounds.lower.events(3)  = q_i(3)-DQ;
    problem.phases(1).bounds.lower.events(4)  = omega_i(0);
    problem.phases(1).bounds.lower.events(5)  = omega_i(1);
    problem.phases(1).bounds.lower.events(6)  = omega_i(2);
    problem.phases(1).bounds.lower.events(7)  = h_i(0);
    problem.phases(1).bounds.lower.events(8)  = h_i(1);
    problem.phases(1).bounds.lower.events(9)  = h_i(2);
    problem.phases(1).bounds.lower.events(10) = q_f(0)-DQ;
    problem.phases(1).bounds.lower.events(11) = q_f(1)-DQ;
    problem.phases(1).bounds.lower.events(12) = q_f(2)-DQ;
    problem.phases(1).bounds.lower.events(13) = q_f(3)-DQ;
    problem.phases(1).bounds.lower.events(14) = omega_f(0)-DWF;
    problem.phases(1).bounds.lower.events(15) = omega_f(1)-DWF;
    problem.phases(1).bounds.lower.events(16) = omega_f(2)-DWF;
    problem.phases(1).bounds.lower.events(17) = h_f(0)-DHF;
    problem.phases(1).bounds.lower.events(18) = h_f(1)-DHF;
    problem.phases(1).bounds.lower.events(19) = h_f(2)-DHF;


    problem.phases(1).bounds.upper.events(0)  = q_i(0)+DQ;
    problem.phases(1).bounds.upper.events(1)  = q_i(1)+DQ;
    problem.phases(1).bounds.upper.events(2)  = q_i(2)+DQ;
    problem.phases(1).bounds.upper.events(3)  = q_i(3)+DQ;
    problem.phases(1).bounds.upper.events(4)  = omega_i(0);
    problem.phases(1).bounds.upper.events(5)  = omega_i(1);
    problem.phases(1).bounds.upper.events(6)  = omega_i(2);
    problem.phases(1).bounds.upper.events(7)  = h_i(0);
    problem.phases(1).bounds.upper.events(8)  = h_i(1);
    problem.phases(1).bounds.upper.events(9)  = h_i(2);
    problem.phases(1).bounds.upper.events(10) = q_f(0)+DQ;
    problem.phases(1).bounds.upper.events(11) = q_f(1)+DQ;
    problem.phases(1).bounds.upper.events(12) = q_f(2)+DQ;
    problem.phases(1).bounds.upper.events(13) = q_f(3)+DQ;
    problem.phases(1).bounds.upper.events(14) = omega_f(0)+DWF;
    problem.phases(1).bounds.upper.events(15) = omega_f(1)+DWF;
    problem.phases(1).bounds.upper.events(16) = omega_f(2)+DWF;
    problem.phases(1).bounds.upper.events(17) = h_f(0)+DHF;
    problem.phases(1).bounds.upper.events(18) = h_f(1)+DHF;
    problem.phases(1).bounds.upper.events(19) = h_f(2)+DHF;

// Path bounds

    double hdotmax = 200.0; // [ ft-lbf ]

    double EQ_TOL = 0.0002;

    problem.phases(1).bounds.lower.path(0) = 1.0-EQ_TOL;
    problem.phases(1).bounds.upper.path(0) = 1.0+EQ_TOL;

    problem.phases(1).bounds.lower.path(1) = 1.0-EQ_TOL;
    problem.phases(1).bounds.upper.path(1) = 1.0+EQ_TOL;


    problem.phases(1).bounds.lower.path(2) = -hmax*hmax;
    problem.phases(1).bounds.upper.path(2) =  0.0;


    problem.phases(1).bounds.lower.path(3) = 0.0;
    problem.phases(1).bounds.upper.path(3) = hdotmax*hdotmax;

    // Time bounds

    double TFINAL;

    if (CASE==1) {
	TFINAL = 6000.0;
    }
    else {
	TFINAL = 7200.0;
    }

    problem.phases(1).bounds.lower.StartTime    = 0.0;
    problem.phases(1).bounds.upper.StartTime    = 0.0;

    problem.phases(1).bounds.lower.EndTime      = TFINAL;
    problem.phases(1).bounds.upper.EndTime      = TFINAL;

////////////////////////////////////////////////////////////////////////////
///////////////////  Register problem functions  ///////////////////////////
////////////////////////////////////////////////////////////////////////////


    problem.integrand_cost      = &integrand_cost;
    problem.endpoint_cost 	     = &endpoint_cost;
    problem.dae             	  = &dae;
    problem.events 		        = &events;
    problem.linkages		        = &linkages;

////////////////////////////////////////////////////////////////////////////
///////////////////  Define & register initial guess ///////////////////////
////////////////////////////////////////////////////////////////////////////

   MatrixXd time_guess;      time_guess      = linspace(0.0, TFINAL, 50 );
   MatrixXd state_guess;     state_guess     = zeros(10,50);
   MatrixXd control_guess;   control_guess   = zeros(4,50);
   MatrixXd parameter_guess; parameter_guess = hmax*hmax*ones(1,1);

   control_guess.row(0) = linspace( q_i(0), q_i(0), 50 );
   control_guess.row(1) = linspace( q_i(1), q_i(1), 50 );
   control_guess.row(2) = linspace( q_i(2), q_i(2), 50 );
   control_guess.row(3) = linspace( q_i(3), q_i(3), 50 );

   state_guess.row(0) = linspace( q_i(0), q_i(0), 50);
   state_guess.row(1) = linspace( q_i(1), q_i(1), 50);
   state_guess.row(2) = linspace( q_i(2), q_i(2), 50);
   state_guess.row(3) = linspace( q_i(3), q_i(3), 50);

   state_guess.row(4) = linspace( omega_i(0), omega_f(0), 50);
   state_guess.row(5) = linspace( omega_i(1), omega_f(1), 50);
   state_guess.row(6) = linspace( omega_i(2), omega_f(2), 50);

   state_guess.row(7) = linspace( h_i(0), h_f(0), 50);
   state_guess.row(8) = linspace( h_i(1), h_f(1), 50);
   state_guess.row(9) = linspace( h_i(2), h_f(2), 50);



   problem.phases(1).guess.controls    = control_guess;
   problem.phases(1).guess.states      = state_guess;
   problem.phases(1).guess.time        = time_guess;
   problem.phases(1).guess.parameters  = parameter_guess;



////////////////////////////////////////////////////////////////////////////
///////////////////  Enter algorithm options  //////////////////////////////
////////////////////////////////////////////////////////////////////////////


    algorithm.nlp_iter_max                = 1000;
    algorithm.nlp_tolerance               = 1.e-5;
    algorithm.nlp_method                  = "IPOPT";
    algorithm.scaling                     = "automatic";
    algorithm.derivatives                 = "automatic";
    algorithm.defect_scaling              = "jacobian-based";
    algorithm.jac_sparsity_ratio          = 0.104;
    

    
    
////////////////////////////////////////////////////////////////////////////
///////////////////  Now call PSOPT to solve the problem   //////////////////
////////////////////////////////////////////////////////////////////////////

    psopt(solution, problem, algorithm);

////////////////////////////////////////////////////////////////////////////
///////////  Extract relevant variables from solution structure   //////////
////////////////////////////////////////////////////////////////////////////


    MatrixXd states, controls, t;

    states      = solution.get_states_in_phase(1);
    controls    = solution.get_controls_in_phase(1);
    t           = solution.get_time_in_phase(1);


////////////////////////////////////////////////////////////////////////////
///////////  Save solution data to files if desired ////////////////////////
////////////////////////////////////////////////////////////////////////////

    Save(states,"states.dat");
    Save(controls,"controls.dat");
    Save(t,"t.dat");

    MatrixXd omega, h, q, phi, theta, psi, qc, euler_angles;

    q      = states.block(0,0,4,length(t)); 
    omega  = states.block(4,0,3,length(t)); 
    h      = states.block(7,0,3,length(t)); 
    qc     = controls;

    quarternion2Euler(phi, theta, psi,  q);


    euler_angles.resize(3,length(t));

    euler_angles << phi , 
                    theta , 
                    psi;

    adouble  qc_ad[4], u_ad[3];

    MatrixXd u(3,length(t));
    MatrixXd hnorm(1,length(t));
    MatrixXd hi;
    MatrixXd hm = hmax*ones(1,length(t));


    int i,j;

  for (i=0; i< length(t); i++ ) {
	for(j=0;j<3;j++) {
	   omega_ad[j] = omega(j,i);
        }
	for(j=0;j<4;j++) {
	   q_ad[j] = q(j,i);
           qc_ad[j]= qc(j,i);
        }

	compute_control_torque(u_ad, q_ad, qc_ad, omega_ad );

	for(j=0;j<3;j++) {
	   u(j,i) = u_ad[j].value();
   }

   hi = h.col(i);  

	hnorm(0,i) = hi.norm();

    }

    omega = omega*(180.0/pi)*1000; // convert to mdeg/s

    phi = phi*180.0/pi; theta=theta*180.0/pi; psi=psi*180.0/pi;


    Save(u,"u.dat");

    Save(euler_angles,"euler_angles.dat");
    
   

///////////////////////////////////////////////////////////////////////////
///////////  Plot some results if desired (requires gnuplot) ///////////////
////////////////////////////////////////////////////////////////////////////


    plot(t,q,problem.name+" quarternion elements: q",  "time (s)", "q", "q");

    plot(t,qc,problem.name+" Control variables: qc",  "time (s)", "qc", "qc");

    plot(t,phi,problem.name+" Euler angles: phi",     "time (s)", "angles (deg)", "phi");

    plot(t,theta,problem.name+" Euler angles: theta", "time (s)", "angles (deg)", "theta");

    plot(t,psi,problem.name+" Euler angle: psi",       "time (s)", "psi (deg)", "psi");

    plot(t,omega.row(0),problem.name+": omega 1","time (s)", "omega1", "omega1");

    plot(t,omega.row(1),problem.name+": omega 2","time (s)", "omega2", "omega2");

    plot(t,omega.row(2),problem.name+": omega 3","time (s)", "omega3", "omega3");

    plot(t,h.row(0),problem.name+": momentum 1","time (s)", "h1", "h1");

    plot(t,h.row(1),problem.name+": momentum 2","time (s)", "h2", "h2");

    plot(t,h.row(2),problem.name+": momentum 3","time (s)", "h3", "h3");

    plot(t,u.row(0),problem.name+": control torque 1","time (s)", "u1", "u1");

    plot(t,u.row(1),problem.name+": control torque 2","time (s)", "u2", "u2");

    plot(t,u.row(2),problem.name+": control torque 3","time (s)", "u3", "u3");

    plot(t,hnorm,t,hm,problem.name+": momentum norm", "time (s)", "h", "h hmax");


    plot(t,phi,problem.name+" Euler angles: phi",     "time (s)", "angles (deg)", "phi",
	       "pdf", "zpm_phi.pdf" );

    plot(t,theta,problem.name+" Euler angles: theta", "time (s)", "angles (deg)", "theta",
	       "pdf", "zpm_theta.pdf");

    plot(t,psi,problem.name+" Euler angle: psi",       "time (s)", "psi (deg)", "psi",
	       "pdf", "zpm_psi.pdf");

    plot(t,omega.row(0),problem.name+": omega 1","time (s)", "omega1", "omega1",
	       "pdf", "zpm_omega1.pdf");

    plot(t,omega.row(1),problem.name+": omega 2","time (s)", "omega2", "omega2",
	       "pdf", "zpm_omega2.pdf");

    plot(t,omega.row(2),problem.name+": omega 3","time (s)", "omega3", "omega3",
	       "pdf", "zpm_omega3.pdf");

    plot(t,h.row(0),problem.name+": momentum 1","time (s)", "h1", "h1",
              "pdf", "zpm_h1.pdf");

    plot(t,h.row(1),problem.name+": momentum 2","time (s)", "h2", "h2",
	      "pdf", "zpm_h2.pdf");

    plot(t,h.row(2),problem.name+": momentum 3","time (s)", "h3", "h3",
	      "pdf", "zpm_h3.pdf");

    plot(t,u.row(0),problem.name+": control torque 1","time (s)", "u1", "u1",
	      "pdf", "zpm_u1.pdf");

    plot(t,u.row(1),problem.name+": control torque 2","time (s)", "u2", "u2",
	      "pdf", "zpm_u2.pdf");

    plot(t,u.row(2),problem.name+": control torque 3","time (s)", "u3", "u3",
	      "pdf", "zpm_u3.pdf");

    plot(t,hnorm,t,hm,problem.name+": momentum norm", "time (s)", "h (ft-lbf-sec)", "h hmax",
	      "pdf", "zpm_hnorm.pdf");

    plot(t,u,problem.name+": control torques","time (s)", "u (ft-lbf)", "u1 u2 u3",
	      "pdf", "zpm_controls.pdf");    
    
}

////////////////////////////////////////////////////////////////////////////
///////////////////////      END OF FILE     ///////////////////////////////
////////////////////////////////////////////////////////////////////////////
