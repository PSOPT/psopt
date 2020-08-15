//////////////////////////////////////////////////////////////////////////
////////////////           manutec.cxx               /////////////////////
//////////////////////////////////////////////////////////////////////////
////////////////           PSOPT example             /////////////////////
//////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////
//////// Title:         Minimum time control of a Manutec R3 robot ///////
//////// Last modified: 25 January 2010                   ////////////////
//////// Reference: Chettibi et al (2007)                 ////////////////
////////                                                  ////////////////
//////////////////////////////////////////////////////////////////////////
////////     Copyright (c) Victor M. Becerra, 2009        ////////////////
//////////////////////////////////////////////////////////////////////////
//////// This is part of the PSOPT software library, which ///////////////
//////// is distributed under the terms of the GNU Lesser ////////////////
//////// General Public License (LGPL)                    ////////////////
//////////////////////////////////////////////////////////////////////////

#include "psopt.h"




// Set option to:  1 for minimum time problem
//                 2 for minimum time with regularization
//                 3 for minimum energy problem
#define OBJ_OPTION 3



typedef struct {
 double   g; // Gravity constant
 double   pl; // Point mass of load;
 double   Ia1; // Moment of inertia of arm 1, element (3,3)
 MatrixXd* Fc; // Voltage-force constant of motor
 MatrixXd* r; // Gear ratio of motor
 MatrixXd* Im; // Moment of inertia of motor
 MatrixXd* m; // Mass of arm 2 and 3
 MatrixXd* L;   // Length of arm 2 and 3 (inc. tool)
 MatrixXd* com; // Center of mass coordinates of arm 2 and 3;
 MatrixXd* Ia; // Moment of inertia arm 2 and 3
} CONSTANTS_;

CONSTANTS_ CONSTANTS;

typedef struct {

adouble*   qdd_ad;
adouble*   F_ad;
adouble*   qd_ad;
adouble*   q_ad;

adouble* kw;
adouble* beta;
adouble* afor;
adouble* wabs;
adouble* zeta;
adouble* cosp;
adouble* ator;
adouble* sinp;
adouble* rhicm;
adouble* genvd ;
adouble* mrel12 ;
adouble* mrel22;
adouble* rhrel2;
adouble* workm3;
adouble* workm6;
adouble* works1;
adouble* works2;
adouble* workv3;
adouble* workv6;
adouble* cmhges;
adouble* ihhges;
adouble* inertc;
adouble* inertg;
adouble* inerth;
adouble* workam;
adouble* workas;
adouble* wworkm;
adouble* wworkv;

} VARS_;

VARS_ VARS;

void r3m2si(double *ml, adouble *aq, adouble *aqd,
	adouble *fgen, adouble *aqdd, VARS_ *vars);

//////////////////////////////////////////////////////////////////////////
///////////////////  Define the end point (Mayer) cost function //////////
//////////////////////////////////////////////////////////////////////////


adouble endpoint_cost(adouble* initial_states, adouble* final_states,
                      adouble* parameters,adouble& t0, adouble& tf,
                      adouble* xad, int iphase, Workspace* workspace)
{
   adouble retval;

   if (OBJ_OPTION==1 || OBJ_OPTION==2 ) retval = tf;
   else retval = 0.0;

   return retval;
}

//////////////////////////////////////////////////////////////////////////
///////////////////  Define the integrand (Lagrange) cost function  //////
//////////////////////////////////////////////////////////////////////////

adouble integrand_cost(adouble* states, adouble* controls, adouble* parameters,
                     adouble& time, adouble* xad, int iphase, Workspace* workspace)
{

   adouble retval;

   adouble* u = controls;

   adouble* qdot = &states[3];

   double rho;


   if (OBJ_OPTION==1)  rho = 0.0;
   else if (OBJ_OPTION==2) rho = 1.e-5;
   else if (OBJ_OPTION==3) rho = 1.0;


   retval = rho*dot( u, u, 3 );


   return (retval);
}


//////////////////////////////////////////////////////////////////////////
///////////////////  Define the DAE's ////////////////////////////////////
//////////////////////////////////////////////////////////////////////////




void dae(adouble* derivatives, adouble* path, adouble* states,
         adouble* controls, adouble* parameters, adouble& time,
         adouble* xad, int iphase, Workspace* workspace)
{

 int i;
 double g    =  CONSTANTS.g;
 MatrixXd& Fc = *CONSTANTS.Fc;
 MatrixXd& r  = *CONSTANTS.r;
 MatrixXd& Im = *CONSTANTS.Im;
 MatrixXd& m  = *CONSTANTS.m;
 double pl   =  CONSTANTS.pl;
 MatrixXd& L  = *CONSTANTS.L;
 MatrixXd& com= *CONSTANTS.com;
 double   Ia1=  CONSTANTS.Ia1;
 MatrixXd& Ia = *CONSTANTS.Ia;


adouble* F = VARS.F_ad;
adouble* qdd = VARS.qdd_ad;
adouble* qd = VARS.qd_ad;
adouble* q = VARS.q_ad;
double ml = pl;

for(i=0;i<3;i++) {
   q[i]  = states[i];
   qd[i] = states[3+i];
   F[i] = Fc(i)*controls[i];
}

r3m2si(&ml, q, qd, F, qdd, &VARS);

derivatives[0] =  qd[0];
derivatives[1] =  qd[1];
derivatives[2] =  qd[2];

derivatives[3]  =  qdd[0];
derivatives[4]  =  qdd[1];
derivatives[5]  =  qdd[2];

}


////////////////////////////////////////////////////////////////////////////
///////////////////  Define the events function ////////////////////////////
////////////////////////////////////////////////////////////////////////////

void events(adouble* e, adouble* initial_states, adouble* final_states,
            adouble* parameters,adouble& t0, adouble& tf, adouble* xad,
            int iphase, Workspace* workspace)

{


    int i;

    for (i=0; i< 6; i++ ) {
       e[ i ]     =  initial_states[ i ];
    }

    for (i=0; i< 6; i++ ) {
       e[6 + i ] =  final_states[ i ];
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
 double g = 9.81;                                     // Gravity constant [m/s2]
 double pl=  0.0;                                     // Point mass of load; [kg]
 double   Ia1 = 1.16;                                 // Moment of inertia of arm 1, element (3,3) [kgm^2]
 MatrixXd Fc(3,1); Fc << -126.0,252.0, 72.0;        // Voltage-force constant of motor [N*m/V]
 MatrixXd r(3,1); r << -105, 210, 60;               // Gear ratio of motor
 MatrixXd Im(3,1); Im << 1.3e-3, 1.3e-3, 1.3e-3;    // Moment of inertia of motor [kg*m^2]
 MatrixXd m(2,1);  m << 56.5, 60.3;                 // Mass of arm 2 and 3 [kg]
 MatrixXd L(2,1);  L << 0.5, 0.98;                  // Length of arm 2 and 3 (inc. tool) [m]
 MatrixXd  com(2,2);                                   // Center of mass coordinates of arm 2 and 3; [m]
 MatrixXd Ia(4,2);                                     // Moment of inertia arm 2 and 3 [kg*m^2]

 com(0,0)=0.172; com(0,1)=0.028; com(1,0)=0.205; com(1,1)=0.202;

 Ia(0,0)=2.58;  Ia(0,1)=11.0;
 Ia(1,0)=2.73;  Ia(1,1)=8.0;
 Ia(2,0)=0.64;  Ia(2,1)=0.80;
 Ia(3,0)=-0.46; Ia(3,1)=0.50;

 Print(Fc,"Fc");
 Print(r,"r");
 Print(Im,"Im");
 Print(m,"m");
 Print(L,"L");
 Print(com,"com");
 Print(Ia,"Ia");


 CONSTANTS.g   = g;
 CONSTANTS.pl  = pl;
 CONSTANTS.Ia1 = Ia1;
 CONSTANTS.Fc = &Fc;
 CONSTANTS.r  = &r;
 CONSTANTS.Im = &Im;
 CONSTANTS.m  = &m;
 CONSTANTS.L  = &L;
 CONSTANTS.com= &com;
 CONSTANTS.Ia = &Ia;


 VARS.qdd_ad = new adouble[3];
 VARS.qd_ad  = new adouble[3];
 VARS.q_ad   = new adouble[3];
 VARS.F_ad   = new adouble[3];



 VARS.kw = new adouble[9];
 VARS.beta = new adouble[18];
 VARS.afor = new adouble[9];
 VARS.wabs = new adouble[9];
 VARS.zeta = new adouble[18];
 VARS.cosp = new adouble[3];
 VARS.ator = new adouble[9];
 VARS.sinp = new adouble[3];
 VARS.rhicm = new adouble[9];
 VARS.genvd = new adouble[3];
 VARS.mrel12 = new adouble[18];
 VARS.mrel22 = new adouble[3];
 VARS.rhrel2 = new adouble[3];
 VARS.workm3 = new adouble[9];
 VARS.workm6 = new adouble[180];
 VARS.works1 = new adouble[12];
 VARS.works2 = new adouble[6];
 VARS.workv3 = new adouble[39];
 VARS.workv6 = new adouble[24];
 VARS.cmhges = new adouble[3];
 VARS.ihhges = new adouble[9];
 VARS.inertc = new adouble[27];
 VARS.inertg = new adouble[108];
 VARS.inerth = new adouble[27];
 VARS.workam = new adouble[144];
 VARS.workas = new adouble[40];
 VARS.wworkm = new adouble[54];
 VARS.wworkv = new adouble[9];


////////////////////////////////////////////////////////////////////////////
///////////////////  Declare key structures ////////////////////////////////
////////////////////////////////////////////////////////////////////////////

    Alg  algorithm;
    Sol  solution;
    Prob problem;

////////////////////////////////////////////////////////////////////////////
///////////////////  Register problem name  ////////////////////////////////
////////////////////////////////////////////////////////////////////////////

    problem.name        		= "Manutec R3 robot problem";
    problem.outfilename                 = "manutec.txt";

////////////////////////////////////////////////////////////////////////////
////////////  Define problem level constants & do level 1 setup ////////////
////////////////////////////////////////////////////////////////////////////

    problem.nphases   			          = 1;
    problem.nlinkages                   = 0;

    psopt_level1_setup(problem);

/////////////////////////////////////////////////////////////////////////////
/////////   Define phase related information & do level 2 setup  ////////////
/////////////////////////////////////////////////////////////////////////////

    problem.phases(1).nstates   		= 6;
    problem.phases(1).ncontrols 		= 3;
    problem.phases(1).nevents   		= 12;
    problem.phases(1).npath     		= 0;
    problem.phases(1).nodes         = (RowVectorXi(5) << 20, 30, 40, 60, 80).finished();   


    psopt_level2_setup(problem, algorithm);



////////////////////////////////////////////////////////////////////////////
///////////////////  Enter problem bounds information //////////////////////
////////////////////////////////////////////////////////////////////////////

    MatrixXd qmax(3,1), qdmax(3,1), qddmax(3,1);

    // Joint position limits in rad

    qmax(0) = 2.97;    qmax(1) = 2.01; qmax(2) = 2.86;

    // Joint angular velocity limits in rad/s

    qdmax(0) = 3.0;   qdmax(1) = 1.5; qdmax(2) = 5.2;


    problem.phases(1).bounds.lower.states << -qmax , -qdmax;
    problem.phases(1).bounds.upper.states <<  qmax ,  qdmax;

    MatrixXd umax(3,1);

    // Control variable limits in V

    umax(0) = 7.5; umax(1) = 7.5; umax(2) = 7.5;


    problem.phases(1).bounds.lower.controls = -umax;

    problem.phases(1).bounds.upper.controls =  umax;



    MatrixXd qi(3,1), qf(3,1), qdi(3,1), qdf(3,1);

    // Initial joint positions in rad
    qi(0) = 0.0;    qi(1) = -1.5; qi(2) = 0.0;

    // Final joint positions in rad
    qf(0) = 1.0; qf(1) = -1.95;  qf(2) = 1.0;

    // Initial joint velocities in rad/s

    qdi = zeros(3,1);

    // Final joint velocities in rad/s

    qdf = zeros(3,1);


    problem.phases(1).bounds.lower.events << qi , qdi , qf , qdf;
    problem.phases(1).bounds.upper.events << qi , qdi , qf , qdf;


    problem.phases(1).bounds.lower.StartTime    = 0.0;
    problem.phases(1).bounds.upper.StartTime    = 0.0;

    if (OBJ_OPTION==1 || OBJ_OPTION==2) {
      problem.phases(1).bounds.lower.EndTime      = 0.0;
      problem.phases(1).bounds.upper.EndTime      = 1.0;
    }
    else {
      problem.phases(1).bounds.lower.EndTime      = 0.53;
      problem.phases(1).bounds.upper.EndTime      = 0.53;
    }




////////////////////////////////////////////////////////////////////////////
///////////////////  Register problem functions  ///////////////////////////
////////////////////////////////////////////////////////////////////////////


    problem.integrand_cost 		= &integrand_cost;
    problem.endpoint_cost 			= &endpoint_cost;
    problem.dae             		= &dae;
    problem.events 					= &events;
    problem.linkages					= &linkages;

////////////////////////////////////////////////////////////////////////////
///////////////////  Define & register initial guess ///////////////////////
////////////////////////////////////////////////////////////////////////////

    int nnodes    			= problem.phases(1).nodes(0);
    int ncontrols          = problem.phases(1).ncontrols;
    int nstates            = problem.phases(1).nstates;

    MatrixXd state_guess    =  zeros(nstates,nnodes);
    MatrixXd control_guess  =  zeros(ncontrols,nnodes);
    MatrixXd time_guess     =  linspace(0.0,0.53,nnodes);


    state_guess << 	linspace( qi(0), qf(0), nnodes ),
     						linspace( qi(1),  qf(1),  nnodes ),
     						linspace( qi(2),  qf(2),  nnodes ),
     						linspace( qdi(0), qdf(0), nnodes ),
     						linspace( qdi(1), qdf(1), nnodes ),
     						linspace( qdi(2), qdf(2), nnodes );


    problem.phases(1).guess.states   = state_guess;
    problem.phases(1).guess.controls = control_guess;
    problem.phases(1).guess.time     = time_guess;



////////////////////////////////////////////////////////////////////////////
///////////////////  Enter algorithm options  //////////////////////////////
////////////////////////////////////////////////////////////////////////////


    algorithm.nlp_iter_max                = 1000;
    algorithm.nlp_tolerance               = 1.e-6;
    algorithm.nlp_method                  = "IPOPT";
    algorithm.scaling                     = "automatic";
    algorithm.derivatives                 = "automatic";
    algorithm.mesh_refinement             = "automatic";
    algorithm.ode_tolerance               = 1.e-5;
    algorithm.mr_max_iterations           = 5;


////////////////////////////////////////////////////////////////////////////
///////////////////  Now call PSOPT to solve the problem   //////////////////
////////////////////////////////////////////////////////////////////////////

    psopt(solution, problem, algorithm);

////////////////////////////////////////////////////////////////////////////
///////////  Extract relevant variables from solution structure   //////////
////////////////////////////////////////////////////////////////////////////


    MatrixXd x, u, t;

    x           = solution.get_states_in_phase(1);
    u           = solution.get_controls_in_phase(1);
    t           = solution.get_time_in_phase(1);




////////////////////////////////////////////////////////////////////////////
///////////  Save solution data to files if desired ////////////////////////
////////////////////////////////////////////////////////////////////////////

    Save(x,"x.dat");
    Save(u,"u.dat");
    Save(t,"t.dat");

    MatrixXd q = x.block(0,0,3,length(t)); 
    MatrixXd qd= x.block(3,0,3,length(t)); 


////////////////////////////////////////////////////////////////////////////
///////////  Plot some results if desired (requires gnuplot) ///////////////
////////////////////////////////////////////////////////////////////////////



    plot(t,q,problem.name+": positions","t", "q (rad)", "q1 q2 q3");

    plot(t,qd,problem.name+": velocities","t", "qdot (rad/s)", "qd1 qd2 qd3");

    plot(t,u.row(0),problem.name+": control u1","time (s)", "control", "u1");

    plot(t,u.row(1),problem.name+": control u2","time (s)", "control", "u2");

    plot(t,u.row(2),problem.name+": control u3","time (s)", "control", "u3");

    plot(t,q,problem.name+": positions","t", "q (rad)", "q1 q2 q3", "pdf",
	                     "positions.pdf");

    plot(t,qd,problem.name+": velocities","t", "qdot (rad/s)", "qd1 qd2 qd3", "pdf",
	                      "velocities.pdf");

    plot(t,u,problem.name+": controls","time (s)", "controls", "u1 u2 u3",
	 "pdf", "controls.pdf");


}


 void r3m2si(double *ml, adouble *aq, adouble *aqd,
             adouble *fgen, adouble *aqdd, VARS_ *vars)
{

adouble mhges;

adouble* kw = vars->kw;
adouble* beta = vars->beta;
adouble* afor = vars->afor;
adouble* wabs = vars->wabs;
adouble* zeta = vars->zeta;
adouble* cosp = vars->cosp;
adouble* ator = vars->ator;
adouble* sinp = vars->sinp;
adouble* rhicm = vars->rhicm;
adouble* genvd = vars->genvd;
adouble* mrel12 = vars->mrel12;
adouble* mrel22 = vars->mrel22;
adouble* rhrel2 = vars->rhrel2;
adouble* workm3 = vars->workm3;
adouble* workm6 = vars->workm6;
adouble* works1 = vars->works1;
adouble* works2 = vars->works2;
adouble* workv3 = vars->workv3;
adouble* workv6 = vars->workv6;
adouble* cmhges = vars->cmhges;
adouble* ihhges = vars->ihhges;
adouble* inertc = vars->inertc;
adouble* inertg = vars->inertg;
adouble* inerth = vars->inerth;
adouble* workam = vars->workam;
adouble* workas = vars->workas;
adouble* wworkm = vars->wworkm;
adouble* wworkv = vars->wworkv;



/* Simulation model equations of the Manutec r3 robot (DFVLR model 2) */
/* This is a modified C++ tranlation of the original Fortran subroutine R3M2SI */
/* by Otter and co-workers, published here with kind persmission from Dr. Martin Otter */
/* DSL, Germany */

/* Procedure purpose: */
/*   This subroutine calculates the generalized accelerations (AQDD) for */
/*   the Manutec r3 robot of the DFVLR model 2. This model is based on */
/*   the following assumptions: */

/*   - The last 3 joints (joints 4,5,6) of the robot hand do not move. */
/*     This corresponds to the case that the brakes in these joints are */
/*     blocked (joints 4,5,6 are taken into account in model 1). */

/*   - The robot consists of a base body (the environment), three arms */
/*     and three rotors. The rotors represent the inertia effects of the */
/*     motors and of the wheels of the gear boxes. The rotors are */
/*     embedded in the preceeding arms, e.g. the rotor 2, which drives */
/*     arm2, is embedded in arm1. */

/*   - Arms and rotors are considered to be rigid bodies. */

/*   - Elasticity, backlash and friction in the gear boxes are neglected. */

/*   - The motors are modelled as ideal voltage-torque converters */
/*     without any dynamic effects. As input arguments of the subroutine */
/*     the torques at the gear output (FGEN in <Nm>) must be given. If the */
/*     voltage (U in <V>) at the input of the current regulator of the */
/*     motor is given, FGEN must be calculated (before calling this */
/*     subroutine) in the following way: */
/*         FGEN(1) = -126.0*U(1) */
/*         FGEN(2) =  252.0*U(2) */
/*         FGEN(3) =   72.0*U(3) */

/*   - At the robot's tip a load mass (ML) is attached. It can range */
/*     between 0 ... 15 kg. */

/*   Given the actual value ML of the load mass, the joint angles AQ(i), */
/*   the derivatives of the joint angles AQD(i), and the driving torques */
/*   in the joints FGEN(i), the second derivatives of the joint angles */
/*   AQDD(i) are calculated by this subroutine. */

/* Usage: */
/*   CALL  R3M2SI (ML, AQ, AQD, FGEN, AQDD) */

/*   ML    :  IN, DOUBLE, <kg>, 0<=ML<=15 */
/*            Load mass. */
/*   AQ    :  IN, DOUBLE(3), <rad>, -2.97 <= AQ(1) <= 2.97  (+/- 170 deg), */
/*                                  -2.01 <= AQ(2) <= 2.01  (+/- 115 deg), */
/*                                  -2.86 <= AQ(3) <= 2.86  (+/- 164 deg) */
/*            (The given limits are hardware constraints by limit switches). */
/*            Vector of generalized coordinates (relative angles between */
/*            two contigous robot arms). */
/*            AQ(1) is not used in this subroutine, because it is not */
/*            needed to calculate AQDD. */
/*   AQD   :  IN, DOUBLE(3), <rad/sec>, */
/*                                -3.0 <= AQD(1) <= 3.0  (+/- 172 deg/sec), */
/*                                -1.5 <= AQD(2) <= 1.5  (+/-  86 deg/sec), */
/*                                -5.2 <= AQD(3) <= 5.2  (+/- 298 deg/sec) */
/*            Derivative of AQ. */
/*   FGEN  :  IN, DOUBLE(3), <Nm>,     -945.0 <= FGEN(1) <=  945.0, */
/*                                    -1890.0 <= FGEN(2) <= 1890.0, */
/*                                     -540.0 <= FGEN(3) <=  540.0 */
/*            Torque at the gear output. */
/*   AQDD  :  OUT, DOUBLE(3), <rad/sec**2> */
/*            Second derivative of AQ. */

/* Bibliography: */
/*   A detailed description of the model, together with the body-data */
/*   (mass, center of mass etc. of the arms and rotors) is given in */
/*     Otter M., Tuerk S., Mathematical Model of the Manutec r3 Robot, */
/*            (DFVLR Model No. 2), DFVLR - Oberpfaffenhofen, Institut fuer */
/*            Dynamik der Flugsysteme, D-8031 Wessling, West Germany, */
/*            corrected version april 1988. */

/*   This subroutine was generated by the program MYROBOT. See */
/*     Otter M., Schlegel S., Symbolic generation of efficient simulation */
/*            codes for robots. 2nd European Simulation Multiconference, */
/*            June 1-3, 1988, Nice. */

/*   The underlying multibody algorithm is a modified version of */
/*     Brandl H., Johanni R., Otter M., A very efficient algorithm for the */
/*            simulation of robots and similar multibody systems without */
/*            inversion of the mass matrix. IFAC/IFIP/IMACS International */
/*            Symposium on Theory of Robots, december 3-5, 1986, Vienna. */

/* Remarks: */
/*   - The limits given for the input variables are not checked in */
/*     this subroutine. The user is responsible for proper data. */

/*   - If a SINGLE PRECISION version of this subroutine is desired, just */
/*     change all strings from DOUBLE PRECISION to REAL in the declaration */
/*     part. */

/* Copyright: */
/*   1988  DFVLR - Institut fuer Dynamik der Flugsysteme */

/* Life cycle: */
/*   1988 APR  M. Otter, S. Tuerk (DFVLR)             : specified. */
/*   1988 APR  M. Otter (DFVLR)                       : generated. */
/*   1988 APR  M. Otter, C. Schlegel, S. Tuerk (DFVLR): tested. */
/* ----------------------------------------------------------------------- */

/*  Statistical information (MySymbol-Version 1.2) */
/*  ============================================== */

/*  1. Number of Operations: */
/*  ------------------------ */

/*                            + - |  *  |  /  | sin | cos | sqrt| */
/*                           -----|-----|-----|-----|-----|-----| */
/*  Without simplifications   3180| 3702|   26|    3|    3|    0| */
/*                           -----|-----|-----|-----|-----|-----| */
/*  Terms simplified           222|  247|   23|    3|    3|    0| */
/*                           -----|-----|-----|-----|-----|-----| */
/*  Unnecessary statements        |     |     |     |     |     | */
/*  removed                    140|  159|   17|    2|    2|    0| */


/*  2. Storage information: */
/*  ----------------------- */

/*                      INTEGER |         | */
/*                       words  |         | */
/*                      --------|---------| */
/*  available storage    100000 | 100.0 % | */
/*                      --------|---------| */
/*  array-elements         1412 |   1.4 % | */
/*                      --------|---------| */
/*  double-numbers          104 |   0.1 % | */
/*                      --------|---------| */
/*  statement buffer       2552 |   2.6 % | */
/*  +++++++++++++++++++++++++++++++++++++++ */
/*  used storage           4068 |   4.1 % | */
/*                      --------|---------| */
/*  free storage          95932 |  95.9 % | */

/* ----------------------------------------------------------------------- */
/* cccccccccccccccccccccccc Procedural section cccccccccccccccccccccccccc */
    /* Parameter adjustments */
    --aqdd;
    --fgen;
    --aqd;
    --aq;

    /* Function Body */
    workv3[14] = *ml * .98;
    workv3[17] = workv3[14] + 12.1806;
    mhges = *ml + 60.3;
    cmhges[0] = 1.6884 / mhges;
    cmhges[2] = workv3[17] / mhges;
    workv3[20] = workv3[14] * .98;
    ihhges[0] = workv3[20] + 7.8704812;
    ihhges[4] = workv3[20] + 8.1077564;
    sinp[1] = sin(aq[2]);
    cosp[1] = cos(aq[2]);
    wworkv[0] = mhges * cmhges[0];
    wworkv[2] = mhges * cmhges[2];
    wworkv[3] = cmhges[0] * wworkv[0];
    wworkv[5] = cmhges[2] * wworkv[2];
    wworkv[6] = wworkv[3] + wworkv[5];
    wworkm[0] = wworkv[6] - wworkv[3];
    wworkm[8] = wworkv[6] - wworkv[5];
    wworkm[6] = -wworkv[0] * cmhges[2];
    inertc[18] = ihhges[0] - wworkm[0];
    inertc[22] = ihhges[4] - wworkv[6];
    inertc[24] = -.0110568 - wworkm[6];
    inertc[26] = .4372752 - wworkm[8];
    sinp[2] = sin(aq[3]);
    cosp[2] = cos(aq[3]);
/* ** The equations of motion have been generated with the algorithm */
/* ** of Brandl/Johanni/Otter (variant 6) */
/* ** Forward recursion -------------------------------------------------- */
/*  -- Quantities of body 1 */
/*  -- Quantities of body 2 */
    wabs[4] = sinp[1] * aqd[1];
    wabs[5] = cosp[1] * aqd[1];
    workv3[1] = aqd[1] * aqd[2];
    zeta[7] = cosp[1] * workv3[1];
    zeta[8] = -sinp[1] * workv3[1];
    kw[3] = aqd[2] * 4.9544125 - wabs[5] * 2.45219;
    kw[4] = wabs[4] * 6.7759085;
    kw[5] = aqd[2] * -2.45219 + wabs[5] * 2.311496;
/*  -- Quantities of body 3 */
    workv3[19] = -sinp[2] * cmhges[2];
    workv3[20] = cosp[2] * cmhges[2];
    workm3[4] = cosp[2] * inertc[22];
    workm3[5] = sinp[2] * inertc[22];
    workm3[7] = -sinp[2] * inertc[26];
    workm3[8] = cosp[2] * inertc[26];
    inertc[21] = -inertc[24] * sinp[2];
    inertc[22] = workm3[4] * cosp[2] - workm3[7] * sinp[2];
    inertc[24] *= cosp[2];
    inertc[25] = workm3[4] * sinp[2] + workm3[7] * cosp[2];
    inertc[26] = workm3[5] * sinp[2] + workm3[8] * cosp[2];
    workv3[21] = mhges * cmhges[0];
    workv3[22] = mhges * workv3[19];
    workv3[23] = mhges * workv3[20];
    workv3[24] = cmhges[0] * workv3[21];
    workv3[25] = workv3[19] * workv3[22];
    workv3[26] = workv3[20] * workv3[23];
    workv3[27] = workv3[24] + workv3[25] + workv3[26];
    inerth[18] = workv3[27] - workv3[24];
    inerth[22] = workv3[27] - workv3[25];
    inerth[26] = workv3[27] - workv3[26];
    inerth[21] = -workv3[21] * workv3[19];
    inerth[24] = -workv3[21] * workv3[20];
    inerth[25] = -workv3[22] * workv3[20];
    inerth[18] += inertc[18];
    inerth[21] += inertc[21];
    inerth[22] += inertc[22];
    inerth[24] += inertc[24];
    inerth[25] += inertc[25];
    inerth[26] += inertc[26];
    wabs[6] = aqd[3] + aqd[2];
    workv3[1] = wabs[5] * aqd[3];
    workv3[2] = -wabs[4] * aqd[3];
    workv3[6] = wabs[4] * .5;
    workv3[7] = -aqd[2] * .5;
    workv3[9] = -wabs[5] * workv3[7];
    workv3[10] = wabs[5] * workv3[6];
    workv3[11] = aqd[2] * workv3[7] - wabs[4] * workv3[6];
    rhicm[6] = mhges * cmhges[0];
    rhicm[7] = mhges * workv3[19];
    rhicm[8] = mhges * workv3[20];
    kw[6] = inerth[18] * wabs[6] + inerth[21] * wabs[4] + inerth[24] * wabs[5]
	    ;
    kw[7] = inerth[21] * wabs[6] + inerth[22] * wabs[4] + inerth[25] * wabs[5]
	    ;
    kw[8] = inerth[24] * wabs[6] + inerth[25] * wabs[4] + inerth[26] * wabs[5]
	    ;
    ator[3] = -wabs[4] * kw[5] + wabs[5] * kw[4];
    ator[4] = -wabs[5] * kw[3] + aqd[2] * kw[5];
    ator[5] = -aqd[2] * kw[4] + wabs[4] * kw[3];
    ator[6] = -wabs[4] * kw[8] + wabs[5] * kw[7];
    ator[7] = -wabs[5] * kw[6] + wabs[6] * kw[8];
    ator[8] = -wabs[6] * kw[7] + wabs[4] * kw[6];
    workv3[15] = wabs[4] * rhicm[8] - wabs[5] * rhicm[7];
    workv3[16] = wabs[5] * rhicm[6] - wabs[6] * rhicm[8];
    workv3[17] = wabs[6] * rhicm[7] - wabs[4] * rhicm[6];
    afor[6] = -wabs[4] * workv3[17] + wabs[5] * workv3[16];
    afor[7] = -wabs[5] * workv3[15] + wabs[6] * workv3[17];
/* ** Backward recursion ------------------------------------------------- */
/*  -- Quantities of body 3 */
    mrel22[2] = inerth[18] + 4.68;
    workv6[0] = ator[6] - inerth[21] * workv3[1] - inerth[24] * workv3[2] +
	    rhicm[8] * workv3[10] - rhicm[7] * workv3[11];
    workv6[1] = ator[7] - inerth[22] * workv3[1] - inerth[25] * workv3[2] -
	    rhicm[8] * workv3[9] + rhicm[6] * workv3[11];
    workv6[2] = ator[8] - inerth[25] * workv3[1] - inerth[26] * workv3[2] +
	    rhicm[7] * workv3[9] - rhicm[6] * workv3[10];
    workv6[3] = afor[6] - rhicm[8] * workv3[1] + rhicm[7] * workv3[2] - mhges
	    * workv3[9];
    workv6[4] = afor[7] - rhicm[6] * workv3[2] - mhges * workv3[10];
    rhrel2[2] = fgen[3] + workv6[0];
    mrel12[12] = inerth[18] + .078 + rhicm[8] * .5;
    workv6[18] = workv6[0] - workv6[4] * .5;
    workv6[19] = workv6[1] + workv6[3] * .5;
    workm6[108] = inerth[18] - rhicm[8] * -.5;
    workm6[112] = -rhicm[8] + mhges * -.5;
    workm6[115] = inerth[22] + rhicm[8] * .5;
    workm6[117] = rhicm[8] + mhges * .5;
    inertg[36] = workm6[108] + 4.9544125 - workm6[112] * .5;
    inertg[43] = workm6[115] + 6.7759085 + workm6[117] * .5;
    inertg[48] = inerth[24] - 2.45219 - rhicm[6] * .5;
    inertg[49] = inerth[25] - rhicm[7] * .5;
    inertg[50] = inerth[26] + 2.311496;
    inertg[60] = -11.5825 - rhicm[8] - mhges * .5;
    inertg[62] = rhicm[6] + 9.718;
    inertg[67] = -9.718 - rhicm[6];
    beta[6] = ator[3] + workv6[18];
    beta[7] = ator[4] + workv6[19];
    beta[8] = ator[5] + workv6[2];
    works2[0] = mrel12[12] / mrel22[2];
    works2[1] = inerth[21] / mrel22[2];
    works2[2] = inerth[24] / mrel22[2];
    works2[4] = -rhicm[8] / mrel22[2];
    works2[5] = rhicm[7] / mrel22[2];
    inertg[36] -= mrel12[12] * works2[0];
    inertg[42] = inerth[21] - mrel12[12] * works2[1];
    inertg[43] -= inerth[21] * works2[1];
    inertg[48] -= mrel12[12] * works2[2];
    inertg[49] -= inerth[21] * works2[2];
    inertg[50] -= inerth[24] * works2[2];
    inertg[60] -= mrel12[12] * works2[4];
    inertg[61] = -inerth[21] * works2[4];
    inertg[62] -= inerth[24] * works2[4];
    inertg[66] = rhicm[7] - mrel12[12] * works2[5];
    inertg[67] -= inerth[21] * works2[5];
    inertg[68] = -inerth[24] * works2[5];
    beta[6] -= works2[0] * rhrel2[2];
    beta[7] -= works2[1] * rhrel2[2];
    beta[8] -= works2[2] * rhrel2[2];
/*  -- Quantities of body 2 */
    mrel22[1] = inertg[36] + 57.33;
    workv6[0] = beta[6] - inertg[42] * zeta[7] - inertg[48] * zeta[8];
    workv6[1] = beta[7] - inertg[43] * zeta[7] - inertg[49] * zeta[8];
    workv6[2] = beta[8] - inertg[49] * zeta[7] - inertg[50] * zeta[8];
    rhrel2[1] = fgen[2] + workv6[0];
    works1[8] = sinp[1] * inertg[42] + cosp[1] * inertg[48];
    works1[11] = sinp[1] * inertg[60] + cosp[1] * inertg[66];
    workv6[14] = sinp[1] * workv6[1] + cosp[1] * workv6[2];
    workas[30] = cosp[1] + sinp[1];
    workas[31] = cosp[1] - sinp[1];
    workas[33] = workas[30] * workas[31];
    workas[32] = cosp[1] * sinp[1];
    workas[34] = workas[32] + workas[32];
    workas[0] = inertg[50] + inertg[43];
    workas[1] = inertg[50] - inertg[43];
    workas[4] = workas[0] / 2.;
    workas[5] = workas[1] / 2.;
    workas[8] = workas[33] * workas[5] + workas[34] * inertg[49];
    workam[115] = workas[4] + workas[8];
    workas[20] = inertg[68] + inertg[61];
    workas[21] = inertg[68] - inertg[61];
    workas[24] = workas[20] / 2.;
    workas[25] = workas[21] / 2.;
    workas[23] = -inertg[62] - inertg[67];
    workas[27] = workas[23] / 2.;
    workas[28] = workas[33] * workas[25] - workas[34] * workas[27];
    workam[133] = workas[24] + workas[28];
    inertg[14] = workam[115] + 1.16;
    works2[2] = works1[8] / mrel22[1];
    works2[5] = works1[11] / mrel22[1];
    inertg[14] -= works1[8] * works2[2];
    inertg[32] = workam[133] - works1[8] * works2[5];
    beta[2] = workv6[14] - works2[2] * rhrel2[1];
/*  -- Quantities of body 1 */
    mrel22[0] = inertg[14] + 14.3325;
    workv6[2] = beta[2] - inertg[32] * 9.81;
    rhrel2[0] = fgen[1] + workv6[2];
/* ** Forward recursion -------------------------------------------------- */
/*  -- Quantities of body 1 */
    genvd[0] = rhrel2[0] / mrel22[0];
/*  -- Quantities of body 2 */
    rhrel2[1] = rhrel2[1] - works1[8] * genvd[0] - works1[11] * 9.81;
    genvd[1] = rhrel2[1] / mrel22[1];
    workv6[7] = zeta[7] + sinp[1] * genvd[0];
    workv6[8] = zeta[8] + cosp[1] * genvd[0];
    workv6[10] = sinp[1] * 9.81;
    workv6[11] = cosp[1] * 9.81;
/*  -- Quantities of body 3 */
    rhrel2[2] = rhrel2[2] - mrel12[12] * genvd[1] - inerth[21] * workv6[7] -
	    inerth[24] * workv6[8] + rhicm[8] * workv6[10] - rhicm[7] *
	    workv6[11];
    genvd[2] = rhrel2[2] / mrel22[2];
    aqdd[1] = genvd[0];
    aqdd[2] = genvd[1];
    aqdd[3] = genvd[2];
    return;
}


////////////////////////////////////////////////////////////////////////////
///////////////////////      END OF FILE     ///////////////////////////////
////////////////////////////////////////////////////////////////////////////
