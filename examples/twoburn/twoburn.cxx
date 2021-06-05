//////////////////////////////////////////////////////////////////////////
////////////////        twoburn.cxx               /////////////////////
//////////////////////////////////////////////////////////////////////////
////////////////           PSOPT  Example             /////////////////////
//////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
//////// Title:   two burn orbit transfer problem         ////////////////
//////// Last modified: 07 February 2010                  ////////////////
//////// Reference:     Betts  (2001)             	  ////////////////
//////// (See PSOPT handbook for full reference)           ////////////////
//////////////////////////////////////////////////////////////////////////
////////     Copyright (c) Victor M. Becerra, 2010        ////////////////
//////////////////////////////////////////////////////////////////////////
//////// This is part of the PSOPT software library, which ////////////////
//////// is distributed under the terms of the GNU Lesser ////////////////
//////// General Public License (LGPL)                    ////////////////
//////////////////////////////////////////////////////////////////////////

#include "psopt.h"

using namespace PSOPT;



//////////////////////////////////////////////////////////////////////////
///////////////////  Define auxiliary functions                 //////////
//////////////////////////////////////////////////////////////////////////

adouble legendre_polynomial( adouble x, int n)
{
// This function computes the value of the legendre polynomials
// for a given value of the argument x and for n=0...5 only

  adouble retval=0.0;

  switch(n) {
     case 0:
         retval=1.0; break;
     case 1:
         retval= x; break;
     case 2:
         retval= 0.5*(3.0*pow(x,2)-1.0); break;
     case 3:
         retval= 0.5*(5.0*pow(x,3)- 3*x); break;
     case 4:
         retval= (1.0/8.0)*(35.0*pow(x,4) - 30.0*pow(x,2) + 3.0); break;
     case 5:
         retval= (1.0/8.0)*(63.0*pow(x,5) - 70.0*pow(x,3) + 15.0*x); break;
     default:
         error_message("legendre_polynomial(x,n) is limited to n=0...5");
  }

  return retval;

}

adouble legendre_polynomial_derivative( adouble x, int n)
{
// This function computes the value of the legendre polynomial derivatives
// for a given value of the argument x and for n=0...5 only.

  adouble retval=0.0;

  switch(n) {
     case 0:
         retval=0.0; break;
     case 1:
         retval= 1.0; break;
     case 2:
         retval= 0.5*(2.0*3.0*x); break;
     case 3:
         retval= 0.5*(3.0*5.0*pow(x,2)-3.0); break;
     case 4:
         retval= (1.0/8.0)*(4.0*35.0*pow(x,3) - 2.0*30.0*x ); break;
     case 5:
         retval= (1.0/8.0)*(5.0*63.0*pow(x,4) - 3.0*70.0*pow(x,2) + 15.0); break;
     default:
         error_message("legendre_polynomial_derivative(x,n) is limited to n=0...5");
  }

  return retval;

}


void compute_cartesian_trajectory(const MatrixXd& x, MatrixXd& xyz )
{

   int npoints = x.cols();
   xyz.resize(3,npoints);

   for(int i=0; i<npoints;i++) {


   double p = x(0,i);
   double f = x(1,i);
   double g = x(2,i);
   double h = x(3,i);
   double k = x(4,i);
   double L = x(5,i);

   double q      =  1.0 + f*cos(L) + g*sin(L);
   double r      =  p/q;
   double alpha2 = h*h - k*k;
   double X      = sqrt( h*h + k*k );
   double s2     = 1 + X*X;

   double r1 = r/s2*( cos(L) + alpha2*cos(L) + 2*h*k*sin(L));
   double r2 = r/s2*( sin(L) - alpha2*sin(L) + 2*h*k*cos(L));
   double r3 = 2*r/s2*( h*sin(L) - k*cos(L) );

   xyz(0,i) = r1;
   xyz(1,i) = r2;
   xyz(2,i) = r3;

   }

}


//////////////////////////////////////////////////////////////////////////
///////////////////  Define the end point (Mayer) cost function //////////
//////////////////////////////////////////////////////////////////////////

adouble endpoint_cost(adouble* initial_states, adouble* final_states,
                      adouble* parameters,adouble& t0, adouble& tf,
                      adouble* xad, int iphase, Workspace* workspace)
{

   if (iphase == 4) {
	   adouble w = final_states[6];
	   return (-w);
   }
   else {
	   return (0);
   }
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
   // Local integers
   int i, j;
   // Define constants:
   double Isp  = 300.0;            // [sec]
   double mu   = 1.407645794e16;   // [f2^2/sec^2]
   double g0   = 32.174;           // [ft/sec^2]
   double Re   = 20925662.73;      // [ft]
   double T    = 1.2;              // [lb]
   double  CM2W = 32.174;

   double J[5];
   J[2] =  1082.639e-6;
   J[3] = -2.565e-6;
   J[4] = -1.608e-6;

   // Extract individual variables

   adouble p = states[ 0 ];
   adouble f = states[ 1 ];
   adouble g = states[ 2 ];
   adouble h = states[ 3 ];
   adouble k = states[ 4 ];
   adouble L = states[ 5 ];
   adouble w;


   adouble* u  = controls;

   // Define some dependent variables

   adouble q      =  1.0 + f*cos(L) + g*sin(L);
   adouble r      =  p/q;
   adouble alpha2 = h*h - k*k;
   adouble X      = sqrt( h*h + k*k );
   adouble s2     = 1 + X*X;

   // r and v

   adouble r1 = r/s2*( cos(L) + alpha2*cos(L) + 2*h*k*sin(L));
   adouble r2 = r/s2*( sin(L) - alpha2*sin(L) + 2*h*k*cos(L));
   adouble r3 = 2*r/s2*( h*sin(L) - k*cos(L) );

   adouble rvec[3];

   rvec[ 0 ] = r1; rvec[ 1] = r2; rvec[ 2 ] = r3;

   adouble v1 = -(1.0/s2)*sqrt(mu/p)*(  sin(L) + alpha2*sin(L) - 2*h*k*cos(L) + g - 2*f*h*k + alpha2*g);
   adouble v2 = -(1.0/s2)*sqrt(mu/p)*( -cos(L) + alpha2*cos(L) + 2*h*k*sin(L) - f + 2*g*h*k + alpha2*f);
   adouble v3 =  (2.0/s2)*sqrt(mu/p)*(h*cos(L) +      k*sin(L) + f*h + g*k);

   adouble vvec[3];

   vvec[ 0 ] = v1; vvec[ 1] = v2; vvec[ 2 ] = v3;

   // compute Qr

   adouble ir[3], ith[3], ih[3];
   adouble rv[3];
   adouble rvr[3];

   cross( rvec, vvec, rv );

   cross( rv, rvec, rvr );

   adouble norm_r = sqrt( dot(rvec, rvec, 3) );

   adouble norm_rv = sqrt( dot(rv, rv, 3) );

   for (i=0; i<3; i++) {

      ir[i]  = rvec[i]/norm_r;

      ith[i] = rvr[i]/( norm_rv*norm_r );

      ih[i]  = rv[i]/norm_rv;

   }

   adouble Qr1[3], Qr2[3], Qr3[3];

   for(i=0; i< 3; i++)
   {
        // Columns of matrix Qr
        Qr1[i] = ir[i];
        Qr2[i] = ith[i];
        Qr3[i] = ih[i];
   }

   adouble Qv1[3], Qv2[3], Qv3[3];

   adouble norm_v = sqrt( dot(vvec, vvec,3) );

   for(i=0;i<3;i++)
   {
     Qv1[i] = vvec[i]/norm_v;
   }

   adouble vr[3];
   cross( vvec, rvec, vr);

   adouble norm_vr = sqrt( dot(vr, vr,3) );

   for(i=0;i<3;i++)
   {
     Qv2[i] = vr[i]/norm_vr;
   }

   cross( Qv1, Qv2, Qv3 );

   // Compute in

   adouble en[3];
   en[ 0 ] = 0.0; en[ 1 ]= 0.0; en[ 2 ] = 1.0;

   adouble enir = dot(en,ir,3);

   adouble in[3];

   for(i=0;i<3;i++) {
      in[i] = en[i] - enir*ir[i];
   }

   adouble norm_in = sqrt( dot( in, in, 3 ) );

   for(i=0;i<3;i++) {
      in[i] = in[i]/norm_in;
   }

   // Geocentric latitude angle:

   adouble sin_phi =  rvec[ 2 ]/ sqrt( dot(rvec,rvec,3) ) ;
   adouble cos_phi =  sqrt(1.0- pow(sin_phi,2.0));

   adouble deltagn = 0.0;
   adouble deltagr = 0.0;
   for (j=2; j<=4;j++) {
     adouble Pdash_j = legendre_polynomial_derivative( sin_phi, j );
     adouble P_j     = legendre_polynomial( sin_phi, j );
     deltagn += -mu*cos_phi/(r*r)*pow(Re/r,j)*Pdash_j*J[j];
     deltagr += -mu/(r*r)* (j+1)*pow( Re/r,j)*P_j*J[j];
   }

   // Compute vector delta_g

   adouble delta_g[3];
   for (i=0;i<3;i++) {
       delta_g[i] = deltagn*in[i] - deltagr*ir[i];
   }

   // Compute vector DELTA_g

   adouble DELTA_g[3];

   DELTA_g[ 0 ] = dot(Qr1, delta_g,3);
   DELTA_g[ 1 ] = dot(Qr2, delta_g,3);
   DELTA_g[ 2 ] = dot(Qr3, delta_g,3);

   // Compute DELTA_T

   adouble DELTA_T[3];

   if (iphase == 1 || iphase==3) {
        for(i=0;i<3;i++) {
              DELTA_T[i] =  0.0;
        }
   }

   else {
       adouble Qr[9], Qrt[9], Qv[9], Tvel[3];

       for(i=0;i<3;i++) {
	  Qr[i]   = Qr1[i];
	  Qr[3+i] = Qr2[i];
	  Qr[6+i] = Qr3[i];

	  Qv[i]   = Qv1[i];
	  Qv[3+i] = Qv2[i];
	  Qv[6+i] = Qv3[i];
       }

       transpose_ad(Qr, 3, 3, Qrt );

       adouble theta = u[ 0 ];
       adouble phi   = u[ 1 ];

       adouble Tvec[3];

       w = states[6];

       adouble mass = w/CM2W;

       adouble Tacc = T/mass;

       Tvec[ 0 ] = Tacc*cos(theta)*cos(phi);
       Tvec[ 1 ] = Tacc*cos(theta)*sin(phi);
       Tvec[ 2 ] = Tacc*sin(theta);

       product_ad( Qv, Tvec, 3, 3, 3, 1, Tvel );

       product_ad(Qrt, Tvel, 3, 3, 3, 1, DELTA_T );


   }

   // Compute DELTA

   adouble DELTA[3];

   for(i=0;i<3;i++) {
      DELTA[i] =  DELTA_g[i] + DELTA_T[i];
   }

   adouble delta1= DELTA[ 0 ];
   adouble delta2= DELTA[ 1 ];
   adouble delta3= DELTA[ 2 ];


  // derivatives

   adouble pdot =                                      2*p/q*sqrt(p/mu)               * delta2;
   adouble fdot =  sqrt(p/mu)*sin(L) * delta1 + sqrt(p/mu)*(1.0/q)*((q+1.0)*cos(L)+f) * delta2
                   -  sqrt(p/mu)*(g/q)*(h*sin(L)-k*cos(L))  * delta3;
   adouble gdot = -sqrt(p/mu)*cos(L) * delta1 + sqrt(p/mu)*(1.0/q)*((q+1.0)*sin(L)+g) * delta2
                   +  sqrt(p/mu)*(f/q)*(h*sin(L)-k*cos(L))  * delta3;
   adouble hdot =  sqrt(p/mu)*s2*cos(L)/(2.0*q)             * delta3;
   adouble kdot =  sqrt(p/mu)*s2*sin(L)/(2.0*q)             * delta3;
   adouble Ldot =  sqrt(p/mu)*(1.0/q)*(h*sin(L)-k*cos(L))* delta3     + sqrt(mu*p)*pow( (q/p),2.);

   adouble wdot;

   if (iphase==2 || iphase==4) {
     wdot = -T/Isp;
   }
   else {
     wdot = 0.0;
   }

   derivatives[ 0 ] = pdot;
   derivatives[ 1 ] = fdot;
   derivatives[ 2 ] = gdot;
   derivatives[ 3 ] = hdot;
   derivatives[ 4 ] = kdot;
   derivatives[ 5 ] = Ldot;
   if (iphase==2 || iphase==4) {
    derivatives[ 6 ] = wdot;
   }


}

////////////////////////////////////////////////////////////////////////////
///////////////////  Define the events function ////////////////////////////
////////////////////////////////////////////////////////////////////////////

void events(adouble* e, adouble* initial_states, adouble* final_states,
            adouble* parameters,adouble& t0, adouble& tf, adouble* xad,
            int iphase, Workspace* workspace)

{


   adouble pti = initial_states[ 0 ];
   adouble fti = initial_states[ 1 ];
   adouble gti = initial_states[ 2 ];
   adouble hti = initial_states[ 3 ];
   adouble kti = initial_states[ 4 ];
   adouble Lti = initial_states[ 5 ];
   adouble wti;
   if (iphase==2) {
      wti = initial_states[ 6 ];
   }

   adouble ptf = final_states[ 0 ];
   adouble ftf = final_states[ 1 ];
   adouble gtf = final_states[ 2 ];
   adouble htf = final_states[ 3 ];
   adouble ktf = final_states[ 4 ];
   adouble Ltf = final_states[ 5 ];


   if (iphase==1) {
   	e[ 0 ]  = pti;
   	e[ 1 ]  = fti;
   	e[ 2 ]  = gti;
   	e[ 3 ]  = hti;
   	e[ 4 ]  = kti;
   	e[ 5 ]  = Lti;
   }

   if (iphase==2) {
     e[ 0 ] = wti;
   }

   if (iphase == 4) {
   	e[ 0 ] = ptf;
   	e[ 1 ] = ftf;
   	e[ 2 ] = gtf;
   	e[ 3 ] = htf;
   	e[ 4 ] = ktf;
   }

}



///////////////////////////////////////////////////////////////////////////
///////////////////  Define the phase linkages function ///////////////////
///////////////////////////////////////////////////////////////////////////

void linkages( adouble* linkages, adouble* xad, Workspace* workspace)
{
    // Numeber of linkages:
    // Boundary of phases 1,2: 6 state continuity + 1 time continuity
    // Boundary of phases 2,3: 6 state continuity + 1 time continuity
    // Boundary of phases 3,4: 6 state continuity + 1 time continuity
    // 1 extra linkage w(tf2) = w(ti4)
    // Total: 22 linkage constraints

    adouble xf[7], xi[7], t0a, t0b, tfa, tfb, wtf2, wti4;

    // Linking phases 1 and 2

    get_final_states(       xf, xad, 1, workspace );
    get_initial_states(     xi, xad, 2, workspace );
    tfa = get_final_time(       xad, 1, workspace );
    t0b = get_initial_time(     xad, 2, workspace );


    linkages[ 0 ] = xf[ 0 ] - xi[ 0 ];
    linkages[ 1 ] = xf[ 1 ] - xi[ 1 ];
    linkages[ 2 ] = xf[ 2 ] - xi[ 2 ];
    linkages[ 3 ] = xf[ 3 ] - xi[ 3 ];
    linkages[ 4 ] = xf[ 4 ] - xi[ 4 ];
    linkages[ 5 ] = xf[ 5 ] - xi[ 5 ];
    linkages[ 6 ] = t0b - tfa;

    // Linking phases 2 and 3

    get_final_states(        xf, xad, 2, workspace );
    get_initial_states(      xi, xad, 3, workspace );
    tfa = get_final_time(        xad, 2, workspace );
    t0b = get_initial_time(      xad, 3, workspace );

    wtf2 = xf[ 6 ];

    linkages[ 7 ]  = xf[ 0 ] - xi[ 0 ];
    linkages[ 8 ]  = xf[ 1 ] - xi[ 1 ];
    linkages[ 9 ] = xf[ 2 ] - xi[ 2 ];
    linkages[ 10 ] = xf[ 3 ] - xi[ 3 ];
    linkages[ 11 ] = xf[ 4 ] - xi[ 4 ];
    linkages[ 12 ] = xf[ 5 ] - xi[ 5 ];
    linkages[ 13 ] = t0b - tfa;

    // Linking phases 3 and 4

    get_final_states(      xf, xad, 3, workspace );
    get_initial_states(    xi, xad, 4, workspace );
    tfa = get_final_time(      xad, 3, workspace );
    t0b = get_initial_time(    xad, 4, workspace );

    wti4 = xi[ 6 ];


    linkages[ 14 ] = xf[ 0 ] - xi[ 0 ];
    linkages[ 15 ] = xf[ 1 ] - xi[ 1 ];
    linkages[ 16 ] = xf[ 2 ] - xi[ 2 ];
    linkages[ 17 ] = xf[ 3 ] - xi[ 3 ];
    linkages[ 18 ] = xf[ 4 ] - xi[ 4 ];
    linkages[ 19 ] = xf[ 5 ] - xi[ 5 ];
    linkages[ 20 ] = t0b - tfa;

    // Linking the weight at the end of phase 2 with the weight at the beginning of phase 4

    linkages[ 21 ] = wtf2 - wti4;

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

    MSdata msdata;

////////////////////////////////////////////////////////////////////////////
///////////////////  Register problem name  ////////////////////////////////
////////////////////////////////////////////////////////////////////////////

    problem.name        		= "Two burn transfer problem";
    problem.outfilename                 = "twoburn.txt";

////////////////////////////////////////////////////////////////////////////
////////////  Define problem level constants & do level 1 setup ////////////
////////////////////////////////////////////////////////////////////////////

    problem.nphases   			          = 4;
    problem.nlinkages                   = 22;

    psopt_level1_setup(problem);

/////////////////////////////////////////////////////////////////////////////
/////////   Define phase related information & do level 2 setup  ////////////
/////////////////////////////////////////////////////////////////////////////

    problem.phases(1).nstates   						= 6;
    problem.phases(1).ncontrols 						= 0;
    problem.phases(1).nparameters               = 0;
    problem.phases(1).nevents   	        			= 6;
    problem.phases(1).npath     						= 0;
    problem.phases(1).nodes                     << 10;

    problem.phases(2).nstates   						= 7;
    problem.phases(2).ncontrols 						= 2;
    problem.phases(2).nparameters               = 0;
    problem.phases(2).nevents   	        			= 1;
    problem.phases(2).npath     						= 0;
    problem.phases(2).nodes                     << 10;

    problem.phases(3).nstates   						= 6;
    problem.phases(3).ncontrols 						= 0;
    problem.phases(3).nparameters               = 0;
    problem.phases(3).nevents   	        			= 0;
    problem.phases(3).npath     						= 0;
    problem.phases(3).nodes                     << 10;

    problem.phases(4).nstates   						= 7;
    problem.phases(4).ncontrols 						= 2;
    problem.phases(4).nparameters               = 0;
    problem.phases(4).nevents   	        			= 5;
    problem.phases(4).npath     						= 0;
    problem.phases(4).nodes                     << 10;

    psopt_level2_setup(problem, algorithm);


////////////////////////////////////////////////////////////////////////////
///////////////////  Enter problem bounds information //////////////////////
////////////////////////////////////////////////////////////////////////////




    double pti = 21837080.052835;
    double fti = 0.0;
    double gti = 0.0;
    double hti = -0.25396764647494;
    double kti = 0.0;
    double Lti = pi;
    double wti = 1.0;

    double  SISP = 300.0;
    double  DELTAV2 = 8000;
    double  DELTAV4 = 3000;
    double  CM2W = 32.174;

    double sigma = 1.0/6076.1154855643;

    double Re = 20925662.73;

    double wtf2  = wti*exp(-DELTAV2/(CM2W*SISP));
    double wtf4  = wtf2*exp(-DELTAV4/(CM2W*SISP));

    double ptf = 19323/sigma + Re;
    double ftf = 0.0;
    double gtf = 0.0;
    double htf = 0.0;
    double ktf = 0.0;

    double D2R = pi/180.0;


   // BOUNDS FOR PHASE 1


    problem.phases(1).bounds.lower.states(0) = 10.e6;
    problem.phases(1).bounds.lower.states(1) = -1;
    problem.phases(1).bounds.lower.states(2) = -1;
    problem.phases(1).bounds.lower.states(3) = -1;
    problem.phases(1).bounds.lower.states(4) = -1;
    problem.phases(1).bounds.lower.states(5) = pi;


    problem.phases(1).bounds.upper.states(0) = 2e8;
    problem.phases(1).bounds.upper.states(1) = 1;
    problem.phases(1).bounds.upper.states(2) = 1;
    problem.phases(1).bounds.upper.states(3) = 1;
    problem.phases(1).bounds.upper.states(4) = 1;
    problem.phases(1).bounds.upper.states(5) = 30*pi;



    problem.phases(1).bounds.lower.events(0)  = pti;
    problem.phases(1).bounds.lower.events(1)  = fti;
    problem.phases(1).bounds.lower.events(2)  = gti;
    problem.phases(1).bounds.lower.events(3)  = hti;
    problem.phases(1).bounds.lower.events(4)  = kti;
    problem.phases(1).bounds.lower.events(5)  = Lti;



    problem.phases(1).bounds.upper.events(0)  = pti;
    problem.phases(1).bounds.upper.events(1)  = fti;
    problem.phases(1).bounds.upper.events(2)  = gti;
    problem.phases(1).bounds.upper.events(3)  = hti;
    problem.phases(1).bounds.upper.events(4)  = kti;
    problem.phases(1).bounds.upper.events(5)  = Lti;


    problem.phases(1).bounds.lower.StartTime    = 0.0;
    problem.phases(1).bounds.upper.StartTime    = 0.0;

    problem.phases(1).bounds.lower.EndTime      = 2000.0;
    problem.phases(1).bounds.upper.EndTime      = 3000.0;

    // BOUNDS FOR PHASE 2

    problem.phases(2).bounds.lower.states(0) = 10.e6;
    problem.phases(2).bounds.lower.states(1) = -1;
    problem.phases(2).bounds.lower.states(2) = -1;
    problem.phases(2).bounds.lower.states(3) = -1;
    problem.phases(2).bounds.lower.states(4) = -1;
    problem.phases(2).bounds.lower.states(5) = pi;
    problem.phases(2).bounds.lower.states(6) = 0.0;

    problem.phases(2).bounds.upper.states(0) = 2.e8;
    problem.phases(2).bounds.upper.states(1) = 1;
    problem.phases(2).bounds.upper.states(2) = 1;
    problem.phases(2).bounds.upper.states(3) = 1;
    problem.phases(2).bounds.upper.states(4) = 1;
    problem.phases(2).bounds.upper.states(5) = 30*pi;
    problem.phases(2).bounds.upper.states(6) = 2.0;


    problem.phases(2).bounds.lower.controls(0) = -pi;
    problem.phases(2).bounds.lower.controls(1) = -pi;
    problem.phases(2).bounds.upper.controls(0) = pi;
    problem.phases(2).bounds.upper.controls(1) = pi;

    problem.phases(2).bounds.lower.events(0)  = wti;
    problem.phases(2).bounds.upper.events(0)  = wti;


    problem.phases(2).bounds.lower.StartTime    = 2000;
    problem.phases(2).bounds.upper.StartTime    = 3000;

    problem.phases(2).bounds.lower.EndTime      = 2100;
    problem.phases(2).bounds.upper.EndTime      = 3100;


    // BOUNDS FOR PHASE 
    problem.phases(3).bounds.lower.states(0) = 10.e6;
    problem.phases(3).bounds.lower.states(1) = -1;
    problem.phases(3).bounds.lower.states(2) = -1;
    problem.phases(3).bounds.lower.states(3) = -1;
    problem.phases(3).bounds.lower.states(4) = -1;
    problem.phases(3).bounds.lower.states(5) = pi;


    problem.phases(3).bounds.upper.states(0) = 2.e8;
    problem.phases(3).bounds.upper.states(1) = 1.0;
    problem.phases(3).bounds.upper.states(2) = 1.0;
    problem.phases(3).bounds.upper.states(3) = 1.0;
    problem.phases(3).bounds.upper.states(4) = 1.0;
    problem.phases(3).bounds.upper.states(5) = 30*pi;


    problem.phases(3).bounds.lower.StartTime    = 2100;
    problem.phases(3).bounds.upper.StartTime    = 3100;

    problem.phases(3).bounds.lower.EndTime      = 21600;
    problem.phases(3).bounds.upper.EndTime      = 21800;


    // BOUNDS FOR PHASE 4

    problem.phases(4).bounds.lower.states(0) = 10.e6;
    problem.phases(4).bounds.lower.states(1) = -1;
    problem.phases(4).bounds.lower.states(2) = -1;
    problem.phases(4).bounds.lower.states(3) = -1;
    problem.phases(4).bounds.lower.states(4) = -1;
    problem.phases(4).bounds.lower.states(5) = pi;
    problem.phases(4).bounds.lower.states(6) = 0.0;

    problem.phases(4).bounds.upper.states(0) = 2.e8;
    problem.phases(4).bounds.upper.states(1) = 1;
    problem.phases(4).bounds.upper.states(2) = 1;
    problem.phases(4).bounds.upper.states(3) = 1;
    problem.phases(4).bounds.upper.states(4) = 1;
    problem.phases(4).bounds.upper.states(5) = 30*pi;
    problem.phases(4).bounds.upper.states(6) = 2.0;


    problem.phases(4).bounds.lower.controls(0) = -pi;
    problem.phases(4).bounds.lower.controls(1) = -pi;
    problem.phases(4).bounds.upper.controls(0) = pi;
    problem.phases(4).bounds.upper.controls(1) = pi;


    problem.phases(4).bounds.lower.events(0)  = ptf;
    problem.phases(4).bounds.lower.events(1)  = ftf;
    problem.phases(4).bounds.lower.events(2)  = gtf;
    problem.phases(4).bounds.lower.events(3)  = htf;
    problem.phases(4).bounds.lower.events(4)  = ktf;

    problem.phases(4).bounds.upper.events(0)  = ptf;
    problem.phases(4).bounds.upper.events(1)  = ftf;
    problem.phases(4).bounds.upper.events(2)  = gtf;
    problem.phases(4).bounds.upper.events(3)  = htf;
    problem.phases(4).bounds.upper.events(4)  = ktf;

    problem.phases(4).bounds.lower.StartTime    = 21600;
    problem.phases(4).bounds.upper.StartTime    = 21800;

    problem.phases(4).bounds.lower.EndTime      = 21650;
    problem.phases(4).bounds.upper.EndTime      = 21900;


////////////////////////////////////////////////////////////////////////////
///////////////////  Register problem functions  ///////////////////////////
////////////////////////////////////////////////////////////////////////////


    problem.integrand_cost 			= &integrand_cost;
    problem.endpoint_cost 				= &endpoint_cost;
    problem.dae             			= &dae;
    problem.events 						= &events;
    problem.linkages						= &linkages;

////////////////////////////////////////////////////////////////////////////
///////////////////  Define & register initial guess ///////////////////////
////////////////////////////////////////////////////////////////////////////

    int nnodes;
    int ncontrols;
    int nstates;
    int iphase;
    MatrixXd x_guess, u_guess, time_guess, param_guess, xini, xfinal;

    // Phase 1

    nnodes    							= problem.phases(1).nodes(0);
    nstates                     	= problem.phases(1).nstates;
    iphase = 1;

    x_guess    =  zeros(nstates,nnodes);
    time_guess =  linspace(0.0,2690,nnodes);

    xini.resize(6,1);
    xini(0)= pti; xini(1)=fti;xini(2)=gti;xini(3)=hti;xini(4)=kti;xini(5)=Lti;



    rk4_propagate( dae, u_guess, time_guess, xini, param_guess, problem, iphase, x_guess, NULL);

//    tra(x_guess).Print("x_guess(iphase=1)");

    xfinal = x_guess.col(nnodes-1); 

    problem.phases(1).guess.states = x_guess;
    problem.phases(1).guess.time   = time_guess;




    // Phase 2

    nnodes    			           = problem.phases(2).nodes(0);
    nstates                     = problem.phases(2).nstates;
    ncontrols                   = problem.phases(2).ncontrols;
    iphase = 2;

    u_guess    =  zeros(ncontrols,nnodes);
    x_guess    =  zeros(nstates,nnodes);
    time_guess =  linspace(2690,2840,nnodes);

    xini.resize(7,1);
    xini(0)= xfinal(0); xini(1)=xfinal(1);xini(2)=xfinal(2);xini(3)=xfinal(3);xini(4)=xfinal(4);xini(5)=xfinal(5);
    xini(6)= wti;

    u_guess.row(0) =  0.148637e-2*D2R*ones(1,nnodes);
    u_guess.row(1) = -9.08446*D2R*ones(1,nnodes);

    rk4_propagate( dae, u_guess, time_guess, xini, param_guess, problem, iphase, x_guess, NULL);

//    tra(x_guess).Print("x_guess(iphase=2)");


    xfinal = x_guess.col(nnodes-1); 
    double wtf2__ = xfinal(6);

    problem.phases(2).guess.states   = x_guess;
    problem.phases(2).guess.controls = u_guess;
    problem.phases(2).guess.time     = time_guess;

    // Phase 3

    nnodes    			           = problem.phases(3).nodes(0);
    nstates                     = problem.phases(3).nstates;
    iphase = 3;


    x_guess    =  zeros(nstates,nnodes);
    time_guess =  linspace(2840,21650,nnodes);


    xini.resize(6,1);
    xini(0)= xfinal(0); xini(1)=xfinal(1);xini(2)=xfinal(2);xini(3)=xfinal(3);xini(4)=xfinal(4);xini(5)=xfinal(5);

    rk4_propagate( dae, u_guess, time_guess, xini, param_guess, problem, iphase, x_guess, NULL);

//    tra(x_guess).Print("x_guess(iphase=3)");

    xfinal = x_guess.col(nnodes-1); 

    problem.phases(3).guess.states = x_guess;
    problem.phases(3).guess.time   = time_guess;

    // Phase 4

    nnodes    			           = problem.phases(4).nodes(0);
    nstates                     = problem.phases(4).nstates;
    ncontrols                   = problem.phases(4).ncontrols;
    iphase = 4;

    u_guess    =  zeros(ncontrols,nnodes);
    x_guess    =  zeros(nstates,nnodes);
    time_guess =  linspace(21650,21700,nnodes);

    u_guess.row(0) =  -0.136658e-2*D2R*ones(1,nnodes);
    u_guess.row(1) =       49.7892*D2R*ones(1,nnodes);

    xini.resize(7,1);
    xini(0)= xfinal(0); xini(1)=xfinal(1);xini(2)=xfinal(2);xini(3)=xfinal(3);xini(4)=xfinal(4);xini(5)=xfinal(5);
    xini(6)= wtf2__;

    rk4_propagate( dae, u_guess, time_guess, xini, param_guess, problem, iphase, x_guess, NULL);

//    tra(x_guess).Print("x_guess(iphase=4)");


    problem.phases(4).guess.states   = x_guess;
    problem.phases(4).guess.controls = u_guess;
    problem.phases(4).guess.time     = time_guess;



////////////////////////////////////////////////////////////////////////////
///////////////////  Enter algorithm options  //////////////////////////////
////////////////////////////////////////////////////////////////////////////


    algorithm.nlp_iter_max                = 1000;
    algorithm.nlp_tolerance               = 1.e-6;
    algorithm.nlp_method                  = "IPOPT";
    algorithm.scaling                     = "automatic";
    algorithm.derivatives                 = "automatic";
    algorithm.defect_scaling              = "jacobian-based";
    algorithm.jac_sparsity_ratio          =  0.11;
    algorithm.collocation_method          = "trapezoidal";
//    algorithm.diff_matrix                 = "central-differences";
    algorithm.mesh_refinement             = "automatic";
    algorithm.mr_max_iterations           = 5;
    algorithm.ode_tolerance               = 1.0e-6;



////////////////////////////////////////////////////////////////////////////
///////////////////  Now call PSOPT to solve the problem   //////////////////
////////////////////////////////////////////////////////////////////////////

    psopt(solution, problem, algorithm);

////////////////////////////////////////////////////////////////////////////
///////////  Extract relevant variables from solution structure   //////////
////////////////////////////////////////////////////////////////////////////



    MatrixXd x, t, xp1, xp2, xp3, xp4, tp1, tp2, tp3, tp4;

    xp1      = solution.get_states_in_phase(1);
    tp1      = solution.get_time_in_phase(1);

    xp2     = solution.get_states_in_phase(2);
    tp2     = solution.get_time_in_phase(2);

    xp3     = solution.get_states_in_phase(3);
    tp3     = solution.get_time_in_phase(3);
    
    xp3     = solution.get_states_in_phase(3);
    tp3     = solution.get_time_in_phase(3);
    
    xp4     = solution.get_states_in_phase(4);
    tp4     = solution.get_time_in_phase(4);

    x.resize(6, length(tp1)+length(tp2)+length(tp3)+length(tp4));
    t.resize(1, length(tp1)+length(tp2)+length(tp3)+length(tp4));

    x  <<  xp1, xp2.block(0,0,6,length(tp2)), xp3, xp4.block(0,0,6,length(tp4));
    t  <<  tp1, tp2, tp3, tp4;


    MatrixXd u_phase2     = solution.get_controls_in_phase(2);
    MatrixXd u_phase4     = solution.get_controls_in_phase(4);





////////////////////////////////////////////////////////////////////////////
///////////  Save solution data to files if desired ////////////////////////
////////////////////////////////////////////////////////////////////////////

//    x.Save("x.dat");
//    u.Save("u.dat");
//    t.Save("t.dat");

////////////////////////////////////////////////////////////////////////////
///////////  Plot some results if desired (requires gnuplot) ///////////////
////////////////////////////////////////////////////////////////////////////

    double R2D = 180/pi;

    MatrixXd x1 = x.row(0)/1.e6;
    MatrixXd x2 = x.row(1); 
    MatrixXd x3 = x.row(2); 
    MatrixXd x4 = x.row(3); 
    MatrixXd x5 = x.row(4); 
    MatrixXd x6 = x.row(5); 
//    MatrixXd x7 = x(7,colon());
    MatrixXd theta_phase2;
    MatrixXd theta_phase4;
    MatrixXd phi_phase2;
    MatrixXd phi_phase4;
    MatrixXd r;

    theta_phase2 =  u_phase2.row(0)*R2D; 
    phi_phase2    = u_phase2.row(1)*R2D; 
    theta_phase4  = u_phase4.row(0)*R2D; 
    phi_phase4    = u_phase4.row(1)*R2D; 

    compute_cartesian_trajectory(x,r);

    Save(r,"r.dat");

    double ft2km = 0.0003048;

    r = r*ft2km;


    plot(tp2,theta_phase2,problem.name+": thrust theta phase 2","time (s)", "theta (deg)", "theta");

    plot(tp2,phi_phase2,problem.name+": thrust phi phase 2","time (s)", "phi (deg)", "phi");

    plot(tp4,theta_phase4,problem.name+": thrust theta phase 4","time (s)", "theta (deg)", "theta");

    plot(tp4,phi_phase4,problem.name+": thrust phi phase 4","time (s)", "phi (deg)", "phi");

    plot(tp2,theta_phase2,problem.name+":  thrust pitch angle phase 2","time (s)", "theta (deg)", "theta",
	 "pdf", "theta2.pdf");

    plot(tp2,phi_phase2,problem.name+": thrust angle phase 2","time (s)", "phi (deg)", "phi",
	 "pdf", "phi2.pdf");

    plot(tp4,theta_phase4,problem.name+": thrust pitch angle phase 4","time (s)", "theta (deg)", "theta",
	 "pdf", "theta4.pdf");

    plot(tp4,phi_phase4,problem.name+": thrust yaw angle phase 4","time (s)", "phi (deg)", "phi",
    	 "pdf", "phi4.pdf");

    plot3(r.row(0) , r.row(1) , r.row(2) , "Two burn trasnfer trajectory", "x (km)", "y (km)", "z (km)",
	   NULL, NULL, "30,110");

    plot3(r.row(0) , r.row(1) , r.row(2) , "Two burn transfer trajectory", "x (km)", "y (km)", "z (km)",
	   "pdf", "trajectory.pdf", "30,110");

    plot(r.row(0) , r.row(1),  "Two burn trajectory - projection on the equatorial plane",
	    "x (km)", "y (km)");



}

////////////////////////////////////////////////////////////////////////////
///////////////////////      END OF FILE     ///////////////////////////////
////////////////////////////////////////////////////////////////////////////
