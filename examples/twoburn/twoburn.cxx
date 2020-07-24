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


void compute_cartesian_trajectory(const DMatrix& x, DMatrix& xyz )
{

   int npoints = x.GetNoCols();
   xyz.Resize(3,npoints);

   for(int i=1; i<=npoints;i++) {


   double p = x(1,i);
   double f = x(2,i);
   double g = x(3,i);
   double h = x(4,i);
   double k = x(5,i);
   double L = x(6,i);

   double q      =  1.0 + f*cos(L) + g*sin(L);
   double r      =  p/q;
   double alpha2 = h*h - k*k;
   double X      = sqrt( h*h + k*k );
   double s2     = 1 + X*X;

   double r1 = r/s2*( cos(L) + alpha2*cos(L) + 2*h*k*sin(L));
   double r2 = r/s2*( sin(L) - alpha2*sin(L) + 2*h*k*cos(L));
   double r3 = 2*r/s2*( h*sin(L) - k*cos(L) );

   xyz(1,i) = r1;
   xyz(2,i) = r2;
   xyz(3,i) = r3;

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
	   adouble w = final_states[CINDEX(7)];
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

   adouble p = states[ CINDEX(1) ];
   adouble f = states[ CINDEX(2) ];
   adouble g = states[ CINDEX(3) ];
   adouble h = states[ CINDEX(4) ];
   adouble k = states[ CINDEX(5) ];
   adouble L = states[ CINDEX(6) ];
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

   rvec[ CINDEX(1) ] = r1; rvec[ CINDEX(2)] = r2; rvec[ CINDEX(3) ] = r3;

   adouble v1 = -(1.0/s2)*sqrt(mu/p)*(  sin(L) + alpha2*sin(L) - 2*h*k*cos(L) + g - 2*f*h*k + alpha2*g);
   adouble v2 = -(1.0/s2)*sqrt(mu/p)*( -cos(L) + alpha2*cos(L) + 2*h*k*sin(L) - f + 2*g*h*k + alpha2*f);
   adouble v3 =  (2.0/s2)*sqrt(mu/p)*(h*cos(L) +      k*sin(L) + f*h + g*k);

   adouble vvec[3];

   vvec[ CINDEX(1) ] = v1; vvec[ CINDEX(2)] = v2; vvec[ CINDEX(3) ] = v3;

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
   en[ CINDEX(1) ] = 0.0; en[ CINDEX(2) ]= 0.0; en[ CINDEX(3) ] = 1.0;

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

   adouble sin_phi =  rvec[ CINDEX(3) ]/ sqrt( dot(rvec,rvec,3) ) ;
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

   DELTA_g[ CINDEX(1) ] = dot(Qr1, delta_g,3);
   DELTA_g[ CINDEX(2) ] = dot(Qr2, delta_g,3);
   DELTA_g[ CINDEX(3) ] = dot(Qr3, delta_g,3);

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

       adouble theta = u[ CINDEX(1) ];
       adouble phi   = u[ CINDEX(2) ];

       adouble Tvec[3];

       w = states[CINDEX(7)];

       adouble mass = w/CM2W;

       adouble Tacc = T/mass;

       Tvec[ CINDEX(1) ] = Tacc*cos(theta)*cos(phi);
       Tvec[ CINDEX(2) ] = Tacc*cos(theta)*sin(phi);
       Tvec[ CINDEX(3) ] = Tacc*sin(theta);

       product_ad( Qv, Tvec, 3, 3, 3, 1, Tvel );

       product_ad(Qrt, Tvel, 3, 3, 3, 1, DELTA_T );


   }

   // Compute DELTA

   adouble DELTA[3];

   for(i=0;i<3;i++) {
      DELTA[i] =  DELTA_g[i] + DELTA_T[i];
   }

   adouble delta1= DELTA[ CINDEX(1) ];
   adouble delta2= DELTA[ CINDEX(2) ];
   adouble delta3= DELTA[ CINDEX(3) ];


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

   derivatives[ CINDEX(1) ] = pdot;
   derivatives[ CINDEX(2) ] = fdot;
   derivatives[ CINDEX(3) ] = gdot;
   derivatives[ CINDEX(4) ] = hdot;
   derivatives[ CINDEX(5) ] = kdot;
   derivatives[ CINDEX(6) ] = Ldot;
   if (iphase==2 || iphase==4) {
    derivatives[ CINDEX(7) ] = wdot;
   }


}

////////////////////////////////////////////////////////////////////////////
///////////////////  Define the events function ////////////////////////////
////////////////////////////////////////////////////////////////////////////

void events(adouble* e, adouble* initial_states, adouble* final_states,
            adouble* parameters,adouble& t0, adouble& tf, adouble* xad,
            int iphase, Workspace* workspace)

{


   adouble pti = initial_states[ CINDEX(1) ];
   adouble fti = initial_states[ CINDEX(2) ];
   adouble gti = initial_states[ CINDEX(3) ];
   adouble hti = initial_states[ CINDEX(4) ];
   adouble kti = initial_states[ CINDEX(5) ];
   adouble Lti = initial_states[ CINDEX(6) ];
   adouble wti;
   if (iphase==2) {
      wti = initial_states[ CINDEX(7) ];
   }

   adouble ptf = final_states[ CINDEX(1) ];
   adouble ftf = final_states[ CINDEX(2) ];
   adouble gtf = final_states[ CINDEX(3) ];
   adouble htf = final_states[ CINDEX(4) ];
   adouble ktf = final_states[ CINDEX(5) ];
   adouble Ltf = final_states[ CINDEX(6) ];


   if (iphase==1) {
   	e[ CINDEX(1) ]  = pti;
   	e[ CINDEX(2) ]  = fti;
   	e[ CINDEX(3) ]  = gti;
   	e[ CINDEX(4) ]  = hti;
   	e[ CINDEX(5) ]  = kti;
   	e[ CINDEX(6) ]  = Lti;
   }

   if (iphase==2) {
     e[ CINDEX(1) ] = wti;
   }

   if (iphase == 4) {
   	e[ CINDEX(1) ] = ptf;
   	e[ CINDEX(2) ] = ftf;
   	e[ CINDEX(3) ] = gtf;
   	e[ CINDEX(4) ] = htf;
   	e[ CINDEX(5) ] = ktf;
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


    linkages[ CINDEX(1) ] = xf[ CINDEX(1) ] - xi[ CINDEX(1) ];
    linkages[ CINDEX(2) ] = xf[ CINDEX(2) ] - xi[ CINDEX(2) ];
    linkages[ CINDEX(3) ] = xf[ CINDEX(3) ] - xi[ CINDEX(3) ];
    linkages[ CINDEX(4) ] = xf[ CINDEX(4) ] - xi[ CINDEX(4) ];
    linkages[ CINDEX(5) ] = xf[ CINDEX(5) ] - xi[ CINDEX(5) ];
    linkages[ CINDEX(6) ] = xf[ CINDEX(6) ] - xi[ CINDEX(6) ];
    linkages[ CINDEX(7) ] = t0b - tfa;

    // Linking phases 2 and 3

    get_final_states(        xf, xad, 2, workspace );
    get_initial_states(      xi, xad, 3, workspace );
    tfa = get_final_time(        xad, 2, workspace );
    t0b = get_initial_time(      xad, 3, workspace );

    wtf2 = xf[ CINDEX(7) ];

    linkages[ CINDEX(8) ]  = xf[ CINDEX(1) ] - xi[ CINDEX(1) ];
    linkages[ CINDEX(9) ]  = xf[ CINDEX(2) ] - xi[ CINDEX(2) ];
    linkages[ CINDEX(10) ] = xf[ CINDEX(3) ] - xi[ CINDEX(3) ];
    linkages[ CINDEX(11) ] = xf[ CINDEX(4) ] - xi[ CINDEX(4) ];
    linkages[ CINDEX(12) ] = xf[ CINDEX(5) ] - xi[ CINDEX(5) ];
    linkages[ CINDEX(13) ] = xf[ CINDEX(6) ] - xi[ CINDEX(6) ];
    linkages[ CINDEX(14) ] = t0b - tfa;

    // Linking phases 3 and 4

    get_final_states(      xf, xad, 3, workspace );
    get_initial_states(    xi, xad, 4, workspace );
    tfa = get_final_time(      xad, 3, workspace );
    t0b = get_initial_time(    xad, 4, workspace );

    wti4 = xi[ CINDEX(7) ];


    linkages[ CINDEX(15) ] = xf[ CINDEX(1) ] - xi[ CINDEX(1) ];
    linkages[ CINDEX(16) ] = xf[ CINDEX(2) ] - xi[ CINDEX(2) ];
    linkages[ CINDEX(17) ] = xf[ CINDEX(3) ] - xi[ CINDEX(3) ];
    linkages[ CINDEX(18) ] = xf[ CINDEX(4) ] - xi[ CINDEX(4) ];
    linkages[ CINDEX(19) ] = xf[ CINDEX(5) ] - xi[ CINDEX(5) ];
    linkages[ CINDEX(20) ] = xf[ CINDEX(6) ] - xi[ CINDEX(6) ];
    linkages[ CINDEX(21) ] = t0b - tfa;

    // Linking the weight at the end of phase 2 with the weight at the beginning of phase 4

    linkages[ CINDEX(22) ] = wtf2 - wti4;

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

    problem.nphases   			= 4;
    problem.nlinkages                   = 22;

    psopt_level1_setup(problem);

/////////////////////////////////////////////////////////////////////////////
/////////   Define phase related information & do level 2 setup  ////////////
/////////////////////////////////////////////////////////////////////////////

    problem.phases(1).nstates   		= 6;
    problem.phases(1).ncontrols 		= 0;
    problem.phases(1).nparameters               = 0;
    problem.phases(1).nevents   	        = 6;
    problem.phases(1).npath     		= 0;
    problem.phases(1).nodes                     = 10;

    problem.phases(2).nstates   		= 7;
    problem.phases(2).ncontrols 		= 2;
    problem.phases(2).nparameters               = 0;
    problem.phases(2).nevents   	        = 1;
    problem.phases(2).npath     		= 0;
    problem.phases(2).nodes                     = 10;

    problem.phases(3).nstates   		= 6;
    problem.phases(3).ncontrols 		= 0;
    problem.phases(3).nparameters               = 0;
    problem.phases(3).nevents   	        = 0;
    problem.phases(3).npath     		= 0;
    problem.phases(3).nodes                     = 10;

    problem.phases(4).nstates   		= 7;
    problem.phases(4).ncontrols 		= 2;
    problem.phases(4).nparameters               = 0;
    problem.phases(4).nevents   	        = 5;
    problem.phases(4).npath     		= 0;
    problem.phases(4).nodes                     = 10;

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


    problem.phases(1).bounds.lower.states(1) = 10.e6;
    problem.phases(1).bounds.lower.states(2) = -1;
    problem.phases(1).bounds.lower.states(3) = -1;
    problem.phases(1).bounds.lower.states(4) = -1;
    problem.phases(1).bounds.lower.states(5) = -1;
    problem.phases(1).bounds.lower.states(6) = pi;


    problem.phases(1).bounds.upper.states(1) = 2e8;
    problem.phases(1).bounds.upper.states(2) = 1;
    problem.phases(1).bounds.upper.states(3) = 1;
    problem.phases(1).bounds.upper.states(4) = 1;
    problem.phases(1).bounds.upper.states(5) = 1;
    problem.phases(1).bounds.upper.states(6) = 30*pi;



    problem.phases(1).bounds.lower.events(1)  = pti;
    problem.phases(1).bounds.lower.events(2)  = fti;
    problem.phases(1).bounds.lower.events(3)  = gti;
    problem.phases(1).bounds.lower.events(4)  = hti;
    problem.phases(1).bounds.lower.events(5)  = kti;
    problem.phases(1).bounds.lower.events(6)  = Lti;



    problem.phases(1).bounds.upper.events(1)  = pti;
    problem.phases(1).bounds.upper.events(2)  = fti;
    problem.phases(1).bounds.upper.events(3)  = gti;
    problem.phases(1).bounds.upper.events(4)  = hti;
    problem.phases(1).bounds.upper.events(5)  = kti;
    problem.phases(1).bounds.upper.events(6)  = Lti;


    problem.phases(1).bounds.lower.StartTime    = 0.0;
    problem.phases(1).bounds.upper.StartTime    = 0.0;

    problem.phases(1).bounds.lower.EndTime      = 2000.0;
    problem.phases(1).bounds.upper.EndTime      = 3000.0;

    // BOUNDS FOR PHASE 2

    problem.phases(2).bounds.lower.states(1) = 10.e6;
    problem.phases(2).bounds.lower.states(2) = -1;
    problem.phases(2).bounds.lower.states(3) = -1;
    problem.phases(2).bounds.lower.states(4) = -1;
    problem.phases(2).bounds.lower.states(5) = -1;
    problem.phases(2).bounds.lower.states(6) = pi;
    problem.phases(2).bounds.lower.states(7) = 0.0;

    problem.phases(2).bounds.upper.states(1) = 2.e8;
    problem.phases(2).bounds.upper.states(2) = 1;
    problem.phases(2).bounds.upper.states(3) = 1;
    problem.phases(2).bounds.upper.states(4) = 1;
    problem.phases(2).bounds.upper.states(5) = 1;
    problem.phases(2).bounds.upper.states(6) = 30*pi;
    problem.phases(2).bounds.upper.states(7) = 2.0;


    problem.phases(2).bounds.lower.controls(1) = -pi;
    problem.phases(2).bounds.lower.controls(2) = -pi;
    problem.phases(2).bounds.upper.controls(1) = pi;
    problem.phases(2).bounds.upper.controls(2) = pi;

    problem.phases(2).bounds.lower.events(1)  = wti;
    problem.phases(2).bounds.upper.events(1)  = wti;


    problem.phases(2).bounds.lower.StartTime    = 2000;
    problem.phases(2).bounds.upper.StartTime    = 3000;

    problem.phases(2).bounds.lower.EndTime      = 2100;
    problem.phases(2).bounds.upper.EndTime      = 3100;


    // BOUNDS FOR PHASE 3

    problem.phases(3).bounds.lower.states(1) = 10.e6;
    problem.phases(3).bounds.lower.states(2) = -1;
    problem.phases(3).bounds.lower.states(3) = -1;
    problem.phases(3).bounds.lower.states(4) = -1;
    problem.phases(3).bounds.lower.states(5) = -1;
    problem.phases(3).bounds.lower.states(6) = pi;


    problem.phases(3).bounds.upper.states(1) = 2.e8;
    problem.phases(3).bounds.upper.states(2) = 1.0;
    problem.phases(3).bounds.upper.states(3) = 1.0;
    problem.phases(3).bounds.upper.states(4) = 1.0;
    problem.phases(3).bounds.upper.states(5) = 1.0;
    problem.phases(3).bounds.upper.states(6) = 30*pi;


    problem.phases(3).bounds.lower.StartTime    = 2100;
    problem.phases(3).bounds.upper.StartTime    = 3100;

    problem.phases(3).bounds.lower.EndTime      = 21600;
    problem.phases(3).bounds.upper.EndTime      = 21800;


    // BOUNDS FOR PHASE 4

    problem.phases(4).bounds.lower.states(1) = 10.e6;
    problem.phases(4).bounds.lower.states(2) = -1;
    problem.phases(4).bounds.lower.states(3) = -1;
    problem.phases(4).bounds.lower.states(4) = -1;
    problem.phases(4).bounds.lower.states(5) = -1;
    problem.phases(4).bounds.lower.states(6) = pi;
    problem.phases(4).bounds.lower.states(7) = 0.0;

    problem.phases(4).bounds.upper.states(1) = 2.e8;
    problem.phases(4).bounds.upper.states(2) = 1;
    problem.phases(4).bounds.upper.states(3) = 1;
    problem.phases(4).bounds.upper.states(4) = 1;
    problem.phases(4).bounds.upper.states(5) = 1;
    problem.phases(4).bounds.upper.states(6) = 30*pi;
    problem.phases(4).bounds.upper.states(7) = 2.0;


    problem.phases(4).bounds.lower.controls(1) = -pi;
    problem.phases(4).bounds.lower.controls(2) = -pi;
    problem.phases(4).bounds.upper.controls(1) = pi;
    problem.phases(4).bounds.upper.controls(2) = pi;


    problem.phases(4).bounds.lower.events(1)  = ptf;
    problem.phases(4).bounds.lower.events(2)  = ftf;
    problem.phases(4).bounds.lower.events(3)  = gtf;
    problem.phases(4).bounds.lower.events(4)  = htf;
    problem.phases(4).bounds.lower.events(5)  = ktf;

    problem.phases(4).bounds.upper.events(1)  = ptf;
    problem.phases(4).bounds.upper.events(2)  = ftf;
    problem.phases(4).bounds.upper.events(3)  = gtf;
    problem.phases(4).bounds.upper.events(4)  = htf;
    problem.phases(4).bounds.upper.events(5)  = ktf;

    problem.phases(4).bounds.lower.StartTime    = 21600;
    problem.phases(4).bounds.upper.StartTime    = 21800;

    problem.phases(4).bounds.lower.EndTime      = 21650;
    problem.phases(4).bounds.upper.EndTime      = 21900;


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

    int nnodes;
    int ncontrols;
    int nstates;
    int iphase;
    DMatrix x_guess, u_guess, time_guess, param_guess, xini, xfinal;

    // Phase 1

    nnodes    			= problem.phases(1).nodes(1);
    nstates                     = problem.phases(1).nstates;
    iphase = 1;

    x_guess    =  zeros(nstates,nnodes);
    time_guess =  linspace(0.0,2690,nnodes);

    xini.Resize(6,1);
    xini(1)= pti; xini(2)=fti;xini(3)=gti;xini(4)=hti;xini(5)=kti;xini(6)=Lti;



    rk4_propagate( dae, u_guess, time_guess, xini, param_guess, problem, iphase, x_guess, NULL);

    tra(x_guess).Print("x_guess(iphase=1)");

    xfinal = x_guess(colon(),nnodes);

    problem.phases(1).guess.states = x_guess;
    problem.phases(1).guess.time   = time_guess;




    // Phase 2

    nnodes    			= problem.phases(2).nodes(1);
    nstates                     = problem.phases(2).nstates;
    ncontrols                   = problem.phases(2).ncontrols;
    iphase = 2;

    u_guess    =  zeros(ncontrols,nnodes);
    x_guess    =  zeros(nstates,nnodes);
    time_guess =  linspace(2690,2840,nnodes);

    xini.Resize(7,1);
    xini(1)= xfinal(1); xini(2)=xfinal(2);xini(3)=xfinal(3);xini(4)=xfinal(4);xini(5)=xfinal(5);xini(6)=xfinal(6);
    xini(7)= wti;

    u_guess(1,colon()) = 0.148637e-2*D2R*ones(1,nnodes);
    u_guess(2,colon()) = -9.08446*D2R*ones(1,nnodes);

    rk4_propagate( dae, u_guess, time_guess, xini, param_guess, problem, iphase, x_guess, NULL);

    tra(x_guess).Print("x_guess(iphase=2)");


    xfinal = x_guess(colon(),nnodes);
    double wtf2__ = xfinal(7);

    problem.phases(2).guess.states   = x_guess;
    problem.phases(2).guess.controls = u_guess;
    problem.phases(2).guess.time     = time_guess;

    // Phase 3

    nnodes    			= problem.phases(3).nodes(1);
    nstates                     = problem.phases(3).nstates;
    iphase = 3;


    x_guess    =  zeros(nstates,nnodes);
    time_guess =  linspace(2840,21650,nnodes);


    xini.Resize(6,1);
    xini(1)= xfinal(1); xini(2)=xfinal(2);xini(3)=xfinal(3);xini(4)=xfinal(4);xini(5)=xfinal(5);xini(6)=xfinal(6);

    rk4_propagate( dae, u_guess, time_guess, xini, param_guess, problem, iphase, x_guess, NULL);

    tra(x_guess).Print("x_guess(iphase=3)");

    xfinal = x_guess(colon(),nnodes);

    problem.phases(3).guess.states = x_guess;
    problem.phases(3).guess.time   = time_guess;

    // Phase 4

    nnodes    			= problem.phases(4).nodes(1);
    nstates                     = problem.phases(4).nstates;
    ncontrols                   = problem.phases(4).ncontrols;
    iphase = 4;

    u_guess    =  zeros(ncontrols,nnodes);
    x_guess    =  zeros(nstates,nnodes);
    time_guess =  linspace(21650,21700,nnodes);

    u_guess(1,colon()) = -0.136658e-2*D2R*ones(1,nnodes);
    u_guess(2,colon()) = 49.7892*D2R*ones(1,nnodes);

    xini.Resize(7,1);
    xini(1)= xfinal(1); xini(2)=xfinal(2);xini(3)=xfinal(3);xini(4)=xfinal(4);xini(5)=xfinal(5);xini(6)=xfinal(6);
    xini(7)= wtf2__;

    rk4_propagate( dae, u_guess, time_guess, xini, param_guess, problem, iphase, x_guess, NULL);

    tra(x_guess).Print("x_guess(iphase=4)");


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



    DMatrix x, t, w2, xi, w4,  ti;

    x      = solution.get_states_in_phase(1);
    t      = solution.get_time_in_phase(1);

    w2     = solution.get_states_in_phase(2);
    w2     = w2(7,colon());

    w4     = solution.get_states_in_phase(4);
    w4     = w4(7,colon());


    for(int i=2;i<=problem.nphases;i++) {
      	     xi      = solution.get_states_in_phase(i);
    	     ti      = solution.get_time_in_phase(i);

	xi = xi(colon(1,6),colon());

	x = x || xi;

        t = t || ti;

    }

    DMatrix u_phase2     = solution.get_controls_in_phase(2);
    DMatrix u_phase4     = solution.get_controls_in_phase(4);
    DMatrix t2           = solution.get_time_in_phase(2);
    DMatrix t4           = solution.get_time_in_phase(4);




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

    DMatrix x1 = x(1,colon())/1.e6;
    DMatrix x2 = x(2,colon());
    DMatrix x3 = x(3,colon());
    DMatrix x4 = x(4,colon());
    DMatrix x5 = x(5,colon());
    DMatrix x6 = x(6,colon());
    DMatrix x7 = x(7,colon());
    DMatrix theta_phase2;
    DMatrix theta_phase4;
    DMatrix phi_phase2;
    DMatrix phi_phase4;
    DMatrix r;

    theta_phase2 = u_phase2(1,colon())*R2D;
    phi_phase2    = u_phase2(2,colon())*R2D;
    theta_phase4  = u_phase4(1,colon())*R2D;
    phi_phase4    = u_phase4(2,colon())*R2D;

    compute_cartesian_trajectory(x,r);

    r.Save("r.dat");

    double ft2km = 0.0003048;

    r = r*ft2km;


    plot(t2,theta_phase2,problem.name+": thrust theta phase 2","time (s)", "theta (deg)", "theta");

    plot(t2,phi_phase2,problem.name+": thrust phi phase 2","time (s)", "phi (deg)", "phi");

    plot(t4,theta_phase4,problem.name+": thrust theta phase 4","time (s)", "theta (deg)", "theta");

    plot(t4,phi_phase4,problem.name+": thrust phi phase 4","time (s)", "phi (deg)", "phi");

    plot(t2,theta_phase2,problem.name+":  thrust pitch angle phase 2","time (s)", "theta (deg)", "theta",
	 "pdf", "theta2.pdf");

    plot(t2,phi_phase2,problem.name+": thrust angle phase 2","time (s)", "phi (deg)", "phi",
	 "pdf", "phi2.pdf");

    plot(t4,theta_phase4,problem.name+": thrust pitch angle phase 4","time (s)", "theta (deg)", "theta",
	 "pdf", "theta4.pdf");

    plot(t4,phi_phase4,problem.name+": thrust yaw angle phase 4","time (s)", "phi (deg)", "phi",
    	 "pdf", "phi4.pdf");

    plot3(r(1,colon()), r(2,colon()), r(3,colon()), "Two burn trasnfer trajectory", "x (km)", "y (km)", "z (km)",
	   NULL, NULL, "30,110");

    plot3(r(1,colon()), r(2,colon()), r(3,colon()), "Two burn transfer trajectory", "x (km)", "y (km)", "z (km)",
	   "pdf", "trajectory.pdf", "30,110");

    plot(r(1,colon()), r(2,colon()),  "Two burn trajectory - projection on the equatorial plane",
	    "x (km)", "y (km)");



}

////////////////////////////////////////////////////////////////////////////
///////////////////////      END OF FILE     ///////////////////////////////
////////////////////////////////////////////////////////////////////////////
