//////////////////////////////////////////////////////////////////////////
////////////////        low_thrust.cxx               /////////////////////
//////////////////////////////////////////////////////////////////////////
////////////////           PSOPT  Example             ////////////////////
//////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////
//////// Title: Low thrust orbit transfer problem         ////////////////
//////// Last modified: 16 February 2009                  ////////////////
//////// Reference:     Betts  (2001)             	  ////////////////
//////// (See PSOPT handbook for full reference)           ///////////////
//////////////////////////////////////////////////////////////////////////
////////     Copyright (c) Victor M. Becerra, 2009        ////////////////
//////////////////////////////////////////////////////////////////////////
//////// This is part of the PSOPT software library, which ///////////////
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

   if (iphase == 1) {
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
   double Isp  = 450.0;            // [sec]
   double mu   = 1.407645794e16;   // [f2^2/sec^2]
   double g0   = 32.174;           // [ft/sec^2]
   double T    = 4.446618e-3;      // [lb]
   double Re   = 20925662.73;      // [ft]
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
   adouble w = states[ CINDEX(7) ];

   adouble* u  = controls;

   adouble tau  = parameters[ CINDEX(1) ];

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



   for(i=0;i<3;i++) {
      DELTA_T[i] =  g0*T*(1.0+0.01*tau)*u[i]/w;
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

   adouble wdot = -T*(1.0+0.01*tau)/Isp;

   derivatives[ CINDEX(1) ] = pdot;
   derivatives[ CINDEX(2) ] = fdot;
   derivatives[ CINDEX(3) ] = gdot;
   derivatives[ CINDEX(4) ] = hdot;
   derivatives[ CINDEX(5) ] = kdot;
   derivatives[ CINDEX(6) ] = Ldot;
   derivatives[ CINDEX(7) ] = wdot;

   path[ CINDEX(1) ] = pow( u[CINDEX(1)] , 2)  + pow( u[CINDEX(2)], 2) + pow( u[CINDEX(3)], 2);

}

////////////////////////////////////////////////////////////////////////////
///////////////////  Define the events function ////////////////////////////
////////////////////////////////////////////////////////////////////////////

void events(adouble* e, adouble* initial_states, adouble* final_states,
            adouble* parameters,adouble& t0, adouble& tf, adouble* xad,
            int iphase, Workspace* workspace)

{


   int offset;

   adouble pti = initial_states[ CINDEX(1) ];
   adouble fti = initial_states[ CINDEX(2) ];
   adouble gti = initial_states[ CINDEX(3) ];
   adouble hti = initial_states[ CINDEX(4) ];
   adouble kti = initial_states[ CINDEX(5) ];
   adouble Lti = initial_states[ CINDEX(6) ];
   adouble wti = initial_states[ CINDEX(7) ];


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
   	e[ CINDEX(7) ]  = wti;
   }

   if (1 == 1) offset = 7;
   else offset = 0;

   if (iphase == 1 ) {
   	e[ offset + CINDEX(1) ]  = ptf;
   	e[ offset + CINDEX(2) ]  = sqrt( ftf*ftf + gtf*gtf );
   	e[ offset + CINDEX(3) ]  = sqrt( htf*htf + ktf*ktf );
   	e[ offset + CINDEX(4) ] = ftf*htf + gtf*ktf;
   	e[ offset + CINDEX(5) ] = gtf*htf - ktf*ftf;
   }

}



///////////////////////////////////////////////////////////////////////////
///////////////////  Define the phase linkages function ///////////////////
///////////////////////////////////////////////////////////////////////////

void linkages( adouble* linkages, adouble* xad, Workspace* workspace)
{
 //  auto_link_multiple(linkages, xad, 1);
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

    problem.name        		= "Low thrust transfer problem";
    problem.outfilename                 = "lowthrust.txt";

////////////////////////////////////////////////////////////////////////////
////////////  Define problem level constants & do level 1 setup ////////////
////////////////////////////////////////////////////////////////////////////

    problem.nphases   			= 1;
    problem.nlinkages                   = 0;

    psopt_level1_setup(problem);

/////////////////////////////////////////////////////////////////////////////
/////////   Define phase related information & do level 2 setup  ////////////
/////////////////////////////////////////////////////////////////////////////

    problem.phases(1).nstates   		= 7;
    problem.phases(1).ncontrols 		= 3;
    problem.phases(1).nparameters               = 1;
    problem.phases(1).nevents   	        = 12;
    problem.phases(1).npath     		= 1;
    problem.phases(1).nodes                     = 80;

    psopt_level2_setup(problem, algorithm);


////////////////////////////////////////////////////////////////////////////
///////////////////  Enter problem bounds information //////////////////////
////////////////////////////////////////////////////////////////////////////

    double tauL = -50.0;
    double tauU =   0.0;


    double pti = 21837080.052835;
    double fti = 0.0;
    double gti = 0.0;
    double hti = -0.25396764647494;
    double kti = 0.0;
    double Lti = pi;
    double wti = 1.0;

    double wtf_guess;

    double  SISP = 450.0;
    double  DELTAV = 22741.1460;
    double  CM2W = 32.174;

    wtf_guess  = wti*exp(-DELTAV/(CM2W*SISP));

    double ptf            = 40007346.015232;
    double event_final_9  = 0.73550320568829;
    double event_final_10  = 0.61761258786099;
    double event_final_11 = 0.0;
    double event_final_12_upper = 0.0;
    double event_final_12_lower = -10.0;



    problem.phases(1).bounds.lower.parameters(1) =  tauL;
    problem.phases(1).bounds.upper.parameters(1) =  tauU;

    problem.phases(1).bounds.lower.states(1) = 10.e6;
    problem.phases(1).bounds.lower.states(2) = -0.20;
    problem.phases(1).bounds.lower.states(3) = -0.10;
    problem.phases(1).bounds.lower.states(4) = -1.0;
    problem.phases(1).bounds.lower.states(5) = -0.20;
    problem.phases(1).bounds.lower.states(6) = pi;
    problem.phases(1).bounds.lower.states(7) = 0.0;

    problem.phases(1).bounds.upper.states(1) = 60.e6;
    problem.phases(1).bounds.upper.states(2) = 0.20;
    problem.phases(1).bounds.upper.states(3) = 1.0;
    problem.phases(1).bounds.upper.states(4) = 1.0;
    problem.phases(1).bounds.upper.states(5) = 0.20;
    problem.phases(1).bounds.upper.states(6) = 20*pi;
    problem.phases(1).bounds.upper.states(7) = 2.0;


    problem.phases(1).bounds.lower.controls(1) = -1.0;
    problem.phases(1).bounds.lower.controls(2) = -1.0;
    problem.phases(1).bounds.lower.controls(3) = -1.0;
    problem.phases(1).bounds.upper.controls(1) = 1.0;
    problem.phases(1).bounds.upper.controls(2) = 1.0;
    problem.phases(1).bounds.upper.controls(3) = 1.0;


    problem.phases(1).bounds.lower.events(1)  = pti;
    problem.phases(1).bounds.lower.events(2)  = fti;
    problem.phases(1).bounds.lower.events(3)  = gti;
    problem.phases(1).bounds.lower.events(4)  = hti;
    problem.phases(1).bounds.lower.events(5)  = kti;
    problem.phases(1).bounds.lower.events(6)  = Lti;
    problem.phases(1).bounds.lower.events(7)  = wti;


    problem.phases(1).bounds.lower.events(8) = ptf;
    problem.phases(1).bounds.lower.events(9) = event_final_9;
    problem.phases(1).bounds.lower.events(10) = event_final_10;
    problem.phases(1).bounds.lower.events(11) = event_final_11;
    problem.phases(1).bounds.lower.events(12) = event_final_12_lower;


    problem.phases(1).bounds.upper.events(1)  = pti;
    problem.phases(1).bounds.upper.events(2)  = fti;
    problem.phases(1).bounds.upper.events(3)  = gti;
    problem.phases(1).bounds.upper.events(4)  = hti;
    problem.phases(1).bounds.upper.events(5)  = kti;
    problem.phases(1).bounds.upper.events(6)  = Lti;
    problem.phases(1).bounds.upper.events(7)  = wti;
    problem.phases(1).bounds.upper.events(8)  = ptf;
    problem.phases(1).bounds.upper.events(9)  = event_final_9;
    problem.phases(1).bounds.upper.events(10) = event_final_10;
    problem.phases(1).bounds.upper.events(11) = event_final_11;
    problem.phases(1).bounds.upper.events(12) = event_final_12_upper;

    double EQ_TOL = 0.001;

    problem.phases(1).bounds.upper.path(1) = 1.0+EQ_TOL;
    problem.phases(1).bounds.lower.path(1) = 1.0-EQ_TOL;


    problem.phases(1).bounds.lower.StartTime    = 0.0;
    problem.phases(1).bounds.upper.StartTime    = 0.0;

    problem.phases(1).bounds.lower.EndTime      = 50000.0;
    problem.phases(1).bounds.upper.EndTime      = 100000.0;


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

    int nnodes    			= 141;
    int ncontrols                       = problem.phases(1).ncontrols;
    int nstates                         = problem.phases(1).nstates;

    DMatrix u_guess    =  zeros(ncontrols,nnodes);
    DMatrix x_guess    =  zeros(nstates,nnodes);
    DMatrix time_guess =  linspace(0.0,86810.0,nnodes);

    DMatrix param_guess = -25.0*ones(1,1);

    u_guess.Load("U0.dat");
    x_guess.Load("X0.dat");
    time_guess.Load("T0.dat");

    auto_phase_guess(problem, u_guess, x_guess, param_guess, time_guess);


////////////////////////////////////////////////////////////////////////////
///////////////////  Enter algorithm options  //////////////////////////////
////////////////////////////////////////////////////////////////////////////


    algorithm.nlp_iter_max                = 1000;
    algorithm.nlp_tolerance               = 1.e-6;
    algorithm.nlp_method                  = "IPOPT";
    algorithm.scaling                     = "automatic";
    algorithm.derivatives                 = "automatic";
    algorithm.defect_scaling              = "jacobian-based";
    algorithm.jac_sparsity_ratio          =  0.11; // 0.05;
    algorithm.collocation_method          = "trapezoidal";
    algorithm.mesh_refinement             = "automatic";
    algorithm.mr_max_increment_factor     = 0.2;





////////////////////////////////////////////////////////////////////////////
///////////////////  Now call PSOPT to solve the problem   //////////////////
////////////////////////////////////////////////////////////////////////////

    psopt(solution, problem, algorithm);

////////////////////////////////////////////////////////////////////////////
///////////  Extract relevant variables from solution structure   //////////
////////////////////////////////////////////////////////////////////////////


    DMatrix x, u, t;

    x      = solution.get_states_in_phase(1);
    u      = solution.get_controls_in_phase(1);
    t      = solution.get_time_in_phase(1);


    t = t/3600.0;


    DMatrix tau = solution.get_parameters_in_phase(1);

    tau.Print("tau");

////////////////////////////////////////////////////////////////////////////
///////////  Save solution data to files if desired ////////////////////////
////////////////////////////////////////////////////////////////////////////

    x.Save("x.dat");
    u.Save("u.dat");
    t.Save("t.dat");

////////////////////////////////////////////////////////////////////////////
///////////  Plot some results if desired (requires gnuplot) ///////////////
////////////////////////////////////////////////////////////////////////////

    DMatrix x1 = x(1,colon())/1.e6;
    DMatrix x2 = x(2,colon());
    DMatrix x3 = x(3,colon());
    DMatrix x4 = x(4,colon());
    DMatrix x5 = x(5,colon());
    DMatrix x6 = x(6,colon());
    DMatrix x7 = x(7,colon());
    DMatrix u1 = u(1,colon());
    DMatrix u2 = u(2,colon());
    DMatrix u3 = u(3,colon());

    plot(t,x1,problem.name+": states", "time (h)", "p (1000000 ft)","p (1000000 ft)");

    plot(t,x2,problem.name+": states", "time (h)", "f","f");

    plot(t,x3,problem.name+": states", "time (h)", "g","g");

    plot(t,x4,problem.name+": states", "time (h)", "h","h");

    plot(t,x5,problem.name+": states", "time (h)", "k","k");

    plot(t,x6,problem.name+": states", "time (h)", "L (rev)","L (rev)");

    plot(t,x7,problem.name+": states", "time (h)", "w (lb)","w (lb)");

    plot(t,u1,problem.name+": controls","time (h)", "ur", "ur");

    plot(t,u2,problem.name+": controls","time (h)", "ut", "ut");

    plot(t,u3,problem.name+": controls","time (h)", "uh", "uh");

    plot(t,x1,problem.name+": states", "time (h)", "p (1000000 ft)","p (1000000 ft)",
	         "pdf","lowthr_x1.pdf");

    plot(t,x2,problem.name+": states", "time (h)", "f","f",
	 	 "pdf","lowthr_x2.pdf");

    plot(t,x3,problem.name+": states", "time (h)", "g","g",
	 	 "pdf","lowthr_x3.pdf");

    plot(t,x4,problem.name+": states", "time (h)", "h","h",
	 	 "pdf","lowthr_x4.pdf");

    plot(t,x5,problem.name+": states", "time (h)", "k","k",
    	        "pdf","lowthr_x5.pdf");

    plot(t,x6,problem.name+": states", "time (h)", "L (rev)","L (rev)",
	 	 "pdf","lowthr_x6.pdf");

    plot(t,x7,problem.name+": states", "time (h)", "w (lb)","w (lb)",
	 	 "pdf","lowthr_x7.pdf");


    plot(t,u1,problem.name+": controls","time (h)", "ur", "ur",
	 	 "pdf","lowthr_u1.pdf");

    plot(t,u2,problem.name+": controls","time (h)", "ut", "ut",
	 	 "pdf","lowthr_u2.pdf");

    plot(t,u3,problem.name+": controls","time (h)", "uh", "uh",
	 	 "pdf","lowthr_u3.pdf");

    DMatrix r;

    compute_cartesian_trajectory(x,r);

    double ft2km = 0.0003048;

    r = r*ft2km;

    DMatrix rnew, tnew;

    tnew = linspace(0.0, t("end"), 1000);

    resample_trajectory( rnew, tnew, r, t );


    plot3(rnew(1,colon()), rnew(2,colon()), rnew(3,colon()),
	        "Low thrust transfer trajectory", "x (km)", "y (km)", "z (km)",
	         NULL, NULL, "30,97");

    plot3(rnew(1,colon()), rnew(2,colon()), rnew(3,colon()),
	       "Low thrust transfer trajectory", "x (km)", "y (km)", "z (km)",
	       "pdf", "trajectory.pdf", "30,97");


}

////////////////////////////////////////////////////////////////////////////
///////////////////////      END OF FILE     ///////////////////////////////
////////////////////////////////////////////////////////////////////////////
