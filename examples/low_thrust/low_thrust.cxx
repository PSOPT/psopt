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

   if (iphase == 1) {
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

   adouble p = states[ 0 ];
   adouble f = states[ 1 ];
   adouble g = states[ 2 ];
   adouble h = states[ 3 ];
   adouble k = states[ 4 ];
   adouble L = states[ 5 ];
   adouble w = states[ 6 ];

   adouble* u  = controls;

   adouble tau  = parameters[ 0 ];

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

   vvec[ 0 ] = v1; vvec[ 1 ] = v2; vvec[ 2 ] = v3;

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



   for(i=0;i<3;i++) {
      DELTA_T[i] =  g0*T*(1.0+0.01*tau)*u[i]/w;
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

   adouble wdot = -T*(1.0+0.01*tau)/Isp;

   derivatives[ 0 ] = pdot;
   derivatives[ 1 ] = fdot;
   derivatives[ 2 ] = gdot;
   derivatives[ 3 ] = hdot;
   derivatives[ 4 ] = kdot;
   derivatives[ 5 ] = Ldot;
   derivatives[ 6 ] = wdot;

   path[ 0 ] = pow( u[0] , 2)  + pow( u[1], 2) + pow( u[2], 2);

}

////////////////////////////////////////////////////////////////////////////
///////////////////  Define the events function ////////////////////////////
////////////////////////////////////////////////////////////////////////////

void events(adouble* e, adouble* initial_states, adouble* final_states,
            adouble* parameters,adouble& t0, adouble& tf, adouble* xad,
            int iphase, Workspace* workspace)

{


   int offset;

   adouble pti = initial_states[ 0 ];
   adouble fti = initial_states[ 1 ];
   adouble gti = initial_states[ 2 ];
   adouble hti = initial_states[ 3 ];
   adouble kti = initial_states[ 4 ];
   adouble Lti = initial_states[ 5 ];
   adouble wti = initial_states[ 6 ];


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
   	e[ 6 ]  = wti;
   }

   if (1 == 1) offset = 7;
   else offset = 0;

   if (iphase == 1 ) {
   	e[ offset + 0 ]  = ptf;
   	e[ offset + 1 ]  = sqrt( ftf*ftf + gtf*gtf );
   	e[ offset + 2 ]  = sqrt( htf*htf + ktf*ktf );
   	e[ offset + 3 ] = ftf*htf + gtf*ktf;
   	e[ offset + 4 ] = gtf*htf - ktf*ftf;
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

    problem.name        					 = "Low thrust transfer problem";
    problem.outfilename                 = "lowthrust.txt";

////////////////////////////////////////////////////////////////////////////
////////////  Define problem level constants & do level 1 setup ////////////
////////////////////////////////////////////////////////////////////////////

    problem.nphases   			          = 1;
    problem.nlinkages                   = 0;

    psopt_level1_setup(problem);

/////////////////////////////////////////////////////////////////////////////
/////////   Define phase related information & do level 2 setup  ////////////
/////////////////////////////////////////////////////////////////////////////

    problem.phases(1).nstates   				= 7;
    problem.phases(1).ncontrols 				= 3;
    problem.phases(1).nparameters         = 1;
    problem.phases(1).nevents   	        	= 12;
    problem.phases(1).npath     				= 1;
    problem.phases(1).nodes               << 80;

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



    problem.phases(1).bounds.lower.parameters <<  tauL;

    problem.phases(1).bounds.upper.parameters <<  tauU;

    problem.phases(1).bounds.lower.states << 10.e6, -0.20, -0.10, -1.0, -0.20, pi, 0.0;

    problem.phases(1).bounds.upper.states << 60.e6, 0.20, 1.0, 1.0, 0.20, 20*pi, 2.0;

    problem.phases(1).bounds.lower.controls << -1.0, -1.0, -1.0;

    problem.phases(1).bounds.upper.controls << 1.0, 1.0, 1.0;

    problem.phases(1).bounds.lower.events << pti, fti, gti, hti, kti, Lti, wti, ptf, event_final_9, event_final_10, event_final_11, event_final_12_lower;

    problem.phases(1).bounds.upper.events << pti, fti, gti, hti, kti, Lti, wti, ptf, event_final_9, event_final_10, event_final_11, event_final_12_upper;

    double EQ_TOL = 0.001;

    problem.phases(1).bounds.upper.path << 1.0+EQ_TOL;
    problem.phases(1).bounds.lower.path << 1.0-EQ_TOL;


    problem.phases(1).bounds.lower.StartTime    = 0.0;
    problem.phases(1).bounds.upper.StartTime    = 0.0;

    problem.phases(1).bounds.lower.EndTime      = 50000.0;
    problem.phases(1).bounds.upper.EndTime      = 100000.0;


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

    int nnodes    			= 141;
    int ncontrols          = problem.phases(1).ncontrols;
    int nstates            = problem.phases(1).nstates;

    MatrixXd u_guess    =  zeros(ncontrols,nnodes);
    MatrixXd x_guess    =  zeros(nstates,nnodes);
    MatrixXd time_guess =  linspace(0.0,86810.0,nnodes);

    MatrixXd param_guess = -25.0*ones(1,1);

    u_guess    = load_data("../../../examples/low_thrust/U0.dat",ncontrols, nnodes );
    x_guess    = load_data("../../../examples/low_thrust/X0.dat",nstates  , nnodes );
    time_guess = load_data("../../../examples/low_thrust/T0.dat",1        , nnodes );

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


    MatrixXd x, u, t;

    x      = solution.get_states_in_phase(1);
    u      = solution.get_controls_in_phase(1);
    t      = solution.get_time_in_phase(1);


    t = t/3600.0;


    MatrixXd tau = solution.get_parameters_in_phase(1);

    Print(tau,"tau");

////////////////////////////////////////////////////////////////////////////
///////////  Save solution data to files if desired ////////////////////////
////////////////////////////////////////////////////////////////////////////

    Save(x,"x.dat");
    Save(u,"u.dat");
    Save(t,"t.dat");

////////////////////////////////////////////////////////////////////////////
///////////  Plot some results if desired (requires gnuplot) ///////////////
////////////////////////////////////////////////////////////////////////////

    MatrixXd x1 = x.row(0)/1.e6; 
    MatrixXd x2 = x.row(1); 
    MatrixXd x3 = x.row(2); 
    MatrixXd x4 = x.row(3); 
    MatrixXd x5 = x.row(4); 
    MatrixXd x6 = x.row(5); 
    MatrixXd x7 = x.row(6); 
    MatrixXd u1 = u.row(0); 
    MatrixXd u2 = u.row(1); 
    MatrixXd u3 = u.row(2); 

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

    MatrixXd r;

    compute_cartesian_trajectory(x,r);

    double ft2km = 0.0003048;

    r = r*ft2km;

    MatrixXd rnew, tnew;

    tnew = linspace(0.0, t(length(t)-1), 1000);

    resample_trajectory( rnew, tnew, r, t );


    plot3(rnew.row(0), rnew.row(1), rnew.row(2),
	        "Low thrust transfer trajectory", "x (km)", "y (km)", "z (km)",
	         NULL, NULL, "30,97");

    plot3(rnew.row(0), rnew.row(1), rnew.row(2),
	       "Low thrust transfer trajectory", "x (km)", "y (km)", "z (km)",
	       "pdf", "trajectory.pdf", "30,97");


}

////////////////////////////////////////////////////////////////////////////
///////////////////////      END OF FILE     ///////////////////////////////
////////////////////////////////////////////////////////////////////////////
