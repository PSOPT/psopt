//////////////////////////////////////////////////////////////////////////
////////////////        geodesic.cxx                 /////////////////////
//////////////////////////////////////////////////////////////////////////
////////////////           PSOPT  Example             ////////////////////
//////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////
//////// Title:         Geodesic calculation problem      ////////////////
//////// Last modified: 22 February 2019                   ////////////////
//////// Reference:  PROPT User's Guide                   ////////////////
//////// (See PSOPT handbook for full reference)          ////////////////
//////////////////////////////////////////////////////////////////////////
////////     Copyright (c) Victor M. Becerra, 2019        ////////////////
//////////////////////////////////////////////////////////////////////////
//////// This is part of the PSOPT software library, which ///////////////
//////// is distributed under the terms of the GNU Lesser ////////////////
//////// General Public License (LGPL)                    ////////////////
//////////////////////////////////////////////////////////////////////////

#include "psopt.h"

using namespace std;

using namespace PSOPT;

typedef struct {
   double V;           
   double a;
   double b; 
} Constants;

//////////////////////////////////////////////////////////////////////////
///////////////////  Define the end point (Mayer) cost function //////////
//////////////////////////////////////////////////////////////////////////

adouble endpoint_cost(adouble* initial_states, adouble* final_states,
                      adouble* parameters,adouble& t0, adouble& tf,
                      adouble* xad, int iphase,Workspace* workspace)
{
   return 0;
}

//////////////////////////////////////////////////////////////////////////
///////////////////  Define the integrand (Lagrange) cost function  //////
//////////////////////////////////////////////////////////////////////////

adouble integrand_cost(adouble* states, adouble* controls,
                       adouble* parameters, adouble& time, adouble* xad,
                       int iphase, Workspace* workspace)
{
     
    Constants* C = (Constants*) workspace->user_data;

    double V = C->V;    
    
    adouble theta = controls[ 0 ];
    adouble phi   = controls[ 1 ];

    // These are the components of the velocity vector in spherical coordinates.
    adouble dxdt = V*sin(theta)*cos(phi);
    adouble dydt = V*sin(theta)*sin(phi);
    adouble dzdt = V*cos(theta);

    // The integrand is the norm of the speed vector
    adouble L =  sqrt( pow(dxdt,2.0) + pow(dydt,2.0)+pow(dzdt,2.0) );

    return  L;
}


//////////////////////////////////////////////////////////////////////////
///////////////////  Define the DAE's ////////////////////////////////////
//////////////////////////////////////////////////////////////////////////

void dae(adouble* derivatives, adouble* path, adouble* states,
         adouble* controls, adouble* parameters, adouble& time,
         adouble* xad, int iphase, Workspace* workspace)
{

   Constants* C = (Constants*) workspace->user_data;

   adouble x    = states[ 0 ];
   adouble y    = states[ 1 ];
   adouble z    = states[ 2 ];

   double V = C->V; // Speed 
   double a = C->a; // Semi-major axis
   double b = C->b; // Semi-minor axis
   
   // These are the angles of the velocity vector in spherical coordinates
   adouble theta = controls[ 0 ];
   adouble phi   = controls[ 1 ];

   // Simple kinematic equations of motion in spherical coordinates
   adouble dxdt = V*sin(theta)*cos(phi);
   adouble dydt = V*sin(theta)*sin(phi);
   adouble dzdt = V*cos(theta);


   derivatives[ 0 ] = dxdt;
   derivatives[ 1 ] = dydt;
   derivatives[ 2 ] = dzdt;

   // This is the geodesic constraint to stay on the surface of the spheroid
   
   path[ 0 ] = x*x/(a*a) + y*y/(a*a) + z*z/(b*b) - 1.0;

}

////////////////////////////////////////////////////////////////////////////
///////////////////  Define the events function ////////////////////////////
////////////////////////////////////////////////////////////////////////////

void events(adouble* e, adouble* initial_states, adouble* final_states,
            adouble* parameters,adouble& t0, adouble& tf, adouble* xad,
            int iphase, Workspace* workspace)

{
   adouble x0 = initial_states[ 0 ];
   adouble y0 = initial_states[ 1 ];
   adouble z0 = initial_states[ 2 ];
   adouble xf = final_states[   0 ];
   adouble yf = final_states[   1 ];
   adouble zf = final_states[   2 ];

   e[ 0 ] = x0;
   e[ 1 ] = y0;
   e[ 2 ] = z0;
   e[ 3 ] = xf;
   e[ 4 ] = yf;
   e[ 5 ] = zf;

}



///////////////////////////////////////////////////////////////////////////
///////////////////  Define the phase linkages function ///////////////////
///////////////////////////////////////////////////////////////////////////

void linkages( adouble* linkages, adouble* xad, Workspace* workspace)
{
  // No linkages as this is a single phase problem
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

    problem.name        		          = "Geodesic problem";
    problem.outfilename                 = "geodesic.txt";

////////////////////////////////////////////////////////////////////////////
////////////  Define problem level constants & do level 1 setup ////////////
////////////////////////////////////////////////////////////////////////////

    problem.nphases   			          = 1;
    problem.nlinkages                   = 0;

    psopt_level1_setup(problem);

/////////////////////////////////////////////////////////////////////////////
/////////   Define phase related information & do level 2 setup  ////////////
/////////////////////////////////////////////////////////////////////////////

    problem.phases(1).nstates   		= 3;
    problem.phases(1).ncontrols 		= 2;
    problem.phases(1).nevents   		= 6;
    problem.phases(1).npath         = 1;
    problem.phases(1).nodes         << 40;

    psopt_level2_setup(problem, algorithm);


////////////////////////////////////////////////////////////////////////////
///////////////////  Enter problem bounds information //////////////////////
////////////////////////////////////////////////////////////////////////////

   Constants* C = new Constants;

   problem.user_data = (void*) C;

   C->V = 900.00; // Speed in km/h
   C->a = 6384.0; // Earth's semi-major axis in km 
   C->b = 6353.0; // Earth's semi-minor axis in km

    double a = C->a;
    double b = C->b;
   
    double xL = -a;
    double yL = -a;
    double zL = -b;
    double xU =  a;
    double yU =  a;
    double zU =  b;

    double thetaL =  0.0;
    double thetaU =  pi;
    double phiL   = 0.0;
    double phiU   =  2.0*pi;    
    
    // Coordinates of LHR: 51.4700° N, 0.4543° W
    double lat_lhr = 51.4700*pi/180.0;
    double lon_lhr = 0.4543*pi/180.0;    
     // Coordinates of JFK: 40.6413° N, 73.7781° W
    double lat_jfk = 40.6413*pi/180.0;
    double lon_jfk = 73.7781*pi/180.0;    
    // Below, theta=0 corresponds to 90 deg latitude north, growing positive towards the south
    // while phi=0 corresponds to 0 longitude, growing positive towards the east.
    
    double theta0 = pi/2.0 - lat_jfk;    // Initial elevation angle, for JFK in New York 
    double phi0   = 2.0*pi - lon_jfk;    // Initial azimuth angle, for JFK in New York
    
    double thetaf = pi/2.0 - lat_lhr; // Final elevation angle, for LHR in London
    double phif   = 2.0*pi - lon_lhr; // Final azimuth angle, for LHR in London

    // Here we calculate initial and final Cartesian coordinates using
    // the parametric equations of the ellipsoid.
    
    double x0 = a*sin(theta0)*cos(phi0);
    double y0 = a*sin(theta0)*sin(phi0);
    double z0 = b*cos(theta0);
    double xf = a*sin(thetaf)*cos(phif);
    double yf = a*sin(thetaf)*sin(phif);
    double zf = b*cos(thetaf);


    problem.phases(1).bounds.lower.states << xL, yL, zL;
    problem.phases(1).bounds.upper.states << xU, yU, zU;



    problem.phases(1).bounds.lower.controls << thetaL, phiL;
    problem.phases(1).bounds.upper.controls << thetaU, phiU;


    problem.phases(1).bounds.lower.events << x0, y0, z0, xf, yf, zf;
    
    problem.phases(1).bounds.upper.events << x0, y0, z0, xf, yf, zf;

    problem.phases(1).bounds.lower.path << 0.0;
    problem.phases(1).bounds.upper.path << 0.0;




    problem.phases(1).bounds.lower.StartTime    = 0.0;
    problem.phases(1).bounds.upper.StartTime    = 0.0;

    problem.phases(1).bounds.lower.EndTime      = 3.0;  // lower bound in hours
    problem.phases(1).bounds.upper.EndTime      = 10.0; // upper bound in hours



////////////////////////////////////////////////////////////////////////////
///////////////////  Register problem functions  ///////////////////////////
////////////////////////////////////////////////////////////////////////////


    problem.integrand_cost 	= &integrand_cost;
    problem.endpoint_cost 	   = &endpoint_cost;
    problem.dae             	= &dae;
    problem.events 		      = &events;
    problem.linkages		      = &linkages;

////////////////////////////////////////////////////////////////////////////
///////////////////  Define & register initial guess ///////////////////////
////////////////////////////////////////////////////////////////////////////

    int nnodes    			             = 30;
    int ncontrols                       = problem.phases(1).ncontrols;
    int nstates                         = problem.phases(1).nstates;

    MatrixXd u_guess    =  zeros(ncontrols,nnodes);
    MatrixXd x_guess    =  zeros(nstates,nnodes);
    MatrixXd time_guess =  linspace(0.0,7.0,nnodes);


    // Position guess: interpolate the parametric angles between the two airports and
    // place the points on the spheroid, so the guess satisfies the path constraint.

    MatrixXd ang_guess(2,nnodes);

    ang_guess << linspace(theta0,thetaf,nnodes),
                 linspace(phi0,phif,nnodes);

    for (int i = 0;i< nnodes;i++) {

      x_guess(0,i) = a*sin(ang_guess(0,i))*cos(ang_guess(1,i));
      x_guess(1,i) = a*sin(ang_guess(0,i))*sin(ang_guess(1,i));
      x_guess(2,i) = b*cos(ang_guess(0,i));

    }

    // Control guess: the controls are the spherical angles of the VELOCITY vector, not
    // of the position, so guess them as the direction of travel along the guessed path.
    // Guessing the position angles here instead leaves the initial guess inconsistent
    // with the dynamics, and the solver is then liable to converge to a path several
    // kilometres longer than the geodesic while reporting success.

    double L_guess = 0.0;

    for (int i = 0;i< nnodes;i++) {

      int  ia = (i < nnodes-1) ? i : nnodes-2;
      double dx = x_guess(0,ia+1)-x_guess(0,ia);
      double dy = x_guess(1,ia+1)-x_guess(1,ia);
      double dz = x_guess(2,ia+1)-x_guess(2,ia);
      double dn = sqrt(dx*dx+dy*dy+dz*dz);

      u_guess(0,i) = acos(dz/dn);                        // polar angle of the velocity
      u_guess(1,i) = atan2(dy,dx);                       // azimuth of the velocity
      if ( u_guess(1,i) < 0.0 ) u_guess(1,i) += 2.0*pi;  // keep phi within its bounds

      if (i < nnodes-1) L_guess += dn;

    }

    // Time guess: the length of the guessed path divided by the speed.

    time_guess = linspace(0.0, L_guess/C->V, nnodes);

    problem.phases(1).guess.controls       = u_guess;
    problem.phases(1).guess.states         = x_guess;
    problem.phases(1).guess.time           = time_guess;


////////////////////////////////////////////////////////////////////////////
///////////////////  Enter algorithm options  //////////////////////////////
////////////////////////////////////////////////////////////////////////////
    algorithm.nlp_iter_max                = 1000;
    algorithm.nlp_tolerance               = 1.e-6;
    algorithm.nlp_method                  = "IPOPT";
    algorithm.scaling                     = "automatic";
    algorithm.derivatives                 = "automatic";
    algorithm.collocation_method          = "trapezoidal";
    algorithm.mesh_refinement             = "automatic";
    // The default ODE tolerance of 1e-3 is met by the first trapezoidal mesh, so the
    // refinement stops immediately and the answer is 0.2 km above the geodesic. At 1e-6
    // the refinement runs to completion and the computed length agrees with the exact
    // geodesic on this spheroid (5549.6819 km) to six significant figures.
    algorithm.ode_tolerance               = 1.e-6;


////////////////////////////////////////////////////////////////////////////
///////////////////  Now call PSOPT to solve the problem   /////////////////
////////////////////////////////////////////////////////////////////////////

    psopt(solution, problem, algorithm);

////////////////////////////////////////////////////////////////////////////
///////////  Extract relevant variables from solution structure   //////////
////////////////////////////////////////////////////////////////////////////


    MatrixXd states    = solution.get_states_in_phase(1);
    MatrixXd controls  = solution.get_controls_in_phase(1);
    MatrixXd t         = solution.get_time_in_phase(1);


    MatrixXd x = states.row(0); 
    MatrixXd y = states.row(1); 
    MatrixXd z = states.row(2); 
    
    MatrixXd theta = controls.row(0); 
    MatrixXd phi   = controls.row(1); 

////////////////////////////////////////////////////////////////////////////
///////////  Save solution data to files if desired ////////////////////////
////////////////////////////////////////////////////////////////////////////

    Save(x,"x.dat");
    Save(y,"y.dat");
    Save(z, "z.dat");
    Save(theta,"theta.dat");
    Save(phi, "phi.dat");
    Save(t,"t.dat");


////////////////////////////////////////////////////////////////////////////
///////////  Plot some results if desired (requires gnuplot) ///////////////
////////////////////////////////////////////////////////////////////////////

    plot3(x, y, z,
	        "Geodesic problem", "x", "y", "z",
	         NULL, NULL, "30,97");

    plot3(x, y, z,
	        "Geodesic problem", "x", "y", "z",
	       "pdf", "trajectory.pdf", "30,97");
	       
	 plot(t,states,problem.name, "time (h)", "states (km)", "x y z");

    plot(t,controls,problem.name, "time (h)", "controls (rad)", "theta phi");
                    
	       
	 plot(t,states,problem.name, "time (h)", "states (km)", "x y z",
                           "pdf", "geodesic_states.pdf");

    plot(t,controls,problem.name, "time (h)", "controls (rad)", "theta phi",
                           "pdf", "geodesic_controls.pdf");

}

////////////////////////////////////////////////////////////////////////////
///////////////////////      END OF FILE     ///////////////////////////////
////////////////////////////////////////////////////////////////////////////
