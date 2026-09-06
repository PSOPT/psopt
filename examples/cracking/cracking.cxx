//////////////////////////////////////////////////////////////////////////
//////////////////         cracking.cxx       ////////////////////////////
//////////////////////////////////////////////////////////////////////////
////////////////           PSOPT  Example             ////////////////////
//////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////
//////// Title:  Catalytic Cracking of Gas Oil            ////////////////
//////// Last modified: 06 September 2026                 ////////////////
//////// Reference:     User's guide for DIRCOL  	  ////////////////
//////// (See PSOPT handbook for full reference)           ///////////////
//////////////////////////////////////////////////////////////////////////
////////     Copyright (c) Victor M. Becerra, 2009        ////////////////
//////////////////////////////////////////////////////////////////////////
//////// This is part of the PSOPT software library, which ///////////////
//////// is distributed under the terms of the GNU Lesser ////////////////
//////// General Public License (LGPL)                    ////////////////
//////////////////////////////////////////////////////////////////////////
//
// Run with no arguments, this estimates the three rate constants and prints
// the 95 per cent confidence intervals PSOPT computes from the linearized
// (Wald) covariance of the estimates.
//
// Run with arguments, it additionally computes a PROFILE LIKELIHOOD for one
// of the parameters:
//
//     cracking <index> [lo] [hi] [ngrid] [nobserved]
//
//     index      1, 2 or 3: which parameter to profile
//     lo, hi     the range of that parameter to scan   (default 0 to 2.2)
//     ngrid      number of grid points                 (default 45)
//     nobserved  1 fits y1 only, 2 fits y1 and y2      (default 2)
//
// The profile is obtained by fixing the chosen parameter at each grid value,
// re-estimating the others, and recording the optimal sum of squares J(theta).
// Fixing a parameter needs no special support: its lower and upper bounds are
// set to the same number. The 95 per cent profile-likelihood interval is the
// set of values for which
//
//     J(theta) <= J* [ 1 + t^2 / (Ns - nf) ],   t = t^{0.975}_{Ns-nf}
//
// which is the standard F-test threshold for one parameter. Unlike the Wald
// interval it needs no assumption that the model is linear near the estimate,
// and it cannot return an interval for a parameter the data do not determine.
//
// Two runs are worth comparing:
//
//     ./cracking 3                 profile theta3 with both variables observed
//     ./cracking 3 0 8 17 1        profile theta3 with only y1 observed
//
// In the second, theta1 and theta3 enter the y1 equation only through their
// sum, so they are structurally non-identifiable: the profile is flat, and
// theta1 moves to hold theta1 + theta3 fixed. The Wald interval, computed at
// the same solution, is finite, narrow, and includes negative rate constants.
// That contrast is the reason for computing a profile at all.
//
//////////////////////////////////////////////////////////////////////////

#include "psopt.h"

#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <vector>

using namespace std;
using namespace PSOPT;

//////////////////////////////////////////////////////////////////////////
///////////////////  Problem configuration ////////////////////////////////
//////////////////////////////////////////////////////////////////////////

// Number of observed variables: 2 fits y1 and y2, 1 fits y1 alone. Set from
// the command line; the observation function has to see it, so it is global.
static int NOBSERVED = 2;

//////////////////////////////////////////////////////////////////////////
///////////////////  Define the observation function //////////
//////////////////////////////////////////////////////////////////////////

void  observation_function( adouble* observations,
                            adouble* states, adouble* controls,
                            adouble* parameters, adouble& time, int k,
                            adouble* xad, int iphase, Workspace* workspace)
{

      observations[ 0 ] = states[ 0 ];
      if ( NOBSERVED > 1 ) observations[ 1 ] = states[ 1 ];
}


//////////////////////////////////////////////////////////////////////////
///////////////////  Define the DAE's ////////////////////////////////////
//////////////////////////////////////////////////////////////////////////

void dae(adouble* derivatives, adouble* path, adouble* states,
         adouble* controls, adouble* parameters, adouble& time,
         adouble* xad, int iphase, Workspace* workspace)
{

   adouble y1 = states[0];
   adouble y2 = states[1];

   adouble theta1 = parameters[ 0 ];
   adouble theta2 = parameters[ 1 ];
   adouble theta3 = parameters[ 2 ];

   derivatives[0] = -(theta1 + theta3)*y1*y1;
   derivatives[1] =  theta1*y1*y1 - theta2*y2;

}

////////////////////////////////////////////////////////////////////////////
///////////////////  Define the events function ////////////////////////////
////////////////////////////////////////////////////////////////////////////

void events(adouble* e, adouble* initial_states, adouble* final_states,
            adouble* parameters,adouble& t0, adouble& tf, adouble* xad,
            int iphase, Workspace* workspace)
{
   // The published problem (COPS gas-oil cracking) fixes y(0) = (1,0);
   // estimating the initial conditions as well makes this a different problem.
   e[0] = initial_states[0];
   e[1] = initial_states[1];
}


///////////////////////////////////////////////////////////////////////////
///////////////////  Define the phase linkages function ///////////////////
///////////////////////////////////////////////////////////////////////////

void linkages( adouble* linkages, adouble* xad, Workspace* workspace)
{
  // No linkages as this is a single phase problem
}



////////////////////////////////////////////////////////////////////////////
///////////////////  Problem setup, done once  /////////////////////////////
////////////////////////////////////////////////////////////////////////////

// Everything that does not change from one solve to the next. The profile
// below then changes exactly two numbers per point -- the lower and the upper
// bound of the parameter being profiled -- and calls psopt() again.

void setup_cracking(Prob& problem, Alg& algorithm)
{
   MatrixXd y1meas(1,21), y2meas(1,21), tmeas(1,21);
   // Measured values of y1
   y1meas << 	   1.0,0.8105,0.6208,0.5258,0.4345,0.3903,0.3342,0.3034, \
                  0.2735,0.2405,0.2283,0.2071,0.1669,0.153,0.1339,0.1265, \
                  0.12,0.099,0.087,0.077,0.069;
   // Measured values of y2
   y2meas << 	0.0,0.2,0.2886,0.301,0.3215,0.3123,0.2716,0.2551,0.2258, \
                 0.1959,0.1789,0.1457,0.1198,0.0909,0.0719,0.0561,0.046, \
                0.028,0.019,0.014,0.01;
   // Sampling instants
   tmeas  <<  0.0,0.025,0.05,0.075,0.1,0.125,0.15,0.175,0.2,0.225,0.25, \
              0.3,0.35,0.4,0.45,0.5,0.55,0.65,0.75,0.85,0.95;

////////////////////////////////////////////////////////////////////////////
///////////////////  Register problem name  ////////////////////////////////
////////////////////////////////////////////////////////////////////////////

    problem.name          		= "Catalytic cracking of gas oil";
    problem.outfilename                 = "cracking.txt";

////////////////////////////////////////////////////////////////////////////
////////////  Define problem level constants & do level 1 setup ////////////
////////////////////////////////////////////////////////////////////////////

    problem.nphases   			= 1;
    problem.nlinkages                   = 0;

    psopt_level1_setup(problem);


/////////////////////////////////////////////////////////////////////////////
/////////   Define phase related information & do level 2 setup /////////////
/////////////////////////////////////////////////////////////////////////////

    problem.phases(1).nstates   		= 2;
    problem.phases(1).ncontrols 		= 0;
    problem.phases(1).nevents   		= 2;
    problem.phases(1).npath     		= 0;
    problem.phases(1).nparameters        	= 3;
    problem.phases(1).nodes    		    	<< 80;
    problem.phases(1).nobserved   = NOBSERVED;
    problem.phases(1).nsamples    = 21;

    psopt_level2_setup(problem, algorithm);

////////////////////////////////////////////////////////////////////////////
////////////  Enter estimation information                      ////////////
////////////////////////////////////////////////////////////////////////////

    MatrixXd observations(NOBSERVED, 21);

    if ( NOBSERVED > 1 ) observations << y1meas, y2meas;
    else                 observations << y1meas;

    problem.phases(1).observation_nodes      = tmeas;
    problem.phases(1).observations           = observations;
    problem.phases(1).residual_weights       = ones(NOBSERVED,21);


////////////////////////////////////////////////////////////////////////////
///////////////////  Enter problem bounds information //////////////////////
////////////////////////////////////////////////////////////////////////////


    problem.phases(1).bounds.lower.states(0) =  0.0;
    problem.phases(1).bounds.lower.states(1) =  0.0;


    problem.phases(1).bounds.upper.states(0) =  2.0;
    problem.phases(1).bounds.upper.states(1) =  2.0;


    problem.phases(1).bounds.lower.parameters(0) = 0.0;
    problem.phases(1).bounds.lower.parameters(1) = 0.0;
    problem.phases(1).bounds.lower.parameters(2) = 0.0;
    problem.phases(1).bounds.upper.parameters(0) = 20.0;
    problem.phases(1).bounds.upper.parameters(1) = 20.0;
    problem.phases(1).bounds.upper.parameters(2) = 20.0;


    problem.phases(1).bounds.lower.events(0) = 1.0;
    problem.phases(1).bounds.upper.events(0) = 1.0;
    problem.phases(1).bounds.lower.events(1) = 0.0;
    problem.phases(1).bounds.upper.events(1) = 0.0;

    problem.phases(1).bounds.lower.StartTime    = 0.0;
    problem.phases(1).bounds.upper.StartTime    = 0.0;

    problem.phases(1).bounds.lower.EndTime      = 0.95;
    problem.phases(1).bounds.upper.EndTime      = 0.95;

////////////////////////////////////////////////////////////////////////////
///////////////////  Register problem functions  ///////////////////////////
////////////////////////////////////////////////////////////////////////////

    problem.dae 		= &dae;
    problem.events 		= &events;
    problem.linkages		= &linkages;
    problem.observation_function = & observation_function;

////////////////////////////////////////////////////////////////////////////
///////////////////  Define & register initial guess ///////////////////////
////////////////////////////////////////////////////////////////////////////

    MatrixXd state_guess(2, 40);

    state_guess.row(0) =  linspace(1.0,0.069, 40);
    state_guess.row(1) =  linspace(0.30,0.01,  40);


    problem.phases(1).guess.states         = state_guess;
    problem.phases(1).guess.time           = linspace(0.0, 0.95, 40);
    problem.phases(1).guess.parameters     = zeros(3,1);


////////////////////////////////////////////////////////////////////////////
///////////////////  Enter algorithm options  //////////////////////////////
////////////////////////////////////////////////////////////////////////////

    algorithm.nlp_method                  = "IPOPT";
    algorithm.scaling                     = "automatic";
    algorithm.derivatives                 = "automatic";
    algorithm.collocation_method          = "Hermite-Simpson";
    algorithm.parameter_statistics        = "yes";
    algorithm.nlp_iter_max                = 1000;
    algorithm.nlp_tolerance               = 1.e-6;
}


////////////////////////////////////////////////////////////////////////////
///////////////////  What one solve reports  ///////////////////////////////
////////////////////////////////////////////////////////////////////////////

struct Fit {
   double   J;              // optimal sum of squared residuals
   double   theta[3];       // estimated parameters
   double   lo[3], hi[3];   // 95 per cent Wald limits, when available
   double   sigma_hat;      // estimated residual standard deviation
   long     dof;            // number of fitted quantities n_f
   long     ns;             // number of scalar observations N_s
   bool     stats_ok;
   bool     solved;
};

Fit read_fit(Sol& solution, int status)
{
    Fit fit;
    fit.solved = ( status == 0 && solution.error_flag == 0 );
    fit.J      = solution.cost;

    MatrixXd p = solution.get_parameters_in_phase(1);
    for (int i = 0; i < 3; i++) fit.theta[i] = p(i,0);

    fit.sigma_hat = solution.sigma_hat;
    fit.dof       = solution.parameter_dof;
    fit.ns        = solution.n_observations;
    fit.stats_ok  = solution.parameter_statistics_ok
                    && solution.parameter_confidence_low.rows() >= 3;
    for (int i = 0; i < 3; i++) {
        fit.lo[i] = fit.stats_ok ? solution.parameter_confidence_low(i,0)  : 0.0;
        fit.hi[i] = fit.stats_ok ? solution.parameter_confidence_high(i,0) : 0.0;
    }
    return fit;
}

// Fix a parameter, or release it. Fixing needs no special support: equal bounds
// remove it from the estimation and leave everything else untouched.
void fix_parameter(Prob& problem, int k, double value)
{
    problem.phases(1).bounds.lower.parameters(k) = value;
    problem.phases(1).bounds.upper.parameters(k) = value;
}

void release_parameter(Prob& problem, int k)
{
    problem.phases(1).bounds.lower.parameters(k) =  0.0;
    problem.phases(1).bounds.upper.parameters(k) = 20.0;
}


////////////////////////////////////////////////////////////////////////////
///////////////////  Define the main routine ///////////////////////////////
////////////////////////////////////////////////////////////////////////////

int main(int argc, char* argv[])
{
    int    index = (argc > 1) ? atoi(argv[1]) : 0;      // 0 means: no profile
    double lo    = (argc > 2) ? atof(argv[2]) : 0.0;
    double hi    = (argc > 3) ? atof(argv[3]) : 2.2;
    int    ngrid = (argc > 4) ? atoi(argv[4]) : 45;
    if (argc > 5) NOBSERVED = atoi(argv[5]);

    if (index < 0 || index > 3 || ngrid < 2 || NOBSERVED < 1 || NOBSERVED > 2) {
        printf("usage: %s <index 1..3> [lo] [hi] [ngrid] [nobserved 1|2]\n",
               argv[0]);
        return 1;
    }

    const bool profiling = (index > 0);

    Alg  algorithm;
    Prob problem;
    setup_cracking(problem, algorithm);
    algorithm.print_level = profiling ? 0 : 1;

    // The unconstrained fit, with all three parameters estimated. Its solution
    // is kept in its own Sol for the whole run, because it is the starting
    // point of both halves of the profile below.
    Sol nominal;
    Fit nom = read_fit(nominal, psopt(nominal, problem, algorithm));

    printf("\n Estimated parameters\n");
    for (int i = 0; i < 3; i++) printf("   theta%d = %12.6f\n", i+1, nom.theta[i]);

    if (!profiling) {
        // Exactly the behaviour this example has always had.
        MatrixXd x = nominal.get_states_in_phase(1);
        MatrixXd t = nominal.get_time_in_phase(1);
        Save(x,"x.dat");
        Save(t,"t.dat");
        plot(t,x,"Catalytic cracking of gas oil", "time (s)", "states", "y1 y2");
        plot(t,x,"Catalytic cracking of gas oil", "time (s)", "states", "y1 y2",
                 "pdf", "cracking_states.pdf");
        return 0;
    }

    ////////////////////////////////////////////////////////////////////////
    ///////////////////  Profile likelihood  ///////////////////////////////
    ////////////////////////////////////////////////////////////////////////

    const int    k   = index - 1;
    const long   nu  = nom.ns - nom.dof;                  // residual degrees of freedom
    const double tq  = inverse_twotailed_t_cdf(0.95, (int) nu);
    const double thr = 1.0 + tq*tq/(double) nu;           // threshold on J/J*

    printf("\n Profile likelihood for theta%d\n", index);
    printf("   observed variables : %d\n", NOBSERVED);
    printf("   J*                 : %.10e\n", nom.J);
    printf("   sigma_hat          : %.6e\n", nom.sigma_hat);
    printf("   Ns, nf, Ns-nf      : %ld, %ld, %ld\n", nom.ns, nom.dof, nu);
    printf("   t^0.975_%ld         : %.4f\n", nu, tq);
    printf("   threshold J/J*     : %.6f\n", thr);
    if (nom.stats_ok)
        printf("   Wald 95%% interval  : [%.6f, %.6f]\n", nom.lo[k], nom.hi[k]);
    else
        printf("   Wald 95%% interval  : not available\n");

    MatrixXd grid(1, ngrid), ratio(1, ngrid);
    vector<Fit> fits(ngrid);
    for (int j = 0; j < ngrid; j++) grid(0,j) = lo + (hi - lo)*j/(double)(ngrid - 1);

    // The grid point nearest the unconstrained estimate, from which the profile
    // is traced outwards in both directions.
    int jstart = 0;
    for (int j = 1; j < ngrid; j++)
        if (fabs(grid(0,j) - nom.theta[k]) < fabs(grid(0,jstart) - nom.theta[k]))
            jstart = j;

    // The profile is traced by CONTINUATION: each solve starts from the
    // solution of the one before it, and each half starts from the
    // unconstrained fit. Solving every point from the same cold guess instead
    // does not fail -- it converges, reports success, and occasionally returns
    // a different local minimum, which appears on a plot of the profile as
    // structure in the likelihood that is not there. Two of seventeen points
    // did exactly that on the y1-only run below, on a profile that is
    // otherwise flat to two parts in 100,000.
    Sol solution;

    set_guess_from_solution(problem, nominal);
    for (int j = jstart; j < ngrid; j++) {
        fix_parameter(problem, k, grid(0,j));
        fits[j] = read_fit(solution, psopt(solution, problem, algorithm));
        set_guess_from_solution(problem, solution);
    }

    set_guess_from_solution(problem, nominal);
    for (int j = jstart - 1; j >= 0; j--) {
        fix_parameter(problem, k, grid(0,j));
        fits[j] = read_fit(solution, psopt(solution, problem, algorithm));
        set_guess_from_solution(problem, solution);
    }
    release_parameter(problem, k);

    char fname[64];
    if (NOBSERVED == 2)
        snprintf(fname, sizeof fname, "cracking_profile_theta%d.dat", index);
    else
        snprintf(fname, sizeof fname, "cracking_profile_theta%d_y%d.dat",
                 index, NOBSERVED);
    FILE* fh = fopen(fname, "w");
    fprintf(fh, "# profile likelihood, catalytic cracking of gas oil\n");
    fprintf(fh, "# nobserved = %d, J* = %.10g, threshold J/J* = %.10g\n",
            NOBSERVED, nom.J, thr);
    fprintf(fh, "# theta%d  J  J/J*  theta1  theta2  theta3\n", index);

    printf("\n%12s %18s %12s %10s %10s %10s\n",
           "theta", "J", "J/J*", "theta1", "theta2", "theta3");
    for (int j = 0; j < ngrid; j++) {
        const Fit& f = fits[j];
        ratio(0,j) = f.J/nom.J;
        printf("%12.6f %18.10e %12.6f %10.5f %10.5f %10.5f%s\n",
               grid(0,j), f.J, f.J/nom.J, f.theta[0], f.theta[1], f.theta[2],
               f.solved ? "" : "   NOT SOLVED");
        fprintf(fh, "%.10g %.10g %.10g %.10g %.10g %.10g\n",
                grid(0,j), f.J, f.J/nom.J, f.theta[0], f.theta[1], f.theta[2]);
    }
    fclose(fh);

    // Where the profile crosses the threshold, by linear interpolation. A side
    // that never crosses is reported as open: that is the signal that the data
    // do not bound the parameter on that side, and it is exactly the case a
    // Wald interval cannot express.
    int jmin = 0;
    for (int j = 1; j < ngrid; j++) if (ratio(0,j) < ratio(0,jmin)) jmin = j;
    double plo = 0.0, phi = 0.0; bool haslo = false, hashi = false;
    for (int j = 0; j < jmin; j++)
        if ((ratio(0,j) - thr)*(ratio(0,j+1) - thr) < 0.0) {
            double s = (thr - ratio(0,j))/(ratio(0,j+1) - ratio(0,j));
            plo = grid(0,j) + s*(grid(0,j+1) - grid(0,j)); haslo = true;
        }
    for (int j = jmin; j < ngrid-1; j++)
        if ((ratio(0,j) - thr)*(ratio(0,j+1) - thr) < 0.0) {
            double s = (thr - ratio(0,j))/(ratio(0,j+1) - ratio(0,j));
            phi = grid(0,j) + s*(grid(0,j+1) - grid(0,j)); hashi = true;
        }

    printf("\n Profile 95%% interval for theta%d: ", index);
    if (haslo) printf("[%.4f, ", plo); else printf("(open below %.4f, ", lo);
    if (hashi) printf("%.4f]\n", phi);  else printf("open above %.4f)\n", hi);
    if (!haslo || !hashi)
        printf("   The profile does not cross the threshold on at least one side\n"
               "   over the range scanned: the data do not bound theta%d there,\n"
               "   whatever the Wald interval above may say.\n", index);
    printf(" wrote %s\n", fname);

    return 0;
}

////////////////////////////////////////////////////////////////////////////
///////////////////////      END OF FILE     ///////////////////////////////
////////////////////////////////////////////////////////////////////////////
