//////////////////////////////////////////////////////////////////////////////
// test_simple_constraints.cpp
//
// Maximum-Principle validation of three "simple / nonstandard constraint"
// double-integrator problems from A. Locatelli, "Optimal Control of a Double
// Integrator", Springer 2017 (Chaps. 3-6). All share xdot1=x2, xdot2=u and
// J = INT_0^1 u^2/2 dt with tf = 1 fixed; they differ in the final-state data.
//////////////////////////////////////////////////////////////////////////////

#include "gtest/gtest.h"
#include <psopt.h>
#include <cmath>

// ---- shared dynamics / cost ------------------------------------------------
namespace di_shared {
adouble endpoint_cost(adouble*, adouble*, adouble*, adouble&, adouble& tf, adouble*, int, Workspace*) { return 0.0; }
adouble integrand_cost(adouble*, adouble* c, adouble*, adouble&, adouble*, int, Workspace*) { adouble u=c[0]; return 0.5*u*u; }
void dae(adouble* d, adouble*, adouble* s, adouble* c, adouble*, adouble&, adouble*, int, Workspace*) { d[0]=s[1]; d[1]=c[0]; }
void linkages(adouble*, adouble*, Workspace*) {}
}

// Build + solve a min-energy double integrator with `ne` events set by `ev`.
static void solve_di(Sol& solution, void (*ev)(adouble*,adouble*,adouble*,adouble*,adouble&,adouble&,adouble*,int,Workspace*),
                     int ne, const MatrixXd& elo, const MatrixXd& eup)
{
    using namespace di_shared;
    const int N = 40;
    Alg algorithm; Prob problem;
    problem.name = "min-energy double integrator"; problem.outfilename = "se.txt";
    problem.nphases = 1; problem.nlinkages = 0;
    psopt_level1_setup(problem);
    problem.phases(1).nstates = 2; problem.phases(1).ncontrols = 1;
    problem.phases(1).nevents = ne; problem.phases(1).npath = 0; problem.phases(1).nodes << N;
    psopt_level2_setup(problem, algorithm);
    problem.phases(1).bounds.lower.states << -5.0,-5.0;  problem.phases(1).bounds.upper.states << 5.0,5.0;
    problem.phases(1).bounds.lower.controls(0) = -50.0;  problem.phases(1).bounds.upper.controls(0) = 50.0;
    problem.phases(1).bounds.lower.events = elo;         problem.phases(1).bounds.upper.events = eup;
    problem.phases(1).bounds.lower.StartTime = 0.0; problem.phases(1).bounds.upper.StartTime = 0.0;
    problem.phases(1).bounds.lower.EndTime   = 1.0; problem.phases(1).bounds.upper.EndTime   = 1.0;
    problem.integrand_cost = &integrand_cost; problem.endpoint_cost = &endpoint_cost;
    problem.dae = &dae; problem.events = ev; problem.linkages = &linkages;
    problem.phases(1).guess.states = zeros(2,N);
    problem.phases(1).guess.states.row(0) = linspace(0.0,1.0,N);
    problem.phases(1).guess.controls = zeros(1,N);
    problem.phases(1).guess.time = linspace(0.0,1.0,N);
    algorithm.nlp_method="IPOPT"; algorithm.scaling="automatic"; algorithm.derivatives="automatic";
    algorithm.collocation_method="Legendre"; algorithm.nlp_iter_max=1000; algorithm.nlp_tolerance=1.e-6;
    psopt(solution, problem, algorithm);
}

// ---- (1) fixed endpoints (Ch 3):  (0,0)->(1,0),  u=6-12t,  J*=6 ------------
namespace se1 { void events(adouble* e,adouble* i,adouble* f,adouble*,adouble&,adouble&,adouble*,int,Workspace*)
{ e[0]=i[0]; e[1]=i[1]; e[2]=f[0]; e[3]=f[1]; } }

TEST(SimpleConstraints, MinEnergyFixedEndpoints)
{
    Sol solution; MatrixXd lo(4,1), up(4,1); lo<<0,0,1,0; up<<0,0,1,0;
    solve_di(solution, &se1::events, 4, lo, up);
    ASSERT_EQ(solution.error_flag, false);
    DMatrix x=solution.get_states_in_phase(1), u=solution.get_controls_in_phase(1), t=solution.get_time_in_phase(1);
    int n=(int)t.cols();
    EXPECT_NEAR(solution.get_cost(), 6.0, 5.0e-3) << "J* = 6";
    EXPECT_NEAR(x(0,n-1), 1.0, 1.0e-4);  EXPECT_NEAR(x(1,n-1), 0.0, 1.0e-4);
    EXPECT_NEAR(u(0,0),  6.0, 1.0e-1) << "u(0) = 6";
    EXPECT_NEAR(u(0,n-1),-6.0, 1.0e-1) << "u(tf) = -6";
    double vmax=0; for(int k=0;k<n;++k) vmax=std::max(vmax,x(1,k));
    EXPECT_NEAR(vmax, 1.5, 2.0e-2) << "peak speed x2(1/2) = 3/2";
}

// ---- (2) free final velocity (Ch 3.3):  x1(1)=1, x2 free -> u(tf)=0 --------
namespace se2 { void events(adouble* e,adouble* i,adouble* f,adouble*,adouble&,adouble&,adouble*,int,Workspace*)
{ e[0]=i[0]; e[1]=i[1]; e[2]=f[0]; } }

TEST(SimpleConstraints, FreeFinalVelocityTransversality)
{
    Sol solution; MatrixXd lo(3,1), up(3,1); lo<<0,0,1; up<<0,0,1;
    solve_di(solution, &se2::events, 3, lo, up);
    ASSERT_EQ(solution.error_flag, false);
    DMatrix x=solution.get_states_in_phase(1), u=solution.get_controls_in_phase(1), t=solution.get_time_in_phase(1);
    int n=(int)t.cols();
    EXPECT_NEAR(solution.get_cost(), 1.5, 5.0e-3) << "J* = 3/2";
    EXPECT_NEAR(x(0,n-1), 1.0, 1.0e-4) << "x1(tf) = 1";
    EXPECT_NEAR(x(1,n-1), 1.5, 5.0e-3) << "free final velocity x2(tf) = 3/2";
    EXPECT_NEAR(u(0,n-1), 0.0, 5.0e-2) << "TRANSVERSALITY: l2(tf)=0 => u(tf)=0";
}

// ---- (3) target line (Ch 6):  x1(1)-2 x2(1)=1  (orthogonality) -------------
namespace se3 { void events(adouble* e,adouble* i,adouble* f,adouble*,adouble&,adouble&,adouble*,int,Workspace*)
{ e[0]=i[0]; e[1]=i[1]; e[2]=f[0]-2.0*f[1]; } }

TEST(SimpleConstraints, TargetLineOrthogonality)
{
    Sol solution; MatrixXd lo(3,1), up(3,1); lo<<0,0,1; up<<0,0,1;
    solve_di(solution, &se3::events, 3, lo, up);
    ASSERT_EQ(solution.error_flag, false);
    DMatrix x=solution.get_states_in_phase(1), t=solution.get_time_in_phase(1);
    int n=(int)t.cols();
    EXPECT_NEAR(solution.get_cost(), 3.0/14.0, 5.0e-3) << "J* = 3/14";
    EXPECT_NEAR(x(0,n-1) - 2.0*x(1,n-1), 1.0, 1.0e-4) << "final state on the line x1-2x2=1";
    EXPECT_NEAR(x(0,n-1), -2.0/7.0, 5.0e-3) << "x1(tf) = -2/7";
    EXPECT_NEAR(x(1,n-1), -9.0/14.0, 5.0e-3) << "x2(tf) = -9/14";
}
