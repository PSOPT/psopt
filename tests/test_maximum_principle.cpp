//////////////////////////////////////////////////////////////////////////////
// test_maximum_principle.cpp
//
// Validates PSOPT's numerical solution of the MINIMUM-TIME double integrator
// against the analytical Maximum-Principle result from
//   A. Locatelli, "Optimal Control of a Double Integrator: A Primer on Maximum
//    Principle", Springer 2017 (Ch. 2.2.4 / Ch. 7).
//
// Problem:  xdot = v,  vdot = u,  |u| <= 1,  (x,v): (1,0) -> (0,0), minimise tf.
// Maximum Principle:  H = lambda0 + l1 v + l2 u,  l1dot=0, l2dot=-l1  => the
// switching function l2(t) is LINEAR, so the extremal control is BANG-BANG with
// AT MOST n-1 = 1 switch (Locatelli Thm 2.4). The closed-form optimum is
//   tf* = 2*sqrt(|x0|) = 2,   u = -1 on [0,1] then +1 on [1,2],  switch at tf*/2.
//////////////////////////////////////////////////////////////////////////////

#include "gtest/gtest.h"
#include <psopt.h>
#include <cmath>

namespace mintime_di {

adouble endpoint_cost(adouble* i, adouble* f, adouble* p, adouble& t0,
                      adouble& tf, adouble* xad, int iphase, Workspace* w)
{ return tf; }                                   // minimum time

adouble integrand_cost(adouble* s, adouble* c, adouble* p, adouble& t,
                       adouble* xad, int iphase, Workspace* w)
{ return 0.0; }

void dae(adouble* d, adouble* path, adouble* states, adouble* controls,
         adouble* p, adouble& time, adouble* xad, int iphase, Workspace* w)
{ d[0] = states[1];  d[1] = controls[0]; }       // xdot=v, vdot=u

void events(adouble* e, adouble* i, adouble* f, adouble* p, adouble& t0,
            adouble& tf, adouble* xad, int iphase, Workspace* w)
{ e[0]=i[0]; e[1]=i[1]; e[2]=f[0]; e[3]=f[1]; }  // (1,0) -> (0,0)

void linkages(adouble* l, adouble* xad, Workspace* w) {}

} // namespace mintime_di

TEST(MaximumPrinciple, MinTimeDoubleIntegrator)
{
    using namespace mintime_di;
    const int N = 40;

    Alg algorithm; Sol solution; Prob problem;
    problem.name = "Minimum-time double integrator (MP validation)";
    problem.outfilename = "mp_mintime.txt";
    problem.nphases = 1; problem.nlinkages = 0;
    psopt_level1_setup(problem);

    problem.phases(1).nstates   = 2;
    problem.phases(1).ncontrols = 1;
    problem.phases(1).nevents   = 4;
    problem.phases(1).npath     = 0;
    problem.phases(1).nodes     << N;
    psopt_level2_setup(problem, algorithm);

    problem.phases(1).bounds.lower.states << -5.0, -5.0;
    problem.phases(1).bounds.upper.states <<  5.0,  5.0;
    problem.phases(1).bounds.lower.controls(0) = -1.0;   // |u| <= 1
    problem.phases(1).bounds.upper.controls(0) =  1.0;
    problem.phases(1).bounds.lower.events << 1.0, 0.0, 0.0, 0.0;
    problem.phases(1).bounds.upper.events << 1.0, 0.0, 0.0, 0.0;
    problem.phases(1).bounds.lower.StartTime = 0.0; problem.phases(1).bounds.upper.StartTime = 0.0;
    problem.phases(1).bounds.lower.EndTime   = 0.5; problem.phases(1).bounds.upper.EndTime   = 5.0;

    problem.integrand_cost = &integrand_cost;
    problem.endpoint_cost  = &endpoint_cost;
    problem.dae            = &dae;
    problem.events         = &events;
    problem.linkages       = &linkages;

    // Warm start with the analytical bang-bang solution (helps the nonsmooth solve).
    MatrixXd xg(2,N), ug(1,N), tg(1,N);
    for (int k = 0; k < N; ++k) {
        double t = 2.0 * k / (N - 1);
        tg(0,k) = t;
        ug(0,k) = (t < 1.0) ? -1.0 : 1.0;
        if (t < 1.0) { xg(0,k) = 1.0 - 0.5*t*t;      xg(1,k) = -t; }
        else         { double s=t-1.0; xg(0,k)=0.5 - s + 0.5*s*s; xg(1,k) = -1.0 + s; }
    }
    problem.phases(1).guess.states   = xg;
    problem.phases(1).guess.controls = ug;
    problem.phases(1).guess.time     = tg;

    algorithm.nlp_method         = "IPOPT";
    algorithm.scaling            = "automatic";
    algorithm.derivatives        = "automatic";
    algorithm.collocation_method = "Legendre";
    algorithm.nlp_iter_max       = 1000;
    algorithm.nlp_tolerance      = 1.e-6;

    psopt(solution, problem, algorithm);

    ASSERT_EQ(solution.error_flag, false) << "PSOPT failed to solve the min-time problem";

    DMatrix x = solution.get_states_in_phase(1);
    DMatrix u = solution.get_controls_in_phase(1);
    DMatrix t = solution.get_time_in_phase(1);
    const int n = (int) t.cols();
    const double tf = t(0, n-1);

    // (1) Minimum time matches the MP closed form  tf* = 2*sqrt(|x0|) = 2.
    EXPECT_NEAR(tf, 2.0, 3.0e-2) << "min time should equal 2*sqrt(|x0|)";

    // (2) Boundary conditions reached.
    EXPECT_NEAR(x(0, n-1), 0.0, 1.0e-3);
    EXPECT_NEAR(x(1, n-1), 0.0, 1.0e-3);

    // (3) Control is BANG-BANG: |u|~1 and exactly ONE sign switch (<= n-1 = 1).
    int switches = 0; double prev = u(0,0);
    for (int k = 1; k < n; ++k) {
        if (prev * u(0,k) < 0.0) switches++;      // sign change
        if (std::fabs(u(0,k)) > 1e-3) prev = u(0,k);
    }
    EXPECT_LE(switches, 1) << "bang-bang extremal must switch at most n-1 = 1 time";
    EXPECT_NEAR(std::fabs(u(0,0)),     1.0, 5.0e-2) << "control saturates the bound";
    EXPECT_NEAR(std::fabs(u(0,n-1)),   1.0, 5.0e-2);
    EXPECT_LT(u(0,0), 0.0)   << "starts with u = -1 (decelerate toward origin)";
    EXPECT_GT(u(0,n-1), 0.0) << "ends with u = +1";

    // (4) Peak speed at the switch equals sqrt(|x0|) = 1 (v is a tent 0 -> -1 -> 0,
    //     peaking at t = tf/2 = 1). The sampled minimum sits slightly above -1
    //     because the Legendre-Gauss-Lobatto nodes cluster at the endpoints and
    //     none lands exactly on the corner, so allow for that node spacing.
    double vmin = 0.0; for (int k = 0; k < n; ++k) vmin = std::min(vmin, x(1,k));
    EXPECT_GE(vmin, -1.02) << "velocity must not exceed the analytical peak |v|=1";
    EXPECT_LE(vmin, -0.90) << "velocity must approach the analytical peak |v|=1";
}
