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
//   tf* = 2*sqrt(|x0|) = 2,   u = -1 on [0,1] then +1 on [1,2],  switch at tf*/2,
//   peak speed |v| = sqrt(|x0|) = 1 at the switch.
//
// The single switch is modelled EXACTLY with TWO phases joined at the switch: each
// arc has constant control so the velocity is piecewise-linear and the peak |v|=1
// at the phase boundary is captured exactly (a single-phase mesh either rounds the
// bang-bang corner or, if refined, chatters). auto_link enforces state/time
// continuity across the switch; the switch time itself is a free variable.
//////////////////////////////////////////////////////////////////////////////

#include "gtest/gtest.h"
#include <psopt.h>
#include <cmath>

namespace mintime_di {

adouble endpoint_cost(adouble* i, adouble* f, adouble* p, adouble& t0,
                      adouble& tf, adouble* xad, int iphase, Workspace* w)
{ return (iphase == 2) ? tf : (adouble) 0.0; }    // minimise the final time

adouble integrand_cost(adouble* s, adouble* c, adouble* p, adouble& t,
                       adouble* xad, int iphase, Workspace* w)
{ return 0.0; }

void dae(adouble* d, adouble* path, adouble* states, adouble* controls,
         adouble* p, adouble& time, adouble* xad, int iphase, Workspace* w)
{ d[0] = states[1];  d[1] = controls[0]; }        // xdot=v, vdot=u (both phases)

void events(adouble* e, adouble* i, adouble* f, adouble* p, adouble& t0,
            adouble& tf, adouble* xad, int iphase, Workspace* w)
{
    if (iphase == 1) { e[0]=i[0]; e[1]=i[1]; }     // phase 1 start: (x,v)=(1,0)
    else             { e[0]=f[0]; e[1]=f[1]; }     // phase 2 end:   (x,v)=(0,0)
}

void linkages(adouble* linkages, adouble* xad, Workspace* w)
{ int index = 0; auto_link(linkages, &index, xad, 1, 2, w); }  // continuity at the switch

} // namespace mintime_di

TEST(MaximumPrinciple, MinTimeDoubleIntegrator)
{
    using namespace mintime_di;
    const int N = 20;                              // nodes per arc

    Alg algorithm; Sol solution; Prob problem;
    problem.name = "Minimum-time double integrator (MP validation)";
    problem.outfilename = "mp_mintime.txt";
    problem.nphases = 2; problem.nlinkages = 3;    // auto_link: 2 states + time
    psopt_level1_setup(problem);

    for (int ph = 1; ph <= 2; ++ph) {
        problem.phases(ph).nstates   = 2;
        problem.phases(ph).ncontrols = 1;
        problem.phases(ph).nevents   = 2;
        problem.phases(ph).npath     = 0;
        problem.phases(ph).nodes     << N;
    }
    psopt_level2_setup(problem, algorithm);

    // Phase 1 is the decelerating arc (u -> -1), phase 2 the accelerating arc
    // (u -> +1); the sign is fixed to structure the two arcs, the magnitude is free.
    problem.phases(1).bounds.lower.states << -5.0, -5.0;  problem.phases(1).bounds.upper.states << 5.0, 5.0;
    problem.phases(2).bounds.lower.states << -5.0, -5.0;  problem.phases(2).bounds.upper.states << 5.0, 5.0;
    problem.phases(1).bounds.lower.controls(0) = -1.0;    problem.phases(1).bounds.upper.controls(0) = 0.0;
    problem.phases(2).bounds.lower.controls(0) =  0.0;    problem.phases(2).bounds.upper.controls(0) = 1.0;

    problem.phases(1).bounds.lower.events << 1.0, 0.0;    problem.phases(1).bounds.upper.events << 1.0, 0.0;
    problem.phases(2).bounds.lower.events << 0.0, 0.0;    problem.phases(2).bounds.upper.events << 0.0, 0.0;

    problem.phases(1).bounds.lower.StartTime = 0.0; problem.phases(1).bounds.upper.StartTime = 0.0;
    problem.phases(1).bounds.lower.EndTime   = 0.2; problem.phases(1).bounds.upper.EndTime   = 4.0; // switch time (free)
    problem.phases(2).bounds.lower.StartTime = 0.2; problem.phases(2).bounds.upper.StartTime = 4.0;
    problem.phases(2).bounds.lower.EndTime   = 0.5; problem.phases(2).bounds.upper.EndTime   = 6.0; // final time (free)

    problem.integrand_cost = &integrand_cost;
    problem.endpoint_cost  = &endpoint_cost;
    problem.dae            = &dae;
    problem.events         = &events;
    problem.linkages       = &linkages;

    // Analytical bang-bang warm start on each arc.
    MatrixXd x1(2,N), u1(1,N), t1(1,N), x2(2,N), u2(1,N), t2(1,N);
    for (int k = 0; k < N; ++k) {
        double a = (double) k / (N - 1);
        double ta = a;             t1(0,k)=ta;  u1(0,k)=-1.0; x1(0,k)=1.0-0.5*ta*ta;      x1(1,k)=-ta;
        double tb = 1.0 + a;       t2(0,k)=tb;  u2(0,k)= 1.0; x2(0,k)=0.5-a+0.5*a*a;      x2(1,k)=-1.0+a;
    }
    problem.phases(1).guess.states = x1; problem.phases(1).guess.controls = u1; problem.phases(1).guess.time = t1;
    problem.phases(2).guess.states = x2; problem.phases(2).guess.controls = u2; problem.phases(2).guess.time = t2;

    algorithm.nlp_method         = "IPOPT";
    algorithm.scaling            = "automatic";
    algorithm.derivatives        = "automatic";
    algorithm.collocation_method = "trapezoidal";  // exact for the piecewise-linear arcs
    algorithm.nlp_iter_max       = 1000;
    algorithm.nlp_tolerance      = 1.e-6;

    psopt(solution, problem, algorithm);

    ASSERT_EQ(solution.error_flag, false) << "PSOPT failed to solve the min-time problem";

    DMatrix x1s = solution.get_states_in_phase(1),  x2s = solution.get_states_in_phase(2);
    DMatrix u1s = solution.get_controls_in_phase(1), u2s = solution.get_controls_in_phase(2);
    DMatrix t1s = solution.get_time_in_phase(1),     t2s = solution.get_time_in_phase(2);
    const int n1 = (int) t1s.cols(), n2 = (int) t2s.cols();
    const double t_switch = t1s(0, n1-1);
    const double tf       = t2s(0, n2-1);

    // (1) Minimum time and switch time match the MP closed form.
    EXPECT_NEAR(tf,       2.0, 1.0e-2) << "min time = 2*sqrt(|x0|)";
    EXPECT_NEAR(t_switch, 1.0, 1.0e-2) << "switch at tf*/2";

    // (2) Boundary conditions reached (phase-2 end).
    EXPECT_NEAR(x2s(0, n2-1), 0.0, 1.0e-4);
    EXPECT_NEAR(x2s(1, n2-1), 0.0, 1.0e-4);

    // (3) Bang-bang: each arc saturates the bound, opposite signs (one switch).
    EXPECT_NEAR(u1s(0, n1/2), -1.0, 1.0e-3) << "arc 1 rides u = -1";
    EXPECT_NEAR(u2s(0, n2/2), +1.0, 1.0e-3) << "arc 2 rides u = +1";

    // (4) Peak speed at the switch equals sqrt(|x0|) = 1 -- exact at the phase join.
    EXPECT_NEAR(x1s(1, n1-1), -1.0, 1.0e-3) << "peak speed at the switch = sqrt(|x0|) = 1";
    EXPECT_NEAR(x2s(1, 0),    -1.0, 1.0e-3) << "velocity continuous across the switch";
}
