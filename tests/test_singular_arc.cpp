//////////////////////////////////////////////////////////////////////////////
// test_singular_arc.cpp
//
// Validates PSOPT on a SINGULAR-ARC double integrator against the analytical
// Maximum-Principle solution from
//   A. Locatelli, "Optimal Control of a Double Integrator", Springer 2017,
//   Problem 11.2-1a.
//
//   xdot1 = x2, xdot2 = u, |u| <= 1,  x(0)=(2,0), x(tf)=(5,0), tf FREE,
//   minimise  J = INT_0^tf (1 + x2^2)/2 dt.
//
// H = (1+x2^2)/2 + l1 x2 + l2 u is LINEAR in u, so u = -sign(l2) is bang-bang
// unless the switching function l2(t) vanishes on an interval -> a SINGULAR arc,
// where dH/du = 0 identically and u follows from d^2/dt^2(dH/du)=0 => u = 0,
// which holds x2 = 1. Closed-form optimum (Locatelli 11.2-1a):
//   u = +1 on [0,1],  u = 0 (singular, x2 = 1) on [1,3],  u = -1 on [3,4];
//   tf* = 4,   J* = 10/3.
//////////////////////////////////////////////////////////////////////////////

#include "gtest/gtest.h"
#include <psopt.h>
#include <cmath>

namespace singular_di {

adouble endpoint_cost(adouble* i, adouble* f, adouble* p, adouble& t0,
                      adouble& tf, adouble* xad, int iphase, Workspace* w)
{ return 0.0; }

adouble integrand_cost(adouble* states, adouble* controls, adouble* p,
                       adouble& time, adouble* xad, int iphase, Workspace* w)
{ adouble x2 = states[1]; return 0.5*(1.0 + x2*x2); }

void dae(adouble* d, adouble* path, adouble* states, adouble* controls,
         adouble* p, adouble& time, adouble* xad, int iphase, Workspace* w)
{ d[0] = states[1];  d[1] = controls[0]; }

void events(adouble* e, adouble* i, adouble* f, adouble* p, adouble& t0,
            adouble& tf, adouble* xad, int iphase, Workspace* w)
{ e[0]=i[0]; e[1]=i[1]; e[2]=f[0]; e[3]=f[1]; }

void linkages(adouble* l, adouble* xad, Workspace* w) {}

} // namespace singular_di

TEST(MaximumPrinciple, SingularArcDoubleIntegrator)
{
    using namespace singular_di;
    const int N = 60;

    Alg algorithm; Sol solution; Prob problem;
    problem.name = "Singular-arc double integrator (MP validation)";
    problem.outfilename = "singular_di.txt";
    problem.nphases = 1; problem.nlinkages = 0;
    psopt_level1_setup(problem);

    problem.phases(1).nstates   = 2;
    problem.phases(1).ncontrols = 1;
    problem.phases(1).nevents   = 4;
    problem.phases(1).npath     = 0;
    problem.phases(1).nodes     << N;
    psopt_level2_setup(problem, algorithm);

    problem.phases(1).bounds.lower.states << -10.0, -5.0;
    problem.phases(1).bounds.upper.states <<  10.0,  5.0;
    problem.phases(1).bounds.lower.controls(0) = -1.0;
    problem.phases(1).bounds.upper.controls(0) =  1.0;
    problem.phases(1).bounds.lower.events << 2.0, 0.0, 5.0, 0.0;
    problem.phases(1).bounds.upper.events << 2.0, 0.0, 5.0, 0.0;
    problem.phases(1).bounds.lower.StartTime = 0.0; problem.phases(1).bounds.upper.StartTime = 0.0;
    problem.phases(1).bounds.lower.EndTime   = 2.0; problem.phases(1).bounds.upper.EndTime   = 8.0;

    problem.integrand_cost = &integrand_cost;
    problem.endpoint_cost  = &endpoint_cost;
    problem.dae            = &dae;
    problem.events         = &events;
    problem.linkages       = &linkages;

    // Analytical bang-singular-bang warm start (tf=4).
    MatrixXd xg(2,N), ug(1,N), tg(1,N);
    for (int k = 0; k < N; ++k) {
        double t = 4.0 * k / (N - 1); tg(0,k) = t;
        if (t < 1.0)      { ug(0,k)= 1.0; xg(0,k)=2.0+0.5*t*t;  xg(1,k)=t; }
        else if (t < 3.0) { ug(0,k)= 0.0; xg(0,k)=2.5+(t-1.0);  xg(1,k)=1.0; }
        else              { double s=t-3.0; ug(0,k)=-1.0; xg(0,k)=4.5+s-0.5*s*s; xg(1,k)=1.0-s; }
    }
    problem.phases(1).guess.states = xg; problem.phases(1).guess.controls = ug; problem.phases(1).guess.time = tg;

    algorithm.nlp_method         = "IPOPT";
    algorithm.scaling            = "automatic";
    algorithm.derivatives        = "automatic";
    algorithm.collocation_method = "Legendre";
    algorithm.nlp_iter_max       = 1000;
    algorithm.nlp_tolerance      = 1.e-6;

    psopt(solution, problem, algorithm);

    ASSERT_EQ(solution.error_flag, false) << "PSOPT failed to solve the singular-arc problem";

    DMatrix x = solution.get_states_in_phase(1);
    DMatrix u = solution.get_controls_in_phase(1);
    DMatrix t = solution.get_time_in_phase(1);
    const int n = (int) t.cols();
    const double tf = t(0, n-1);

    // (1) Optimal cost and final time match the analytical MP result.
    EXPECT_NEAR(solution.get_cost(), 10.0/3.0, 5.0e-3) << "J* = 10/3";
    EXPECT_NEAR(tf, 4.0, 3.0e-2) << "tf* = 4";

    // (2) Boundary conditions reached.
    EXPECT_NEAR(x(0, n-1), 5.0, 1.0e-3);
    EXPECT_NEAR(x(1, n-1), 0.0, 1.0e-3);

    // (3) Bang structure: accelerate first (u=+1), decelerate last (u=-1).
    EXPECT_GT(u(0,0),   0.5) << "starts on u = +1";
    EXPECT_LT(u(0,n-1),-0.5) << "ends on u = -1";

    // (4) SINGULAR ARC over t in [1,3]: velocity held at x2 = 1, control u ~ 0.
    //     Sample strictly inside the arc to avoid the bang transitions.
    int cnt = 0; double x2err_max = 0.0, uabs_sum = 0.0;
    for (int k = 0; k < n; ++k) {
        double tk = t(0,k);
        if (tk > 1.3 && tk < 2.7) {                       // interior of the singular arc
            x2err_max = std::max(x2err_max, std::fabs(x(1,k) - 1.0));
            uabs_sum += std::fabs(u(0,k));
            cnt++;
        }
    }
    ASSERT_GT(cnt, 3) << "expected several nodes inside the singular arc";
    EXPECT_LT(x2err_max, 2.0e-2)      << "singular arc holds x2 = 1 (peak speed = sqrt of the arc value)";
    EXPECT_LT(uabs_sum / cnt, 1.0e-1) << "singular control u ~ 0 on the arc";
}
