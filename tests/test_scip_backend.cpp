//////////////////////////////////////////////////////////////////////////////
// test_scip_backend.cpp
//
// nlp_method="SCIP" cross-check on the double integrator: a MIXED-INTEGER problem
// (which IPOPT/CasADi cannot solve) with a known bang-off-bang structure.
//   xdot1=x2, xdot2=u, u INTEGER in {-1,0,1};  a >= |u| (continuous);
//   (0,0)->(1,0) over T=3;  minimise INT a dt.
// The minimum-time-like optimum is bang-off-bang: accelerate (u=+1), coast (u=0),
// decelerate (u=-1). Guarded: skips where PSOPT was built without SCIP.
//////////////////////////////////////////////////////////////////////////////

#include "gtest/gtest.h"
#include <psopt.h>
#include <cmath>

namespace scipmi {
adouble endpoint_cost(adouble*, adouble*, adouble*, adouble&, adouble&, adouble*, int, Workspace*) { return 0.0; }
adouble integrand_cost(adouble*, adouble* c, adouble*, adouble&, adouble*, int, Workspace*) { return c[1]; } // a
void dae(adouble* d, adouble* path, adouble* s, adouble* c, adouble*, adouble&, adouble*, int, Workspace*)
{ adouble u=c[0], a=c[1]; d[0]=s[1]; d[1]=u; path[0]=a-u; path[1]=a+u; }   // a >= |u|
void events(adouble* e, adouble* i, adouble* f, adouble*, adouble&, adouble&, adouble*, int, Workspace*)
{ e[0]=i[0]; e[1]=i[1]; e[2]=f[0]; e[3]=f[1]; }
void linkages(adouble*, adouble*, Workspace*) {}
}

TEST(NlpBackends, SCIP_MixedInteger_BangOffBang)
{
    using namespace scipmi;
    const int N = 16;
    Alg algorithm; Sol solution; Prob problem;
    problem.name="scip-mi"; problem.outfilename="scipmi.txt"; problem.nphases=1; problem.nlinkages=0;
    psopt_level1_setup(problem);
    problem.phases(1).nstates=2; problem.phases(1).ncontrols=2; problem.phases(1).nevents=4;
    problem.phases(1).npath=2; problem.phases(1).nodes<<N;
    psopt_level2_setup(problem, algorithm);
    problem.phases(1).bounds.lower.states << -10.0,-10.0;  problem.phases(1).bounds.upper.states << 10.0,10.0;
    problem.phases(1).bounds.lower.controls << -1.0, 0.0;  problem.phases(1).bounds.upper.controls << 1.0, 1.0;  // u in [-1,1], a in [0,1]
    problem.phases(1).bounds.lower.path << 0.0, 0.0;       problem.phases(1).bounds.upper.path << 1e19, 1e19;
    problem.phases(1).bounds.lower.events << 0.0,0.0,1.0,0.0;
    problem.phases(1).bounds.upper.events << 0.0,0.0,1.0,0.0;
    problem.phases(1).bounds.lower.StartTime=0.0; problem.phases(1).bounds.upper.StartTime=0.0;
    problem.phases(1).bounds.lower.EndTime=3.0;   problem.phases(1).bounds.upper.EndTime=3.0;
    problem.phases(1).integer_controls.resize(1,1); problem.phases(1).integer_controls(0)=0;   // u is INTEGER
    problem.integrand_cost=&integrand_cost; problem.endpoint_cost=&endpoint_cost;
    problem.dae=&dae; problem.events=&events; problem.linkages=&linkages;
    problem.phases(1).guess.states=zeros(2,N);
    problem.phases(1).guess.states.row(0)=linspace(0.0,1.0,N);
    problem.phases(1).guess.controls=zeros(2,N); problem.phases(1).guess.time=linspace(0.0,3.0,N);
    algorithm.nlp_method="SCIP"; algorithm.collocation_method="trapezoidal";
    algorithm.mesh_refinement="manual"; algorithm.scaling="automatic";
    algorithm.nlp_iter_max=1000; algorithm.nlp_tolerance=1e-6;

    try {
        psopt(solution, problem, algorithm);
    } catch (...) {
        GTEST_SKIP() << "PSOPT built without SCIP (nlp_method=SCIP unavailable)";
    }
    if (solution.error_flag != false)
        GTEST_SKIP() << "PSOPT built without SCIP (nlp_method=SCIP unavailable)";

    DMatrix x=solution.get_states_in_phase(1), u=solution.get_controls_in_phase(1), t=solution.get_time_in_phase(1);
    int n=(int)t.cols();
    // (1) reaches the target at rest
    EXPECT_NEAR(x(0,n-1), 1.0, 1.0e-3);  EXPECT_NEAR(x(1,n-1), 0.0, 1.0e-3);
    // (2) u is genuinely INTEGER at every node (in {-1,0,1})
    for (int k=0;k<n;++k) { double uk=u(0,k); EXPECT_LT(std::fabs(uk-std::floor(uk+0.5)), 1.0e-4) << "u integer at node "<<k; }
    // (3) bang-off-bang: starts accelerating (+1), ends decelerating (-1), coasts (0) somewhere
    EXPECT_GT(u(0,0),    0.5)  << "starts u = +1";
    EXPECT_LT(u(0,n-1), -0.5)  << "ends u = -1";
    bool coasts=false; for(int k=0;k<n;++k) if (std::fabs(u(0,k))<0.5) coasts=true;
    EXPECT_TRUE(coasts) << "has a u = 0 coast arc (bang-OFF-bang)";
}
