//////////////////////////////////////////////////////////////////////////////
// test_bound_infinity.cpp
//
// What PSOPT means by an absent bound, and that scaling leaves it absent.
//
// Every bound PSOPT takes is two-sided, so a one-sided constraint is written by
// putting the unused side out of reach. Ipopt's nlp_lower_bound_inf defaults to
// -1e19 and it reads anything at or past that as -infinity, so a model written
// for Ipopt -- which is every model in examples/ -- says "no lower bound" by
// writing -1.0e19. Five shipped examples do exactly that on a path constraint:
// chance_constraint, chance_covariance, conic_sdp, conic_soc and path_window.
//
// Two things then went wrong, and each alone was enough to stop the SQP dead at
// its first subproblem on all five.
//
//   1. The SQP read an absent bound at 1e20 rather than 1e19, so -1e19 was a
//      finite bound. It reached the QP backend as a constraint row with a
//      right-hand side of ten million million million. Clarabel returned "dual
//      infeasible" after one iteration, with an objective of NaN, on a model
//      whose Hessian was the identity -- which cannot be unbounded -- and the
//      SQP reported that the backend had declined a non-convex subproblem.
//
//   2. Scaling multiplies every bound by a factor between 1e-7 and 1e7, and the
//      sentinel is a number like any other to a multiplication. On conic_sdp the
//      automatic constraint scaling took -1e19 to -3.3e+18, which is finite
//      under any threshold at all. Fixing the threshold alone left four of the
//      five failing.
//
// The remedy is PSOPT::bound_inf, PSOPT::no_lower_bound / no_upper_bound and
// PSOPT::scaled_lower_bound / scaled_upper_bound in psopt.h, the last pair
// mapping an absent bound to an IEEE infinity so that it stays absent under any
// further scaling.
//////////////////////////////////////////////////////////////////////////////

#include "gtest/gtest.h"
#include <psopt.h>

#include <cmath>
#include <limits>
#include <string>

bool psopt_qp_plugin_available(const std::string& backend, std::string& message);

namespace boundinf {

// A QP backend this build can actually load. GALAHAD is left out deliberately: it
// needs OMP_CANCELLATION and OMP_PROC_BIND set before the process starts, which is
// a condition on the run rather than on the build, and this file has nothing to say
// about backends. Empty means there is none and the SQP tests below skip.
static std::string usable_backend()
{
    const char* candidates[] = { "Clarabel", "PIQP", "ProxQP", "QPALM", "OSQP" };
    std::string ignored;
    for (size_t k = 0; k < sizeof candidates / sizeof candidates[0]; k++)
        if (psopt_qp_plugin_available(candidates[k], ignored)) return candidates[k];
    return std::string();
}

// ---------------------------------------------------------------------------
// The convention itself.
// ---------------------------------------------------------------------------

TEST(BoundInfinity, TheSentinelIsIpoptsAndTheTestsAreInclusive)
{
    // Ipopt's own default. If this number ever moves, every model in examples/
    // means a different problem to the two solvers, so it is worth pinning.
    EXPECT_EQ(PSOPT::bound_inf, 1.0e19);

    // At the sentinel, not merely past it -- Ipopt's test is "less or equal".
    EXPECT_TRUE (PSOPT::no_lower_bound(-1.0e19));
    EXPECT_TRUE (PSOPT::no_upper_bound( 1.0e19));
    EXPECT_TRUE (PSOPT::no_lower_bound(-1.0e20));
    EXPECT_TRUE (PSOPT::no_lower_bound(-PSOPT::inf));
    EXPECT_TRUE (PSOPT::no_upper_bound( PSOPT::inf));

    // And a bound a user could mean is not absent.
    EXPECT_FALSE(PSOPT::no_lower_bound(-1.0e18));
    EXPECT_FALSE(PSOPT::no_upper_bound( 1.0e18));
    EXPECT_FALSE(PSOPT::no_lower_bound(0.0));
}

TEST(BoundInfinity, ScalingLeavesAnAbsentBoundAbsent)
{
    // The conic_sdp case exactly: a constraint scale factor of one third.
    const double sc = 1.0/3.0;

    EXPECT_TRUE(std::isinf(PSOPT::scaled_lower_bound(-1.0e19, sc)));
    EXPECT_LT  (PSOPT::scaled_lower_bound(-1.0e19, sc), 0.0);
    EXPECT_TRUE(std::isinf(PSOPT::scaled_upper_bound( 1.0e19, sc)));
    EXPECT_GT  (PSOPT::scaled_upper_bound( 1.0e19, sc), 0.0);

    // An infinity stays one whatever it is multiplied by, which is the point of
    // returning one rather than the sentinel.
    EXPECT_TRUE(PSOPT::no_lower_bound(PSOPT::scaled_lower_bound(-1.0e19, sc)*1.0e-7));
    EXPECT_TRUE(PSOPT::no_upper_bound(PSOPT::scaled_upper_bound( 1.0e19, sc)*1.0e-7));

    // A bound the user meant is scaled, as it must be.
    EXPECT_DOUBLE_EQ(PSOPT::scaled_lower_bound(-6.0, sc), -2.0);
    EXPECT_DOUBLE_EQ(PSOPT::scaled_upper_bound( 6.0, sc),  2.0);
    EXPECT_DOUBLE_EQ(PSOPT::scaled_lower_bound( 0.0, sc),  0.0);
}


// ---------------------------------------------------------------------------
// And the same thing end to end.
//
// Minimum-energy double integrator with one one-sided path constraint. The
// constraint is written with a factor of a hundred in it so that the automatic
// constraint scaling gives its rows a factor far from one: that is what turned
// -1e19 into -3.3e+18 on conic_sdp, and a problem whose scale factors were all
// 1.0 would exercise only half of the defect.
//
//   min (1/2) INT_0^1 u^2 dt   s.t.  x1' = x2, x2' = u,
//   (x1,x2)(0) = (0,0), (x1,x2)(1) = (1,0),  100*x2 <= 120.
//
// The exact unconstrained answer is u* = 6 - 12t, J* = 6, whose largest x2 is
// 1.5; the path bound at 1.2 binds, so the constraint is doing something and a
// run that ignored it would give a different cost.
// ---------------------------------------------------------------------------

adouble endpoint_cost(adouble*, adouble*, adouble*, adouble&, adouble&, adouble*,
                      int, Workspace*)
{ return 0.0; }

adouble integrand_cost(adouble*, adouble* c, adouble*, adouble&, adouble*, int,
                       Workspace*)
{ adouble u = c[0]; return 0.5*u*u; }

void dae(adouble* d, adouble* path, adouble* s, adouble* c, adouble*, adouble&,
         adouble*, int, Workspace*)
{
    d[0] = s[1];
    d[1] = c[0];
    path[0] = 100.0*s[1];
}

void events(adouble* e, adouble* i, adouble* f, adouble*, adouble&, adouble&,
            adouble*, int, Workspace*)
{ e[0] = i[0]; e[1] = i[1]; e[2] = f[0]; e[3] = f[1]; }

void linkages(adouble*, adouble*, Workspace*) {}

struct Run {
    int         flag;
    double      cost;
    std::string message;
};

// how_absent: the value written for the unused side of the path constraint.
static Run solve(const std::string& nlp_method, double how_absent,
                 const std::string& backend = std::string())
{
    Alg algorithm; Sol solution; Prob problem;
    Run out; out.flag = -1; out.cost = 0.0;

    problem.name        = "one-sided path constraint";
    problem.outfilename = "test_bound_infinity.txt";
    problem.nphases     = 1;
    problem.nlinkages   = 0;
    psopt_level1_setup(problem);

    problem.phases(1).nstates   = 2;
    problem.phases(1).ncontrols = 1;
    problem.phases(1).nevents   = 4;
    problem.phases(1).npath     = 1;
    problem.phases(1).nodes     << 40;
    psopt_level2_setup(problem, algorithm);

    problem.phases(1).bounds.lower.states   << -5.0, -10.0;
    problem.phases(1).bounds.upper.states   <<  5.0,  10.0;
    problem.phases(1).bounds.lower.controls(0) = -60.0;
    problem.phases(1).bounds.upper.controls(0) =  60.0;

    problem.phases(1).bounds.lower.path(0) = how_absent;
    problem.phases(1).bounds.upper.path(0) = 120.0;

    problem.phases(1).bounds.lower.events   << 0.0, 0.0, 1.0, 0.0;
    problem.phases(1).bounds.upper.events   << 0.0, 0.0, 1.0, 0.0;
    problem.phases(1).bounds.lower.StartTime = 0.0;
    problem.phases(1).bounds.upper.StartTime = 0.0;
    problem.phases(1).bounds.lower.EndTime   = 1.0;
    problem.phases(1).bounds.upper.EndTime   = 1.0;

    problem.integrand_cost = &integrand_cost;
    problem.endpoint_cost  = &endpoint_cost;
    problem.dae            = &dae;
    problem.events         = &events;
    problem.linkages       = &linkages;

    problem.phases(1).guess.states   = zeros(2, 40);
    problem.phases(1).guess.controls = zeros(1, 40);
    problem.phases(1).guess.time     = linspace(0.0, 1.0, 40);

    algorithm.nlp_method         = nlp_method;
    algorithm.scaling            = "automatic";
    algorithm.derivatives        = "automatic";
    algorithm.collocation_method = "Legendre";
    algorithm.nlp_iter_max       = 1000;
    algorithm.nlp_tolerance      = 1.0e-6;
    algorithm.mesh_refinement    = "manual";
    algorithm.print_level        = 0;
    if (nlp_method == "SQP") {
        algorithm.hessian   = "exact";
        algorithm.qp_solver = backend;
    }

    out.flag    = psopt(solution, problem, algorithm);
    out.message = solution.error_msg;
    if (out.flag == 0) out.cost = solution.cost;
    return out;
}

} // namespace boundinf


// The two spellings of "no lower bound" must pose the same problem. Under Ipopt
// they always did, because Ipopt applies its own 1e19 rule to whatever it is
// handed -- so this half of the test would have passed throughout and is here to
// say what the answer is.
TEST(BoundInfinity, TheTwoSpellingsAgreeUnderIpopt)
{
    const boundinf::Run sentinel = boundinf::solve("IPOPT", -1.0e19);
    const boundinf::Run infinite = boundinf::solve("IPOPT", -PSOPT::inf);

    ASSERT_EQ(sentinel.flag, 0) << sentinel.message;
    ASSERT_EQ(infinite.flag, 0) << infinite.message;
    EXPECT_NEAR(sentinel.cost, infinite.cost, 1.0e-6*std::fabs(infinite.cost));

    // The path constraint binds, so this is not the unconstrained answer of 6.
    EXPECT_GT(sentinel.cost, 6.0);
}


// And under the SQP, where they did not. Written -1.0e19, the bound reached the
// subproblem as a finite row of 1e19 -- or, once the threshold was right but the
// scaling was not, of 1e17 -- and the run stopped at iteration zero.
#ifdef USE_SQP
TEST(BoundInfinity, TheTwoSpellingsAgreeUnderTheSqp)
{
    const std::string backend = boundinf::usable_backend();
    if (backend.empty()) GTEST_SKIP() << "no QP backend this build can load";

    const boundinf::Run infinite = boundinf::solve("SQP", -PSOPT::inf, backend);
    if (infinite.flag != 0)
        GTEST_SKIP() << "the SQP does not solve this problem through " << backend
                     << ": " << infinite.message;

    const boundinf::Run sentinel = boundinf::solve("SQP", -1.0e19, backend);

    ASSERT_EQ(sentinel.flag, 0)
        << "a lower bound written -1.0e19, which is how a model written for "
           "Ipopt says there is none, was not read as absent: "
        << sentinel.message;

    EXPECT_NEAR(sentinel.cost, infinite.cost, 1.0e-4*std::fabs(infinite.cost));
}
#endif


// The same question asked of the shipped example that first showed it. conic_sdp
// is a two-by-two linear matrix inequality enforced through its determinant,
// c >= b^2, written as b^2 - c <= 0 with the lower side at -1.0e19. Its
// constraint rows scale by a third, which is what defeated the threshold fix on
// its own.
#ifdef USE_SQP
namespace boundinf_sdp {

adouble endpoint_cost(adouble*, adouble* f, adouble*, adouble&, adouble&, adouble*,
                      int, Workspace*)
{ return f[1]; }

adouble integrand_cost(adouble*, adouble* c, adouble*, adouble&, adouble*, int,
                       Workspace*)
{ adouble ub = c[0], uc = c[1]; return 0.01*(ub*ub + uc*uc); }

void dae(adouble* d, adouble* path, adouble* s, adouble* c, adouble*, adouble&,
         adouble*, int, Workspace*)
{
    adouble b = s[0], cc = s[1];
    d[0] = c[0];  d[1] = c[1];
    path[0] = b*b - cc;
}

void events(adouble* e, adouble* i, adouble* f, adouble*, adouble&, adouble&,
            adouble*, int, Workspace*)
{ e[0] = i[0]; e[1] = i[1]; e[2] = f[0]; }

void linkages(adouble*, adouble*, Workspace*) {}

} // namespace boundinf_sdp

TEST(BoundInfinity, TheSemidefiniteExampleSolvesUnderTheSqp)
{
    const std::string backend = boundinf::usable_backend();
    if (backend.empty()) GTEST_SKIP() << "no QP backend this build can load";

    Alg algorithm; Sol solution; Prob problem;

    problem.name        = "conic_sdp, as shipped";
    problem.outfilename = "test_bound_infinity_sdp.txt";
    problem.nphases     = 1;
    problem.nlinkages   = 0;
    psopt_level1_setup(problem);

    problem.phases(1).nstates   = 2;
    problem.phases(1).ncontrols = 2;
    problem.phases(1).nevents   = 3;
    problem.phases(1).npath     = 1;
    problem.phases(1).nodes     << 30;
    psopt_level2_setup(problem, algorithm);

    problem.phases(1).bounds.lower.states   << -2.0, -0.1;
    problem.phases(1).bounds.upper.states   <<  2.0,  3.0;
    problem.phases(1).bounds.lower.controls << -5.0, -5.0;
    problem.phases(1).bounds.upper.controls <<  5.0,  5.0;
    problem.phases(1).bounds.lower.path(0)  = -1.0e19;
    problem.phases(1).bounds.upper.path(0)  =  0.0;
    problem.phases(1).bounds.lower.events   << 0.0, 1.0, 0.8;
    problem.phases(1).bounds.upper.events   << 0.0, 1.0, 0.8;
    problem.phases(1).bounds.lower.StartTime = 0.0;
    problem.phases(1).bounds.upper.StartTime = 0.0;
    problem.phases(1).bounds.lower.EndTime   = 1.0;
    problem.phases(1).bounds.upper.EndTime   = 1.0;

    problem.integrand_cost = &boundinf_sdp::integrand_cost;
    problem.endpoint_cost  = &boundinf_sdp::endpoint_cost;
    problem.dae            = &boundinf_sdp::dae;
    problem.events         = &boundinf_sdp::events;
    problem.linkages       = &boundinf_sdp::linkages;

    problem.phases(1).guess.states           = zeros(2,30);
    problem.phases(1).guess.states.row(0)    = linspace(0.0,0.8,30);
    problem.phases(1).guess.states.row(1)    = linspace(1.0,1.0,30);
    problem.phases(1).guess.controls         = zeros(2,30);
    problem.phases(1).guess.time             = linspace(0.0,1.0,30);

    algorithm.nlp_method         = "SQP";
    algorithm.hessian            = "exact";
    algorithm.qp_solver          = backend;
    algorithm.scaling            = "automatic";
    algorithm.derivatives        = "automatic";
    algorithm.collocation_method = "Legendre";
    algorithm.nlp_iter_max       = 1000;
    algorithm.nlp_tolerance      = 1.0e-6;
    algorithm.mesh_refinement    = "manual";
    algorithm.print_level        = 0;

    const int flag = psopt(solution, problem, algorithm);
    ASSERT_EQ(flag, 0) << solution.error_msg;

    // Ipopt's answer for the same mesh. c(T) is driven down to the positive
    // semidefinite boundary, c = b^2 = 0.64, and the rest is the control energy.
    EXPECT_NEAR(solution.cost, 0.647696, 1.0e-4);
}
#endif
