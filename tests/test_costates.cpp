//////////////////////////////////////////////////////////////////////////////
// test_costates.cpp
//
// The covector mapping, against a problem whose adjoint is known in closed form.
//
//     minimize  J = int_0^1 (1/2)( x^2 + u^2 ) dt   subject to  xdot = u,
//                   x(0) = 1,  x(1) = 0.75.
//
// With H = L + lambda^T f -- the convention of the book and of PSOPT's own
// Hamiltonian, see solution.dual.Hamiltonian -- stationarity gives u = -lambda and
// the adjoint equation lambdadot = -dH/dx = -x, so xddot = x and
//
//     x(t) = cosh t + B sinh t,   B = (0.75 - cosh 1)/sinh 1,
//     lambda(t) = -xdot(t) = -( sinh t + B cosh t ).
//
// Two things are under test, and the first is what makes this file worth having.
//
// The magnitude. The local-collocation branch of the covector mapping scales the
// defect multipliers by (tf - t0), and tf - t0 was computed from the *second* node
// rather than the first: every costate from a trapezoidal or Hermite-Simpson run
// came out short by exactly one interval's worth, a factor (M-1)/M on a mesh of M
// intervals. It is 2.6 per cent at the 39 intervals used here and a third of one
// per cent by 300, which is why it survived -- it looks like discretization error
// until one notices that two schemes of different order return bit-identical
// costates. The pseudospectral branches divide by quadrature weights instead and
// were never affected, which is the other half of the diagnosis.
//
// The sign. PSOPT's defect is xdot - f while the continuous theory adjoins
// lambda^T (f - xdot), so the relation between the two is worth pinning rather
// than trusting.
//////////////////////////////////////////////////////////////////////////////

#include "gtest/gtest.h"
#include <psopt.h>
#include <cmath>
#include <string>

namespace costate_test {

adouble endpoint_cost(adouble* i0, adouble* xf, adouble* p, adouble& t0, adouble& tf,
                      adouble* xad, int iphase, Workspace* w) { return 0.0; }
adouble integrand_cost(adouble* s, adouble* u, adouble* p, adouble& t,
                       adouble* xad, int iphase, Workspace* w)
{ return 0.5*(s[0]*s[0] + u[0]*u[0]); }
void dae(adouble* d, adouble* path, adouble* st, adouble* u, adouble* p,
         adouble& t, adouble* xad, int iphase, Workspace* w) { d[0] = u[0]; }
void events(adouble* e, adouble* i0, adouble* xf, adouble* p, adouble& t0, adouble& tf,
            adouble* xad, int iphase, Workspace* w) { e[0] = i0[0]; e[1] = xf[0]; }
void linkages(adouble* l, adouble* xad, Workspace* w) { }

// Solve on `nodes` nodes with the given collocation method and return the costate
// and the time it belongs to, in `lam` and `tt`.
static bool solve(const std::string& method, int nodes, MatrixXd& lam, MatrixXd& tt)
{
    Alg algorithm; Sol solution; Prob problem;
    problem.name        = "costate covector mapping";
    problem.outfilename = "test_costates.txt";
    problem.nphases     = 1;
    problem.nlinkages   = 0;
    psopt_level1_setup(problem);

    problem.phases(1).nstates   = 1;
    problem.phases(1).ncontrols = 1;
    problem.phases(1).nevents   = 2;
    problem.phases(1).npath     = 0;
    problem.phases(1).nodes     << nodes;
    psopt_level2_setup(problem, algorithm);

    problem.phases(1).bounds.lower.states   << -10.0;
    problem.phases(1).bounds.upper.states   <<  10.0;
    problem.phases(1).bounds.lower.controls << -10.0;
    problem.phases(1).bounds.upper.controls <<  10.0;
    problem.phases(1).bounds.lower.events   << 1.0, 0.75;
    problem.phases(1).bounds.upper.events   << 1.0, 0.75;
    problem.phases(1).bounds.lower.StartTime = 0.0;
    problem.phases(1).bounds.upper.StartTime = 0.0;
    problem.phases(1).bounds.lower.EndTime   = 1.0;
    problem.phases(1).bounds.upper.EndTime   = 1.0;

    problem.integrand_cost = &integrand_cost;
    problem.endpoint_cost  = &endpoint_cost;
    problem.dae            = &dae;
    problem.events         = &events;
    problem.linkages       = &linkages;

    problem.phases(1).guess.states   = zeros(1, nodes);
    problem.phases(1).guess.controls = zeros(1, nodes);
    problem.phases(1).guess.time     = linspace(0.0, 1.0, nodes);

    algorithm.nlp_method         = "IPOPT";
    algorithm.scaling            = "automatic";
    algorithm.derivatives        = "automatic";
    algorithm.nlp_tolerance      = 1.0e-10;
    algorithm.nlp_iter_max       = 500;
    algorithm.collocation_method = method;
    algorithm.mesh_refinement    = "manual";
    algorithm.print_level        = 0;

    if (psopt(solution, problem, algorithm) != 0) return false;
    lam = solution.get_dual_costates_in_phase(1);
    tt  = solution.get_time_in_phase(1);
    return true;
}

static double exact_costate(double t)
{
    const double B = (0.75 - std::cosh(1.0))/std::sinh(1.0);
    return -(std::sinh(t) + B*std::cosh(t));
}

} // namespace costate_test


// The local schemes. A relative error of 2 per cent here is the (M-1)/M scaling
// returning; the tolerance is deliberately far tighter than that and far looser than
// the schemes' own accuracy, so this fails on the bug and on nothing else.
TEST(Costates, LocalCollocationRecoversTheAdjointToTheRightScale)
{
    for (const std::string method : {std::string("trapezoidal"), std::string("Hermite-Simpson")}) {
        MatrixXd lam, tt;
        ASSERT_TRUE(costate_test::solve(method, 40, lam, tt)) << method;
        ASSERT_GT(lam.cols(), 10) << method;

        // Measured against the size of the costate rather than against its value at the
        // point: this adjoint crosses zero near t = 0.82, where a relative test says
        // nothing. The scale is the costate's own maximum, |lambda(0)|.
        const double scale = std::fabs(costate_test::exact_costate(0.0));

        // Interior points only: the ends of a local mesh carry their own one-sided error.
        for (int k = 4; k < lam.cols() - 4; k++) {
            const double ex = costate_test::exact_costate(tt(0,k));
            EXPECT_NEAR(lam(0,k), ex, 3.0e-3*scale)
                << method << " costate at t = " << tt(0,k)
                << " is " << lam(0,k) << " against an exact " << ex;
        }
    }
}

// The pseudospectral branch divides by the quadrature weights and takes a different
// route through the mapping; it is checked here so that a change to one branch cannot
// silently be made to the other.
TEST(Costates, PseudospectralRecoversTheAdjoint)
{
    MatrixXd lam, tt;
    ASSERT_TRUE(costate_test::solve("Legendre", 40, lam, tt));
    const double scale = std::fabs(costate_test::exact_costate(0.0));
    for (int k = 2; k < lam.cols() - 2; k++) {
        const double ex = costate_test::exact_costate(tt(0,k));
        EXPECT_NEAR(lam(0,k), ex, 5.0e-3*scale)
            << "Legendre costate at t = " << tt(0,k);
    }
}

// The sign, stated once. PSOPT's defect is xdot - f and the continuous theory adjoins
// lambda^T (f - xdot); the costate PSOPT reports is the one belonging to
// H = L + lambda^T f, so on this problem it is positive over most of the interval.
TEST(Costates, SignMatchesTheHamiltonianConvention)
{
    MatrixXd lam, tt;
    ASSERT_TRUE(costate_test::solve("Legendre", 40, lam, tt));
    EXPECT_GT(lam(0,0), 0.0) << "the initial costate should be positive on this problem";
    EXPECT_NEAR(lam(0,0), costate_test::exact_costate(0.0), 5.0e-3);
}
