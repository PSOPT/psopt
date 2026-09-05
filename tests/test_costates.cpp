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

// ---------------------------------------------------------------------------
// A second problem, for the Legendre smoothing filter: the minimum-energy
// rest-to-rest double integrator of examples/mineng_di,
//
//     minimize  J = int_0^1 (1/2) u^2 dt   s.t.  xdot1 = x2, xdot2 = u,
//                   (x1,x2)(0) = (0,0),  (x1,x2)(1) = (1,0),
//
// whose adjoint is lambda1 = -12 and lambda2 = 12t - 6 -- one constant, one
// linear. That is the point of it. The LGL covector map returns both to nine
// figures, and the smoothing filter applied afterwards used the fixed stencil
// (1/4, 1/2, 1/4), which reproduces a linear function only on a uniform grid.
// On the non-uniform LGL nodes it scaled lambda2 by 0.9935 at 20 nodes while
// leaving the constant lambda1 alone, an error falling like 1/N^2 that reads as
// discretization error until one notices the states and the cost are exact.
// ---------------------------------------------------------------------------
namespace di {
adouble endpoint_cost(adouble* i0, adouble* xf, adouble* p, adouble& t0, adouble& tf,
                      adouble* xad, int iphase, Workspace* w) { return 0.0; }
adouble integrand_cost(adouble* s, adouble* u, adouble* p, adouble& t,
                       adouble* xad, int iphase, Workspace* w) { return 0.5*u[0]*u[0]; }
void dae(adouble* d, adouble* path, adouble* st, adouble* u, adouble* p,
         adouble& t, adouble* xad, int iphase, Workspace* w) { d[0] = st[1]; d[1] = u[0]; }
void events(adouble* e, adouble* i0, adouble* xf, adouble* p, adouble& t0, adouble& tf,
            adouble* xad, int iphase, Workspace* w)
{ e[0] = i0[0]; e[1] = i0[1]; e[2] = xf[0]; e[3] = xf[1]; }
void linkages(adouble* l, adouble* xad, Workspace* w) { }
}

static bool solve_di(const std::string& method, int nodes, MatrixXd& lam, MatrixXd& tt,
                     MatrixXd& uu, double& J, const std::string& scaling = "automatic")
{
    Alg algorithm; Sol solution; Prob problem;
    problem.name        = "minimum-energy double integrator";
    problem.outfilename = "test_costates_di.txt";
    problem.nphases     = 1;
    problem.nlinkages   = 0;
    psopt_level1_setup(problem);

    problem.phases(1).nstates   = 2;
    problem.phases(1).ncontrols = 1;
    problem.phases(1).nevents   = 4;
    problem.phases(1).npath     = 0;
    problem.phases(1).nodes     << nodes;
    psopt_level2_setup(problem, algorithm);

    problem.phases(1).bounds.lower.states   << -5.0, -5.0;
    problem.phases(1).bounds.upper.states   <<  5.0,  5.0;
    problem.phases(1).bounds.lower.controls << -50.0;
    problem.phases(1).bounds.upper.controls <<  50.0;
    problem.phases(1).bounds.lower.events   << 0.0, 0.0, 1.0, 0.0;
    problem.phases(1).bounds.upper.events   << 0.0, 0.0, 1.0, 0.0;
    problem.phases(1).bounds.lower.StartTime = 0.0;
    problem.phases(1).bounds.upper.StartTime = 0.0;
    problem.phases(1).bounds.lower.EndTime   = 1.0;
    problem.phases(1).bounds.upper.EndTime   = 1.0;

    problem.integrand_cost = &di::integrand_cost;
    problem.endpoint_cost  = &di::endpoint_cost;
    problem.dae            = &di::dae;
    problem.events         = &di::events;
    problem.linkages       = &di::linkages;

    problem.phases(1).guess.states           = zeros(2, nodes);
    problem.phases(1).guess.states.row(0)    = linspace(0.0, 1.0, nodes);
    problem.phases(1).guess.controls         = zeros(1, nodes);
    problem.phases(1).guess.time             = linspace(0.0, 1.0, nodes);

    algorithm.nlp_method         = "IPOPT";
    algorithm.scaling            = scaling;
    algorithm.derivatives        = "automatic";
    algorithm.nlp_tolerance      = 1.0e-10;
    algorithm.nlp_iter_max       = 500;
    algorithm.collocation_method = method;
    algorithm.mesh_refinement    = "manual";
    algorithm.print_level        = 0;

    if (psopt(solution, problem, algorithm) != 0) return false;
    lam = solution.get_dual_costates_in_phase(1);
    tt  = solution.get_time_in_phase(1);
    uu  = solution.get_controls_in_phase(1);
    J   = solution.cost;
    return true;
}

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


// The Legendre smoothing filter, against an adjoint that is exactly linear.
//
// The LGL scheme solves this problem exactly -- the states are polynomials of degree
// at most three and the cost integrand a quadratic, both inside what 20 nodes
// represent and integrate without error -- so any error left in the reported costate
// is the covector map's or the filter's, and nothing else's. That is what makes the
// tolerance below meaningful at 1e-6 on a costate of size 12.
//
// The fixed (1/4, 1/2, 1/4) stencil this replaces returns lambda2 short by a factor
// 0.9935 at 20 nodes -- an error of 5.8e-2, four orders outside the tolerance -- while
// returning the constant lambda1 exactly, which is the signature of a filter that is
// inconsistent on a non-uniform grid rather than of a discretization error.
TEST(Costates, LegendreSmoothingPreservesALinearAdjoint)
{
    MatrixXd lam, tt, uu;  double J = 0.0;
    ASSERT_TRUE(costate_test::solve_di("Legendre", 20, lam, tt, uu, J));
    ASSERT_EQ(lam.rows(), 2);
    ASSERT_GT(lam.cols(), 10);

    // The discretization really is exact here; if this fails the rest means nothing.
    EXPECT_NEAR(J, 6.0, 1.0e-8) << "the LGL solution of this problem is exact";

    for (int k = 0; k < lam.cols(); k++) {
        const double t = tt(0,k);
        EXPECT_NEAR(lam(0,k), -12.0, 1.0e-4)
            << "lambda1 at t = " << t << " (constant adjoint)";
        EXPECT_NEAR(lam(1,k), 12.0*t - 6.0, 1.0e-4)
            << "lambda2 at t = " << t << " (linear adjoint): a uniform shortfall here is"
               " the smoothing filter, not the covector map";
    }

    // dH/du = u + lambda2 = 0. Stated separately because it is the condition a user
    // checks, and it is the one the filter used to break on every Legendre run.
    for (int k = 0; k < lam.cols(); k++)
        EXPECT_NEAR(uu(0,k) + lam(1,k), 0.0, 1.0e-4)
            << "stationarity dH/du at t = " << tt(0,k);
}

// The same problem under Hermite-Simpson, where the filter does not apply: the local
// branch has to agree with the pseudospectral one on a problem both solve exactly, or
// one of the two mappings is wrong.
TEST(Costates, LocalAndPseudospectralAgreeOnALinearAdjoint)
{
    MatrixXd lam, tt, uu;  double J = 0.0;
    ASSERT_TRUE(costate_test::solve_di("Hermite-Simpson", 20, lam, tt, uu, J));
    EXPECT_NEAR(J, 6.0, 1.0e-8);
    for (int k = 0; k < lam.cols(); k++) {
        EXPECT_NEAR(lam(0,k), -12.0, 1.0e-4) << "lambda1 at t = " << tt(0,k);
        EXPECT_NEAR(lam(1,k), 12.0*tt(0,k) - 6.0, 1.0e-4) << "lambda2 at t = " << tt(0,k);
    }
}


// algorithm.scaling = "user", which with no factors set is unit scaling.
//
// This is a memory test before it is a numerical one. The t0 <= tf row of each phase is
// written at index phase_offset + ncons_phase_i - 1, and under user scaling its scale
// factor was stored at phase_offset + ncons_phase_i: one past. On the last phase that is
// one element past the end of constraint_scaling, which is sized nlp_ncons, so selecting
// a documented option corrupted the heap on a plain single-phase problem -- an abort here,
// a segmentation fault on older builds. On any earlier phase it silently set the FIRST row
// of the next phase to the previous phase's time scaling.
//
// What it asserts afterwards is worth having on its own: the costate PSOPT reports should
// not depend on how the constraints were scaled, since the recovery undoes the scaling it
// applied. Unit scaling and automatic scaling must therefore agree with the closed form,
// and they do.
TEST(Costates, UserScalingIsUnitScalingAndRecoversTheSameCostate)
{
    MatrixXd lam, tt, uu;  double J = 0.0;
    ASSERT_TRUE(costate_test::solve_di("Legendre", 40, lam, tt, uu, J, "user"));
    EXPECT_NEAR(J, 6.0, 1.0e-6);
    for (int k = 0; k < lam.cols(); k++) {
        EXPECT_NEAR(lam(0,k), -12.0, 1.0e-3) << "lambda1 at t = " << tt(0,k);
        EXPECT_NEAR(lam(1,k), 12.0*tt(0,k) - 6.0, 1.0e-3) << "lambda2 at t = " << tt(0,k);
    }
}
