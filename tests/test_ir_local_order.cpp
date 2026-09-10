//////////////////////////////////////////////////////////////////////////////
// test_ir_local_order.cpp
//
// The Nie-Kerrigan local representation (algorithm.ir_local_order = d >= 2),
// and in particular that it converges at d = 2.
//
// It did not. On an element of local order d the state is degree d, so its
// derivative is degree d-1, and the residual box drives x' - f to zero at every
// sample point on the element. For a chain of integrators f carries the next
// state, of degree d, and then the control, also of degree d, so each link of
// the chain annihilates one degree: at d = 2 the control was left constant on
// each element. Elements used to share their end controls, which made the
// control continuous, and a continuous piecewise-constant function is a single
// constant. The answer converged to the best chained-constant control instead of
// the true one, and refining the mesh lengthened the chain rather than helping.
//
// Measured on the minimum-energy double integrator below, whose exact answer is
// J* = 6 with the linear control u* = 6 - 12t, at ir_residual_bound = 1e-8:
//
//     nodes            21      41      81     161
//     before          7.143   7.508   7.744   7.868      u reached +/-4, not +/-6
//     after           6.061   6.015   6.004   6.001      u reaches +/-5.93
//
// where 8 and u = +/-4 are exactly the best piecewise-constant control for that
// problem. The remedy was to stop sharing: a state must be C0 because it is a
// state, a control need not be, so each element now carries its own control at
// its left end. See ir_extra_control_vars and get_element_controls.
//
// The problem here has no path constraint. The defect was found through
// tests/test_ir_path_inequality.cpp, but it never had anything to do with path
// constraints -- that was simply the first test to run this representation
// against a problem whose optimal control is not constant.
//////////////////////////////////////////////////////////////////////////////

#include "gtest/gtest.h"
#include <psopt.h>

#include <cmath>
#include <string>

namespace irorder {

// min (1/2) INT_0^1 u^2 dt  s.t.  x1' = x2, x2' = u,  (x1,x2)(0) = (0,0),
// (x1,x2)(1) = (1,0).  Exactly u* = 6 - 12t and J* = 6.
const double J_EXACT = 6.0;
const double U0_EXACT = 6.0;

adouble endpoint_cost(adouble*, adouble*, adouble*, adouble&, adouble&, adouble*,
                      int, Workspace*)
{ return 0.0; }

adouble integrand_cost(adouble*, adouble* c, adouble*, adouble&, adouble*, int,
                       Workspace*)
{ adouble u = c[0]; return 0.5*u*u; }

void dae(adouble* d, adouble*, adouble* s, adouble* c, adouble*, adouble&,
         adouble*, int, Workspace*)
{ d[0] = s[1]; d[1] = c[0]; }

void events(adouble* e, adouble* i, adouble* f, adouble*, adouble&, adouble&,
            adouble*, int, Workspace*)
{ e[0] = i[0]; e[1] = i[1]; e[2] = f[0]; e[3] = f[1]; }

void linkages(adouble*, adouble*, Workspace*) {}

struct Run {
    int         flag;
    double      cost;
    double      umax;      // largest |u| over the nodes
    std::string message;
};

static Run solve(int order, int nodes, int residual_nodes,
                 bool element_local_controls = true)
{
    Alg algorithm; Sol solution; Prob problem;
    Run out; out.flag = -1; out.cost = 0.0; out.umax = 0.0;

    problem.name        = "minimum-energy double integrator";
    problem.outfilename = "test_ir_local_order.txt";
    problem.nphases     = 1;
    problem.nlinkages   = 0;
    psopt_level1_setup(problem);

    problem.phases(1).nstates   = 2;
    problem.phases(1).ncontrols = 1;
    problem.phases(1).nevents   = 4;
    problem.phases(1).npath     = 0;
    problem.phases(1).nodes     << nodes;
    psopt_level2_setup(problem, algorithm);

    problem.phases(1).bounds.lower.states   << -5.0, -10.0;
    problem.phases(1).bounds.upper.states   <<  5.0,  10.0;
    problem.phases(1).bounds.lower.controls(0) = -60.0;
    problem.phases(1).bounds.upper.controls(0) =  60.0;
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

    problem.phases(1).guess.states   = zeros(2, nodes);
    problem.phases(1).guess.controls = zeros(1, nodes);
    problem.phases(1).guess.time     = linspace(0.0, 1.0, nodes);

    algorithm.nlp_method           = "IPOPT";
    algorithm.scaling              = "automatic";
    algorithm.derivatives          = "automatic";
    algorithm.collocation_method   = "Hermite-Simpson";
    algorithm.transcription_method = "integrated-residual";
    algorithm.ir_local_order       = order;
    algorithm.ir_objective         = "cost";
    algorithm.ir_residual_bound    = 1.0e-8;
    algorithm.ir_residual_nodes    = residual_nodes;
    algorithm.ir_element_local_controls = element_local_controls;
    algorithm.nlp_iter_max         = 3000;
    algorithm.nlp_tolerance        = 1.0e-8;
    algorithm.print_level          = 0;

    out.flag    = psopt(solution, problem, algorithm);
    out.message = solution.error_msg;
    if (out.flag != 0) return out;

    out.cost = solution.cost;
    MatrixXd u = solution.get_controls_in_phase(1);
    for (int k = 0; k < (int) u.cols(); k++)
        if (std::fabs(u(0,k)) > out.umax) out.umax = std::fabs(u(0,k));
    return out;
}

} // namespace irorder


// --------------------------------------------------------------------------
// The test that would have caught it: refine, and require the answer to improve.
//
// A wrong answer that is merely inaccurate is caught by a tolerance; a wrong
// answer a method converges to is not, and this one passed through the right
// value on the way past it as the residual bound was tightened. Only refinement
// distinguishes the two.
// --------------------------------------------------------------------------

TEST(IntegratedResidualLocalOrder, OrderTwoConvergesUnderRefinement)
{
    const irorder::Run coarse = irorder::solve(2,  41, 4);
    const irorder::Run fine   = irorder::solve(2, 161, 4);

    ASSERT_EQ(coarse.flag, 0) << coarse.message;
    ASSERT_EQ(fine.flag,   0) << fine.message;

    const double ec = std::fabs(coarse.cost - irorder::J_EXACT);
    const double ef = std::fabs(fine.cost   - irorder::J_EXACT);

    EXPECT_LT(ef, ec)
        << "refining from 41 to 161 nodes moved the cost from " << coarse.cost
        << " to " << fine.cost << ", away from the exact " << irorder::J_EXACT;

    // Four times the elements. Second order would give sixteen; eight leaves room
    // for the NLP tolerance without admitting a method that merely wanders.
    EXPECT_GT(ec, 8.0*ef)
        << "the error fell only from " << ec << " to " << ef
        << " for a fourfold refinement, which is not convergence";

    EXPECT_NEAR(fine.cost, irorder::J_EXACT, 5.0e-3);

    // And the control reaches the amplitude the exact solution has. While the
    // element controls were chained the whole first arc was one constant and |u|
    // stopped at 4; the cost alone would not have said so.
    EXPECT_GT(fine.umax, 5.5)
        << "the largest |u| was " << fine.umax << ", against an exact 6";
    EXPECT_LT(fine.umax, 6.5) << "the largest |u| was " << fine.umax;
    EXPECT_GT(fine.umax, coarse.umax)
        << "refining did not bring |u| closer to its exact value";
}


// --------------------------------------------------------------------------
// The orders that were already right stay right, and to the same accuracy. The
// element-boundary controls give every order a degree of freedom it did not
// have, so this is not a formality: it is the check that the extra freedom is
// harmless where it was not needed.
// --------------------------------------------------------------------------

TEST(IntegratedResidualLocalOrder, HigherOrdersAreUnchanged)
{
    struct { int order, nodes, rnodes; } cases[] = { {3, 79, 5}, {4, 81, 6}, {5, 81, 7} };

    for (size_t k = 0; k < sizeof cases / sizeof cases[0]; k++) {
        const irorder::Run r = irorder::solve(cases[k].order, cases[k].nodes,
                                              cases[k].rnodes);
        ASSERT_EQ(r.flag, 0) << "ir_local_order = " << cases[k].order
                             << ": " << r.message;
        EXPECT_NEAR(r.cost, irorder::J_EXACT, 1.0e-4)
            << "ir_local_order = " << cases[k].order;
        EXPECT_NEAR(r.umax, irorder::U0_EXACT, 5.0e-2)
            << "ir_local_order = " << cases[k].order;
    }
}


// --------------------------------------------------------------------------
// The cubic-Hermite representation, ir_local_order = 0, shares no elements and
// was never affected. Here so that a change to the element machinery that
// reached the legacy path would be seen.
// --------------------------------------------------------------------------

TEST(IntegratedResidualLocalOrder, TheCubicHermiteRepresentationIsUntouched)
{
    const irorder::Run r = irorder::solve(0, 41, 4);
    ASSERT_EQ(r.flag, 0) << r.message;
    EXPECT_NEAR(r.cost, irorder::J_EXACT, 1.0e-4);
    EXPECT_NEAR(r.umax, irorder::U0_EXACT, 5.0e-2);
}


// --------------------------------------------------------------------------
// The option that turns the remedy off, and what it is for.
//
// algorithm.ir_element_local_controls = false restores the shared control, and
// with it the defect. It exists because what PSOPT carries as a control is not
// always a control: the algebraic variable of a DAE is a function of the state
// and is continuous, so giving it a jump at every element boundary is a
// modelling error rather than a freedom. examples/dae_i3, whose "control" is the
// multiplier of a holonomic constraint, sets it false for exactly that reason
// and stops converging without it.
//
// So the option has to do something, and this is the test that it does: with the
// controls shared, refining the mesh must fail to help, which is the behaviour
// the default was introduced to remove.
// --------------------------------------------------------------------------

TEST(IntegratedResidualLocalOrder, SharedControlsReproduceTheOldNonConvergence)
{
    const irorder::Run coarse = irorder::solve(2,  41, 4, false);
    const irorder::Run fine   = irorder::solve(2, 161, 4, false);

    ASSERT_EQ(coarse.flag, 0) << coarse.message;
    ASSERT_EQ(fine.flag,   0) << fine.message;

    const double ec = std::fabs(coarse.cost - irorder::J_EXACT);
    const double ef = std::fabs(fine.cost   - irorder::J_EXACT);

    EXPECT_GT(ef, ec)
        << "with the controls shared, refining from 41 to 161 nodes moved the cost "
           "from " << coarse.cost << " to " << fine.cost
        << "; it used to move away from the exact " << irorder::J_EXACT
        << ", and if it no longer does then either the option or the defect it "
           "describes has changed";

    // The control cannot reach its exact amplitude while the element constants are
    // chained: it stops at the best two-level step, which for this problem is 4.
    EXPECT_LT(fine.umax, 4.5)
        << "the largest |u| was " << fine.umax
        << ", which is more than a chained-constant control can produce here";
}
