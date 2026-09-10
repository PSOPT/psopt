//////////////////////////////////////////////////////////////////////////////
// test_ir_path_inequality.cpp
//
// The integrated-residual transcription with an INEQUALITY path constraint.
// No shipped example exercises that combination, which matters because the two
// kinds of path constraint are handled by different machinery: an equality
// constraint (lower bound == upper bound) is folded into the integrated
// residual as an algebraic row and its pointwise rows are freed to +/- infinity,
// while an inequality constraint is not folded and must be imposed pointwise, at
// the nodes and at the interval midpoints. ir_algebraic_rows decides which is
// which, and an error there is silent: a freed row is a constraint that is
// simply not imposed.
//
// The test problem is Bryson-Denham posed with its state limit as a path
// constraint rather than as a state bound:
//
//     min  (1/2) INT_0^1 u^2 dt   s.t.  x1' = x2,  x2' = u,
//          x1(0) = x1(1) = 0,  x2(0) = 1,  x2(1) = -1,   x1 <= l
//
// with l = 1/9. The exact answer is J* = 4/(9l) = 4, with a boundary arc on
// [1/3, 2/3] where x1 == l exactly and u == 0, and u(0) = -2/(3l) = -6. Drop the
// path constraint and the same problem has J = 2 with u == -2. A factor of two
// between "imposed" and "not imposed" is what makes this a test rather than an
// illustration.
//////////////////////////////////////////////////////////////////////////////

#include "gtest/gtest.h"
#include <psopt.h>

#include <cmath>
#include <string>

namespace irpath {

const double LIMIT       = 1.0/9.0;   // the state limit l
const double J_CONSTRAINED   = 4.0;   // exact, with the path constraint imposed
const double J_UNCONSTRAINED = 2.0;   // exact, with it absent or not imposed

int NPATH = 1;                        // the dae must not write path rows it has not got

adouble endpoint_cost(adouble*, adouble*, adouble*, adouble&, adouble&, adouble*,
                      int, Workspace*)
{ return 0.0; }

adouble integrand_cost(adouble*, adouble* c, adouble*, adouble&, adouble*, int,
                       Workspace*)
{ adouble u = c[0]; return 0.5*u*u; }

// x1' = x2, x2' = u. path[0] = x1, the inequality constraint. path[1], when it is
// declared, is u2 - u1: an equality path constraint with no effect on the dynamics
// or the cost, which exists only so that a folded and an unfolded path constraint
// are present at once.
//
// It ties one control to another, and deliberately not a control to a state. A
// folded algebraic row is driven to the residual box at every sample point on an
// element, so it imposes an identity between two polynomials; tie objects of
// different representation degree and the higher one is annihilated. u2 = x1 was
// tried first and moved the cost from 4.000 to 4.405 at ir_local_order 0, because
// the control is the quadratic through node, midpoint and node while the state is
// the cubic Hermite, and forcing them equal removes the state's cubic term. Two
// controls share a representation, so nothing is lost.
void dae(adouble* d, adouble* path, adouble* s, adouble* c, adouble*, adouble&,
         adouble*, int, Workspace*)
{
    d[0] = s[1];
    d[1] = c[0];
    if (NPATH >= 1) path[0] = s[0];
    if (NPATH >= 2) path[1] = c[1] - c[0];
}

void events(adouble* e, adouble* i, adouble* f, adouble*, adouble&, adouble&,
            adouble*, int, Workspace*)
{ e[0] = i[0]; e[1] = i[1]; e[2] = f[0]; e[3] = f[1]; }

void linkages(adouble*, adouble*, Workspace*) {}

struct Run {
    int      flag;
    double   cost;
    double   max_path;      // max over the nodes of x1
    double   max_alg_error; // max over the nodes of |u2 - u1|, when declared
    double   umin;
    std::string message;
};

// One solve. `order` is algorithm.ir_local_order; `nodes` must satisfy
// (nodes-1) % order == 0 when order >= 2, and `residual_nodes` must be at least
// order+2, both of which validate_user_input enforces.
static Run solve(int order, int nodes, int residual_nodes, double upper_path_bound,
                 int npath = 1)
{
    Alg algorithm; Sol solution; Prob problem;
    Run out; out.flag = -1; out.cost = 0.0; out.max_path = -1e30;
    out.max_alg_error = 0.0; out.umin = 0.0;

    NPATH = npath;

    const int ncontrols = (npath >= 2) ? 2 : 1;

    problem.name        = "Bryson-Denham with a path constraint";
    problem.outfilename = "test_ir_path_inequality.txt";
    problem.nphases     = 1;
    problem.nlinkages   = 0;
    psopt_level1_setup(problem);

    problem.phases(1).nstates   = 2;
    problem.phases(1).ncontrols = ncontrols;
    problem.phases(1).nevents   = 4;
    problem.phases(1).npath     = npath;
    problem.phases(1).nodes     << nodes;
    psopt_level2_setup(problem, algorithm);

    problem.phases(1).bounds.lower.states   << -1.0, -5.0;
    problem.phases(1).bounds.upper.states   <<  1.0,  5.0;
    for (int c = 0; c < ncontrols; c++) {
        problem.phases(1).bounds.lower.controls(c) = -60.0;
        problem.phases(1).bounds.upper.controls(c) =  60.0;
    }
    problem.phases(1).bounds.lower.path(0) = -1.0;
    problem.phases(1).bounds.upper.path(0) = upper_path_bound;
    if (npath >= 2) {   // an equality: this is the one ir_algebraic_rows folds
        problem.phases(1).bounds.lower.path(1) = 0.0;
        problem.phases(1).bounds.upper.path(1) = 0.0;
    }
    problem.phases(1).bounds.lower.events   << 0.0, 1.0, 0.0, -1.0;
    problem.phases(1).bounds.upper.events   << 0.0, 1.0, 0.0, -1.0;
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
    problem.phases(1).guess.controls = zeros(ncontrols, nodes);
    problem.phases(1).guess.time     = linspace(0.0, 1.0, nodes);

    algorithm.nlp_method           = "IPOPT";
    algorithm.scaling              = "automatic";
    algorithm.derivatives          = "automatic";
    algorithm.collocation_method   = "Hermite-Simpson";
    algorithm.transcription_method = "integrated-residual";
    algorithm.ir_local_order       = order;
    algorithm.ir_objective         = "cost";      // optimality, with the residual boxed
    algorithm.ir_residual_bound    = 1.0e-8;
    algorithm.ir_residual_nodes    = residual_nodes;
    algorithm.nlp_iter_max         = 3000;
    algorithm.nlp_tolerance        = 1.0e-8;
    algorithm.print_level          = 0;

    out.flag    = psopt(solution, problem, algorithm);
    out.message = solution.error_msg;
    if (out.flag != 0) return out;

    out.cost = solution.cost;
    MatrixXd x = solution.get_states_in_phase(1);
    MatrixXd u = solution.get_controls_in_phase(1);
    const int n = (int) x.cols();
    for (int k = 0; k < n; k++) {
        if (x(0,k) > out.max_path) out.max_path = x(0,k);
        if (u(0,k) < out.umin)     out.umin     = u(0,k);
        if (npath >= 2) {
            const double e = std::fabs(u(1,k) - u(0,k));
            if (e > out.max_alg_error) out.max_alg_error = e;
        }
    }
    return out;
}

} // namespace irpath


// --------------------------------------------------------------------------
// The constraint is imposed, and imposed as an inequality.
// --------------------------------------------------------------------------

TEST(IntegratedResidualPath, AnInequalityPathConstraintIsImposed)
{
    // ir_local_order 0 is the cubic-Hermite local representation; 3 and 4 are the
    // Nie-Kerrigan flexible-order one. Order 2 is deliberately absent; see the
    // disabled test at the foot of this file.
    struct { int order, nodes, rnodes; } cases[] = { {0, 41, 4}, {3, 40, 5}, {4, 41, 6} };

    for (size_t k = 0; k < sizeof cases / sizeof cases[0]; k++) {
        const irpath::Run r = irpath::solve(cases[k].order, cases[k].nodes,
                                            cases[k].rnodes, irpath::LIMIT);
        const int d = cases[k].order;

        ASSERT_EQ(r.flag, 0) << "ir_local_order = " << d << ": " << r.message;

        // The constraint holds at the nodes, to the NLP's own constraint tolerance.
        EXPECT_LE(r.max_path, irpath::LIMIT + 1.0e-6)
            << "ir_local_order = " << d << ": x1 reached " << r.max_path
            << ", above the limit " << irpath::LIMIT;

        // And it was actually binding: the cost is the constrained optimum and not
        // the unconstrained one. This is the check that a freed or dropped row
        // cannot pass.
        EXPECT_NEAR(r.cost, irpath::J_CONSTRAINED, 5.0e-3)
            << "ir_local_order = " << d;
        EXPECT_GT(std::fabs(r.cost - irpath::J_UNCONSTRAINED), 1.0)
            << "ir_local_order = " << d
            << ": the cost is the unconstrained optimum, so the path constraint was "
               "not imposed";

        // The exact solution starts at u(0) = -2/(3l) = -6. A representation that
        // cannot reach it is not solving this problem, whatever its cost says.
        EXPECT_NEAR(r.umin, -6.0, 0.05) << "ir_local_order = " << d;
    }
}


// --------------------------------------------------------------------------
// And it costs nothing when it does not bind. The pointwise path rows exist in
// both cases; if their presence perturbed the solve, this is where it would show.
// --------------------------------------------------------------------------

TEST(IntegratedResidualPath, AnInactiveInequalityPathConstraintChangesNothing)
{
    const irpath::Run r = irpath::solve(0, 41, 4, 10.0);   // a bound far above the arc
    ASSERT_EQ(r.flag, 0) << r.message;
    EXPECT_NEAR(r.cost, irpath::J_UNCONSTRAINED, 5.0e-3)
        << "an inactive path constraint moved the answer away from the "
           "unconstrained optimum";
    EXPECT_NEAR(r.umin, -2.0, 0.05);
}


// --------------------------------------------------------------------------
// An equality and an inequality path constraint at once, which is where the
// folding has to pick the right one. path[0] is the inequality x1 <= l and must
// stay pointwise; path[1] is the equality u2 - x1 = 0 and is folded into the
// residual. Fold the wrong one and path[0]'s rows are freed to +/- infinity, so
// the cost falls to the unconstrained 2.
// --------------------------------------------------------------------------

TEST(IntegratedResidualPath, AnEqualityAndAnInequalityPathConstraintTogether)
{
    const irpath::Run r = irpath::solve(0, 41, 4, irpath::LIMIT, 2);
    ASSERT_EQ(r.flag, 0) << r.message;

    EXPECT_LE(r.max_path, irpath::LIMIT + 1.0e-6)
        << "the inequality path constraint was not imposed";
    EXPECT_NEAR(r.cost, irpath::J_CONSTRAINED, 5.0e-3)
        << "the cost is not the constrained optimum, so the wrong path constraint "
           "was folded into the residual";

    // The folded equality is enforced through the residual rather than pointwise,
    // and the residual box is 1e-8, so u2 should track u1 to about that. A loose
    // tolerance here would pass on an equality that was not enforced at all.
    EXPECT_LT(r.max_alg_error, 1.0e-4)
        << "the folded equality path constraint u2 = u1 is not being enforced: "
           "the largest departure at a node was " << r.max_alg_error;
}


// --------------------------------------------------------------------------
// ir_local_order = 2 is covered by tests/test_ir_local_order.cpp.
//
// Writing this file turned up a defect that had nothing to do with path
// constraints: the Nie-Kerrigan elements shared their end controls, which made
// the control continuous, and at d = 2 the residual box leaves the control
// constant on each element -- so a continuous piecewise-constant control is one
// constant, and refining the mesh lengthened the chain instead of helping. It
// was found here only because this was the first test to run that
// representation against a problem whose optimal control is not constant. The
// remedy and the tests for it live in the file named after it.
//
// What belongs here is that the constraint machinery works at d = 2 as well.
// --------------------------------------------------------------------------

TEST(IntegratedResidualPath, AnInequalityPathConstraintIsImposedAtLocalOrderTwo)
{
    const irpath::Run coarse = irpath::solve(2,  41, 4, irpath::LIMIT);
    const irpath::Run fine   = irpath::solve(2, 161, 4, irpath::LIMIT);
    ASSERT_EQ(coarse.flag, 0) << coarse.message;
    ASSERT_EQ(fine.flag,   0) << fine.message;

    EXPECT_LE(fine.max_path, irpath::LIMIT + 1.0e-6);
    EXPECT_NEAR(fine.cost, irpath::J_CONSTRAINED, 5.0e-3);
    EXPECT_LT(std::fabs(fine.cost   - irpath::J_CONSTRAINED),
              std::fabs(coarse.cost - irpath::J_CONSTRAINED))
        << "refining the mesh moved the answer away from the exact one";
}
