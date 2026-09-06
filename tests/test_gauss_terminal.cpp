//////////////////////////////////////////////////////////////////////////////
// test_gauss_terminal.cpp
//
// The Gauss scheme collocates strictly interior points, so x(+1) is an appended
// NLP variable rather than one of the norder+1 stored nodes. The event
// constraints are imposed on it correctly -- that was never in doubt -- but it
// was in nothing the solution accessors returned, so the reported trajectory
// stopped at the last Gauss node.
//
// On the linear tangent steering problem of the book, at 40 nodes, that put the
// reported final time 0.39 s short of t_f and left the last reported state at
// y = 407.999044, vx = 7.653981, vy = 4.919e-3 against required terminal values
// of 408, 7.66 and 0. It reads as a converged-to-the-wrong-answer failure and is
// nothing of the kind: solution.cost agreed with the Radau run to eight decimals.
//
// The problem below is the minimum-energy rest-to-rest double integrator of
// examples/mineng_di,
//
//     minimize  J = int_0^1 (1/2) u^2 dt   s.t.  xdot1 = x2, xdot2 = u,
//                   (x1,x2)(0) = (0,0),  (x1,x2)(1) = (1,0),   t_f = 1,
//
// whose exact solution is u* = 6 - 12t, x1* = 3t^2 - 2t^3, x2* = 6t - 6t^2 and
// J* = 6. Legendre and Gauss both solve it exactly -- the states are cubic, the
// control linear and the integrand quadratic -- so the terminal state is (1, 0)
// to machine precision and a tolerance of 1e-8 on it is meaningful.
//////////////////////////////////////////////////////////////////////////////

#include "gtest/gtest.h"
#include <psopt.h>
#include <cmath>
#include <string>

namespace gauss_terminal_test {

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

struct Result {
    MatrixXd x, u, t, lam, terminal, H;
    double   J = 0.0;
    bool     ok = false;
};

// hp = true splits the mesh into three intervals of unequal order, so that the
// terminal point is the right end of an interval carrying fewer nodes than the
// mesh as a whole -- where the gap the fix closes is largest.
static Result solve(const std::string& method, int nodes, bool hp = false)
{
    Result r;
    Alg algorithm; Sol solution; Prob problem;
    problem.name        = "Gauss terminal point";
    problem.outfilename = "test_gauss_terminal.txt";
    problem.nphases     = 1;
    problem.nlinkages   = 0;
    psopt_level1_setup(problem);

    problem.phases(1).nstates   = 2;
    problem.phases(1).ncontrols = 1;
    problem.phases(1).nevents   = 4;
    problem.phases(1).npath     = 0;
    problem.phases(1).nodes     << nodes;
    if (hp) {
        problem.phases(1).hp_breakpoints.resize(2); problem.phases(1).hp_breakpoints << 0.3, 0.7;
        problem.phases(1).hp_orders.resize(3);      problem.phases(1).hp_orders      << 8, 12, 8;
    }
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

    problem.integrand_cost = &integrand_cost;
    problem.endpoint_cost  = &endpoint_cost;
    problem.dae            = &dae;
    problem.events         = &events;
    problem.linkages       = &linkages;

    problem.phases(1).guess.states        = zeros(2, nodes);
    problem.phases(1).guess.states.row(0) = linspace(0.0, 1.0, nodes);
    problem.phases(1).guess.controls      = zeros(1, nodes);
    problem.phases(1).guess.time          = linspace(0.0, 1.0, nodes);

    algorithm.nlp_method         = "IPOPT";
    algorithm.scaling            = "automatic";
    algorithm.derivatives        = "automatic";
    algorithm.nlp_tolerance      = 1.0e-10;
    algorithm.nlp_iter_max       = 500;
    algorithm.collocation_method = method;
    algorithm.mesh_refinement    = "manual";
    algorithm.print_level        = 0;

    if (psopt(solution, problem, algorithm) != 0) return r;
    r.x        = solution.get_states_in_phase(1);
    r.u        = solution.get_controls_in_phase(1);
    r.t        = solution.get_time_in_phase(1);
    r.lam      = solution.get_dual_costates_in_phase(1);
    r.terminal = solution.get_terminal_state_in_phase(1);
    r.H        = solution.get_dual_hamiltonian_in_phase(1);
    r.J        = solution.cost;
    r.ok       = true;
    return r;
}

} // namespace gauss_terminal_test


// The reported trajectory must reach t_f and satisfy the terminal events there.
//
// Before the fix the last stored node was at t = 0.9990736915 -- tau = 0.99814738,
// the largest Legendre-Gauss node on 40 points -- carrying x1 = 0.9999974274 and
// x2 = 0.0055527025. Both are the EXACT solution at that time, to ten figures; the
// trajectory was right everywhere it was reported and simply stopped early.
TEST(GaussTerminal, ReportedTrajectoryReachesTheFinalTime)
{
    const auto r = gauss_terminal_test::solve("Gauss", 40);
    ASSERT_TRUE(r.ok);
    const long M = r.t.cols();
    ASSERT_GT(M, 10);

    EXPECT_NEAR(r.t(0,0),   0.0, 1.0e-12) << "the phase starts at t0";
    EXPECT_NEAR(r.t(0,M-1), 1.0, 1.0e-12)
        << "the last reported node is t_f, not the last Gauss point";

    EXPECT_NEAR(r.x(0,M-1), 1.0, 1.0e-8) << "x1(t_f) is the event value";
    EXPECT_NEAR(r.x(1,M-1), 0.0, 1.0e-8) << "x2(t_f) is the event value";
    EXPECT_NEAR(r.J,        6.0, 1.0e-8);

    // Every array a caller may pair with the trajectory has to be the same width,
    // or a column loop over one runs off the end of another. solution_diagnostics
    // does exactly that.
    EXPECT_EQ(r.x.cols(),   M);
    EXPECT_EQ(r.u.cols(),   M);
    EXPECT_EQ(r.lam.cols(), M);

    // The costate at t_f is lambda(+1) and this problem's is the constant -12.
    EXPECT_NEAR(r.lam(0,M-1), -12.0, 1.0e-3);
    EXPECT_NEAR(r.lam(1,M-1),   6.0, 1.0e-3);

    // The dedicated accessor and the last column of the trajectory are the same point.
    ASSERT_EQ(r.terminal.rows(), 2);
    EXPECT_NEAR(r.terminal(0,0), r.x(0,M-1), 1.0e-14);
    EXPECT_NEAR(r.terminal(1,0), r.x(1,M-1), 1.0e-14);
}

// The same, on a three-interval hp mesh. The gap is larger there because the last
// interval carries fewer Gauss points: before the fix the reported trajectory ended
// at t = 0.9940434785 with x2 = 0.0355262482 against a required 0.
TEST(GaussTerminal, ReachesTheFinalTimeOnAnHpMesh)
{
    const auto r = gauss_terminal_test::solve("Gauss", 28, /*hp=*/true);
    ASSERT_TRUE(r.ok);
    const long M = r.t.cols();
    ASSERT_GT(M, 10);

    EXPECT_NEAR(r.t(0,M-1), 1.0, 1.0e-12);
    EXPECT_NEAR(r.x(0,M-1), 1.0, 1.0e-8);
    EXPECT_NEAR(r.x(1,M-1), 0.0, 1.0e-8);
    EXPECT_EQ(r.u.cols(),   M);
    EXPECT_EQ(r.lam.cols(), M);
}

// Radau has the same non-collocated terminal point and has always reported it. The
// two siblings must agree about what the trajectory covers, since a caller should
// not have to know which of them produced a solution in order to read its end.
TEST(GaussTerminal, GaussAndRadauAgreeOnWhatTheTrajectoryCovers)
{
    const auto g = gauss_terminal_test::solve("Gauss", 40);
    const auto d = gauss_terminal_test::solve("Radau", 40);
    ASSERT_TRUE(g.ok);
    ASSERT_TRUE(d.ok);

    EXPECT_NEAR(g.t(0,0),              d.t(0,0),              1.0e-12);
    EXPECT_NEAR(g.t(0,g.t.cols()-1),   d.t(0,d.t.cols()-1),   1.0e-12);
    EXPECT_NEAR(g.x(0,g.x.cols()-1),   d.x(0,d.x.cols()-1),   1.0e-8);
    EXPECT_NEAR(g.x(1,g.x.cols()-1),   d.x(1,d.x.cols()-1),   1.0e-8);
    EXPECT_NEAR(g.J,                   d.J,                   1.0e-8);

    // Radau stores its terminal node, so the accessor reports nothing extra for it.
    const auto dr = d;
    EXPECT_EQ(dr.terminal.size(), 0)
        << "get_terminal_state_in_phase is empty when the last stored node IS the terminal";
}


// The Hamiltonian at each Gauss interval's left breakpoint.
//
// H = L + lambda^T f is formed from solution.integrand_cost, and that array is filled
// by the objective evaluation from the raw NLP controls. At a Gauss breakpoint the
// control is a variable that appears in no defect and no quadrature weight, so the
// barrier alone decides it and it comes back at the midpoint of its bounds. The
// reported control is corrected for that -- it is replaced by the interval's own
// interpolant -- but the running cost was left holding the artefact's value, because
// ff_num re-evaluates the objective from the decision vector immediately afterwards.
//
// This problem is autonomous with a fixed final time, so H is constant, and its exact
// value is -18. The control is bounded in [-50,50], so the artefact is 0 and the
// running cost u^2/2 with it: the breakpoint reported H = -36 while every other node
// reported -18, and a Hamiltonian that is constant, and is, read as though it were not.
//
// The tolerance is 1e-6 against a defect of 18.
TEST(GaussTerminal, HamiltonianIsConstantIncludingAtTheBreakpoints)
{
    const auto r = gauss_terminal_test::solve("Gauss", 40);
    ASSERT_TRUE(r.ok);
    const long M = r.t.cols();
    ASSERT_EQ(r.H.cols(), M) << "the Hamiltonian must span the same nodes as the trajectory";

    for (long k = 0; k < M; k++)
        EXPECT_NEAR(r.H(0,k), -18.0, 1.0e-6)
            << "H at t = " << r.t(0,k) << " (node " << k << " of " << M << ")";
}

// The same on a three-interval hp mesh, where there are three breakpoints rather than
// one and each carries its own artefact.
TEST(GaussTerminal, HamiltonianIsConstantOnAnHpMesh)
{
    const auto r = gauss_terminal_test::solve("Gauss", 28, /*hp=*/true);
    ASSERT_TRUE(r.ok);
    const long M = r.t.cols();
    ASSERT_EQ(r.H.cols(), M);

    for (long k = 0; k < M; k++)
        EXPECT_NEAR(r.H(0,k), -18.0, 1.0e-6)
            << "H at t = " << r.t(0,k) << " (node " << k << " of " << M << ")";
}

// The schemes that collocate their first node never had the problem, and must not
// acquire one.
TEST(GaussTerminal, HamiltonianIsConstantUnderTheOtherSchemes)
{
    for (const std::string m : { std::string("Legendre"), std::string("Radau"),
                                 std::string("Hermite-Simpson") }) {
        const auto r = gauss_terminal_test::solve(m, 40);
        ASSERT_TRUE(r.ok) << m;
        for (long k = 0; k < r.t.cols(); k++)
            EXPECT_NEAR(r.H(0,k), -18.0, 1.0e-5) << m << " at t = " << r.t(0,k);
    }
}
