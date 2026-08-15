//////////////////////////////////////////////////////////////////////////////
// test_sqp.cpp
//
// PSOPT's own sequential quadratic programming solver, checked against a closed-form
// optimum and against IPOPT on the same transcription.
//
// The solver is dense at this stage, so the problems here are small on purpose. What
// is under test is the algorithm and its wiring into PSOPT rather than its speed:
// that the constraint bounds are handed over correctly, that the multiplier sign
// convention is right, and that the merit function and the quasi-Newton update
// converge to the same point PSOPT already reaches by an entirely different method.
//
// The problem is the linear-quadratic regulator
//
//     minimise  int_0^1 ( x^2 + u^2 ) dt      subject to  x' = u,  x(0) = 1,  x(1) = 0.75.
//
// Its Euler-Lagrange equation is x'' = x, so on a segment of length T joining x = a to
// x = b the optimal cost is [(a^2+b^2) cosh T - 2ab]/sinh T, here 0.775240441234. It is
// solved twice: once with the control bound inactive, where the closed form applies,
// and once with |u| <= 0.4, which the unconstrained solution violates near t = 0, so
// that the active-set machinery of the QP subproblem is exercised as well.
//
// The tests are skipped when PSOPT is built without WITH_SQP.
//////////////////////////////////////////////////////////////////////////////

#include "gtest/gtest.h"
#include <psopt.h>
#include <cmath>
#include <string>

#ifdef USE_SQP
#include <qpOASES.hpp>
#endif

#ifdef USE_PROXQP
#include <proxsuite/proxqp/sparse/sparse.hpp>
#endif

#ifdef USE_QPALM
#include <qpalm.hpp>
#include <qpalm/constants.h>
#endif

namespace sqp_test {

adouble endpoint_cost(adouble* i0, adouble* xf, adouble* p, adouble& t0, adouble& tf,
                      adouble* xad, int iphase, Workspace* w)
{ return 0.0; }

adouble integrand_cost(adouble* s, adouble* u, adouble* p, adouble& t,
                       adouble* xad, int iphase, Workspace* w)
{ return s[0]*s[0] + u[0]*u[0]; }

void dae(adouble* d, adouble* path, adouble* states, adouble* controls, adouble* p,
         adouble& time, adouble* xad, int iphase, Workspace* w)
{ d[0] = controls[0]; }

void events(adouble* e, adouble* i0, adouble* xf, adouble* p, adouble& t0,
            adouble& tf, adouble* xad, int iphase, Workspace* w)
{ e[0] = i0[0]; e[1] = xf[0]; }

void linkages(adouble* linkages, adouble* xad, Workspace* w) { }

// The Jacobian evaluations of the last solve, which the SQP performs once per accepted
// step and so counts its iterations.
static int last_nlp_iterations = 0;

// Solve with the requested NLP method and control bound; return the objective.
static double solve_lq(const std::string& nlp_method, double u_bound, int& error_flag,
                       const std::string& hessian   = "limited-memory",
                       const std::string& qp_solver = "qpOASES")
{
    Alg algorithm; Sol solution; Prob problem;
    const int nodes = 20;

    problem.name        = "Linear-quadratic regulator (SQP test)";
    problem.outfilename = "test_sqp.txt";
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
    problem.phases(1).bounds.lower.controls << -u_bound;
    problem.phases(1).bounds.upper.controls <<  u_bound;
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

    problem.phases(1).guess.states   = ones(1, nodes);
    problem.phases(1).guess.controls = zeros(1, nodes);
    problem.phases(1).guess.time     = linspace(0.0, 1.0, nodes);

    algorithm.nlp_method         = nlp_method;
    algorithm.hessian            = hessian;
    algorithm.qp_solver          = qp_solver;
    algorithm.scaling            = "automatic";
    algorithm.derivatives        = "automatic";
    algorithm.nlp_iter_max       = 500;
    algorithm.nlp_tolerance      = 1.e-8;
    algorithm.collocation_method = "Legendre";
    algorithm.print_level        = 0;

    psopt(solution, problem, algorithm);

    error_flag = solution.error_flag;
    last_nlp_iterations = solution.mesh_stats ? solution.mesh_stats[0].n_jacobian_evals : 0;
    return solution.cost;
}

} // namespace sqp_test


#ifdef USE_SQP

// The sign convention of qpOASES's multipliers, which SQP_interface converts to
// PSOPT's. It is the one thing in the interface that the documentation does not settle
// unambiguously, so it is pinned here against a QP whose answer is known: minimising
// 1/2 x'x subject to x1 + x2 = 2 gives x = (1,1), and in the convention
// grad f + A' lambda = 0 the multiplier is -1. qpOASES returns +1, which is why
// SQP_interface negates it.
TEST(SQPSolver, QpOasesDualSignConvention)
{
    USING_NAMESPACE_QPOASES

    real_t H[4]   = {1,0, 0,1};
    real_t gv[2]  = {0,0};
    real_t A[2]   = {1,1};
    real_t lbA[1] = {2}, ubA[1] = {2};
    real_t lb[2]  = {-1e20,-1e20}, ub[2] = {1e20,1e20};

    QProblem qp(2,1);
    Options o; o.printLevel = PL_NONE; qp.setOptions(o);
    int_t nWSR = 100;
    ASSERT_EQ(qp.init(H, gv, A, lb, ub, lbA, ubA, nWSR), SUCCESSFUL_RETURN);

    real_t x[2], y[3];
    qp.getPrimalSolution(x);
    qp.getDualSolution(y);

    EXPECT_NEAR(x[0], 1.0, 1e-10);
    EXPECT_NEAR(x[1], 1.0, 1e-10);
    EXPECT_NEAR(y[2], 1.0, 1e-8);

    const double lambda_psopt = -y[2];
    EXPECT_NEAR(x[0] + 1.0*lambda_psopt, 0.0, 1e-8);   // stationarity in PSOPT's sign
}

// qpOASES's own BLAS and LAPACK replacements must not be linked in. They are stand-ins
// for use when no BLAS is available, and if they are present they capture dgemm_ and
// dpotrf_ for the whole program, including for MUMPS inside IPOPT, which then fails
// inside its linear solver on problems it otherwise solves without difficulty. The
// build guards against this in CMakeLists.txt; this checks the guard held, because the
// symptom appears nowhere near the SQP code.
TEST(SQPSolver, RealBlasIsNotOverriddenByQpOases)
{
    int flag = -1;
    const double J = sqp_test::solve_lq("IPOPT", 10.0, flag);
    ASSERT_EQ(flag, 0) << "IPOPT failed; qpOASES's BLAS replacements are the first suspect";
    EXPECT_NEAR(J, 0.775240441234, 1.0e-9);
}

// The solver against the closed-form optimum, with no bound active.
TEST(SQPSolver, LinearQuadraticAgainstClosedForm)
{
    int flag = -1;
    const double J = sqp_test::solve_lq("SQP", 10.0, flag);

    ASSERT_EQ(flag, 0) << "the SQP solver failed";
    EXPECT_NEAR(J, 0.775240441234, 1.0e-9);
}

// The same problem with |u| <= 0.4, which the unconstrained optimal control violates
// near t = 0 (it starts at -0.675), so the bound is active over part of the horizon and
// the QP subproblem has a non-trivial active set. There is no closed form once the
// bound bites, so the comparison is against IPOPT: the two share nothing but the
// problem functions, one being a primal-dual interior-point method and the other an
// active-set SQP with a BFGS model.
TEST(SQPSolver, AgreesWithIpoptWithAnActiveControlBound)
{
    int flag_ipopt = -1, flag_sqp = -1;

    const double J_ipopt = sqp_test::solve_lq("IPOPT", 0.4, flag_ipopt);
    const double J_sqp   = sqp_test::solve_lq("SQP",   0.4, flag_sqp);

    ASSERT_EQ(flag_ipopt, 0) << "IPOPT failed on the reference problem";
    ASSERT_EQ(flag_sqp,   0) << "the SQP solver failed";

    // The bound must actually bite, or the test would prove nothing about the active set.
    EXPECT_GT(J_ipopt, 0.775240441234 + 1.0e-4);

    EXPECT_NEAR(J_sqp, J_ipopt, 1.0e-7*std::fabs(J_ipopt));
}

// The same problem with the exact Hessian of the Lagrangian in place of the
// quasi-Newton model. The closed form does not care which model was used to reach it,
// so this is a direct check that the sparse second derivatives, the convexification
// and the trust region that goes with them all agree on the same answer.
TEST(SQPSolver, ExactHessianReachesTheClosedForm)
{
    int flag = -1;
    const double J = sqp_test::solve_lq("SQP", 10.0, flag, "exact");

    ASSERT_EQ(flag, 0) << "the SQP solver failed with an exact Hessian";
    EXPECT_NEAR(J, 0.775240441234, 1.0e-9);
}

// With the control bound active the subproblem has a non-trivial working set, and the
// exact Hessian has to agree with IPOPT there too.
TEST(SQPSolver, ExactHessianAgreesWithIpoptWithAnActiveControlBound)
{
    int flag_ipopt = -1, flag_sqp = -1;

    const double J_ipopt = sqp_test::solve_lq("IPOPT", 0.4, flag_ipopt);
    const double J_sqp   = sqp_test::solve_lq("SQP",   0.4, flag_sqp, "exact");

    ASSERT_EQ(flag_ipopt, 0) << "IPOPT failed on the reference problem";
    ASSERT_EQ(flag_sqp,   0) << "the SQP solver failed with an exact Hessian";

    EXPECT_GT(J_ipopt, 0.775240441234 + 1.0e-4);      // the bound must actually bite
    EXPECT_NEAR(J_sqp, J_ipopt, 1.0e-7*std::fabs(J_ipopt));
}

// The exact Hessian is second-order information and should show it: on a problem whose
// Lagrangian has real curvature, it must reach the same answer in materially fewer
// iterations than the quasi-Newton model, which learns that curvature one step at a
// time. The margin is deliberately loose -- what is under test is that the second
// derivatives are being used at all, not the precise rate.
TEST(SQPSolver, ExactHessianCostsFewerIterationsThanBfgs)
{
    int flag = -1;
    (void) sqp_test::solve_lq("SQP", 0.4, flag, "limited-memory");
    ASSERT_EQ(flag, 0);
    const int it_bfgs = sqp_test::last_nlp_iterations;

    (void) sqp_test::solve_lq("SQP", 0.4, flag, "exact");
    ASSERT_EQ(flag, 0);
    const int it_exact = sqp_test::last_nlp_iterations;

    EXPECT_GT(it_bfgs, 0);
    EXPECT_GT(it_exact, 0);
    EXPECT_LT(it_exact, it_bfgs) << "exact " << it_exact << " vs BFGS " << it_bfgs;
}

#ifdef USE_PROXQP

// ProxQP's dual sign convention, pinned the same way qpOASES's is above and against the
// same QP: minimising 1/2 x'x subject to x1 + x2 = 2 gives x = (1,1) and, in the
// convention grad f + A' lambda = 0, a multiplier of -1. ProxQP states its own
// stationarity as H x + g + A' y + C' z = 0, so its y is already PSOPT's lambda and,
// unlike qpOASES's, must not be negated. Getting this backwards costs nothing in the
// primal solution and everything in the costates.
TEST(SQPSolver, ProxQpDualSignConvention)
{
    typedef long long Idx;
    typedef Eigen::SparseMatrix<double, Eigen::ColMajor, Idx> SpMat;

    SpMat H(2,2); H.insert(0,0) = 1.0; H.insert(1,1) = 1.0; H.makeCompressed();
    SpMat A(1,2); A.insert(0,0) = 1.0; A.insert(0,1) = 1.0; A.makeCompressed();
    SpMat C(0,2);

    Eigen::VectorXd g = Eigen::VectorXd::Zero(2), b(1), l(0), u(0);
    b << 2.0;

    proxsuite::proxqp::sparse::QP<double, Idx> qp(2, 1, 0);
    qp.settings.verbose = false;
    qp.settings.eps_abs = 1.0e-12;
    qp.init(H, g, A, b, C, l, u);
    qp.solve();

    EXPECT_NEAR(qp.results.x(0), 1.0, 1.0e-9);
    EXPECT_NEAR(qp.results.x(1), 1.0, 1.0e-9);

    const double lambda_psopt = qp.results.y(0);          // no negation
    EXPECT_NEAR(lambda_psopt, -1.0, 1.0e-8);
    EXPECT_NEAR(qp.results.x(0) + 1.0*lambda_psopt, 0.0, 1.0e-8);
}

// The two QP backends are different algorithms -- an active-set method and a proximal
// augmented-Lagrangian one -- reached through the same SQP. They must agree, with the
// quasi-Newton model and with the exact Hessian, and with a bound active so that the
// working set is not trivial.
TEST(SQPSolver, ProxQpAgreesWithQpOases)
{
    int f1 = -1, f2 = -1, f3 = -1, f4 = -1;

    const double a = sqp_test::solve_lq("SQP", 10.0, f1, "limited-memory", "qpOASES");
    const double b = sqp_test::solve_lq("SQP", 10.0, f2, "limited-memory", "ProxQP");
    const double c = sqp_test::solve_lq("SQP",  0.4, f3, "exact",          "qpOASES");
    const double d = sqp_test::solve_lq("SQP",  0.4, f4, "exact",          "ProxQP");

    ASSERT_EQ(f1, 0); ASSERT_EQ(f2, 0); ASSERT_EQ(f3, 0); ASSERT_EQ(f4, 0);

    EXPECT_NEAR(a, 0.775240441234, 1.0e-9);
    EXPECT_NEAR(b, 0.775240441234, 1.0e-9);
    EXPECT_NEAR(d, c, 1.0e-7*std::fabs(c));
}

#endif // USE_PROXQP

#ifdef USE_QPALM

// QPALM's dual sign convention, pinned against the same QP as the other two backends:
// minimising 1/2 x'x subject to x1 + x2 = 2 gives x = (1,1) and, in the convention
// grad f + A' lambda = 0, a multiplier of -1. QPALM states stationarity as
// Q x + q + A' y = 0, so its y is PSOPT's lambda unnegated, as ProxQP's is and as
// qpOASES's is not.
TEST(SQPSolver, QpalmDualSignConvention)
{
    using namespace qpalm;

    Data data(2, 1);
    sparse_mat_t Q(2,2); Q.insert(0,0) = 1.0; Q.insert(1,1) = 1.0; Q.makeCompressed();
    sparse_mat_t A(1,2); A.insert(0,0) = 1.0; A.insert(0,1) = 1.0; A.makeCompressed();
    data.set_Q(Q);
    data.set_A(A);
    data.q    = vec_t::Zero(2);
    data.c    = 0.0;
    data.bmin = vec_t::Constant(1, 2.0);
    data.bmax = vec_t::Constant(1, 2.0);

    Settings settings;
    settings.verbose = 0;
    settings.eps_abs = 1.0e-12;
    settings.eps_rel = 0.0;

    Solver solver(data, settings);
    solver.solve();

    ASSERT_EQ(solver.get_info().status_val, QPALM_SOLVED);
    EXPECT_NEAR(solver.get_solution().x(0), 1.0, 1.0e-8);
    EXPECT_NEAR(solver.get_solution().x(1), 1.0, 1.0e-8);

    const double lambda_psopt = solver.get_solution().y(0);      // no negation
    EXPECT_NEAR(lambda_psopt, -1.0, 1.0e-7);
    EXPECT_NEAR(solver.get_solution().x(0) + 1.0*lambda_psopt, 0.0, 1.0e-7);
}

// QPALM reached through the SQP, against the closed form and against qpOASES with a
// bound active. A third algorithm again -- proximal augmented Lagrangian with a sparse
// LDL underneath -- so agreement here is agreement between three different methods.
TEST(SQPSolver, QpalmAgreesWithQpOases)
{
    int f1 = -1, f2 = -1, f3 = -1;

    const double a = sqp_test::solve_lq("SQP", 10.0, f1, "limited-memory", "QPALM");
    const double b = sqp_test::solve_lq("SQP",  0.4, f2, "exact",          "QPALM");
    const double c = sqp_test::solve_lq("SQP",  0.4, f3, "exact",          "qpOASES");

    ASSERT_EQ(f1, 0); ASSERT_EQ(f2, 0); ASSERT_EQ(f3, 0);

    EXPECT_NEAR(a, 0.775240441234, 1.0e-9);
    EXPECT_NEAR(b, c, 1.0e-7*std::fabs(c));
}

#endif // USE_QPALM

#else

TEST(SQPSolver, DISABLED_NotBuilt)
{
    GTEST_SKIP() << "PSOPT was built without WITH_SQP";
}

#endif // USE_SQP
