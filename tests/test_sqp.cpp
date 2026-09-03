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
#include <vector>
#include <cstdlib>

#include "psopt_qp_plugin.h"

bool psopt_qp_plugin_available(const std::string& backend, std::string& message);
bool psopt_qp_plugin_solve(const std::string& backend,
                           const psopt_qp_problem* problem,
                           psopt_qp_solution*      solution,
                           std::string&            message);

#ifdef USE_SQP
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

// GALAHAD's default linear solver requires OMP_CANCELLATION and OMP_PROC_BIND to be
// TRUE in the environment, and cannot be made to work without them from inside the
// process: the OpenMP runtime reads them when it first initialises, long before a
// plugin is loaded. So whether GALAHAD can be exercised at all is a property of how
// this binary was started, and is asked at run time rather than assumed. CTest runs the
// suite a second time with those variables set, so the backend is covered.
bool galahad_environment_ok()
{
    const char* c = getenv("OMP_CANCELLATION");
    const char* b = getenv("OMP_PROC_BIND");
    return c != NULL && b != NULL && (c[0] == 'T' || c[0] == 't')
                                  && (b[0] != 'F' && b[0] != 'f');
}

// Which backend the tests use when they do not name one. qpOASES was this default while
// it existed, and asked nothing of its environment; every remaining backend is a plugin,
// and GALAHAD additionally needs OMP_CANCELLATION set before the process started. So the
// default is whichever backend this binary was built with, preferring one that carries
// no such condition, and tests that rely on it skip when there is none.
const char* default_backend()
{
#if   defined(USE_PROXQP)
    return "ProxQP";
#elif defined(USE_QPALM)
    return "QPALM";
#elif defined(USE_OSQP)
    return "OSQP";
#else
    return "GALAHAD";
#endif
}

bool default_backend_usable()
{
#if defined(USE_PROXQP) || defined(USE_QPALM) || defined(USE_OSQP)
    return true;
#else
    return galahad_environment_ok();
#endif
}

// Solve with the requested NLP method and control bound; return the objective.
static double solve_lq(const std::string& nlp_method, double u_bound, int& error_flag,
                       const std::string& hessian   = "limited-memory",
                       const std::string& qp_solver = default_backend())
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

// The sign convention the QP backends use for their multipliers, which SQP_interface
// converts to PSOPT's, is pinned against a QP whose answer is known by
// SQPSolver.GalahadPluginSolvesAndUsesTheRightDualSign below -- through the plugin ABI,
// which is now the only way a subproblem reaches a solver. Two earlier tests here did
// the same thing by linking qpOASES directly and checked that its BLAS replacements had
// not captured the real BLAS; both went with qpOASES itself.

// The solver against the closed-form optimum, with no bound active.
TEST(SQPSolver, LinearQuadraticAgainstClosedForm)
{
    // These use whichever backend solve_lq defaults to. When that is GALAHAD -- which
    // it is whenever GALAHAD is the only backend built -- it cannot run unless the
    // process was started with OMP_CANCELLATION set, so the test would be measuring the
    // environment rather than the solver. CTest's second pass sets it and covers this.
    if (!sqp_test::default_backend_usable())
        GTEST_SKIP() << "the default QP backend is GALAHAD, which needs "
                        "OMP_CANCELLATION=TRUE and OMP_PROC_BIND=TRUE in the environment";

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
    // These use whichever backend solve_lq defaults to. When that is GALAHAD -- which
    // it is whenever GALAHAD is the only backend built -- it cannot run unless the
    // process was started with OMP_CANCELLATION set, so the test would be measuring the
    // environment rather than the solver. CTest's second pass sets it and covers this.
    if (!sqp_test::default_backend_usable())
        GTEST_SKIP() << "the default QP backend is GALAHAD, which needs "
                        "OMP_CANCELLATION=TRUE and OMP_PROC_BIND=TRUE in the environment";

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
    // These use whichever backend solve_lq defaults to. When that is GALAHAD -- which
    // it is whenever GALAHAD is the only backend built -- it cannot run unless the
    // process was started with OMP_CANCELLATION set, so the test would be measuring the
    // environment rather than the solver. CTest's second pass sets it and covers this.
    if (!sqp_test::default_backend_usable())
        GTEST_SKIP() << "the default QP backend is GALAHAD, which needs "
                        "OMP_CANCELLATION=TRUE and OMP_PROC_BIND=TRUE in the environment";

    int flag = -1;
    const double J = sqp_test::solve_lq("SQP", 10.0, flag, "exact");

    ASSERT_EQ(flag, 0) << "the SQP solver failed with an exact Hessian";
    EXPECT_NEAR(J, 0.775240441234, 1.0e-9);
}

// With the control bound active the subproblem has a non-trivial working set, and the
// exact Hessian has to agree with IPOPT there too.
TEST(SQPSolver, ExactHessianAgreesWithIpoptWithAnActiveControlBound)
{
    // These use whichever backend solve_lq defaults to. When that is GALAHAD -- which
    // it is whenever GALAHAD is the only backend built -- it cannot run unless the
    // process was started with OMP_CANCELLATION set, so the test would be measuring the
    // environment rather than the solver. CTest's second pass sets it and covers this.
    if (!sqp_test::default_backend_usable())
        GTEST_SKIP() << "the default QP backend is GALAHAD, which needs "
                        "OMP_CANCELLATION=TRUE and OMP_PROC_BIND=TRUE in the environment";

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
    // These use whichever backend solve_lq defaults to. When that is GALAHAD -- which
    // it is whenever GALAHAD is the only backend built -- it cannot run unless the
    // process was started with OMP_CANCELLATION set, so the test would be measuring the
    // environment rather than the solver. CTest's second pass sets it and covers this.
    if (!sqp_test::default_backend_usable())
        GTEST_SKIP() << "the default QP backend is GALAHAD, which needs "
                        "OMP_CANCELLATION=TRUE and OMP_PROC_BIND=TRUE in the environment";

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

// The backends, reached through the SQP, all in one process.
//
// This is the test the plugin architecture exists for. Every one of these QP libraries
// carries a sparse factorisation, every factorisation carries an ordering, and the
// orderings share C symbol names: linked into one program, QPALM's LADEL exports
// amd_order built for 64-bit indices, another backend's copy is built for 32, the
// linker picks one, and the loser writes past the end of something. Before the backends
// were separated into plugins this test file crashed -- inside ProxQP, which had not
// changed, because QPALM had been linked beside it. Now each backend is a shared object
// loaded with RTLD_LOCAL, and what is inside one is invisible to the next.
//
// So the value of solving the same problem through each in turn is not that the answers
// agree, though they must; it is that the process survives doing so.
namespace {

struct Backend { const char* name; bool enabled; };

// Built at run time rather than declared as an array: with GALAHAD the only backend
// compiled in -- which is now a perfectly ordinary way to build PSOPT -- every entry
// below is conditioned out and an empty array is not valid C++.
std::vector<Backend> plugin_backends()
{
    std::vector<Backend> v;
#ifdef USE_PROXQP
    { Backend b = { "ProxQP", true }; v.push_back(b); }
#endif
#ifdef USE_QPALM
    { Backend b = { "QPALM",  true }; v.push_back(b); }
#endif
#ifdef USE_OSQP
    { Backend b = { "OSQP",   true }; v.push_back(b); }
#endif
    return v;
}

std::vector<Backend> enabled_backends()
{
    std::vector<Backend> out = plugin_backends();
#ifdef USE_GALAHAD
    if (sqp_test::galahad_environment_ok()) { Backend g = { "GALAHAD", true }; out.push_back(g); }
#endif
    return out;
};

} // namespace

TEST(SQPSolver, EveryBackendReachesTheClosedFormInOneProcess)
{
    const std::vector<Backend> use = enabled_backends();
    for (size_t k = 0; k < use.size(); k++) {
        int flag = -1;
        const double J = sqp_test::solve_lq("SQP", 10.0, flag, "limited-memory",
                                            use[k].name);
        EXPECT_EQ(flag, 0)          << "backend " << use[k].name;
        EXPECT_NEAR(J, 0.775240441234, 1.0e-9) << "backend " << use[k].name;
    }
}

// The same, with the exact sparse Hessian and a bound active, so that the subproblem
// has a non-trivial working set and an indefinite model to cope with. The backends
// differ in kind -- an active-set method and two proximal augmented-Lagrangian methods
// over different factorisations -- so agreement here is agreement between three
// algorithms sharing nothing but the problem.
TEST(SQPSolver, EveryBackendAgreesWithAnActiveBoundInOneProcess)
{
    int flag_ref = -1;
    const double J_ref = sqp_test::solve_lq("IPOPT", 0.4, flag_ref);
    ASSERT_EQ(flag_ref, 0);
    ASSERT_GT(J_ref, 0.775240441234 + 1.0e-4);      // the bound must actually bite

    const std::vector<Backend> use = enabled_backends();
    for (size_t k = 0; k < use.size(); k++) {
        int flag = -1;
        const double J = sqp_test::solve_lq("SQP", 0.4, flag, "exact", use[k].name);
        EXPECT_EQ(flag, 0) << "backend " << use[k].name;
        EXPECT_NEAR(J, J_ref, 1.0e-6*std::fabs(J_ref)) << "backend " << use[k].name;
    }
}

#ifdef USE_GALAHAD

// GALAHAD's QPA through the plugin, on the QP every backend's dual convention is pinned
// against: minimising 1/2 x'x subject to x1 + x2 = 2 gives x = (1,1) and, in the
// convention grad f + A' lambda = 0, a multiplier of -1. QPA states its own
// stationarity as H x + g - A' y - z = 0, so its y carries the opposite sign to PSOPT's
// and the plugin negates it, as the removed active-set backend also required and
// ProxQP and QPALM do not.
//
// This also exercises the whole plugin path end to end: the loader, the ABI, the
// coordinate-format conversion and the bound handling. Since the plugin ABI is now the
// only route from the driver to a solver, this is where that convention is pinned.
TEST(SQPSolver, GalahadPluginSolvesAndUsesTheRightDualSign)
{
    if (!sqp_test::galahad_environment_ok())
        GTEST_SKIP() << "GALAHAD needs OMP_CANCELLATION=TRUE and OMP_PROC_BIND=TRUE "
                        "in the environment; CTest runs a second pass with them set";

    const long long H_p[3] = {0, 1, 2};
    const long long H_i[2] = {0, 1};
    const double    H_x[2] = {1.0, 1.0};

    const long long A_p[3] = {0, 1, 2};
    const long long A_i[2] = {0, 0};
    const double    A_x[2] = {1.0, 1.0};

    const double g[2]   = {0.0, 0.0};
    const double lbA[1] = {2.0}, ubA[1] = {2.0};
    const double lb[2]  = {-PSOPT_QP_INFINITY, -PSOPT_QP_INFINITY};
    const double ub[2]  = { PSOPT_QP_INFINITY,  PSOPT_QP_INFINITY};

    psopt_qp_problem q;
    psopt_qp_problem_init(&q);
    q.n = 2; q.m = 1;
    q.H_p = H_p; q.H_i = H_i; q.H_x = H_x; q.H_dense = NULL;
    q.g = g;
    q.A_p = A_p; q.A_i = A_i; q.A_x = A_x;
    q.lbA = lbA; q.ubA = ubA; q.lb = lb; q.ub = ub;
    q.tolerance = 1.0e-10; q.max_iter = 200; q.nonconvex = 0;

    double d[2] = {0,0}, lambda[1] = {0}, z[2] = {0,0};
    psopt_qp_solution r;
    r.d = d; r.lambda = lambda; r.z = z; r.iterations = 0; r.status = -1;

    std::string message;
    ASSERT_TRUE(psopt_qp_plugin_solve("GALAHAD", &q, &r, message)) << message;
    ASSERT_EQ(r.status, PSOPT_QP_SOLVED);

    EXPECT_NEAR(d[0], 1.0, 1.0e-8);
    EXPECT_NEAR(d[1], 1.0, 1.0e-8);
    EXPECT_NEAR(lambda[0], -1.0, 1.0e-7);
    EXPECT_NEAR(d[0] + 1.0*lambda[0] - z[0], 0.0, 1.0e-7);   // PSOPT's stationarity
}

#endif // USE_GALAHAD

#ifdef USE_PIQP

// PIQP through the plugin, on the same QP, and then on one whose answer is set by an
// active inequality row. PIQP is the first backend to be given the equalities as
// equalities rather than as two-sided inequalities, so the row split is exercised here
// as well as the sign: the first problem's single row is an equality and the second's is
// not, and the multiplier has to come back in PSOPT's convention from either block.
TEST(SQPSolver, PiqpPluginSolvesAndUsesTheRightDualSign)
{
    const long long H_p[3] = {0, 1, 2};
    const long long H_i[2] = {0, 1};
    const double    H_x[2] = {1.0, 1.0};

    const long long A_p[3] = {0, 1, 2};
    const long long A_i[2] = {0, 0};
    const double    A_x[2] = {1.0, 1.0};

    const double g[2]   = {0.0, 0.0};
    const double lbA[1] = {2.0}, ubA[1] = {2.0};
    const double lb[2]  = {-PSOPT_QP_INFINITY, -PSOPT_QP_INFINITY};
    const double ub[2]  = { PSOPT_QP_INFINITY,  PSOPT_QP_INFINITY};

    psopt_qp_problem q;
    psopt_qp_problem_init(&q);
    q.n = 2; q.m = 1;
    q.H_p = H_p; q.H_i = H_i; q.H_x = H_x; q.H_dense = NULL;
    q.g = g;
    q.A_p = A_p; q.A_i = A_i; q.A_x = A_x;
    q.lbA = lbA; q.ubA = ubA; q.lb = lb; q.ub = ub;
    q.tolerance = 1.0e-10; q.max_iter = 200; q.nonconvex = 0;

    double d[2] = {0,0}, lambda[1] = {0}, z[2] = {0,0};
    psopt_qp_solution r;
    r.d = d; r.lambda = lambda; r.z = z; r.iterations = 0; r.status = -1;

    std::string message;
    ASSERT_TRUE(psopt_qp_plugin_solve("PIQP", &q, &r, message)) << message;
    ASSERT_EQ(r.status, PSOPT_QP_SOLVED);

    EXPECT_NEAR(d[0], 1.0, 1.0e-7);
    EXPECT_NEAR(d[1], 1.0, 1.0e-7);
    EXPECT_NEAR(lambda[0], -1.0, 1.0e-6);
    EXPECT_NEAR(d[0] + 1.0*lambda[0] - z[0], 0.0, 1.0e-6);   // PSOPT's stationarity

    // min 1/2 x'x - 3x1 - x2  s.t.  x1 + x2 <= 2. The unconstrained optimum is (3,1);
    // the row is active, the answer is (2,0) and the multiplier is +1.
    const double g2[2]   = {-3.0, -1.0};
    const double lbA2[1] = {-PSOPT_QP_INFINITY}, ubA2[1] = {2.0};
    q.g = g2; q.lbA = lbA2; q.ubA = ubA2;

    double d2[2] = {0,0}, lambda2[1] = {0}, z2[2] = {0,0};
    r.d = d2; r.lambda = lambda2; r.z = z2; r.iterations = 0; r.status = -1;

    ASSERT_TRUE(psopt_qp_plugin_solve("PIQP", &q, &r, message)) << message;
    ASSERT_EQ(r.status, PSOPT_QP_SOLVED);

    EXPECT_NEAR(d2[0], 2.0, 1.0e-6);
    EXPECT_NEAR(d2[1], 0.0, 1.0e-6);
    EXPECT_NEAR(lambda2[0], 1.0, 1.0e-6);
    EXPECT_NEAR(d2[0] + g2[0] + lambda2[0] - z2[0], 0.0, 1.0e-6);
    EXPECT_NEAR(d2[1] + g2[1] + lambda2[0] - z2[1], 0.0, 1.0e-6);
}

#endif // USE_PIQP

#ifdef USE_CLARABEL

// Clarabel through the plugin, on the same two QPs as PIQP. Clarabel is conic, so a
// two-sided row is split into one row per side and a simple bound into two more; the
// multiplier of an original row is therefore assembled here rather than passed through,
// and both halves of that assembly are exercised: the first problem's row is an equality
// and the second's is an active inequality.
TEST(SQPSolver, ClarabelPluginSolvesAndUsesTheRightDualSign)
{
    const long long H_p[3] = {0, 1, 2};
    const long long H_i[2] = {0, 1};
    const double    H_x[2] = {1.0, 1.0};

    const long long A_p[3] = {0, 1, 2};
    const long long A_i[2] = {0, 0};
    const double    A_x[2] = {1.0, 1.0};

    const double g[2]   = {0.0, 0.0};
    const double lbA[1] = {2.0}, ubA[1] = {2.0};
    const double lb[2]  = {-PSOPT_QP_INFINITY, -PSOPT_QP_INFINITY};
    const double ub[2]  = { PSOPT_QP_INFINITY,  PSOPT_QP_INFINITY};

    psopt_qp_problem q;
    psopt_qp_problem_init(&q);
    q.n = 2; q.m = 1;
    q.H_p = H_p; q.H_i = H_i; q.H_x = H_x; q.H_dense = NULL;
    q.g = g;
    q.A_p = A_p; q.A_i = A_i; q.A_x = A_x;
    q.lbA = lbA; q.ubA = ubA; q.lb = lb; q.ub = ub;
    q.tolerance = 1.0e-10; q.max_iter = 200; q.nonconvex = 0;

    double d[2] = {0,0}, lambda[1] = {0}, z[2] = {0,0};
    psopt_qp_solution r;
    r.d = d; r.lambda = lambda; r.z = z; r.iterations = 0; r.status = -1;

    std::string message;
    ASSERT_TRUE(psopt_qp_plugin_solve("Clarabel", &q, &r, message)) << message;
    ASSERT_EQ(r.status, PSOPT_QP_SOLVED);

    EXPECT_NEAR(d[0], 1.0, 1.0e-6);
    EXPECT_NEAR(d[1], 1.0, 1.0e-6);
    EXPECT_NEAR(lambda[0], -1.0, 1.0e-6);
    EXPECT_NEAR(d[0] + 1.0*lambda[0] - z[0], 0.0, 1.0e-6);   // PSOPT's stationarity

    // min 1/2 x'x - 3x1 - x2  s.t.  x1 + x2 <= 2, answer (2,0) and a multiplier of +1.
    const double g2[2]   = {-3.0, -1.0};
    const double lbA2[1] = {-PSOPT_QP_INFINITY}, ubA2[1] = {2.0};
    q.g = g2; q.lbA = lbA2; q.ubA = ubA2;

    double d2[2] = {0,0}, lambda2[1] = {0}, z2[2] = {0,0};
    r.d = d2; r.lambda = lambda2; r.z = z2; r.iterations = 0; r.status = -1;

    ASSERT_TRUE(psopt_qp_plugin_solve("Clarabel", &q, &r, message)) << message;
    ASSERT_EQ(r.status, PSOPT_QP_SOLVED);

    EXPECT_NEAR(d2[0], 2.0, 1.0e-5);
    EXPECT_NEAR(d2[1], 0.0, 1.0e-5);
    EXPECT_NEAR(lambda2[0], 1.0, 1.0e-5);
    EXPECT_NEAR(d2[0] + g2[0] + lambda2[0] - z2[0], 0.0, 1.0e-5);
    EXPECT_NEAR(d2[1] + g2[1] + lambda2[0] - z2[1], 0.0, 1.0e-5);
}

// The Euclidean trust region, on a subproblem whose answer is known in closed form.
//
//   min 1/2 d'd - 3 d1 - 4 d2   s.t.  ||d||_2 <= 1
//
// The unconstrained minimiser is (3,4), of length 5, so the region is active and the
// answer is the point of the unit circle in that direction: (0.6, 0.8), exactly. The
// same problem under a box of half-width 1 answers (1,1), of length sqrt(2) -- so this
// test would fail if the region were quietly being imposed componentwise, which is the
// mistake it exists to catch.
TEST(SQPSolver, ClarabelImposesTheEuclideanTrustRegion)
{
    const long long H_p[3] = {0, 1, 2};
    const long long H_i[2] = {0, 1};
    const double    H_x[2] = {1.0, 1.0};
    const double    g[2]   = {-3.0, -4.0};
    const double    lb[2]  = {-PSOPT_QP_INFINITY, -PSOPT_QP_INFINITY};
    const double    ub[2]  = { PSOPT_QP_INFINITY,  PSOPT_QP_INFINITY};

    psopt_qp_problem q;
    psopt_qp_problem_init(&q);
    q.n = 2; q.m = 0;
    q.H_p = H_p; q.H_i = H_i; q.H_x = H_x;
    q.g = g; q.lb = lb; q.ub = ub;
    q.tolerance = 1.0e-10; q.max_iter = 200;
    q.trust_radius = 1.0; q.trust_dim = 2;

    double d[2] = {0,0}, lambda[1] = {0}, z[2] = {0,0};
    psopt_qp_solution r;
    r.d = d; r.lambda = lambda; r.z = z; r.iterations = 0; r.status = -1;

    std::string message;
    ASSERT_TRUE(psopt_qp_plugin_solve("Clarabel", &q, &r, message)) << message;
    ASSERT_EQ(r.status, PSOPT_QP_SOLVED);

    EXPECT_NEAR(d[0], 0.6, 1.0e-6);
    EXPECT_NEAR(d[1], 0.8, 1.0e-6);
    EXPECT_NEAR(sqrt(d[0]*d[0] + d[1]*d[1]), 1.0, 1.0e-6);

    // The region's own multiplier is not reported as a bound multiplier. Were it added
    // to z, the stationarity residual below would vanish and the SQP would read a step
    // held back by the region as an optimal one.
    EXPECT_NEAR(z[0], 0.0, 1.0e-9);
    EXPECT_NEAR(z[1], 0.0, 1.0e-9);

    // Widen the region past the unconstrained minimiser and it stops binding.
    q.trust_radius = 10.0;
    r.status = -1;
    ASSERT_TRUE(psopt_qp_plugin_solve("Clarabel", &q, &r, message)) << message;
    ASSERT_EQ(r.status, PSOPT_QP_SOLVED);
    EXPECT_NEAR(d[0], 3.0, 1.0e-5);
    EXPECT_NEAR(d[1], 4.0, 1.0e-5);

    // A region over the leading variable only leaves the other one free.
    q.trust_radius = 1.0; q.trust_dim = 1;
    r.status = -1;
    ASSERT_TRUE(psopt_qp_plugin_solve("Clarabel", &q, &r, message)) << message;
    ASSERT_EQ(r.status, PSOPT_QP_SOLVED);
    EXPECT_NEAR(d[0], 1.0, 1.0e-5);
    EXPECT_NEAR(d[1], 4.0, 1.0e-5);
}

#endif // USE_CLARABEL

// A backend with no cones must refuse a subproblem carrying a Euclidean trust region
// rather than solve the one it can see and return the step as though it were restricted.
// This is the contract of psopt_qp_plugin.h and the reason the ABI version was raised;
// it is pinned here on whichever non-conic backend this build has.
TEST(SQPSolver, ANonConicBackendRefusesAEuclideanTrustRegion)
{
    const char* backend = NULL;
#if defined(USE_PIQP)
    backend = "PIQP";
#elif defined(USE_GALAHAD)
    backend = "GALAHAD";
#elif defined(USE_OSQP)
    backend = "OSQP";
#elif defined(USE_PROXQP)
    backend = "ProxQP";
#elif defined(USE_QPALM)
    backend = "QPALM";
#endif
    if (backend == NULL) GTEST_SKIP() << "no non-conic backend in this build";

    const long long H_p[3] = {0, 1, 2};
    const long long H_i[2] = {0, 1};
    const double    H_x[2] = {1.0, 1.0};
    const double    g[2]   = {-3.0, -4.0};
    const double    lb[2]  = {-1.0e3, -1.0e3};
    const double    ub[2]  = { 1.0e3,  1.0e3};

    psopt_qp_problem q;
    psopt_qp_problem_init(&q);
    q.n = 2; q.m = 0;
    q.H_p = H_p; q.H_i = H_i; q.H_x = H_x;
    q.g = g; q.lb = lb; q.ub = ub;
    q.tolerance = 1.0e-10; q.max_iter = 200;

    double d[2] = {0,0}, lambda[1] = {0}, z[2] = {0,0};
    psopt_qp_solution r;
    r.d = d; r.lambda = lambda; r.z = z; r.iterations = 0; r.status = -1;

    std::string message;

    // Without a region it solves this perfectly well, so a refusal below is about the
    // region and not about the problem.
    ASSERT_TRUE(psopt_qp_plugin_solve(backend, &q, &r, message)) << message;
    EXPECT_NE(r.status, PSOPT_QP_FAILED);

    q.trust_radius = 1.0; q.trust_dim = 2;
    r.status = -1;
    psopt_qp_plugin_solve(backend, &q, &r, message);
    EXPECT_EQ(r.status, PSOPT_QP_FAILED) << backend << " did not refuse a trust region it cannot impose";
}

// A plugin that cannot be found must be reported as such, at once and in words, rather
// than surfacing as a failed subproblem partway through a solve.
TEST(SQPSolver, AMissingPluginIsReportedClearly)
{
    std::string message;
    EXPECT_FALSE(psopt_qp_plugin_available("NoSuchBackend", message));
    EXPECT_NE(message.find("could not load"), std::string::npos) << message;
}

#else

TEST(SQPSolver, DISABLED_NotBuilt)
{
    GTEST_SKIP() << "PSOPT was built without WITH_SQP";
}

#endif // USE_SQP
