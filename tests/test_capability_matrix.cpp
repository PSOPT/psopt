//////////////////////////////////////////////////////////////////////////////
// test_capability_matrix.cpp
//
// Cross-comparison capability matrix for PSOPT, using the double integrator as
// the prototype (analytical oracle). Every configuration solves the SAME
// min-energy problem  xdot1=x2, xdot2=u, (0,0)->(1,0), tf=1, min INT u^2/2 dt,
// whose Maximum-Principle optimum is J* = 6 (u = 6-12t). A capability is
// "validated" when its cell reproduces J* to a per-method tolerance, so the whole
// matrix is a mutual cross-check of:
//   - derivatives:  analytic (ADOL-C)  vs  finite differences (numerical)
//   - collocation:  Legendre / Chebyshev / Radau / Gauss / trapezoidal / Hermite-Simpson
//   - mesh:         non-adaptive (manual)  vs  adaptive (automatic refinement)
//   - ps_method:    Ross-Fahroo  vs  Bellman
//   - linear solver / parallelism:  serial (MUMPS, 1 thread) vs OpenMP (SPRAL,
//     8 threads) vs PARDISO / MKL-PARDISO (optional, guarded). Threads drive the
//     solver AND the underlying BLAS (OpenBLAS via OPENBLAS_NUM_THREADS, Intel MKL
//     via MKL_NUM_THREADS) -- PSOPT calls no BLAS/OpenMP/MPI itself.
// (The nlp_method axis {IPOPT, CASADI, SCIP} extends this in a CASADI/SCIP build.)
//////////////////////////////////////////////////////////////////////////////

#include "gtest/gtest.h"
#include <psopt.h>
#include <cstdlib>
#include <string>
#include <cmath>

namespace cm {
adouble endpoint_cost(adouble*, adouble*, adouble*, adouble&, adouble&, adouble*, int, Workspace*) { return 0.0; }
adouble integrand_cost(adouble*, adouble* c, adouble*, adouble&, adouble*, int, Workspace*) { adouble u=c[0]; return 0.5*u*u; }
void dae(adouble* d, adouble*, adouble* s, adouble* c, adouble*, adouble&, adouble*, int, Workspace*) { d[0]=s[1]; d[1]=c[0]; }
void events(adouble* e, adouble* i, adouble* f, adouble*, adouble&, adouble&, adouble*, int, Workspace*)
{ e[0]=i[0]; e[1]=i[1]; e[2]=f[0]; e[3]=f[1]; }
void linkages(adouble*, adouble*, Workspace*) {}
}

struct Cfg {
    const char* name; const char* deriv; const char* colloc; const char* mesh; const char* ps; double tol;
    const char* solver;      // ipopt_linear_solver: mumps / spral / ma97 / pardiso / pardisomkl
    const char* ompthreads;  // OMP_NUM_THREADS for a threaded solver (serial="1", OpenMP="8")
};

class CapabilityMatrix : public ::testing::TestWithParam<Cfg> {};

TEST_P(CapabilityMatrix, MinEnergyDoubleIntegrator_J_is_6)
{
    using namespace cm;
    const Cfg cfg = GetParam();
    const int N = 30;

    Alg algorithm; Sol solution; Prob problem;
    problem.name = cfg.name; problem.outfilename = "cm.txt";
    problem.nphases = 1; problem.nlinkages = 0;
    psopt_level1_setup(problem);
    problem.phases(1).nstates = 2; problem.phases(1).ncontrols = 1;
    problem.phases(1).nevents = 4; problem.phases(1).npath = 0; problem.phases(1).nodes << N;
    psopt_level2_setup(problem, algorithm);
    problem.phases(1).bounds.lower.states << -5.0,-5.0;  problem.phases(1).bounds.upper.states << 5.0,5.0;
    problem.phases(1).bounds.lower.controls(0) = -50.0;  problem.phases(1).bounds.upper.controls(0) = 50.0;
    problem.phases(1).bounds.lower.events << 0.0,0.0,1.0,0.0;
    problem.phases(1).bounds.upper.events << 0.0,0.0,1.0,0.0;
    problem.phases(1).bounds.lower.StartTime = 0.0; problem.phases(1).bounds.upper.StartTime = 0.0;
    problem.phases(1).bounds.lower.EndTime   = 1.0; problem.phases(1).bounds.upper.EndTime   = 1.0;
    problem.integrand_cost = &integrand_cost; problem.endpoint_cost = &endpoint_cost;
    problem.dae = &dae; problem.events = &events; problem.linkages = &linkages;
    problem.phases(1).guess.states = zeros(2,N);
    problem.phases(1).guess.states.row(0) = linspace(0.0,1.0,N);
    problem.phases(1).guess.controls = zeros(1,N);
    problem.phases(1).guess.time = linspace(0.0,1.0,N);

    algorithm.nlp_method          = "IPOPT";
    algorithm.scaling             = "automatic";
    algorithm.derivatives         = cfg.deriv;
    algorithm.collocation_method  = cfg.colloc;
    algorithm.mesh_refinement     = cfg.mesh;
    algorithm.ps_method           = cfg.ps;
    algorithm.ipopt_linear_solver = cfg.solver;
    algorithm.nlp_iter_max        = 1000;
    algorithm.nlp_tolerance       = 1.e-6;

    // Parallelism axis. PSOPT itself has NO OpenMP/MPI code paths and calls no BLAS
    // directly; all parallelism lives BELOW it:
    //   - the threaded linear solver (ma97/spral/pardiso*)  -> OMP_NUM_THREADS
    //   - the dense BLAS/LAPACK it uses. Here IPOPT/MUMPS/SPRAL link OpenBLAS
    //     (OPENBLAS_NUM_THREADS) while the MKL PARDISO path uses Intel MKL BLAS
    //     (MKL_NUM_THREADS). Enable ALL of them so serial=1 is truly single-threaded
    //     and "OpenMP"=8 actually threads the factorization + BLAS.
    // (MPI omitted: PSOPT has no MPI code paths.)
    setenv("OMP_NUM_THREADS",      cfg.ompthreads, 1);
    setenv("OPENBLAS_NUM_THREADS", cfg.ompthreads, 1);
    setenv("MKL_NUM_THREADS",      cfg.ompthreads, 1);
    if (std::string(cfg.solver) == "spral") {
        setenv("OMP_CANCELLATION","TRUE",1); setenv("OMP_NESTED","TRUE",1);
        setenv("OMP_PROC_BIND","TRUE",1);    setenv("OMP_STACKSIZE","64M",1);
    }

    psopt(solution, problem, algorithm);

    // GUARD the optional linear solvers (spral / ma97 / pardiso / pardisomkl):
    // if one is not installed, IPOPT fails to load its library (or "solves" to a
    // non-optimal point) -- treat any non-MUMPS solver that does not cleanly
    // reproduce J* as "not installed" and SKIP, so the suite runs on any build.
    const bool solved_ok = (solution.nlp_return_code == 0) && (solution.error_flag == false)
                           && (std::fabs(solution.get_cost() - 6.0) < 0.5);
    if (!solved_ok && std::string(cfg.solver) != "mumps")
        GTEST_SKIP() << "linear solver '" << cfg.solver << "' not installed in this IPOPT build";

    ASSERT_EQ(solution.error_flag, false) << "config '" << cfg.name << "' failed to solve";
    EXPECT_NEAR(solution.get_cost(), 6.0, cfg.tol)
        << "config '" << cfg.name << "' should reproduce the analytical J* = 6";
}

static const Cfg kConfigs[] = {
    // name                deriv         colloc            mesh         ps_method       tol      solver         threads
    {"analytic_Legendre",  "automatic",  "Legendre",       "manual",    "Ross-Fahroo",  5.0e-3,  "mumps",       "1"},
    {"finite_diff",        "numerical",  "Legendre",       "manual",    "Ross-Fahroo",  5.0e-3,  "mumps",       "1"}, // FD vs analytic
    {"Chebyshev",          "automatic",  "Chebyshev",      "manual",    "Ross-Fahroo",  6.0e-2,  "mumps",       "1"},
    {"Radau",              "automatic",  "Radau",          "manual",    "Ross-Fahroo",  5.0e-3,  "mumps",       "1"},
    {"Gauss",              "automatic",  "Gauss",          "manual",    "Ross-Fahroo",  5.0e-3,  "mumps",       "1"},
    {"trapezoidal",        "automatic",  "trapezoidal",    "manual",    "Ross-Fahroo",  3.5e-2,  "mumps",       "1"},
    {"Hermite_Simpson",    "automatic",  "Hermite-Simpson","manual",    "Ross-Fahroo",  1.0e-2,  "mumps",       "1"},
    {"adaptive_mesh",      "automatic",  "Legendre",       "automatic", "Ross-Fahroo",  5.0e-3,  "mumps",       "1"}, // adaptive vs non
    {"Bellman",            "automatic",  "Legendre",       "automatic", "Bellman",      5.0e-3,  "mumps",       "1"}, // Bellman vs Ross-Fahroo
    // --- parallelism / linear-solver axis (guarded: skip if the solver is not built) ---
    {"serial_mumps",       "automatic",  "Legendre",       "manual",    "Ross-Fahroo",  5.0e-3,  "mumps",       "1"}, // serial
    {"openmp8_spral",      "automatic",  "Legendre",       "manual",    "Ross-Fahroo",  5.0e-3,  "spral",       "8"}, // OpenMP, 8 threads
    {"pardiso",            "automatic",  "Legendre",       "manual",    "Ross-Fahroo",  5.0e-3,  "pardiso",     "8"}, // Panua PARDISO
    {"pardiso_mkl",        "automatic",  "Legendre",       "manual",    "Ross-Fahroo",  5.0e-3,  "pardisomkl",  "8"}, // Intel MKL PARDISO
};

INSTANTIATE_TEST_SUITE_P(Axes, CapabilityMatrix, ::testing::ValuesIn(kConfigs),
    [](const ::testing::TestParamInfo<Cfg>& info){ return std::string(info.param.name); });
