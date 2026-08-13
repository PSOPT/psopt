//////////////////////////////////////////////////////////////////////////////
// test_parameter_covariance.cpp
//
// The constrained least-squares covariance of the estimated parameters,
//
//     C_z = sigma_hat^2 Z ( Z^T Jr^T Jr Z )^{-1} Z^T,
//
// with Z a basis for the null space of the active constraints (Kostina et al.,
// 2003). Every part of that expression has been wrong at some point in this
// code's history -- the inverse omitted, a constraint row dropped, the null space
// built without the active bounds, the covariance left in scaled units, the
// variance factor absent -- and each error is invisible in the sense that the
// numbers still look like a covariance. This test pins the result of the whole
// calculation against values obtained independently of PSOPT.
//
// The problem is the COPS catalytic cracking of gas oil, as distributed in
// examples/cracking:
//
//     y1' = -(theta1 + theta3) y1^2,   y2' = theta1 y1^2 - theta2 y2,
//     y(0) = (1,0),   t in [0,1],
//
// fitted to 21 samples of both states. Solving it outside PSOPT by nonlinear
// least squares on the reduced parameter space, with an accurate ODE solver,
// gives
//
//     theta1 = 11.847 in [11.186, 12.507]
//     theta2 =  8.344 in [ 7.722,  8.967]
//     theta3 =  1.001 in [ 0.295,  1.708]
//     correlations 0.786, -0.844, -0.870,    sigma_hat = 1.1588e-02
//
// The tolerances below are loose enough to absorb the transcription error of the
// 80 node mesh, which is what separates PSOPT's objective from the COPS reference
// value, and tight enough that any of the five defects above would fail: each of
// them moved the half-widths by a factor of 5 or more, or the covariance by 30.
//////////////////////////////////////////////////////////////////////////////

#include "gtest/gtest.h"
#include <psopt.h>
#include <cmath>

namespace cracking_cov {

static const double t_meas[21] = {0.0,0.025,0.05,0.075,0.1,0.125,0.15,0.175,0.2,
                                  0.225,0.25,0.3,0.35,0.4,0.45,0.5,0.55,0.65,0.75,
                                  0.85,0.95};
static const double y1_meas[21] = {1.0,0.8105,0.6208,0.5258,0.4345,0.3903,0.3342,
                                   0.3034,0.2735,0.2405,0.2283,0.2071,0.1669,0.153,
                                   0.1339,0.1265,0.12,0.099,0.087,0.077,0.069};
static const double y2_meas[21] = {0.0,0.2,0.2886,0.301,0.3215,0.3123,0.2716,0.2551,
                                   0.2258,0.1959,0.1789,0.1457,0.1198,0.0909,0.0719,
                                   0.0561,0.046,0.028,0.019,0.014,0.01};

void observation_function(adouble* observations, adouble* states, adouble* controls,
                          adouble* parameters, adouble& time, int k, adouble* xad,
                          int iphase, Workspace* w)
{ observations[0] = states[0]; observations[1] = states[1]; }

adouble endpoint_cost(adouble* i, adouble* f, adouble* p, adouble& t0, adouble& tf,
                      adouble* xad, int iphase, Workspace* w)
{ return 0.0; }

adouble integrand_cost(adouble* s, adouble* c, adouble* p, adouble& t, adouble* xad,
                       int iphase, Workspace* w)
{ return 0.0; }

void dae(adouble* d, adouble* path, adouble* states, adouble* controls, adouble* p,
         adouble& time, adouble* xad, int iphase, Workspace* w)
{
    adouble y1 = states[0], y2 = states[1];
    adouble th1 = p[0], th2 = p[1], th3 = p[2];

    d[0] = -(th1 + th3)*y1*y1;
    d[1] =  th1*y1*y1 - th2*y2;
}

void events(adouble* e, adouble* i, adouble* f, adouble* p, adouble& t0, adouble& tf,
            adouble* xad, int iphase, Workspace* w)
{ e[0] = i[0]; e[1] = i[1]; }          // the published problem fixes y(0) = (1,0)

void linkages(adouble* linkages, adouble* xad, Workspace* w) { }

} // namespace cracking_cov

TEST(ParameterStatistics, CatalyticCrackingConfidenceLimits)
{
    using namespace cracking_cov;

    Alg algorithm; Sol solution; Prob problem;

    problem.name        = "Catalytic cracking (covariance validation)";
    problem.outfilename = "test_parameter_covariance.txt";
    problem.nphases     = 1;
    problem.nlinkages   = 0;
    psopt_level1_setup(problem);

    problem.phases(1).nstates     = 2;
    problem.phases(1).ncontrols   = 0;
    problem.phases(1).nevents     = 2;
    problem.phases(1).npath       = 0;
    problem.phases(1).nparameters = 3;
    problem.phases(1).nodes       << 80;
    problem.phases(1).nobserved   = 2;
    problem.phases(1).nsamples    = 21;
    psopt_level2_setup(problem, algorithm);

    MatrixXd tm(1,21), obs(2,21);
    for (int k = 0; k < 21; k++) {
        tm(0,k)  = t_meas[k];
        obs(0,k) = y1_meas[k];
        obs(1,k) = y2_meas[k];
    }
    problem.phases(1).observation_nodes = tm;
    problem.phases(1).observations      = obs;
    problem.phases(1).residual_weights  = ones(2,21);

    problem.phases(1).bounds.lower.states     << 0.0, 0.0;
    problem.phases(1).bounds.upper.states     << 2.0, 2.0;
    problem.phases(1).bounds.lower.parameters << 0.0, 0.0, 0.0;
    problem.phases(1).bounds.upper.parameters << 20.0, 20.0, 20.0;
    problem.phases(1).bounds.lower.events     << 1.0, 0.0;
    problem.phases(1).bounds.upper.events     << 1.0, 0.0;
    problem.phases(1).bounds.lower.StartTime  = 0.0;
    problem.phases(1).bounds.upper.StartTime  = 0.0;
    problem.phases(1).bounds.lower.EndTime    = 1.0;
    problem.phases(1).bounds.upper.EndTime    = 1.0;

    problem.integrand_cost      = &integrand_cost;
    problem.endpoint_cost       = &endpoint_cost;
    problem.dae                 = &dae;
    problem.events              = &events;
    problem.linkages            = &linkages;
    problem.observation_function = &observation_function;

    MatrixXd x_guess = zeros(2,80);
    x_guess.row(0)   = ones(1,80);
    problem.phases(1).guess.states     = x_guess;
    problem.phases(1).guess.time       = linspace(0.0,1.0,80);
    problem.phases(1).guess.parameters = zeros(3,1);
    problem.phases(1).guess.parameters(0) = 12.0;
    problem.phases(1).guess.parameters(1) =  8.0;
    problem.phases(1).guess.parameters(2) =  1.0;

    algorithm.nlp_method           = "IPOPT";
    algorithm.scaling              = "automatic";
    algorithm.derivatives          = "automatic";
    algorithm.nlp_iter_max         = 1000;
    algorithm.nlp_tolerance        = 1.e-6;
    algorithm.collocation_method   = "Hermite-Simpson";
    algorithm.parameter_statistics = "yes";
    algorithm.print_level          = 0;

    psopt(solution, problem, algorithm);

    ASSERT_EQ(solution.error_flag, 0);
    ASSERT_TRUE(solution.parameter_statistics_ok);

    // The estimates themselves, against the independent least-squares solution.
    MatrixXd p = solution.get_parameters_in_phase(1);
    EXPECT_NEAR(p(0), 11.847, 0.02);
    EXPECT_NEAR(p(1),  8.344, 0.02);
    EXPECT_NEAR(p(2),  1.001, 0.02);

    // Degrees of freedom and the number of scalar observations: N_s counts
    // measurements, not sampling instants, and n_f is n_z - n_c, not the number of
    // parameters -- though here the two coincide.
    EXPECT_EQ(solution.n_observations, 42);
    EXPECT_EQ(solution.parameter_dof,   3);

    // The residual variance factor. Without it the intervals are wrong by a factor
    // of 1/sigma_hat = 86 for this problem.
    EXPECT_NEAR(solution.sigma_hat, 1.1588e-02, 5.0e-4);

    // The confidence limits. The half-widths below are the quantity that every one
    // of the historical defects got wrong.
    const MatrixXd& lo = solution.parameter_confidence_low;
    const MatrixXd& hi = solution.parameter_confidence_high;
    ASSERT_EQ(lo.size(), 3);
    ASSERT_EQ(hi.size(), 3);

    const double half_ref[3] = {0.660, 0.622, 0.707};
    for (int i = 0; i < 3; i++) {
        double half = 0.5*(hi(i) - lo(i));
        EXPECT_GT(half, 0.0) << "parameter " << i;
        EXPECT_NEAR(half, half_ref[i], 0.06) << "parameter " << i;   // within ~10 per cent
        EXPECT_LT(lo(i), p(i));
        EXPECT_GT(hi(i), p(i));
    }

    // The covariance is in physical units, is symmetric, and has full rank; a
    // covariance computed in scaled variables would be out by the square of the
    // parameter scale factors, which are 1/20 here.
    const MatrixXd& C = solution.parameter_covariance;
    ASSERT_EQ(C.rows(), 3);
    ASSERT_EQ(C.cols(), 3);
    for (int i = 0; i < 3; i++) {
        EXPECT_GT(C(i,i), 0.0);
        for (int j = 0; j < 3; j++) EXPECT_NEAR(C(i,j), C(j,i), 1.0e-12);
    }
    EXPECT_EQ(C.fullPivLu().rank(), 3);

    // The correlations, which are what tell a reader whether the parameters are
    // separately identifiable. These are strong, and they are the reason the book
    // warns against reading the diagonal on its own.
    double r12 = C(0,1)/std::sqrt(C(0,0)*C(1,1));
    double r13 = C(0,2)/std::sqrt(C(0,0)*C(2,2));
    double r23 = C(1,2)/std::sqrt(C(1,1)*C(2,2));
    EXPECT_NEAR(r12,  0.786, 0.02);
    EXPECT_NEAR(r13, -0.844, 0.02);
    EXPECT_NEAR(r23, -0.870, 0.02);
}
