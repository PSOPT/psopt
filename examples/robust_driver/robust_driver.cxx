//////////////////////////////////////////////////////////////////////////
//////////////////      robust_driver.cxx       //////////////////////////
//////////////////////////////////////////////////////////////////////////
////////////////           PSOPT  Example             ////////////////////
//////////////////////////////////////////////////////////////////////////
//////// Title:   The two-link arm with an uncertain payload,        //////
////////          through psopt_solve_robust                         //////
//////// Last modified: 27 September 2026                            //////
//////// Reference:     the arm is the PROPT user's guide problem, as //////
////////                in examples/twolinkarm and examples/robust_arm //////
//////////////////////////////////////////////////////////////////////////
////////     Copyright (c) Victor M. Becerra, 2026         ///////////////
//////////////////////////////////////////////////////////////////////////
//////// This is part of the PSOPT software library, which ///////////////
//////// is distributed under the terms of the GNU Lesser ////////////////
//////// General Public License (LGPL)                    ////////////////
//////////////////////////////////////////////////////////////////////////
//
//  WHAT THIS EXAMPLE IS FOR
//
//  examples/robust_arm writes the whole robust design loop out by hand: the
//  augmented problem, the verification integrator, the worst-case scan, the
//  tightening, the warm start. It is the study, and it is where the method is
//  explained.
//
//  This is the same problem through psopt_solve_robust, and it exists to show
//  what the library now owns and what is still the user's. That division is not
//  the same as in the Python interface, and the difference is worth stating.
//
//    THE DRIVER'S:  the scenario rule (an unscented set that reproduces the mean
//                   and the covariance, and for one Gaussian parameter is exactly
//                   the three-point Gauss-Hermite rule); the search for the
//                   parameter a design serves worst; the generation loop with its
//                   warm start; and the diagnosis when the transcription runs out
//                   of degrees of freedom.
//
//    THE USER'S:    the augmented dae and events, and the independent integrator
//                   that says how badly a design serves a given plant.
//
//  In Python the driver builds the augmented problem itself, because the maths
//  arrives as CasADi expressions it can replicate. In C++ the maths arrives as
//  dae() and events() with fixed signatures, taped by CppAD, and nothing can
//  rewrite them -- so the augmentation stays with the user, where it is a loop of
//  three lines, and the library takes the parts that are the same for every
//  problem.
//
//  The verification integrator is the user's for a better reason than convenience:
//  it has to be INDEPENDENT of the transcription, and only the user knows what
//  that means for their problem. Here it shares arm_rhs with the dae and nothing
//  else -- not PSOPT's integrator, not its mesh, not its solution.
//
//////////////////////////////////////////////////////////////////////////

#include "psopt.h"
#include "robust.h"

#include <cmath>

using namespace PSOPT;

//////////////////////////////////////////////////////////////////////////
///////////////////  The problem  ////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////

static const double MU    = 0.50;     // mean payload
static const double SIGMA = 0.15;     // its standard deviation
static const double DELTA = 0.03;     // terminal tolerance
static const double SLACK = 1.0e-3;   // violation of THAT a design may leave
static const int    NODES = 25;

static const double X0[4] = { 0.0, 0.0, 0.500, 0.000 };
static const double XF[4] = { 0.0, 0.0, 0.500, 0.522 };

// The scenario list, read by the augmented dae and events through
// problem.user_data. This is the channel examples/climb uses for its aerodynamic
// tables and is the one any driver would use.
struct Scenarios {
    std::vector<double> mp;
};

// One scenario's dynamics. Templated so that the same source serves the adouble
// tape and the plain-double verification integrator below: a design checked
// against a DIFFERENT implementation of the plant is checking the implementation
// and not the design, and this is the compromise that keeps the physics in one
// place while keeping the integrators apart.
template <class T>
static void arm_rhs(const T* x, const T* u, double mp, T* dx)
{
    const T x1 = x[0], x2 = x[1], x3 = x[2];
    const double M11 = 7.0/3.0 + mp, M22 = 4.0/3.0 + mp, S2 = 3.0/2.0 + mp;
    const T c = cos(x3), s = sin(x3);
    const T M12 = S2*c;
    const T det = M11*M22 - M12*M12;
    const T b1 = (u[0] - u[1]) + S2*s*x2*x2;
    const T b2 =  u[1]         - S2*s*x1*x1;
    dx[0] = ( M22*b1 - M12*b2)/det;
    dx[1] = (-M12*b1 + M11*b2)/det;
    dx[2] = x2 - x1;
    dx[3] = x1;
}

adouble endpoint_cost(adouble* /*x0*/, adouble* /*xf*/, adouble* /*p*/,
                      adouble& /*t0*/, adouble& tf, adouble* /*xad*/,
                      int /*iphase*/, Workspace* /*workspace*/)
{
    return tf;                        // minimum time, shared by every scenario
}

adouble integrand_cost(adouble* /*x*/, adouble* /*u*/, adouble* /*p*/,
                       adouble& /*t*/, adouble* /*xad*/, int /*iphase*/,
                       Workspace* /*workspace*/)
{
    return 0.0;
}

// The augmentation, and it is three lines: M copies of the arm, decoupled from
// each other, meeting only in the single control vector above them and in t_f.
void dae(adouble* derivatives, adouble* /*path*/, adouble* states,
         adouble* controls, adouble* /*parameters*/, adouble& /*time*/,
         adouble* /*xad*/, int /*iphase*/, Workspace* workspace)
{
    const Scenarios& S = *((Scenarios*) workspace->problem->user_data);
    for (size_t i = 0; i < S.mp.size(); i++)
        arm_rhs<adouble>(states + 4*i, controls, S.mp[i], derivatives + 4*i);
}

// Eight per scenario: four initial states, known exactly and shared, and four
// terminal states bounded to a tolerance ball rather than pinned. The relaxation
// is what makes the problem well posed for more than one scenario -- one
// open-loop torque history cannot steer several different plants to the same
// point -- and it is also what keeps the scenario budget open, because an
// inequality event costs no degrees of freedom.
void events(adouble* e, adouble* initial_states, adouble* final_states,
            adouble* /*parameters*/, adouble& /*t0*/, adouble& /*tf*/,
            adouble* /*xad*/, int /*iphase*/, Workspace* workspace)
{
    const Scenarios& S = *((Scenarios*) workspace->problem->user_data);
    const int M = (int) S.mp.size();
    int k = 0;
    for (int i = 0; i < M; i++)
        for (int j = 0; j < 4; j++) e[k++] = initial_states[4*i + j];
    for (int i = 0; i < M; i++)
        for (int j = 0; j < 4; j++) e[k++] = final_states[4*i + j];
}

void linkages(adouble* /*linkages*/, adouble* /*xad*/, Workspace* /*workspace*/) {}

//////////////////////////////////////////////////////////////////////////
///////////////////  The two callbacks the driver needs  /////////////////
//////////////////////////////////////////////////////////////////////////

// Everything the callbacks share. Held in one place and passed as user_data, so
// that nothing here is a global and the example can be read from the top.
struct Context {
    Scenarios scenarios;      // what the dae and events read
    double    margin;         // inward tightening of the terminal ball
};

// 1. BUILD the augmented problem for a scenario list.
//
// Called before every solve. The tightened terminal ball is the detail that makes
// the loop terminate: the scenarios are satisfied to the tolerance EXACTLY -- an
// active constraint is active -- so the miss between two neighbouring scenarios is
// necessarily a little larger, and a loop demanding the full tolerance over the
// whole set from a design only ever asked for it at finitely many points could not
// finish.
static void setup(Prob& problem, Alg& algorithm,
                  const std::vector<RowVectorXd>& scenarios,
                  const RobustDesign& previous, void* user_data)
{
    Context& C = *((Context*) user_data);
    const int M  = (int) scenarios.size();
    const int nx = 4*M;

    C.scenarios.mp.clear();
    for (int i = 0; i < M; i++) C.scenarios.mp.push_back(scenarios[i](0));

    problem.name        = "Two-link arm with an uncertain payload (driver)";
    problem.outfilename = "robust_driver.txt";
    problem.nphases     = 1;
    problem.nlinkages   = 0;
    psopt_level1_setup(problem);

    problem.phases(1).nstates   = nx;
    problem.phases(1).ncontrols = 2;
    problem.phases(1).nevents   = 8*M;
    problem.phases(1).npath     = 0;
    problem.phases(1).nodes     << NODES;
    psopt_level2_setup(problem, algorithm);

    problem.user_data = (void*) &C.scenarios;

    for (int j = 0; j < nx; j++) {
        problem.phases(1).bounds.lower.states(j) = -2.0;
        problem.phases(1).bounds.upper.states(j) =  2.0;
    }
    for (int j = 0; j < 2; j++) {
        problem.phases(1).bounds.lower.controls(j) = -1.0;
        problem.phases(1).bounds.upper.controls(j) =  1.0;
    }

    int k = 0;
    for (int i = 0; i < M; i++)
        for (int j = 0; j < 4; j++) {
            problem.phases(1).bounds.lower.events(k) = X0[j];
            problem.phases(1).bounds.upper.events(k) = X0[j];
            k++;
        }
    for (int i = 0; i < M; i++)
        for (int j = 0; j < 4; j++) {
            problem.phases(1).bounds.lower.events(k) = XF[j] - C.margin;
            problem.phases(1).bounds.upper.events(k) = XF[j] + C.margin;
            k++;
        }

    problem.phases(1).bounds.lower.StartTime = 0.0;
    problem.phases(1).bounds.upper.StartTime = 0.0;
    problem.phases(1).bounds.lower.EndTime   = 1.0;
    problem.phases(1).bounds.upper.EndTime   = 15.0;

    problem.integrand_cost = &integrand_cost;
    problem.endpoint_cost  = &endpoint_cost;
    problem.dae            = &dae;
    problem.events         = &events;
    problem.linkages       = &linkages;

    // The warm start. Each scenario's own plant integrated through the previous
    // control, which is feasible scenario by scenario rather than merely the right
    // shape -- without it the augmented problem starts from M copies of an
    // infeasible arc and the solver is free to wander to a distant local minimum.
    //
    // WARM_SUB is the part that is easy to get wrong, and this example got it
    // wrong first. The horizon being integrated over is a DECISION VARIABLE: the
    // robust final time grows as scenarios are added, here from about three
    // seconds to nearly ten. One RK4 step per node interval is 0.12 s at the
    // start and 0.4 s by the end, and at 0.4 s this plant's integration diverges.
    // The guess then holds infinities, the constraint evaluation returns NaN, and
    // PSOPT reports rows its coverage guard found unwritten -- a message that
    // names a defect in the library when the fault is a guess built with a step
    // size chosen for a shorter horizon. A step size validated on the nominal
    // problem is not validated on the robust one.
    const int WARM_SUB = 8;
    MatrixXd x_guess(nx, NODES);
    bool warm_ok = (previous.valid && previous.time.cols() == NODES);
    if (warm_ok) {
        const MatrixXd& tg = previous.time;
        const MatrixXd& ug = previous.controls;
        for (int i = 0; i < M && warm_ok; i++) {
            double x[4] = { X0[0], X0[1], X0[2], X0[3] };
            for (int j = 0; j < 4; j++) x_guess(4*i + j, 0) = x[j];
            for (int c = 0; c < NODES - 1 && warm_ok; c++) {
                const double h = (tg(0, c+1) - tg(0, c))/WARM_SUB;
                for (int sub = 0; sub < WARM_SUB; sub++) {
                    const double w0 =  sub       /(double) WARM_SUB;
                    const double wh = (sub + 0.5)/(double) WARM_SUB;
                    const double w1 = (sub + 1.0)/(double) WARM_SUB;
                    double uA[2], uH[2], uB[2], k1[4], k2[4], k3[4], k4[4], y[4];
                    for (int j = 0; j < 2; j++) {
                        const double a = ug(j, c), b = ug(j, c+1);
                        uA[j] = a + w0*(b - a);
                        uH[j] = a + wh*(b - a);
                        uB[j] = a + w1*(b - a);
                    }
                    arm_rhs<double>(x, uA, C.scenarios.mp[i], k1);
                    for (int j=0;j<4;j++) y[j] = x[j] + 0.5*h*k1[j];
                    arm_rhs<double>(y, uH, C.scenarios.mp[i], k2);
                    for (int j=0;j<4;j++) y[j] = x[j] + 0.5*h*k2[j];
                    arm_rhs<double>(y, uH, C.scenarios.mp[i], k3);
                    for (int j=0;j<4;j++) y[j] = x[j] + h*k3[j];
                    arm_rhs<double>(y, uB, C.scenarios.mp[i], k4);
                    for (int j=0;j<4;j++)
                        x[j] += (h/6.0)*(k1[j] + 2*k2[j] + 2*k3[j] + k4[j]);
                }
                for (int j = 0; j < 4; j++) {
                    if (!std::isfinite(x[j])) { warm_ok = false; break; }
                    x_guess(4*i + j, c+1) = x[j];
                }
            }
        }
        if (warm_ok) {
            problem.phases(1).guess.controls = ug;
            problem.phases(1).guess.time     = tg;
        }
    }
    if (!warm_ok) {
        for (int i = 0; i < M; i++)
            for (int j = 0; j < 4; j++)
                x_guess.row(4*i + j) = X0[j]*ones(1, NODES);
        problem.phases(1).guess.controls = zeros(2, NODES);
        problem.phases(1).guess.time     = linspace(0.0, 3.0, NODES);
    }
    problem.phases(1).guess.states = x_guess;

    algorithm.nlp_method                  = "IPOPT";
    algorithm.scaling                     = "automatic";
    algorithm.derivatives                 = "automatic";
    algorithm.nlp_iter_max                = 1500;
    algorithm.nlp_tolerance               = 1.e-6;
    algorithm.transcription_method        = "multiple-shooting";
    algorithm.ms_integrator               = "RK4";
    algorithm.ms_steps_per_segment        = 12;
    algorithm.ms_control_parameterisation = "linear";
    algorithm.print_level                 = 0;
}

// A bridge so that the worst-case oracle can be pointed at ANY design, not only at
// the one the driver is currently holding. Used below to ask how badly the
// deterministic design serves the set it was never shown -- which is the number
// the whole comparison rests on, and which cannot be had by giving the driver a
// one-point uncertainty, because the worst case over one point is zero by
// construction.
struct Probe { const RobustDesign* design; void* user_data; };

static double probe_violation(const RowVectorXd& theta, void* p);

// 2. VERIFY a design at one payload, with an integrator that is not PSOPT's.
//
// Fixed-step RK4 through the plant at this payload, driven by the designed control
// by linear interpolation -- which reproduces the designed control exactly, because
// multiple shooting with a linear parameterisation makes it piecewise linear
// through the node table. Returns how far outside the declared terminal ball the
// trajectory arrives, and zero when it arrives inside: the same measure the driver
// compares against its slack.
static double violation(const RowVectorXd& theta, const RobustDesign& design,
                        void* /*user_data*/)
{
    const double   mp = theta(0);
    const MatrixXd& t = design.time;
    const MatrixXd& u = design.controls;
    const int       N = (int) t.cols();
    const int    nsub = 16;

    double x[4] = { X0[0], X0[1], X0[2], X0[3] };
    for (int i = 0; i < N - 1; i++) {
        const double h = (t(0, i+1) - t(0, i))/nsub;
        for (int s = 0; s < nsub; s++) {
            const double w0 =  s       /(double) nsub;
            const double wh = (s + 0.5)/(double) nsub;
            const double w1 = (s + 1.0)/(double) nsub;
            double uA[2], uH[2], uB[2], k1[4], k2[4], k3[4], k4[4], y[4];
            for (int j = 0; j < 2; j++) {
                const double a = u(j, i), b = u(j, i+1);
                uA[j] = a + w0*(b - a);
                uH[j] = a + wh*(b - a);
                uB[j] = a + w1*(b - a);
            }
            arm_rhs<double>(x, uA, mp, k1);
            for (int j=0;j<4;j++) y[j] = x[j] + 0.5*h*k1[j];
            arm_rhs<double>(y, uH, mp, k2);
            for (int j=0;j<4;j++) y[j] = x[j] + 0.5*h*k2[j];
            arm_rhs<double>(y, uH, mp, k3);
            for (int j=0;j<4;j++) y[j] = x[j] + h*k3[j];
            arm_rhs<double>(y, uB, mp, k4);
            for (int j=0;j<4;j++)
                x[j] += (h/6.0)*(k1[j] + 2*k2[j] + 2*k3[j] + k4[j]);
        }
    }

    // The violation is the excess beyond the DECLARED ball, not the miss itself.
    // Those differ by exactly the ball's radius, and confusing them is how a design
    // that misses by 0.053 gets reported as having met a tolerance of 0.03.
    double v = 0.0;
    for (int j = 0; j < 4; j++) {
        const double e = fabs(x[j] - XF[j]) - DELTA;
        if (e > v) v = e;
    }
    return v;
}

static double probe_violation(const RowVectorXd& theta, void* p)
{
    Probe* q = (Probe*) p;
    return violation(theta, *q->design, q->user_data);
}

//////////////////////////////////////////////////////////////////////////

int main(void)
{
    Context C;
    C.margin = 0.9*DELTA;          // the design is given 90% of the ball

    RowVectorXd mean(1);  mean  << MU;
    MatrixXd    cov(1,1); cov   << SIGMA*SIGMA;

    RobustSpec spec;
    spec.uncertainty    = robust_gaussian(mean, cov, 3.0);
    spec.slack          = SLACK;
    spec.max_iterations = 12;
    spec.n_seed         = 128;
    spec.n_refine       = 3;
    spec.setup          = &setup;
    spec.violation      = &violation;
    spec.user_data      = (void*) &C;
    spec.verbose        = true;

    printf("\nTwo-link arm with an uncertain payload, through psopt_solve_robust\n");
    printf("=================================================================\n");
    printf("  payload m_p ~ N(mu = %.3f, sigma = %.3f), set = mu +- 3 sigma\n",
           MU, SIGMA);
    printf("  terminal ball %.3f, slack beyond it %.4f\n", DELTA, SLACK);
    printf("  the %.1f sigma set holds %.2f%% of the distribution in %d dimension\n",
           spec.uncertainty.truncate,
           100.0*robust_set_coverage(spec.uncertainty.truncate, 1), 1);

    {
        std::vector<RowVectorXd> pts;
        RowVectorXd              w;
        (void) robust_sigma_points(spec.uncertainty, pts, w);
        printf("  starting scenarios %.4f %.4f %.4f, weights %.4f %.4f %.4f\n",
               pts[0](0), pts[1](0), pts[2](0), w(0), w(1), w(2));
        printf("  -- which is the three-point Gauss-Hermite rule, mu and "
               "mu +- sqrt(3) sigma\n");
        printf("  the worst-case search will spend %ld integrations per iteration\n",
               robust_worst_case_evaluations(spec.uncertainty, spec.n_seed,
                                             spec.n_refine));
    }

    Sol  solution;
    Prob problem;
    Alg  algorithm;
    const int rc = psopt_solve_robust(solution, spec, problem, algorithm);

    if (!spec.design.valid) {
        printf("\n  the robust design failed (return code %d)\n\n", rc);
        return 1;
    }

    // ---- what it cost, against the deterministic design ---------------------
    //
    // One scenario at the mean, solved through the same driver with generation
    // switched off, so that nothing but the scenario set differs between the two.
    RobustSpec nom;
    nom.uncertainty = robust_explicit(std::vector<RowVectorXd>(1, mean));
    nom.slack       = 1.0e9;                       // accept whatever it gives
    nom.max_iterations = 1;
    nom.setup       = &setup;
    nom.violation   = &violation;
    nom.user_data   = (void*) &C;
    nom.verbose     = false;

    Sol  nsol;
    Prob nprob;
    Alg  nalg;
    (void) psopt_solve_robust(nsol, nom, nprob, nalg);

    // Its worst case over the REAL set. nom.certificate is the worst over the
    // one-point set it was designed against, which is zero and means nothing.
    Probe probe;
    probe.design    = &nom.design;
    probe.user_data = (void*) &C;
    RowVectorXd nom_at;
    const double nom_worst = robust_worst_case(spec.uncertainty, &probe_violation,
                                               &probe, 128, 3, 11u, nom_at);

    printf("\n  %-34s %10s %14s\n", "design", "t_f", "worst in set");
    printf("  %-34s %10.4f %14.3e\n", "nominal (the mean payload alone)",
           nom.design.objective, nom_worst);
    printf("  %-34s %10.4f %14.3e\n", "robust (generated scenarios)",
           spec.design.objective, spec.certificate);
    printf("\n  price of robustness: t_f %.4f -> %.4f, a factor of %.2f\n",
           nom.design.objective, spec.design.objective,
           spec.design.objective/nom.design.objective);

    // ---- an independent check of the certificate ----------------------------
    //
    // The driver reports the worst violation its SEARCH found. For one parameter a
    // dense scan is exhaustive to its resolution, so it costs little to confirm the
    // certificate that way rather than take it on trust.
    {
        const int NS = 4001;
        double worst = 0.0, at = 0.0;
        for (int i = 0; i < NS; i++) {
            RowVectorXd th(1);
            th(0) = (MU - 3.0*SIGMA)
                  + (6.0*SIGMA)*((double) i/(double) (NS - 1));
            const double v = violation(th, spec.design, (void*) &C);
            if (v > worst) { worst = v; at = th(0); }
        }
        printf("\n  the certificate against a dense scan of %d payloads:\n", NS);
        printf("    the search found %.3e at m_p = %.4f\n",
               spec.certificate, spec.certificate_at(0));
        printf("    the scan found   %.3e at m_p = %.4f\n", worst, at);
        printf("    they agree to    %.2e\n", fabs(worst - spec.certificate));
    }

    printf("\n  %d scenarios, %d calls to psopt(), %ld integrations.\n",
           (int) spec.scenarios.size(), spec.n_solves, spec.evaluations);
    printf("  Converged: %s.\n\n", spec.converged ? "yes" : "no");
    return spec.converged ? 0 : 1;
}

////////////////////////////////////////////////////////////////////////////
///////////////////////      END OF FILE     ///////////////////////////////
////////////////////////////////////////////////////////////////////////////
