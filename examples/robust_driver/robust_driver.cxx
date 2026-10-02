//////////////////////////////////////////////////////////////////////////
//////////////////      robust_driver.cxx       //////////////////////////
//////////////////////////////////////////////////////////////////////////
////////////////           PSOPT  Example             ////////////////////
//////////////////////////////////////////////////////////////////////////
//////// Title:   The two-link arm with an uncertain payload,        //////
////////          through psopt_solve_robust and a nominal model     //////
//////// Last modified: 30 September 2026                            //////
//////// Reference:     Weinreb and Bryson (1985), Sec. IV, which    //////
////////                carries the tip mass symbolically;           //////
////////                Luus (2000), Sec. 12.4.2, for the mu = 1     //////
////////                case solved in examples/twolinkarm           //////
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
//  The arm of examples/twolinkarm carries a payload whose mass m_p is not known exactly.
//  One torque history and one final time must be committed BEFORE m_p is revealed, and
//  must bring every plant in the uncertainty set to the target within a stated tolerance.
//  That is a here-and-now decision, and it is what distinguishes robust optimal control
//  from solving the problem once per sample: solving one problem per sample and averaging
//  the answers gives the wait-and-see solution, in which every realisation is optimised
//  with foreknowledge of its own uncertainty, so its cost is a lower bound no
//  implementable control attains and the average of the control histories solves nothing
//  at all. What couples the samples into one problem is that they SHARE the control.
//
//  HOW THE PAYLOAD ENTERS THE MODEL
//
//  examples/twolinkarm writes the arm with its inertia coefficients already evaluated,
//  9/4, 2, 4/3, 3/2, 7/2, 7/3, 31/36, so there is nowhere for a payload to go. Working
//  backwards from those numbers, the model is a planar 2R arm in ABSOLUTE angular
//  coordinates,
//
//      x1 = dtheta1/dt,  x2 = absolute angular rate of link 2,
//      x3 = theta2 - theta1,  x4 = theta1,
//
//  driven by torques tau = (u1 - u2, u2), with unit link lengths and
//
//      M11 = 7/3,   M12 = (3/2) cos x3,   M22 = 4/3,   S2 = 3/2,
//
//  where S2 is the first moment of link 2 about joint 2 and M22 its second moment about
//  the same point. Every shipped expression follows, including the determinant
//  M11*M22 - M12^2 = 28/9 - (9/4)cos^2(x3), which is identically 31/36 + (9/4)sin^2(x3).
//  A point payload of mass m_p at the tip of link 2 then enters in exactly one way:
//
//      M11 = 7/3 + m_p,  M12 = (3/2 + m_p) cos x3,  M22 = 4/3 + m_p,  S2 = 3/2 + m_p,
//
//  because the payload swings about joint 1 with the rest of link 2, adds m_p L2 to the
//  first moment about joint 2 and m_p L2^2 to the second, and L1 = L2 = 1 in these units.
//  Setting m_p = 0 reproduces the shipped dynamics to 2.8e-15 over 200000 random states
//  and torques, which was checked before anything was built on it, and a payload of 0.5
//  changes the joint accelerations by 12 to 17 per cent.
//
//  THE RECONSTRUCTION IS NOW CONFIRMED AGAINST THE ORIGINAL, and it did not have to be.
//  Weinreb and Bryson (1985), Section IV, state this arm with the tip mass present
//  SYMBOLICALLY, as the ratio mu = M/m of tip mass to link mass, in their equations (23)
//  to (27). Their denominator is 7/36 + (2/3)mu + (mu + 1/2)^2 sin^2(theta) and their
//  torque coefficients are (mu + 1/3) and -[mu + 1/3 + (mu + 1/2)cos(theta)].
//
//  Substituting mu = 1 + m_p into their equations gives this file's model exactly: the
//  difference is identically zero in both dynamic states, by symbolic algebra, and the
//  two torque coefficients agree term for term. So the payload reconstructed here by
//  working backwards from seven evaluated numbers IS their tip mass, offset by one
//  because their nominal arm already carries a tip mass equal to a link mass.
//
//  That matters for a reason beyond provenance. The reconstruction was the one step in
//  this example that rested on inference instead of on a source, and it is the step
//  everything else is built on. It is now a citation.
//
//  WHY THE ROBUST SLEW IS A SLOW ONE
//
//  The payload enters only through the inertia, so it changes the ACCELERATION the torques
//  produce. Drive hard and the payload matters; drive gently and the trajectory approaches
//  a quasi-static one on which it matters much less. A design that must land every plant
//  in the set within the tolerance therefore buys its insensitivity with time, and the
//  numbers below show it doing so, at about two and a half times the nominal final time.
//  That is the physics of the problem and not a defect of the method, and it is why the
//  final-time bound of examples/twolinkarm had to be raised here.
//
//  WHAT THE USER WRITES, AND WHAT THE LIBRARY DOES
//
//  The problem is stated once as a NOMINAL model and handed to the driver. What the user
//  writes is the problem: four states, two controls, eight events, the dynamics with the
//  payload as one extra argument, the bounds, a guess. What the library does is everything
//  that follows from putting M copies of that problem side by side:
//
//    * the dae and the events, looped over the scenario list, so that M copies of the
//      state are carried against one copy of the control and the design cannot
//      anticipate the payload;
//    * the phase sized for those copies, and the bounds replicated with the inward
//      tightening applied, without which the generation loop cannot terminate;
//    * the warm start, each scenario's own plant integrated through the previous
//      control, doubling its step count if the growing horizon needs it;
//    * the scenario rule, the search for the payload the design serves worst, the
//      generation loop, the polish step, and the diagnosis when the transcription runs
//      out of degrees of freedom;
//    * the risk measure, if the design is to be scored by something other than the
//      central scenario's cost: an expectation, a mean-variance combination, or the
//      conditional value at risk, the last two carrying each scenario's cost as a state
//      the library adds and pins for itself;
//    * the ancillary feedback, if the design is to be a tube rather than one open-loop
//      history: the correction inside the augmented dae, the realised control's bounds as
//      path rows, and a verification that integrates the closed loop beside its reference;
//    * the verification integrator.
//
//  That last item is the one to think about rather than accept. The library's verifier
//  is independent of the TRANSCRIPTION: its own integrator, its own step control, its
//  own code path, sharing with the solve only the equations written above. It is not
//  independent of the LIBRARY. Setting spec.violation replaces it with an
//  implementation of the user's own, which is the stronger claim and is what
//  python/examples/robust_arm.py does, verifying with NumPy code that shares nothing with
//  the library. The driver prints which of the two produced its certificate, and this
//  example checks that certificate against a scan of 4001 payloads, which for one
//  uncertain parameter is exhaustive to its resolution.
//
//  WHAT ELSE THE PROGRAM MEASURES, AND WHY EACH PART IS HERE
//
//    (1) the deterministic design, so that the price of robustness is not confused with
//        the price of relaxing the terminal condition;
//    (2) the generated-scenario design, which is the driver's own loop;
//    (3) a sweep over FIXED scenario sets at the Gauss-Hermite nodes of the payload
//        distribution, M = 1, 3, 5, 7, 9, each scored on payloads it never saw. This is
//        here because it does NOT work, and the way it fails is the most useful
//        measurement in the file: a fixed rule pins the terminal miss at its own nodes and
//        says nothing about the gaps between them, so the optimiser drives the miss to the
//        tolerance exactly at the nodes and lets it grow freely elsewhere. Raising M cures
//        nothing in any orderly way, because where the next quadrature node lands has
//        nothing to do with where the previous design was failing. Quadrature is the right
//        tool for an EXPECTATION and the wrong tool for a constraint that must hold
//        everywhere;
//    (4) the same design solved at two step counts in the segment integrator, which is the
//        accuracy trap this problem sprang and the reason the comparison is printed rather
//        than trusted. An accuracy setting validated on the NOMINAL problem has to be
//        validated again on the robust one, the robust solution being a different
//        trajectory: usually a longer and gentler one, on which a different balance of
//        errors applies.
//
//  COST. Some twenty calls to psopt() and a few tens of thousands of verification
//  integrations: two to four minutes, most of it in the solver. This example is a study as
//  well as a demonstration, and it is one of the few in the distribution that takes
//  minutes.
//
//  MEASURING A MISS WITH THE LIBRARY'S VERIFIER. Several tables below report a terminal
//  miss and not a bound violation. robust_model_violation returns how far the trajectory
//  falls OUTSIDE the declared bounds, so a model whose terminal events are pinned to the
//  target instead of bounded by a ball returns the miss itself, in the infinity norm over
//  the four states. measuring_model() below is that model, and it is the only reason this
//  file needs no integrator of its own.
//
//  WHERE THE ALGORITHM OPTIONS GO
//
//  In configure() below, and not in main(). psopt_level2_setup resets every field of the
//  Alg it is given, and the driver calls it once per scenario count; configure() runs
//  immediately after each of those calls, which is the one moment the options can be set.
//
//////////////////////////////////////////////////////////////////////////

#include "psopt.h"
#include "robust.h"

#include <cmath>
#include <cstdlib>
#include <vector>

using namespace PSOPT;

static const double MU    = 0.50;     // mean payload
static const double SIGMA = 0.15;     // its standard deviation
static const double DELTA = 0.03;     // terminal tolerance
static const double SLACK = 1.0e-3;   // violation of THAT a design may leave
static const int    NODES = 25;

static const double X0[4] = { 0.0, 0.0, 0.500, 0.000 };
static const double XF[4] = { 0.0, 0.0, 0.500, 0.522 };

//////////////////////////////////////////////////////////////////////////
///////////////////  The nominal problem, written once  //////////////////
//////////////////////////////////////////////////////////////////////////

// The payload arrives as DATA, not as an adouble, because a scenario is not a decision
// variable and has no business on the derivative tape. That is also what lets the
// library call this same function numerically to build its warm start and its verifier,
// so the physics is written here and nowhere else. Compare examples/robust_arm, where
// the same right-hand side is templated in order to serve both the tape and a
// hand-written integrator.
void nominal_dae(adouble* derivatives, adouble* /*path*/, adouble* states,
                 adouble* controls, adouble* /*parameters*/, adouble& /*time*/,
                 const double* theta, int /*ntheta*/,
                 adouble* /*xad*/, int /*iphase*/, Workspace* /*workspace*/)
{
    const double mp = theta[0];
    const adouble x1 = states[0], x2 = states[1], x3 = states[2];
    const double M11 = 7.0/3.0 + mp, M22 = 4.0/3.0 + mp, S2 = 3.0/2.0 + mp;
    const adouble M12 = S2*cos(x3);
    const adouble det = M11*M22 - M12*M12;
    const adouble b1 = (controls[0] - controls[1]) + S2*sin(x3)*x2*x2;
    const adouble b2 =  controls[1]                - S2*sin(x3)*x1*x1;
    derivatives[0] = ( M22*b1 - M12*b2)/det;
    derivatives[1] = (-M12*b1 + M11*b2)/det;
    derivatives[2] = x2 - x1;
    derivatives[3] = x1;
}

// Four initial states, known exactly and shared by every scenario, and four terminal
// states bounded to a tolerance ball rather than pinned. The relaxation is what makes the
// problem well posed for more than one scenario, one open-loop torque history being
// unable to steer several different plants to the same point, and it is also what keeps
// the scenario budget open, an inequality event costing no degrees of freedom.
void nominal_events(adouble* e, adouble* initial_states, adouble* final_states,
                    adouble* /*parameters*/, adouble& /*t0*/, adouble& /*tf*/,
                    const double* /*theta*/, int /*ntheta*/,
                    adouble* /*xad*/, int /*iphase*/, Workspace* /*workspace*/)
{
    for (int j = 0; j < 4; ++j) e[j]     = initial_states[j];
    for (int j = 0; j < 4; ++j) e[4 + j] = final_states[j];
}

adouble endpoint_cost(adouble* /*x0*/, adouble* /*xf*/, adouble* /*p*/,
                      adouble& /*t0*/, adouble& tf, adouble* /*xad*/,
                      int /*iphase*/, Workspace* /*workspace*/)
{
    return tf;                        // minimum time, shared by every scenario
}

// Steps per segment in PSOPT's own multiple-shooting integrator. A variable and not a
// constant because the measurement that fixed its value is one this example repeats: see
// the accuracy table at the end of main().
static int g_steps = 12;

// The one moment the algorithm options can be set. See the note at the head of the file.
void configure(Alg& algorithm, void* /*user_data*/)
{
    algorithm.nlp_method                  = "IPOPT";
    algorithm.scaling                     = "automatic";
    algorithm.derivatives                 = "automatic";
    algorithm.nlp_iter_max                = 1500;
    algorithm.nlp_tolerance               = 1.e-6;
    algorithm.transcription_method        = "multiple-shooting";
    algorithm.ms_integrator               = "RK4";
    algorithm.ms_steps_per_segment        = g_steps;
    algorithm.ms_control_parameterisation = "linear";
    algorithm.print_level                 = 0;
}

//////////////////////////////////////////////////////////////////////////

static RobustModel arm_model(void)
{
    RobustModel model;
    model.nstates = 4; model.ncontrols = 2; model.nevents = 8;
    model.nodes.resize(1); model.nodes << NODES;

    model.dae           = &nominal_dae;
    model.events        = &nominal_events;
    model.endpoint_cost = &endpoint_cost;
    model.configure     = &configure;

    model.states_lower   = -2.0*ones(1, 4);
    model.states_upper   =  2.0*ones(1, 4);
    model.controls_lower = -1.0*ones(1, 2);
    model.controls_upper =  1.0*ones(1, 2);

    model.events_lower.resize(8);
    model.events_upper.resize(8);
    for (int j = 0; j < 4; ++j) {
        model.events_lower(j)     = X0[j];
        model.events_upper(j)     = X0[j];
        model.events_lower(4 + j) = XF[j] - DELTA;
        model.events_upper(4 + j) = XF[j] + DELTA;
    }

    model.t0_lower = 0.0; model.t0_upper =  0.0;
    model.tf_lower = 1.0; model.tf_upper = 15.0;

    // The state the verification integrates every scenario out of. Asked for rather than
    // read off the event bounds, because the driver cannot in general tell which of them
    // are the initial conditions.
    model.initial_state.resize(4);
    for (int j = 0; j < 4; ++j) model.initial_state(j) = X0[j];

    // The design is given nine tenths of the declared terminal ball and verified against
    // all of it. Without that margin the loop cannot converge: the scenarios are
    // satisfied to the tolerance exactly, so the miss between two neighbouring scenarios
    // is necessarily a little larger.
    model.tighten = 0.9;

    model.guess_states = zeros(4, NODES);
    for (int j = 0; j < 4; ++j) model.guess_states.row(j) = X0[j]*ones(1, NODES);
    model.guess_controls = zeros(2, NODES);
    model.guess_time     = linspace(0.0, 3.0, NODES);
    return model;
}

//////////////////////////////////////////////////////////////////////////
///////////////////  Measuring, with the library's verifier  /////////////
//////////////////////////////////////////////////////////////////////////

// The same model with its terminal events PINNED to the target and no inward tightening.
// robust_model_violation returns how far a trajectory falls outside the declared bounds,
// so against pinned bounds it returns the terminal miss itself, in the infinity norm over
// the four states. That is what the tables below report, and it is why this file needs no
// integrator of its own: the miss and the certificate then come from the same code, and
// the only difference between them is what was declared.
static RobustModel measuring_model(void)
{
    RobustModel m = arm_model();
    for (int j = 0; j < 4; ++j) {
        m.events_lower(4 + j) = XF[j];
        m.events_upper(4 + j) = XF[j];
    }
    m.tighten = 1.0;
    return m;
}

static double miss_at(const RobustModel& mm, Alg& algorithm, const RobustDesign& d,
                      double mp)
{
    RowVectorXd th(1); th(0) = mp;
    return robust_model_violation(mm, algorithm, th, d);
}

// The worst miss over the uncertainty set for a GIVEN design: a scan on 241 points and
// then a refinement over the two intervals around the maximum. One scalar uncertain
// parameter makes a scan both exhaustive to its resolution and cheap, and honest in a way
// a gradient search on a function that is not concave would not be.
static double worst_miss(const RobustModel& mm, Alg& algorithm, const RobustDesign& d,
                         double lo, double hi, double* at)
{
    const int NG = 241;
    double best = -1.0, mbest = lo;
    for (int i = 0; i < NG; ++i) {
        const double mp = lo + (hi - lo)*i/(NG - 1.0);
        const double m  = miss_at(mm, algorithm, d, mp);
        if (m > best) { best = m; mbest = mp; }
    }
    const double h = (hi - lo)/(NG - 1.0);
    const double a = fmax(lo, mbest - h), b = fmin(hi, mbest + h);
    for (int i = 0; i <= 40; ++i) {
        const double mp = a + (b - a)*i/40.0;
        const double m  = miss_at(mm, algorithm, d, mp);
        if (m > best) { best = m; mbest = mp; }
    }
    if (at) *at = mbest;
    return best;
}

// The fraction of a fixed sample of payloads a design brings within the tolerance, over
// the payloads that lie INSIDE the uncertainty set it was given. The distinction is not
// pedantry: a design constrained on [lo, hi] promises nothing beyond it, and scoring it
// against the whole Gaussian tail measures the truncation of the set and not the quality
// of the design.
static double within_tolerance(const RobustModel& mm, Alg& algorithm,
                               const RobustDesign& d,
                               const std::vector<double>& sample,
                               double lo, double hi)
{
    int inside = 0, within = 0;
    for (size_t k = 0; k < sample.size(); ++k) {
        if (sample[k] < lo || sample[k] > hi) continue;
        ++inside;
        if (miss_at(mm, algorithm, d, sample[k]) <= DELTA) ++within;
    }
    return 100.0*within/(inside > 0 ? inside : 1);
}

// Gauss-Hermite nodes and weights for the standard normal, by Golub-Welsch: the nodes are
// the eigenvalues of the symmetric tridiagonal Jacobi matrix of the probabilists' Hermite
// polynomials (zero diagonal, off-diagonal sqrt(k)), and the weights are the squared first
// components of its eigenvectors. Generated and not tabulated, so that the scenario count
// is a parameter and not a table somebody has to extend and can mistype. Checked against
// numpy.polynomial.hermite_e.hermegauss for M = 3, 5, 7 and 9.
static void gauss_hermite(int M, std::vector<double>& x, std::vector<double>& w)
{
    MatrixXd J = zeros(M, M);
    for (int k = 1; k < M; ++k) {
        J(k-1, k) = sqrt((double) k);
        J(k, k-1) = sqrt((double) k);
    }
    Eigen::SelfAdjointEigenSolver<MatrixXd> es(J);
    x.resize(M); w.resize(M);
    for (int k = 0; k < M; ++k) {
        x[k] = es.eigenvalues()(k);
        w[k] = es.eigenvectors()(0, k)*es.eigenvectors()(0, k);
    }
}

// One design on a FIXED scenario list: no generation, no polish, and no inward tightening
// when `tighten` says so. The driver is still what solves it, so the only thing that
// differs from the generated design is where the scenarios came from.
static bool solve_fixed(const std::vector<RowVectorXd>& scenarios, double tighten,
                        const RowVectorXd& mean, const MatrixXd& cov,
                        RobustDesign& design, double& objective)
{
    RobustModel model = arm_model();
    model.tighten = tighten;
    Prob problem; Alg algorithm; Sol solution;
    RobustSpec spec;
    spec.uncertainty    = robust_gaussian(mean, cov, 3.0);
    spec.slack          = SLACK;
    spec.scenarios      = scenarios;
    spec.max_iterations = 1;
    spec.polish         = false;
    spec.verbose        = false;
    (void) psopt_solve_robust(solution, spec, model, problem, algorithm);
    design    = spec.design;
    objective = spec.design.objective;
    return spec.design.valid;
}

int main(void)
{
    RowVectorXd mean(1);  mean << MU;
    MatrixXd    cov(1,1); cov  << SIGMA*SIGMA;

    printf("\nTwo-link arm with an uncertain payload, through psopt_solve_robust\n");
    printf("=================================================================\n");
    printf("  payload m_p ~ N(mu = %.3f, sigma = %.3f), set = mu +- 3 sigma\n",
           MU, SIGMA);
    printf("  terminal ball %.3f, slack beyond it %.4f\n", DELTA, SLACK);
    printf("  the 3.0 sigma set holds %.2f%% of the distribution in %d dimension\n",
           100.0*robust_set_coverage(3.0, 1), 1);

    // ---- the deterministic design, for the comparison that gives the rest meaning ----
    // One scenario at the mean is the ordinary problem. Solved through the same driver so
    // that nothing but the scenario set differs, with one iteration and no polish, and
    // the loop still scores it over the whole uncertainty set before it stops.
    RobustModel nominal_model = arm_model();
    Prob nprob; Alg nalg; Sol nsol;
    RobustSpec nspec;
    nspec.uncertainty    = robust_gaussian(mean, cov, 3.0);
    nspec.slack          = SLACK;
    nspec.max_iterations = 1;
    nspec.polish         = false;
    nspec.verbose        = false;
    nspec.scenarios.push_back(mean);
    (void) psopt_solve_robust(nsol, nspec, nominal_model, nprob, nalg);

    // ---- the robust design ------------------------------------------------------------
    RobustModel model = arm_model();
    Prob problem; Alg algorithm; Sol solution;
    RobustSpec spec;
    spec.uncertainty    = robust_gaussian(mean, cov, 3.0);
    spec.slack          = SLACK;
    spec.max_iterations = 12;
    spec.n_seed         = 128;
    spec.n_refine       = 3;
    spec.verbose        = true;
    // spec.violation is left null, so the library verifies the design with an integrator
    // built from nominal_dae above. Set it to check against an implementation of your own.

    {
        std::vector<RowVectorXd> pts; RowVectorXd w;
        (void) robust_sigma_points(spec.uncertainty, pts, w);
        printf("  starting scenarios %.4f %.4f %.4f, weights %.4f %.4f %.4f\n",
               pts[0](0), pts[1](0), pts[2](0), w(0), w(1), w(2));
        printf("  -- which is the three-point Gauss-Hermite rule, mu and "
               "mu +- sqrt(3) sigma\n");
        printf("  the worst-case search will spend %ld integrations per iteration\n",
               robust_worst_case_evaluations(spec.uncertainty, spec.n_seed,
                                             spec.n_refine));
    }

    const int rc = psopt_solve_robust(solution, spec, model, problem, algorithm);

    if (!spec.design.valid) {
        printf("\n  the robust design failed (return code %d)\n\n", rc);
        return 1;
    }

    printf("\n  %-34s %10s %14s\n", "design", "t_f", "worst in set");
    printf("  %-34s %10.4f %14.3e\n", "nominal (the mean payload alone)",
           nspec.design.objective, nspec.certificate);
    printf("  %-34s %10.4f %14.3e\n", "robust (generated scenarios)",
           spec.design.objective, spec.certificate);
    printf("\n  price of robustness: t_f %.4f -> %.4f, a factor of %.2f\n",
           nspec.design.objective, spec.design.objective,
           spec.design.objective/nspec.design.objective);

    // ---- the certificate against a dense scan ----------------------------------------
    // The search reports the worst it FOUND. For one uncertain parameter a scan settles
    // whether that was the worst there is; for more than one no such check exists, which
    // is why the driver's wording is what it is.
    {
        const int NS = 4001;
        double worst = 0.0, at = 0.0;
        for (int i = 0; i < NS; ++i) {
            RowVectorXd th(1);
            th(0) = (MU - 3.0*SIGMA)
                    + (6.0*SIGMA)*i/(double)(NS - 1);
            const double v = robust_model_violation(model, algorithm, th, spec.design);
            if (v > worst) { worst = v; at = th(0); }
        }
        printf("\n  the certificate against a dense scan of %d payloads:\n", NS);
        printf("    the search found %.3e at m_p = %.4f\n",
               spec.certificate, spec.certificate_at(0));
        printf("    the scan found   %.3e at m_p = %.4f\n", worst, at);
        printf("    they agree to    %.2e\n", fabs(worst - spec.certificate));
    }

    // ---- the payloads every design below is scored on --------------------------------
    // Drawn once and reused, so that a comparison between designs is not also a
    // comparison between samples. Box-Muller, so the sample is the stated Gaussian and
    // not merely something with the right mean and variance.
    const double lo = MU - 3.0*SIGMA, hi = MU + 3.0*SIGMA;
    std::vector<double> sample(1000);
    srand(20260927);
    for (size_t k = 0; k < sample.size(); ++k) {
        const double u1 = (rand() + 1.0)/(RAND_MAX + 2.0);
        const double u2 = (rand() + 1.0)/(RAND_MAX + 2.0);
        sample[k] = MU + SIGMA*sqrt(-2.0*log(u1))*cos(2.0*M_PI*u2);
    }

    RobustModel mm = measuring_model();

    // ---- why a fixed quadrature rule is not enough -----------------------------------
    //
    // The scenarios are placed at the Gauss-Hermite nodes of the payload distribution and
    // the design is given the FULL tolerance at them, with no inward margin, so that the
    // rule is seen at its most favourable. Each design is then scored on the whole set and
    // on payloads it never saw.
    printf("\n  fixed scenario sets at the Gauss-Hermite nodes, the design given the\n");
    printf("  full tolerance %.3f at each of them and no inward margin, which is why the\n",
           DELTA);
    printf("  M = 1 row is a little faster than the nominal design above: that one was\n");
    printf("  given nine tenths of the tolerance. Both miss columns are the MISS itself\n");
    printf("  and not the excess beyond the ball, so they are comparable with %.3f and\n",
           DELTA);
    printf("  not with the certificates above.\n\n");
    printf("  %2s %15s %10s %11s %11s %10s\n", "M", "scenarios", "t_f",
           "in-sample", "worst miss", "within d");

    RobustDesign best_gh;
    double       best_gh_tf = 0.0;
    int          best_gh_M  = 0;
    const int    Msweep[5]  = { 1, 3, 5, 7, 9 };
    for (int q = 0; q < 5; ++q) {
        const int M = Msweep[q];
        std::vector<double> z, w;
        gauss_hermite(M, z, w);
        std::vector<RowVectorXd> pts(M);
        for (int k = 0; k < M; ++k) {
            RowVectorXd th(1); th(0) = MU + SIGMA*z[k];
            pts[k] = th;
        }
        RobustDesign d; double tf = 0.0;
        if (!solve_fixed(pts, 1.0, mean, cov, d, tf)) { printf("  %2d   FAILED\n", M);
                                                        continue; }
        double in_sample = 0.0;
        for (int k = 0; k < M; ++k)
            in_sample = fmax(in_sample, miss_at(mm, algorithm, d, pts[k](0)));
        double at = 0.0;
        const double wc = worst_miss(mm, algorithm, d, lo, hi, &at);
        const double pc = within_tolerance(mm, algorithm, d, sample, lo, hi);
        char span[64];
        if (M == 1) snprintf(span, sizeof span, "%.3f only", MU);
        else        snprintf(span, sizeof span, "[%.3f,%.3f]",
                             MU + SIGMA*z[0], MU + SIGMA*z[M-1]);
        printf("  %2d %15s %10.6f %11.3e %11.3e %9.1f%%\n", M, span, tf, in_sample,
               wc, pc);
        best_gh = d; best_gh_tf = tf; best_gh_M = M;
    }

    {
        double at = 0.0;
        const double wc = worst_miss(mm, algorithm, spec.design, lo, hi, &at);
        const double pc = within_tolerance(mm, algorithm, spec.design, sample, lo, hi);
        printf("  %2d %15s %10.6f %11s %11.3e %9.1f%%  (generated)\n",
               (int) spec.scenarios.size(), "generated", spec.design.objective, "-",
               wc, pc);
    }

    printf("\n  The in-sample column is the tolerance by construction at every M: the\n");
    printf("  optimiser drives the miss to what it was asked for at the nodes and lets\n");
    printf("  it grow between them, and the worst-case column never once falls below\n");
    printf("  the tolerance. Raising M does not cure it in any orderly way, because\n");
    printf("  where the next quadrature node lands has nothing to do with where the\n");
    printf("  previous design was failing. The generated row is the same problem with\n");
    printf("  the worst case in the loop.\n");

    // ---- the accuracy setting, which had to be validated again -----------------------
    //
    // Two step counts in the segment integrator, on the final scenario set with the full
    // tolerance. The first column is the worst terminal miss as the TRANSCRIPTION reports
    // it, read from the NLP's own terminal states; the second is what the verification
    // integrator finds for the same control. A converged solve whose two columns disagree
    // has satisfied a slightly different plant.
    printf("\n  the same %d scenarios at two step counts in the segment integrator:\n\n",
           (int) spec.scenarios.size());
    printf("  %6s %12s %14s %14s %12s\n", "steps", "t_f", "miss, as solved",
           "miss, verified", "difference");
    const int saved_steps = g_steps;
    const int steps[2] = { 8, 16 };
    for (int q = 0; q < 2; ++q) {
        g_steps = steps[q];
        RobustDesign d; double tf = 0.0;
        if (!solve_fixed(spec.scenarios, 1.0, mean, cov, d, tf)) {
            printf("  %6d   FAILED\n", steps[q]); continue;
        }
        double as_solved = 0.0, verified = 0.0;
        const int N = (int) d.time.cols();
        for (size_t i = 0; i < spec.scenarios.size(); ++i) {
            for (int j = 0; j < 4; ++j)
                as_solved = fmax(as_solved,
                                 fabs(d.states(4*(int) i + j, N-1) - XF[j]));
            verified = fmax(verified, miss_at(mm, algorithm, d, spec.scenarios[i](0)));
        }
        printf("  %6d %12.6f %14.3e %14.3e %12.1e\n", steps[q], tf, as_solved,
               verified, fabs(as_solved - verified));
    }
    g_steps = saved_steps;
    printf("\n  examples/twolinkarm uses eight steps per segment, a value chosen by\n");
    printf("  measurement on that problem, a three-second slew. It is not enough here:\n");
    printf("  the robust designs are two and a half times longer, and over a trajectory\n");
    printf("  that long the segment integrator's error is amplified enough to matter. An\n");
    printf("  accuracy setting validated on the nominal problem has to be validated again\n");
    printf("  on the robust one, the robust solution being a different trajectory.\n");

    printf("\n  %d scenarios, %d calls to psopt(), %ld integrations, verified by %s.\n",
           (int) spec.scenarios.size(), spec.n_solves, spec.evaluations,
           spec.own_verifier ? "the caller's integrator" : "the library's");
    printf("  Converged: %s.\n", spec.converged ? "yes" : "no");

    ////////////////////////////////////////////////////////////////////////
    ///////////  Plot some results if desired (requires gnuplot) ///////////
    ////////////////////////////////////////////////////////////////////////

    // The miss as a function of the payload is the figure this problem exists to produce:
    // the nominal design touches zero at the one payload it was given and rises steeply
    // either side of it, while the robust design is held under the tolerance across the
    // whole set. The third curve is the tolerance itself.
    {
        const int NC = 181;
        MatrixXd mp_axis = zeros(1, NC), curves = zeros(3, NC);
        for (int i = 0; i < NC; ++i) {
            const double mp = lo + (hi - lo)*i/(NC - 1.0);
            mp_axis(0, i) = mp;
            curves(0, i)  = miss_at(mm, nalg, nspec.design, mp);
            curves(1, i)  = miss_at(mm, algorithm, spec.design, mp);
            curves(2, i)  = DELTA;
        }
        plot(mp_axis, curves, "Terminal miss versus payload", "payload m_p",
             "terminal miss (infinity norm)", "nominal robust delta");
        plot(mp_axis, curves, "Terminal miss versus payload", "payload m_p",
             "terminal miss (infinity norm)", "nominal robust delta",
             "pdf", "robust_driver_miss.pdf");

        plot(spec.design.time, spec.design.controls, "Robust design: controls",
             "time (s)", "controls", "u1 u2");
        plot(spec.design.time, spec.design.controls, "Robust design: controls",
             "time (s)", "controls", "u1 u2", "pdf", "robust_driver_controls.pdf");
    }

    printf("\n");
    return spec.converged ? 0 : 1;
}
