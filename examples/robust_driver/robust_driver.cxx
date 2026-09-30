//////////////////////////////////////////////////////////////////////////
//////////////////      robust_driver.cxx       //////////////////////////
//////////////////////////////////////////////////////////////////////////
////////////////           PSOPT  Example             ////////////////////
//////////////////////////////////////////////////////////////////////////
//////// Title:   The two-link arm with an uncertain payload,        //////
////////          through psopt_solve_robust and a nominal model     //////
//////// Last modified: 30 September 2026                            //////
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
//  examples/robust_arm writes the whole robust design loop out by hand: the augmented
//  problem, the verification integrator, the worst-case scan, the tightening, the warm
//  start. It is the study, and it is where the method is explained.
//
//  This is the same problem stated as a NOMINAL model and handed to the driver, which is
//  the short way to write one. What the user writes is the problem: four states, two
//  controls, eight events, the dynamics with the payload as one extra argument, the
//  bounds, a guess. What the library does is everything that follows from putting M
//  copies of that problem side by side:
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
//  examples/robust_arm does. The driver prints which of the two produced its
//  certificate, and this example ends by checking that certificate against a scan of
//  4001 payloads, which for one uncertain parameter is exhaustive to its resolution.
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
    algorithm.ms_steps_per_segment        = 12;
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

    printf("\n  %d scenarios, %d calls to psopt(), %ld integrations, verified by %s.\n",
           (int) spec.scenarios.size(), spec.n_solves, spec.evaluations,
           spec.own_verifier ? "the caller's integrator" : "the library's");
    printf("  Converged: %s.\n\n", spec.converged ? "yes" : "no");
    return spec.converged ? 0 : 1;
}
