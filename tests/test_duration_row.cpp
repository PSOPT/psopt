//////////////////////////////////////////////////////////////////////////////
// test_duration_row.cpp
//
// The t0 <= tf row, and the lower bound it used to carry.
//
// PSOPT appends one row per phase for t0 - tf <= 0. It used to be stated
// two-sidedly,
//
//     t0MIN - tfMAX  <=  t0 - tf  <=  0,
//
// and the lower bound is the smallest value t0 - tf can take anywhere in the
// variable box: t0 >= t0MIN and tf <= tfMAX give t0 - tf >= t0MIN - tfMAX at
// every point the NLP can reach. It excluded nothing.
//
// What it did instead was manufacture a degenerate active constraint. On a phase
// whose horizon is fixed -- which is most of them -- t0 and tf are pinned by
// coincident bounds, so t0 - tf *equals* t0MIN - tfMAX identically and the row
// sits exactly on that bound at every iterate. Ipopt removes fixed variables
// (fixed_variable_treatment defaults to make_parameter), so the row's gradient in
// what remains is identically zero: an inequality that is always active and has
// no gradient, which is a point at which LICQ fails and whose multiplier the
// barrier drives towards infinity. PSOPT handed the NLP one per phase on nearly
// every example it ships.
//
// The tests below pin the two halves of the claim. The first is that the bound is
// redundant, which is a statement about the problem and is checked by solving with
// and without a horizon so tight that the row would bind if it ever could. The
// second is that the remaining upper bound still does its job when the times are
// free, which is the case the row exists for.
//////////////////////////////////////////////////////////////////////////////

#include "gtest/gtest.h"
#include <psopt.h>

#include <cmath>
#include <string>

namespace durrow {

// Minimum-time double integrator: xdot = v, vdot = u, |u| <= 1, from (0,0) to
// (1,0). The exact answer is tf* = 2, bang-bang with one switch at t = 1.
const double TF_EXACT = 2.0;

adouble endpoint_cost(adouble*, adouble*, adouble*, adouble&, adouble& tf, adouble*,
                      int, Workspace*)
{ return tf; }

adouble integrand_cost(adouble*, adouble*, adouble*, adouble&, adouble*, int,
                       Workspace*)
{ return 0.0; }

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
    double      tf;
    std::string message;
};

// t0 is always fixed at zero. tf_lo == tf_hi is the fixed-horizon case, in which
// the row's old lower bound was exactly attained at every iterate.
static Run solve(double tf_lo, double tf_hi)
{
    Alg algorithm; Sol solution; Prob problem;
    Run out; out.flag = -1; out.cost = 0.0; out.tf = 0.0;

    problem.name        = "minimum-time double integrator";
    problem.outfilename = "test_duration_row.txt";
    problem.nphases     = 1;
    problem.nlinkages   = 0;
    psopt_level1_setup(problem);

    problem.phases(1).nstates   = 2;
    problem.phases(1).ncontrols = 1;
    problem.phases(1).nevents   = 4;
    problem.phases(1).npath     = 0;
    problem.phases(1).nodes     << 40;
    psopt_level2_setup(problem, algorithm);

    problem.phases(1).bounds.lower.states   << -2.0, -2.0;
    problem.phases(1).bounds.upper.states   <<  2.0,  2.0;
    problem.phases(1).bounds.lower.controls(0) = -1.0;
    problem.phases(1).bounds.upper.controls(0) =  1.0;
    problem.phases(1).bounds.lower.events   << 0.0, 0.0, 1.0, 0.0;
    problem.phases(1).bounds.upper.events   << 0.0, 0.0, 1.0, 0.0;
    problem.phases(1).bounds.lower.StartTime = 0.0;
    problem.phases(1).bounds.upper.StartTime = 0.0;
    problem.phases(1).bounds.lower.EndTime   = tf_lo;
    problem.phases(1).bounds.upper.EndTime   = tf_hi;

    problem.integrand_cost = &integrand_cost;
    problem.endpoint_cost  = &endpoint_cost;
    problem.dae            = &dae;
    problem.events         = &events;
    problem.linkages       = &linkages;

    problem.phases(1).guess.states   = zeros(2, 40);
    problem.phases(1).guess.states.row(0) = linspace(0.0, 1.0, 40);
    problem.phases(1).guess.controls = zeros(1, 40);
    problem.phases(1).guess.time     = linspace(0.0, tf_hi, 40);

    algorithm.nlp_method         = "IPOPT";
    algorithm.scaling            = "automatic";
    algorithm.derivatives        = "automatic";
    algorithm.collocation_method = "Legendre";
    algorithm.nlp_iter_max       = 1000;
    algorithm.nlp_tolerance      = 1.0e-6;
    algorithm.mesh_refinement    = "manual";
    algorithm.print_level        = 0;

    out.flag    = psopt(solution, problem, algorithm);
    out.message = solution.error_msg;
    if (out.flag == 0) {
        out.cost = solution.cost;
        MatrixXd t = solution.get_time_in_phase(1);
        out.tf = t(0, (int) t.cols() - 1);
    }
    return out;
}

} // namespace durrow


// ---------------------------------------------------------------------------
// A fixed horizon. This is the configuration in which the row's old lower bound
// was exactly attained at every iterate and its gradient was empty, and it is by
// far the commonest one: of the examples PSOPT ships, nearly all fix the horizon.
//
// The problem is posed with the horizon pinned at the known optimum, so that the
// answer is decided and the only thing under test is whether the transcription
// gets in its own way.
// ---------------------------------------------------------------------------

TEST(DurationRow, AFixedHorizonSolves)
{
    const durrow::Run r = durrow::solve(durrow::TF_EXACT, durrow::TF_EXACT);
    ASSERT_EQ(r.flag, 0) << r.message;
    EXPECT_NEAR(r.cost, durrow::TF_EXACT, 1.0e-6);
    EXPECT_NEAR(r.tf,   durrow::TF_EXACT, 1.0e-9);
}


// ---------------------------------------------------------------------------
// A free final time, which is the case the row exists for: with tf a genuine
// variable, t0 - tf <= 0 is a constraint and the upper bound has to survive.
//
// The window is wide and straddles the answer, so the solver has to find tf = 2
// rather than be given it.
//
// The tolerance is 5e-3 and that is the discretisation, not slack in the test. The
// minimum-time control here is bang-bang with a switch at t = 1, and a single
// degree-39 Legendre polynomial cannot represent a jump: the answer comes out at
// 2.00264, high by 1.3e-3 relative, on this mesh. Tightening the test would be
// testing the mesh rather than the row.
// ---------------------------------------------------------------------------

TEST(DurationRow, AFreeFinalTimeStillFindsTheMinimumTime)
{
    const durrow::Run r = durrow::solve(0.5, 6.0);
    ASSERT_EQ(r.flag, 0) << r.message;
    EXPECT_NEAR(r.cost, durrow::TF_EXACT, 5.0e-3);
    EXPECT_NEAR(r.tf,   durrow::TF_EXACT, 5.0e-3);
    EXPECT_GT(r.tf, 0.5 + 1.0e-6);      // not pinned at either end of the window,
    EXPECT_LT(r.tf, 6.0 - 1.0e-6);      // which is what says the row did not decide it
}


// ---------------------------------------------------------------------------
// And the statement the change rests on, which is about the problem rather than
// about the solver: t0 - tf >= t0MIN - tfMAX holds at every point of the variable
// box, so a lower bound at that value excludes nothing. Written out for the three
// shapes a phase can have.
// ---------------------------------------------------------------------------

TEST(DurationRow, TheOldLowerBoundWasImpliedByTheVariableBox)
{
    struct { double t0lo, t0hi, tflo, tfhi; const char* what; } cases[] = {
        { 0.0, 0.0, 2.0, 2.0, "fixed start, fixed horizon" },
        { 0.0, 0.0, 0.5, 6.0, "fixed start, free final time" },
        { 1.0, 3.0, 2.0, 9.0, "both free, overlapping" },
    };

    for (size_t k = 0; k < sizeof cases / sizeof cases[0]; k++) {
        const double floor_of_t0_minus_tf = cases[k].t0lo - cases[k].tfhi;

        // The corners of the box are where t0 - tf is extreme, so checking them
        // checks the box.
        const double corners[4] = {
            cases[k].t0lo - cases[k].tflo, cases[k].t0lo - cases[k].tfhi,
            cases[k].t0hi - cases[k].tflo, cases[k].t0hi - cases[k].tfhi
        };
        for (int c = 0; c < 4; c++)
            EXPECT_GE(corners[c], floor_of_t0_minus_tf) << cases[k].what;

        // And it is attained, which is why the row was always active rather than
        // merely redundant: the corner (t0MIN, tfMAX) is admissible.
        EXPECT_DOUBLE_EQ(corners[1], floor_of_t0_minus_tf) << cases[k].what;
    }
}
