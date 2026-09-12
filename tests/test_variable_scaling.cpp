//////////////////////////////////////////////////////////////////////////////
// test_variable_scaling.cpp
//
// The variable map, and the case a multiplicative factor alone cannot fix.
//
// PSOPT scales every variable by a factor chosen from its bounds, so that the
// scaled variable has magnitude near one. That is the right thing to do to a
// variable bounded between -1000 and 1000. It is not much use on a variable
// bounded between 1000 and 1200, which maps onto [0.833, 1.0]: the magnitude is
// order one and the *variation* is 0.17, so the solver still has to resolve a
// step of 1e-3 in a quantity whose value is near one, which is the difficulty
// scaling exists to remove. algorithm.scaling = "affine" maps such a variable
// onto [-1/2, 1/2], where the magnitude and the variation are both order one.
//
// The tests below pin three things. The map itself must be a bijection and must
// reduce to the historical one when the shift is zero. The rule that chooses the
// factors must produce [-1/2, 1/2] where it can and must fall back where it
// cannot, since a variable without two finite bounds has no centre to move to.
// And -- the test that matters -- a problem posed on an offset box must give the
// same answer under both maps.
//
// That last one is not a formality. The shift has to be applied by every reader
// of the decision vector and by everything that writes bounds or a guess into
// it, and there are a dozen such places. A site that was missed is a site that
// reads a variable under a different map from the one it was written under, and
// the symptom is not a small error: it is a different problem. Solving the same
// problem both ways and requiring the same answer finds any such omission at
// once, which bit-comparing the default against itself cannot do, because with
// the shift at zero a missed site is still correct.
//////////////////////////////////////////////////////////////////////////////

#include "gtest/gtest.h"
#include <psopt.h>

#include <cmath>
#include <string>

namespace varscale {

// Rest-to-rest minimum-energy move of a double integrator: xdot = v, vdot = u,
// from x = X0 to x = X0 + D in fixed time T, both ends at rest, minimising the
// integral of u^2. The answer does not depend on X0 and is
//
//     J* = 12 D^2 / T^3 ,
//
// so putting X0 far from the origin changes nothing about the problem and
// everything about how well a multiplicative factor describes it.
const double X0    = 1000.0;
const double D     =    1.0;
const double T     =    1.0;
const double J_EXACT = 12.0*D*D/(T*T*T);

adouble endpoint_cost(adouble*, adouble*, adouble*, adouble&, adouble&, adouble*,
                      int, Workspace*)
{ return 0.0; }

adouble integrand_cost(adouble*, adouble* u, adouble*, adouble&, adouble*, int,
                       Workspace*)
{ return u[0]*u[0]; }

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
    std::string message;
};

static Run solve(const std::string& scaling)
{
    Alg algorithm; Sol solution; Prob problem;
    Run out; out.flag = -1; out.cost = 0.0;

    problem.name        = "offset-box double integrator";
    problem.outfilename = "test_variable_scaling_" + scaling + ".txt";
    problem.nphases     = 1;
    problem.nlinkages   = 0;
    psopt_level1_setup(problem);

    problem.phases(1).nstates   = 2;
    problem.phases(1).ncontrols = 1;
    problem.phases(1).nevents   = 4;
    problem.phases(1).npath     = 0;
    problem.phases(1).nodes     << 30;
    psopt_level2_setup(problem, algorithm);

    // The offset box. The position never leaves a window of width three placed a
    // thousand units from the origin; the velocity and the control are ordinary.
    problem.phases(1).bounds.lower.states   << X0 - 1.0, -10.0;
    problem.phases(1).bounds.upper.states   << X0 + 2.0,  10.0;
    problem.phases(1).bounds.lower.controls(0) = -100.0;
    problem.phases(1).bounds.upper.controls(0) =  100.0;
    problem.phases(1).bounds.lower.events   << X0, 0.0, X0 + D, 0.0;
    problem.phases(1).bounds.upper.events   << X0, 0.0, X0 + D, 0.0;
    problem.phases(1).bounds.lower.StartTime = 0.0;
    problem.phases(1).bounds.upper.StartTime = 0.0;
    problem.phases(1).bounds.lower.EndTime   = T;
    problem.phases(1).bounds.upper.EndTime   = T;

    problem.integrand_cost = &integrand_cost;
    problem.endpoint_cost  = &endpoint_cost;
    problem.dae            = &dae;
    problem.events         = &events;
    problem.linkages       = &linkages;

    problem.phases(1).guess.states   = zeros(2, 30);
    problem.phases(1).guess.states.row(0) = linspace(X0, X0 + D, 30);
    problem.phases(1).guess.controls = zeros(1, 30);
    problem.phases(1).guess.time     = linspace(0.0, T, 30);

    algorithm.nlp_method         = "IPOPT";
    algorithm.scaling            = scaling;
    algorithm.derivatives        = "automatic";
    algorithm.collocation_method = "Legendre";
    algorithm.nlp_iter_max       = 1000;
    algorithm.nlp_tolerance      = 1.0e-8;
    algorithm.mesh_refinement    = "manual";
    algorithm.print_level        = 0;

    out.flag    = psopt(solution, problem, algorithm);
    out.message = solution.error_msg;
    if (out.flag == 0) out.cost = solution.cost;
    return out;
}

} // namespace varscale


// ---------------------------------------------------------------------------
// The map is a bijection, and the historical map is the shift = 0 case of it.
// ---------------------------------------------------------------------------

TEST(VariableScaling, TheMapRoundTripsAndReducesToTheOldOne)
{
    const double xs[] = { -3.25, 0.0, 1.0, 1000.5, -1e6 };
    const double scs[] = { 1.0, 1e-3, 7.5, 1e4 };
    const double shs[] = { 0.0, 1.0, -2.5, 1100.0 };

    for (double x : xs)
        for (double sc : scs) {
            // Shift zero is exactly the product and the quotient it replaced.
            EXPECT_DOUBLE_EQ(PSOPT::scale_variable(x, sc, 0.0),   x*sc);
            EXPECT_DOUBLE_EQ(PSOPT::unscale_variable(x, sc, 0.0), x/sc);

            for (double sh : shs) {
                const double back =
                    PSOPT::unscale_variable(PSOPT::scale_variable(x, sc, sh), sc, sh);
                EXPECT_NEAR(back, x, 1.0e-9*std::max(1.0, std::fabs(x)));
            }
        }
}


// ---------------------------------------------------------------------------
// A bound that is absent stays absent under either map. This is the convention
// of PSOPT::bound_inf, and a variable bound is as entitled to it as a
// constraint row is -- a state declared unbounded by the documented means used
// to arrive here as a finite bound of 1e19 and be scaled by 1e-19.
// ---------------------------------------------------------------------------

TEST(VariableScaling, AnAbsentBoundSurvivesTheMap)
{
    EXPECT_EQ(PSOPT::scaled_lower_bound(-PSOPT::bound_inf, 3.0, 5.0), -PSOPT::inf);
    EXPECT_EQ(PSOPT::scaled_upper_bound( PSOPT::bound_inf, 3.0, 5.0),  PSOPT::inf);
    EXPECT_EQ(PSOPT::scaled_lower_bound(-1.0e20, 3.0, 5.0), -PSOPT::inf);
    EXPECT_EQ(PSOPT::scaled_upper_bound( 1.0e20, 3.0, 5.0),  PSOPT::inf);

    // A bound that is present moves with the origin.
    EXPECT_DOUBLE_EQ(PSOPT::scaled_lower_bound(1000.0, 0.005, 1100.0), -0.5);
    EXPECT_DOUBLE_EQ(PSOPT::scaled_upper_bound(1200.0, 0.005, 1100.0),  0.5);
}


// ---------------------------------------------------------------------------
// The same problem under both maps. The answer is known in closed form, so this
// is not merely a comparison of two runs against each other: each is checked
// against 12 D^2 / T^3 as well.
// ---------------------------------------------------------------------------

TEST(VariableScaling, AnOffsetBoxGivesTheSameAnswerUnderBothMaps)
{
    const varscale::Run mult = varscale::solve("automatic");
    const varscale::Run aff  = varscale::solve("affine");

    ASSERT_EQ(mult.flag, 0) << "multiplicative: " << mult.message;
    ASSERT_EQ(aff.flag,  0) << "affine: "         << aff.message;

    EXPECT_NEAR(mult.cost, varscale::J_EXACT, 1.0e-5*varscale::J_EXACT);
    EXPECT_NEAR(aff.cost,  varscale::J_EXACT, 1.0e-5*varscale::J_EXACT);

    // And to each other, more tightly than either is to the exact answer: the
    // two runs discretise the same problem on the same mesh and differ only in
    // the coordinates the NLP sees.
    EXPECT_NEAR(aff.cost, mult.cost, 1.0e-6*std::fabs(mult.cost));
}
