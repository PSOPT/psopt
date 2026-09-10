//////////////////////////////////////////////////////////////////////////////
// test_constraint_coverage.cpp
//
// get_ncons_phase_i counts the NLP constraint rows; gg_ad writes them. They are
// two independent statements of one layout and they have disagreed three times:
// the integrated-residual midpoint path rows (160 rows counted and never
// written), the Gauss terminal point, and the row past the end of the scaling
// vector under user scaling. None announced itself, because the constraint
// buffer is zero-filled and zero is a plausible constraint value: an equality
// row nobody writes is enforced as 0 = 0 and the solve converges to a problem
// that was not posed.
//
// check_constraint_coverage runs once per mesh iteration: it poisons the
// constraint buffer, evaluates the constraints, and requires that nothing is
// left as it found it. These tests cover the two halves of that -- that the
// poison can be told from any value a constraint could legitimately take, and
// that every transcription writes every row it counts.
//////////////////////////////////////////////////////////////////////////////

#include "gtest/gtest.h"
#include <psopt.h>

#include <cmath>
#include <limits>
#include <string>
#include <vector>

// --------------------------------------------------------------------------
// The discrimination the check rests on.
// --------------------------------------------------------------------------

TEST(ConstraintCoverage, PoisonIsDistinguishedFromAnyPlausibleConstraintValue)
{
    const double poison = psopt_constraint_poison_value();

    EXPECT_TRUE(psopt_constraint_row_unwritten(poison))
        << "the poison itself must be recognised";

    // The scaling passes multiply every row by its scale factor before anything
    // reads it, so an unwritten row does not arrive as the poison but as the
    // poison times that factor. Both a very small factor and one that overflows
    // the product to infinity have to stay recognisable.
    EXPECT_TRUE(psopt_constraint_row_unwritten(poison * 1.0e-200))
        << "a heavily down-scaled poisoned row must still be recognised";
    EXPECT_TRUE(psopt_constraint_row_unwritten(poison * 1.0e10))
        << "a poisoned row whose product overflowed must still be recognised";
    EXPECT_TRUE(psopt_constraint_row_unwritten(-std::numeric_limits<double>::infinity()));

    // And nothing a constraint can legitimately be may be mistaken for it. The
    // largest magnitude PSOPT itself treats as a finite constraint value is
    // psopt_inf, 1e19.
    const double ordinary[] = { 0.0, -0.0, 1.0, -1.0, 1.0e-30, -1.0e-30,
                                1.0e19, -1.0e19, 1.0e20, -1.0e20 };
    for (size_t k = 0; k < sizeof ordinary / sizeof ordinary[0]; k++)
        EXPECT_FALSE(psopt_constraint_row_unwritten(ordinary[k]))
            << "a constraint value of " << ordinary[k] << " was read as unwritten";

    // NaN is a different fault, reported differently, and must not be swept into
    // this one: an unwritten row is a defect in PSOPT, whereas a NaN can come out
    // of the user's own model at a bad initial guess.
    EXPECT_FALSE(psopt_constraint_row_unwritten(
                     std::numeric_limits<double>::quiet_NaN()))
        << "NaN must not be reported as an unwritten row";
}


// --------------------------------------------------------------------------
// Every transcription writes every row it counts.
//
// A minimum-energy double integrator with fixed endpoints, four events and one
// path constraint, so that the phase block carries defect rows, event rows and
// path rows at once. The path constraint is inactive, which is deliberate: an
// active one would be held at its bound and a row left at zero could not be
// told from a row held at zero.
// --------------------------------------------------------------------------

namespace cov {

adouble endpoint_cost(adouble*, adouble*, adouble*, adouble&, adouble&, adouble*,
                      int, Workspace*)
{ return 0.0; }

adouble integrand_cost(adouble*, adouble* c, adouble*, adouble&, adouble*, int,
                       Workspace*)
{ adouble u = c[0]; return 0.5*u*u; }

void dae(adouble* d, adouble* path, adouble* s, adouble* c, adouble*, adouble&,
         adouble*, int, Workspace*)
{
    d[0] = s[1];
    d[1] = c[0];
    path[0] = s[1];            // inactive: bounded well away from the solution
}

void events(adouble* e, adouble* i, adouble* f, adouble*, adouble&, adouble&,
            adouble*, int, Workspace*)
{ e[0] = i[0]; e[1] = i[1]; e[2] = f[0]; e[3] = f[1]; }

void linkages(adouble*, adouble*, Workspace*) {}

// Solve under one collocation method and return what psopt() reported. The guard
// is fatal, so a counted-but-unwritten row arrives as a thrown ErrorHandler and
// therefore as error_flag = 1 with the diagnostic in solution.error_msg.
static int solve_under(const std::string& method, int nodes, std::string& message,
                       double& cost)
{
    Alg algorithm; Sol solution; Prob problem;

    problem.name        = "coverage: min-energy double integrator";
    problem.outfilename = "test_constraint_coverage.txt";
    problem.nphases     = 1;
    problem.nlinkages   = 0;
    psopt_level1_setup(problem);

    problem.phases(1).nstates   = 2;
    problem.phases(1).ncontrols = 1;
    problem.phases(1).nevents   = 4;
    problem.phases(1).npath     = 1;
    problem.phases(1).nodes     << nodes;
    psopt_level2_setup(problem, algorithm);

    problem.phases(1).bounds.lower.states   << -5.0, -5.0;
    problem.phases(1).bounds.upper.states   <<  5.0,  5.0;
    problem.phases(1).bounds.lower.controls(0) = -50.0;
    problem.phases(1).bounds.upper.controls(0) =  50.0;
    problem.phases(1).bounds.lower.path(0)     = -4.0;
    problem.phases(1).bounds.upper.path(0)     =  4.0;
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

    problem.phases(1).guess.states          = zeros(2, nodes);
    problem.phases(1).guess.states.row(0)   = linspace(0.0, 1.0, nodes);
    problem.phases(1).guess.controls        = zeros(1, nodes);
    problem.phases(1).guess.time            = linspace(0.0, 1.0, nodes);

    algorithm.nlp_method         = "IPOPT";
    algorithm.scaling            = "automatic";
    algorithm.derivatives        = "automatic";
    algorithm.collocation_method = method;
    algorithm.nlp_iter_max       = 1000;
    algorithm.nlp_tolerance      = 1.0e-6;
    algorithm.print_level        = 0;

    const int flag = psopt(solution, problem, algorithm);
    message = solution.error_msg;
    cost    = solution.cost;
    return flag;
}

} // namespace cov

TEST(ConstraintCoverage, EveryRowIsWrittenUnderEveryCollocationMethod)
{
    // Every transcription PSOPT offers. Each lays the phase block out differently
    // -- Radau appends a terminal-control pin, Gauss appends one quadrature
    // defining constraint per interval, the LGL methods append interface defects
    // on an hp mesh -- and each of those appendices is a place where the count
    // and the assignments can part company.
    const char* methods[] = { "Legendre", "Chebyshev", "Radau", "Gauss",
                              "trapezoidal", "Hermite-Simpson" };

    for (size_t k = 0; k < sizeof methods / sizeof methods[0]; k++) {
        std::string message;
        double cost = 0.0;
        const int flag = cov::solve_under(methods[k], 30, message, cost);

        EXPECT_EQ(flag, 0)
            << "collocation_method = \"" << methods[k]
            << "\" did not complete; diagnostic: " << message;
        EXPECT_EQ(message.find("poison"), std::string::npos)
            << "collocation_method = \"" << methods[k]
            << "\" counted a constraint row that nothing writes: " << message;

        // The transcription is exact on this problem -- u* = 6 - 12t, states
        // cubic, integrand quadratic -- so a run that reached the answer is a run
        // that had the constraints it was supposed to have. J* = 6.
        if (flag == 0)
            EXPECT_NEAR(cost, 6.0, 5.0e-2)
                << "collocation_method = \"" << methods[k] << "\"";
    }
}
