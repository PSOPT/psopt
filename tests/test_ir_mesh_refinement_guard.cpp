//////////////////////////////////////////////////////////////////////////////
// test_ir_mesh_refinement_guard.cpp
//
// Automatic mesh refinement and the Nie-Kerrigan element basis do not compose,
// and the combination is refused rather than run.
//
// Betts local refinement inserts equally spaced points inside the intervals it
// flags. The element basis is evaluated through matrices built once for the
// reference LGL abscissae of the unit element, so after a refinement the nodes
// of an element are no longer where the basis expects them and the residual is
// a polynomial evaluated off its own nodes. Nothing downstream checks where a
// node sits, so the failure is silent in principle. In practice it crashed --
// in sort_vector(), because this branch built snodes as a column where every
// other producer builds a row -- which is luck rather than a diagnostic.
//
// Both halves are fixed: the shape, and the combination. This test pins the
// second. The first is pinned by the fact that the refusal happens at all: a
// column-shaped snodes array would have died before validate could speak.
//////////////////////////////////////////////////////////////////////////////

#include "gtest/gtest.h"
#include <psopt.h>

#include <string>

namespace irguard {

adouble endpoint_cost(adouble*, adouble*, adouble*, adouble&, adouble& tf, adouble*,
                      int, Workspace*)
{ return tf; }

adouble integrand_cost(adouble*, adouble*, adouble*, adouble&, adouble*, int, Workspace*)
{ return 0.0; }

void dae(adouble* d, adouble*, adouble* s, adouble* c, adouble*, adouble&,
         adouble*, int, Workspace*)
{ d[0] = s[1]; d[1] = c[0]; }

void events(adouble* e, adouble* i, adouble* f, adouble*, adouble&, adouble&,
            adouble*, int, Workspace*)
{ e[0] = i[0]; e[1] = i[1]; e[2] = f[0]; e[3] = f[1]; }

void linkages(adouble*, adouble*, Workspace*) {}

// Returns the value of error_flag: nonzero means psopt refused or failed.
static int run(const std::string& refinement, int ir_order)
{
    Alg algorithm; Sol solution; Prob problem;

    problem.name        = "ir refinement guard";
    problem.outfilename = "test_ir_guard.txt";
    problem.nphases     = 1;
    problem.nlinkages   = 0;
    psopt_level1_setup(problem);

    problem.phases(1).nstates   = 2;
    problem.phases(1).ncontrols = 1;
    problem.phases(1).nevents   = 4;
    problem.phases(1).npath     = 0;
    problem.phases(1).nodes     << 25;
    psopt_level2_setup(problem, algorithm);

    problem.phases(1).bounds.lower.states   << -5.0, -5.0;
    problem.phases(1).bounds.upper.states   <<  5.0,  5.0;
    problem.phases(1).bounds.lower.controls(0) = -1.0;
    problem.phases(1).bounds.upper.controls(0) =  1.0;
    problem.phases(1).bounds.lower.events   << 0.0, 0.0, 1.0, 0.0;
    problem.phases(1).bounds.upper.events   << 0.0, 0.0, 1.0, 0.0;
    problem.phases(1).bounds.lower.StartTime = 0.0;
    problem.phases(1).bounds.upper.StartTime = 0.0;
    problem.phases(1).bounds.lower.EndTime   = 0.5;
    problem.phases(1).bounds.upper.EndTime   = 6.0;

    problem.integrand_cost = &integrand_cost;
    problem.endpoint_cost  = &endpoint_cost;
    problem.dae            = &dae;
    problem.events         = &events;
    problem.linkages       = &linkages;

    problem.phases(1).guess.states   = zeros(2, 25);
    problem.phases(1).guess.states.row(0) = linspace(0.0, 1.0, 25);
    problem.phases(1).guess.controls = zeros(1, 25);
    problem.phases(1).guess.time     = linspace(0.0, 2.0, 25);

    algorithm.nlp_method            = "IPOPT";
    algorithm.scaling               = "automatic";
    algorithm.derivatives           = "automatic";
    algorithm.nlp_iter_max          = 500;
    algorithm.nlp_tolerance         = 1.0e-6;
    algorithm.print_level           = 0;
    algorithm.mesh_refinement       = refinement;
    algorithm.collocation_method    = "Legendre";
    if (ir_order >= 2) {
        algorithm.transcription_method = "integrated-residual";
        algorithm.ir_local_order       = ir_order;
        // Patch 156's rule, unrelated to the guard under test but enforced by the
        // same validate pass, so the test has to satisfy it to reach the guard.
        algorithm.ir_residual_nodes    = ir_order + 2;
        algorithm.collocation_method   = "Hermite-Simpson";
    }

    return psopt(solution, problem, algorithm);
}

} // namespace irguard


// ---------------------------------------------------------------------------
// The combination is refused. psopt returns nonzero rather than producing a
// number computed on a mesh the basis does not describe.
// ---------------------------------------------------------------------------

TEST(IRMeshRefinementGuard, AutomaticRefinementWithAnElementBasisIsRefused)
{
    EXPECT_NE(irguard::run("automatic", 4), 0);
}


// ---------------------------------------------------------------------------
// And the two things it must not refuse: automatic refinement without the
// element basis, and the element basis without automatic refinement. A guard
// that is too wide is as much a defect as one that is missing.
// ---------------------------------------------------------------------------

TEST(IRMeshRefinementGuard, TheGuardIsNoWiderThanItNeedsToBe)
{
    EXPECT_EQ(irguard::run("manual",    4), 0) << "element basis, manual refinement";
    EXPECT_EQ(irguard::run("automatic", 0), 0) << "automatic refinement, no element basis";
}
