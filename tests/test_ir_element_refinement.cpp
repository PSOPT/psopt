//////////////////////////////////////////////////////////////////////////////
// test_ir_element_refinement.cpp
//
// Automatic mesh refinement for the integrated-residual transcription, which
// refines ELEMENTS rather than inserting nodes.
//
// This file replaces a test that pinned the opposite. Betts refinement inserts
// equally spaced points inside the intervals it flags, and the element basis is
// evaluated through matrices built once for the reference LGL abscissae of the
// unit element, so after such a refinement the residual is a polynomial
// evaluated off its own nodes -- silently, since nothing downstream checks where
// a node sits -- and the divisibility rule norder % d == 0 stops holding. The
// combination was refused. It is now supported, by refining in the currency the
// transcription is written in: an element is kept or split into equal
// sub-elements, so the node count moves in multiples of the stride and every
// element that results is a proper one.
//
// The problem is the minimum-time double integrator with u in [-1, 1.5], whose
// switch falls at 0.4 t_f. The more usual u in [-1, 2] puts it at t_f/3, which
// bisecting the initial two-element mesh into three lands on EXACTLY -- so a
// fixed mesh passes that test by arithmetic luck rather than by refining well.
// At 0.4 no sequence of equal splits of a uniform mesh reaches it.
//
//     tf* = 1.825741858,   switch at 0.730296743
//////////////////////////////////////////////////////////////////////////////

#include "gtest/gtest.h"
#include <psopt.h>

#include <cmath>
#include <string>

namespace irrefine {

const double TF_EXACT = 1.825741858350554;      // u in [-1, 1.5]

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

struct Run { int flag; int nodes_out; double tf; double rel_err; };

static Run solve(int ir_order, bool flexible, const std::string& refinement, bool dair)
{
    Alg algorithm; Sol solution; Prob problem;
    Run out; out.flag = -1; out.nodes_out = 0; out.tf = 0.0; out.rel_err = 1.0;

    const int nodes = 9;                       // at d=4: two elements of degree four

    problem.name        = "ir element refinement";
    problem.outfilename = "test_ir_element_refinement.txt";
    problem.nphases     = 1;
    problem.nlinkages   = 0;
    psopt_level1_setup(problem);

    problem.phases(1).nstates   = 2;
    problem.phases(1).ncontrols = 1;
    problem.phases(1).nevents   = 4;
    problem.phases(1).npath     = 0;
    problem.phases(1).nodes     << nodes;
    psopt_level2_setup(problem, algorithm);

    problem.phases(1).bounds.lower.states   << -2.0, -2.0;
    problem.phases(1).bounds.upper.states   <<  2.0,  2.0;
    problem.phases(1).bounds.lower.controls(0) = -1.0;
    problem.phases(1).bounds.upper.controls(0) =  1.5;   // switch at 0.4 tf
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

    problem.phases(1).guess.states   = zeros(2, nodes);
    problem.phases(1).guess.states.row(0) = linspace(0.0, 1.0, nodes);
    problem.phases(1).guess.controls = zeros(1, nodes);
    problem.phases(1).guess.time     = linspace(0.0, 1.73, nodes);

    algorithm.nlp_method            = "IPOPT";
    algorithm.scaling               = "automatic";
    algorithm.derivatives           = "automatic";
    algorithm.nlp_iter_max          = 2000;
    algorithm.nlp_tolerance         = 1.0e-8;
    algorithm.print_level           = 0;
    algorithm.mesh_refinement       = refinement;
    algorithm.mr_max_iterations     = 5;
    algorithm.mr_max_growth_factor  = 0.5;
    // The element refinement is driven by the local ODE error, and on this problem the
    // coarse mesh already meets a loose tolerance, so a loose one would stop at the first
    // iteration and the test would measure nothing.
    algorithm.ode_tolerance         = 1.0e-12;
    algorithm.collocation_method    = "Hermite-Simpson";
    algorithm.transcription_method  = "integrated-residual";
    algorithm.ir_local_order        = ir_order;
    algorithm.ir_residual_nodes     = ( ir_order >= 2 ) ? ir_order + 2 : 4;
    algorithm.ir_objective          = "cost";
    algorithm.ir_residual_bound     = 1.0e-6;
    algorithm.ir_flexible_mesh      = flexible;
    if (dair) { algorithm.ir_dair = true; algorithm.ir_dair_delta_factor = 1.0; }

    out.flag = psopt(solution, problem, algorithm);
    if (out.flag == 0) {
        MatrixXd t = solution.get_time_in_phase(1);
        out.nodes_out = (int) t.cols();
        out.tf        = t(0, out.nodes_out - 1);
        out.rel_err   = std::fabs(out.tf - TF_EXACT)/TF_EXACT;
    }
    return out;
}

} // namespace irrefine


// ---------------------------------------------------------------------------
// The combination that used to be refused now runs, and refining is worth two
// and a half orders of magnitude on a problem the coarse mesh cannot state.
//
// The node count is the invariant the old route broke: an element basis of
// degree d needs the interval count to stay a multiple of d, and splitting
// elements keeps it one by construction where inserting nodes did not.
// ---------------------------------------------------------------------------

TEST(IRElementRefinement, AnElementBasisCanNowBeRefinedAutomatically)
{
    const irrefine::Run coarse = irrefine::solve(4, false, "manual",    false);
    const irrefine::Run refine = irrefine::solve(4, false, "automatic", false);

    ASSERT_EQ(coarse.flag, 0);
    ASSERT_EQ(refine.flag, 0) << "automatic refinement with the element basis failed";

    EXPECT_GT(refine.nodes_out, coarse.nodes_out) << "the mesh did not grow";
    EXPECT_EQ((refine.nodes_out - 1) % 4, 0)
        << "the interval count is no longer a multiple of the local order: "
        << refine.nodes_out - 1;
    EXPECT_LT(refine.rel_err, 0.01*coarse.rel_err)
        << "refined " << refine.rel_err << " against coarse " << coarse.rel_err;
}


// ---------------------------------------------------------------------------
// And it composes with the flexible mesh, which is the combination it was built
// for. The two are given disjoint jobs: the error estimator says how many
// elements the phase needs, and the flexible mesh says where they go. Paired
// with DAIR, so that the residual box tightens as the mesh does, this is the
// most accurate of everything measured on this problem -- 4.1e-8, against
// 1.6e-7 for the flexible mesh alone and 1.5e-2 for the coarse fixed mesh.
//
// Letting the estimator choose WHERE as well makes it worse rather than better,
// which is why it does not: once a discontinuity is isolated inside a thin
// element, that element's local error stays large however thin it is, so the
// estimator flags it every iteration and the refinement starves the rest of the
// trajectory to feed a point that was already handled.
// ---------------------------------------------------------------------------

TEST(IRElementRefinement, ItComposesWithTheFlexibleMesh)
{
    const irrefine::Run flex_only = irrefine::solve(4, true, "manual",    false);
    const irrefine::Run both      = irrefine::solve(4, true, "automatic", true);

    ASSERT_EQ(flex_only.flag, 0);
    ASSERT_EQ(both.flag,      0) << "automatic refinement with the flexible mesh failed";

    EXPECT_LT(flex_only.rel_err, 1.0e-5) << "tf = " << flex_only.tf;
    EXPECT_LT(both.rel_err,      1.0e-6) << "tf = " << both.tf;
    EXPECT_EQ((both.nodes_out - 1) % 4, 0)
        << "the interval count is no longer a multiple of the local order";
}


// ---------------------------------------------------------------------------
// The two things the new route must not take over: manual refinement, and the
// legacy cubic-Hermite form on a fixed mesh, which keeps the Betts refinement it
// has always used. A route that is too wide changes answers nobody asked to
// have changed.
// ---------------------------------------------------------------------------

TEST(IRElementRefinement, TheRouteIsNoWiderThanItNeedsToBe)
{
    Alg a;
    a.mesh_refinement      = "automatic";
    a.transcription_method = "integrated-residual";
    a.ir_local_order       = 0;
    a.ir_flexible_mesh     = false;
    EXPECT_FALSE(ir_element_refinement_active(a)) << "cubic Hermite, fixed mesh";

    a.ir_flexible_mesh = true;
    EXPECT_TRUE(ir_element_refinement_active(a))  << "cubic Hermite, flexible mesh";

    a.ir_flexible_mesh = false;
    a.ir_local_order   = 4;
    EXPECT_TRUE(ir_element_refinement_active(a))  << "element basis, fixed mesh";

    a.mesh_refinement = "manual";
    EXPECT_FALSE(ir_element_refinement_active(a)) << "manual refinement";

    a.mesh_refinement      = "automatic";
    a.transcription_method = "collocation";
    EXPECT_FALSE(ir_element_refinement_active(a)) << "collocation";
}
