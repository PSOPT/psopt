//////////////////////////////////////////////////////////////////////////////
// test_ir_flexible_mesh.cpp
//
// The flexible mesh, and the limitation it exists to remove.
//
// A degree-d element polynomial cannot represent a discontinuity in its
// interior. On a fixed mesh a switching time that falls inside an element is
// therefore a wall: the transcription converges, and it converges to the wrong
// answer, and tightening the residual tolerance does not move it -- because the
// residual it is driving down is the residual of a problem the mesh cannot
// state. The flexible mesh makes the element boundaries decision variables so
// the optimisation can put a boundary ON the switching time (Nie & Kerrigan,
// 2022).
//
// The test problem is chosen so that the fixed mesh cannot accidentally win. A
// minimum-time double integrator from (0,0) to (1,0) switches at exactly tf/2 by
// symmetry, which on a uniform mesh of an even number of elements is already a
// boundary -- a test the fixed mesh passes for a reason that has nothing to do
// with the mesh being good. Making the control bounds asymmetric, u in [-1,2],
// moves the switch to tf/3:
//
//     tf* = sqrt(3),   switch at t = sqrt(3)/3,
//
// and no uniform partition into four elements has a boundary at one third.
//////////////////////////////////////////////////////////////////////////////

#include "gtest/gtest.h"
#include <psopt.h>

#include <cmath>

namespace irflex {

const double TF_EXACT = 1.7320508075688772;   // sqrt(3)

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

struct Run { int flag; double tf; double rel_err; };

static Run solve(bool flexible, double residual_bound)
{
    Alg algorithm; Sol solution; Prob problem;
    Run out; out.flag = -1; out.tf = 0.0; out.rel_err = 1.0;

    const int nodes = 17;                      // 16 intervals = 4 elements of degree 4

    problem.name        = "flexible mesh";
    problem.outfilename = flexible ? "test_ir_flex_on.txt" : "test_ir_flex_off.txt";
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
    problem.phases(1).bounds.lower.controls(0) = -1.0;   // asymmetric: the switch
    problem.phases(1).bounds.upper.controls(0) =  2.0;   // lands at tf/3, not tf/2
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
    algorithm.mesh_refinement       = "manual";
    algorithm.collocation_method    = "Hermite-Simpson";
    algorithm.transcription_method  = "integrated-residual";
    algorithm.ir_local_order        = 4;
    algorithm.ir_residual_nodes     = 6;
    // Minimise the cost subject to a residual box. The default minimises the
    // residual itself, and on a free horizon the smoothest trajectory is the
    // longest one, so a minimum-time problem posed that way simply runs tf to
    // its upper bound -- which is correct behaviour for a question nobody meant
    // to ask.
    algorithm.ir_objective          = "cost";
    algorithm.ir_residual_bound     = residual_bound;
    algorithm.ir_flexible_mesh      = flexible;

    out.flag = psopt(solution, problem, algorithm);
    if (out.flag == 0) {
        MatrixXd t = solution.get_time_in_phase(1);
        out.tf      = t(0, (int) t.cols() - 1);
        out.rel_err = std::fabs(out.tf - TF_EXACT)/TF_EXACT;
    }
    return out;
}

} // namespace irflex


// ---------------------------------------------------------------------------
// The wall. A fixed mesh converges to an answer that is wrong by most of a per
// cent, and tightening the residual box by three orders of magnitude does not
// improve it -- which is the signature of a discretisation that cannot state the
// problem, as opposed to one that has not yet been solved accurately enough.
// ---------------------------------------------------------------------------

TEST(IRFlexibleMesh, AFixedMeshCannotResolveASwitchInsideAnElement)
{
    const irflex::Run loose  = irflex::solve(false, 1.0e-4);
    const irflex::Run tight  = irflex::solve(false, 1.0e-7);

    ASSERT_EQ(loose.flag, 0);
    ASSERT_EQ(tight.flag, 0);

    EXPECT_GT(loose.rel_err, 1.0e-3) << "tf = " << loose.tf;
    EXPECT_GT(tight.rel_err, 1.0e-3) << "tf = " << tight.tf;

    // And the error is essentially the same at both tolerances: the residual box
    // is not what is limiting it.
    EXPECT_LT(std::fabs(tight.rel_err - loose.rel_err), 0.25*loose.rel_err);
}


// ---------------------------------------------------------------------------
// And the mesh that can move gets past it, by three orders of magnitude, on the
// same nodes, the same element degree and the same residual box.
// ---------------------------------------------------------------------------

TEST(IRFlexibleMesh, AFlexibleMeshResolvesIt)
{
    const irflex::Run fixed = irflex::solve(false, 1.0e-4);
    const irflex::Run flex  = irflex::solve(true,  1.0e-4);

    ASSERT_EQ(fixed.flag, 0);
    ASSERT_EQ(flex.flag,  0) << "the flexible mesh failed to solve";

    EXPECT_LT(flex.rel_err, 1.0e-4) << "tf = " << flex.tf;
    EXPECT_LT(flex.rel_err, 0.01*fixed.rel_err)
        << "flexible " << flex.rel_err << " against fixed " << fixed.rel_err;
}
