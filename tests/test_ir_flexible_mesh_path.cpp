//////////////////////////////////////////////////////////////////////////////
// test_ir_flexible_mesh_path.cpp
//
// A path constraint has to be imposed at the time the SOLVED mesh puts its node
// at, not at the time the uniform mesh the phase started from would have put it.
//
// On a fixed mesh the distinction does not exist, which is why it went unnoticed
// while the mesh was fixed: the node positions the path rows read are the node
// positions the trajectory has. Once the element boundaries become decision
// variables the two part company, and a path row evaluated at the stale time
// constrains a state that belongs somewhere else -- an extra constraint the user
// never wrote, silently, and with no symptom but a worse answer.
//
// The problem below is built so that symptom is unmissable. Minimum time for a
// double integrator from (0,0) to (1,0) with u in [-2,1]: accelerate at +1 to
// t1 = 2tf/3, then decelerate at -2, so
//
//     tf* = sqrt(3),   switch at t = 2 sqrt(3)/3,   and v(t) = t on the first arc.
//
// The path constraint is v - t <= 0. It holds with EQUALITY on the whole first
// arc and strictly on the second, so it does not bind: the answer with it must
// be the answer without it.
//
// Nine nodes at degree four is TWO elements, so there is exactly one interior
// boundary and only one useful place for it -- the switch, at tau = +1/3, which
// is to the RIGHT of the uniform boundary at tau = 0. The first element must
// therefore widen, every interior node of it moves to a LATER time, and a path
// row at the stale time demands v_k <= t_uniform,k when the truth is v_k = t_k >
// t_uniform,k. With one interior boundary the solver cannot dodge by putting the
// switch on some other boundary. It has to give up the mesh, and the flexible
// mesh's entire advantage with it.
//
// Measured, with the path rows reading the stored mesh:
//
//     fixed mesh,    no path     rel err 3.28e-02
//     fixed mesh,    with path   rel err 3.28e-02      (tangent: costs nothing)
//     flexible mesh, no path     rel err 1.40e-07
//     flexible mesh, with path   rel err 3.28e-02      <-- the whole gain, gone
//
// and with them reading the solved one, the last line is 1.09e-07.
//////////////////////////////////////////////////////////////////////////////

#include "gtest/gtest.h"
#include <psopt.h>

#include <cmath>

namespace irflexpath {

const double TF_EXACT = 1.7320508075688772;   // sqrt(3)

adouble endpoint_cost(adouble*, adouble*, adouble*, adouble&, adouble& tf, adouble*,
                      int, Workspace*)
{ return tf; }

adouble integrand_cost(adouble*, adouble*, adouble*, adouble&, adouble*, int, Workspace*)
{ return 0.0; }

void dae(adouble* d, adouble* path, adouble* s, adouble* c, adouble*, adouble& t,
         adouble*, int, Workspace*)
{ d[0] = s[1]; d[1] = c[0]; path[0] = s[1] - t; }

void events(adouble* e, adouble* i, adouble* f, adouble*, adouble&, adouble&,
            adouble*, int, Workspace*)
{ e[0] = i[0]; e[1] = i[1]; e[2] = f[0]; e[3] = f[1]; }

void linkages(adouble*, adouble*, Workspace*) {}

struct Run { int flag; double tf; double rel_err; };

// The path constraint is present in both variants -- same rows, same layout, same
// Jacobian -- and only its upper bound changes. `active` gives it the tangent bound
// 0; otherwise the bound is slack by miles and the row cannot do anything. That
// keeps the comparison to the one thing being tested.
static Run solve(bool flexible, bool active)
{
    Alg algorithm; Sol solution; Prob problem;
    Run out; out.flag = -1; out.tf = 0.0; out.rel_err = 1.0;

    const int nodes = 9;                       // 8 intervals = 2 elements of degree 4

    problem.name        = "flexible mesh path";
    problem.outfilename = "test_ir_flex_path.txt";
    problem.nphases     = 1;
    problem.nlinkages   = 0;
    psopt_level1_setup(problem);

    problem.phases(1).nstates   = 2;
    problem.phases(1).ncontrols = 1;
    problem.phases(1).nevents   = 4;
    problem.phases(1).npath     = 1;
    problem.phases(1).nodes     << nodes;
    psopt_level2_setup(problem, algorithm);

    problem.phases(1).bounds.lower.states   << -2.0, -2.0;
    problem.phases(1).bounds.upper.states   <<  2.0,  2.0;
    problem.phases(1).bounds.lower.controls(0) = -2.0;   // the switch lands at 2tf/3,
    problem.phases(1).bounds.upper.controls(0) =  1.0;   // right of the uniform boundary
    problem.phases(1).bounds.lower.path(0)     = -6.0;
    problem.phases(1).bounds.upper.path(0)     = active ? 0.0 : 6.0;
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
    algorithm.ir_objective          = "cost";
    algorithm.ir_residual_bound     = 1.0e-8;
    algorithm.ir_flexible_mesh      = flexible;

    out.flag = psopt(solution, problem, algorithm);
    if (out.flag == 0) {
        MatrixXd t = solution.get_time_in_phase(1);
        out.tf      = t(0, (int) t.cols() - 1);
        out.rel_err = std::fabs(out.tf - TF_EXACT)/TF_EXACT;
    }
    return out;
}

} // namespace irflexpath


// ---------------------------------------------------------------------------
// The control. On a fixed mesh the tangent constraint costs nothing, because
// there the node times the path rows read and the node times the trajectory has
// are the same times. If this test ever fails, the constraint is not tangent and
// the test below is measuring something else.
// ---------------------------------------------------------------------------

TEST(IRFlexibleMeshPath, ATangentPathConstraintCostsNothingOnAFixedMesh)
{
    const irflexpath::Run without = irflexpath::solve(false, false);
    const irflexpath::Run with    = irflexpath::solve(false, true);

    ASSERT_EQ(without.flag, 0);
    ASSERT_EQ(with.flag,    0);

    EXPECT_LT(std::fabs(with.tf - without.tf), 1.0e-6)
        << "with " << with.tf << " against without " << without.tf;
}


// ---------------------------------------------------------------------------
// And it must cost nothing on a flexible mesh either. It does only if the rows
// are evaluated at the solved node times; read from the stored mesh they throw
// away the whole of the flexible mesh's advantage, five orders of magnitude on
// this problem.
// ---------------------------------------------------------------------------

TEST(IRFlexibleMeshPath, ATangentPathConstraintCostsNothingOnAFlexibleMeshEither)
{
    const irflexpath::Run without = irflexpath::solve(true, false);
    const irflexpath::Run with    = irflexpath::solve(true, true);

    ASSERT_EQ(without.flag, 0);
    ASSERT_EQ(with.flag,    0) << "the flexible mesh failed to solve with the path row";

    EXPECT_LT(without.rel_err, 1.0e-6) << "tf = " << without.tf;
    EXPECT_LT(with.rel_err,    1.0e-6) << "tf = " << with.tf
        << " -- the path rows are being imposed at the wrong times";
}
