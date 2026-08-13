// ---------------------------------------------------------------------------
// examples/lotka2 -- Lotka-Volterra fishing with TWO independent binary controls
//
// The classical fishing problem lets a single fleet act on both species at once.
// Here the two species are fished independently: w1 switches the effort applied to
// the prey and w2 the effort applied to the predator, each in {0,1}. The phase
// therefore declares two integer controls, and PSOPT convexifies over the product
// of their admissible sets: four weight-controls and one SOS1 constraint, all
// added internally, with the dynamics written once in terms of w1 and w2 as
// ordinary controls.
//
// The point of the example is that independent switching is strictly more
// capable than joint switching. The single-control problem of examples/lotka is
// the special case w1 == w2, so its optimum is an upper bound for this one; the
// run prints both so the reader can see the gap.
//
// Reconstruction uses reconstruct_integer_controls, which rounds over the product
// modes and then decodes each chosen mode into a value for each declared control.
// ---------------------------------------------------------------------------


#include "psopt.h"
#include <vector>
#include "integer_controls.h"
#include "sum_up_rounding.h"

using namespace std;

static const double C0 = 0.4;
static const double C1 = 0.2;
static const double PHI_STAR = 1.34408;   // mintOC relaxed optimum (lower bound)

//////////////////////////////////////////////////////////////////////////
///////////////////  Lotka-Volterra RHS at fishing value w ////////////////
//////////////////////////////////////////////////////////////////////////
// Templated on the scalar so the same expression serves the adouble dynamics
// (with w an ordinary control) and the double forward simulation.

template<class T>
inline void lv_rhs(const T& x0, const T& x1, const T& /*x2*/,
                   const T& w1, const T& w2, T d[3])
{
    d[0] =  x0 - x0*x1 - C0*x0*w1;   // prey fished by w1
    d[1] = -x1 + x0*x1 - C1*x1*w2;   // predator fished by w2
    d[2] = (x0 - 1.0)*(x0 - 1.0) + (x1 - 1.0)*(x1 - 1.0);
}

//////////////////////////////////////////////////////////////////////////
///////////////////  Cost, DAE, events, linkages //////////////////////////
//////////////////////////////////////////////////////////////////////////

adouble endpoint_cost(adouble* initial_states, adouble* final_states,
                      adouble* parameters, adouble& t0, adouble& tf,
                      adouble* xad, int iphase, Workspace* workspace)
{
    return final_states[2];
}

adouble integrand_cost(adouble* states, adouble* controls,
                       adouble* parameters, adouble& time, adouble* xad,
                       int iphase, Workspace* workspace)
{
    return 0.0;
}

// Dynamics written once, with both fishing controls as ordinary controls.
void dae(adouble* derivatives, adouble* path, adouble* states,
         adouble* controls, adouble* parameters, adouble& time,
         adouble* xad, int iphase, Workspace* workspace)
{
    adouble x0 = states[0];
    adouble x1 = states[1];
    adouble x2 = states[2];
    adouble w1 = controls[0];   // integer control: fishing effort on the prey
    adouble w2 = controls[1];   // integer control: fishing effort on the predator

    lv_rhs<adouble>(x0, x1, x2, w1, w2, derivatives);
}

void events(adouble* e, adouble* initial_states, adouble* final_states,
            adouble* parameters, adouble& t0, adouble& tf, adouble* xad,
            int iphase, Workspace* workspace)
{
    e[0] = initial_states[0];
    e[1] = initial_states[1];
    e[2] = initial_states[2];
}

void linkages(adouble* linkages, adouble* xad, Workspace* workspace)
{
}

//////////////////////////////////////////////////////////////////////////
///////////////////  Forward simulation of a fixed binary control /////////
//////////////////////////////////////////////////////////////////////////

static double simulate_integer_objective(const RowVectorXd& hwidth,
                                         const RowVectorXd& w1bin,
                                         const RowVectorXd& w2bin)
{
    double x[3] = {0.5, 0.7, 0.0};
    const int sub = 40;
    for (int i = 0; i < hwidth.size(); ++i) {
        const double a1 = w1bin(i);
        const double a2 = w2bin(i);
        const double dt = hwidth(i) / sub;
        for (int s = 0; s < sub; ++s) {
            double k1[3], k2[3], k3[3], k4[3], xt[3];
            lv_rhs<double>(x[0], x[1], x[2], a1, a2, k1);
            for (int j = 0; j < 3; ++j) xt[j] = x[j] + 0.5*dt*k1[j];
            lv_rhs<double>(xt[0], xt[1], xt[2], a1, a2, k2);
            for (int j = 0; j < 3; ++j) xt[j] = x[j] + 0.5*dt*k2[j];
            lv_rhs<double>(xt[0], xt[1], xt[2], a1, a2, k3);
            for (int j = 0; j < 3; ++j) xt[j] = x[j] + dt*k3[j];
            lv_rhs<double>(xt[0], xt[1], xt[2], a1, a2, k4);
            for (int j = 0; j < 3; ++j)
                x[j] += (dt/6.0)*(k1[j] + 2.0*k2[j] + 2.0*k3[j] + k4[j]);
        }
    }
    return x[2];
}

//////////////////////////////////////////////////////////////////////////
///////////////////  Solve at one mesh and round //////////////////////////
//////////////////////////////////////////////////////////////////////////

struct SweepResult {
    double relaxed;
    double integer;
    double gap;
    int    n_switches;
    int    n_split;      // intervals on which w1 and w2 differ
    int    status;
};

static SweepResult solve_and_round(int nnodes)
{
    Alg  algorithm;
    Sol  solution;
    Prob problem;

    problem.name        = "Lotka-Volterra fishing with two independent binary controls";
    problem.outfilename = "lotka2.txt";

    problem.nphases   = 1;
    problem.nlinkages = 0;
    psopt_level1_setup(problem);

    problem.phases(1).nstates   = 3;
    problem.phases(1).ncontrols = 2;   // two controls: prey effort, predator effort
    problem.phases(1).nevents   = 3;
    problem.phases(1).npath     = 0;   // SOS1 added internally by the declaration
    problem.phases(1).nodes     << nnodes;

    psopt_level2_setup(problem, algorithm);

    // Two integer controls. PSOPT convexifies over the product of the two
    // admissible sets: 2 x 2 = 4 weight-controls and one SOS1 constraint.
    RowVectorXd ivalues(2); ivalues << 0.0, 1.0;
    declare_integer_control(problem, 1, 0, ivalues);   // w1, prey
    declare_integer_control(problem, 1, 1, ivalues);   // w2, predator

    problem.phases(1).bounds.lower.states   << 0.0, 0.0, 0.0;
    problem.phases(1).bounds.upper.states   << 2.0, 2.0, 100.0;

    // Integer-control bounds are ignored; value range written for readability.
    problem.phases(1).bounds.lower.controls << 0.0, 0.0;
    problem.phases(1).bounds.upper.controls << 1.0, 1.0;

    problem.phases(1).bounds.lower.events   << 0.5, 0.7, 0.0;
    problem.phases(1).bounds.upper.events   << 0.5, 0.7, 0.0;

    problem.phases(1).bounds.lower.StartTime = 0.0;
    problem.phases(1).bounds.upper.StartTime = 0.0;
    problem.phases(1).bounds.lower.EndTime   = 12.0;
    problem.phases(1).bounds.upper.EndTime   = 12.0;

    problem.integrand_cost = &integrand_cost;
    problem.endpoint_cost  = &endpoint_cost;
    problem.dae            = &dae;
    problem.events         = &events;
    problem.linkages       = &linkages;

    RowVectorXd tg = linspace(0.0, 12.0, nnodes);

    // State guess by forward-simulating the relaxed midpoint fishing rate; the
    // integer-control guess is given in value terms (0.5 -> nearest admissible).
    DMatrix control_guess(2, nnodes);
    control_guess << 0.5 * ones(1, nnodes), 0.5 * ones(1, nnodes);

    DMatrix state_guess(3, nnodes);
    {
        double xs[3] = {0.5, 0.7, 0.0};
        const int sub = 20;
        for (int i = 0; i < nnodes; ++i) {
            state_guess(0, i) = xs[0];
            state_guess(1, i) = xs[1];
            state_guess(2, i) = xs[2];
            if (i == nnodes - 1) break;
            double dt = (tg(i+1) - tg(i)) / sub;
            for (int s = 0; s < sub; ++s) {
                double k1[3], k2[3], k3[3], k4[3], xt[3];
                lv_rhs<double>(xs[0], xs[1], xs[2], 0.5, 0.5, k1);
                for (int j = 0; j < 3; ++j) xt[j] = xs[j] + 0.5*dt*k1[j];
                lv_rhs<double>(xt[0], xt[1], xt[2], 0.5, 0.5, k2);
                for (int j = 0; j < 3; ++j) xt[j] = xs[j] + 0.5*dt*k2[j];
                lv_rhs<double>(xt[0], xt[1], xt[2], 0.5, 0.5, k3);
                for (int j = 0; j < 3; ++j) xt[j] = xs[j] + dt*k3[j];
                lv_rhs<double>(xt[0], xt[1], xt[2], 0.5, 0.5, k4);
                for (int j = 0; j < 3; ++j)
                    xs[j] += (dt/6.0)*(k1[j] + 2.0*k2[j] + 2.0*k3[j] + k4[j]);
            }
        }
    }

    problem.phases(1).guess.states   = state_guess;
    problem.phases(1).guess.controls = control_guess;
    problem.phases(1).guess.time     = tg;

    algorithm.nlp_method         = "IPOPT";
    algorithm.scaling            = "automatic";
    algorithm.derivatives        = "automatic";
    algorithm.nlp_iter_max       = 3000;
    algorithm.nlp_tolerance      = 1.e-7;
    algorithm.collocation_method = "trapezoidal";
    algorithm.mesh_refinement    = "manual";

    int status = psopt(solution, problem, algorithm);

    DMatrix x = solution.get_states_in_phase(1);
    DMatrix t = solution.get_time_in_phase(1);
    const int N = (int)t.cols();

    // One reconstruction per declared integer control, in declaration order. Rounding
    // is over the product modes; each chosen mode is then decoded into a value for w1
    // and a value for w2.
    std::vector<IntegerControlReconstruction> rec =
        reconstruct_integer_controls(solution, problem, 1);

    SweepResult r;
    r.relaxed    = x(2, N-1);
    r.integer    = simulate_integer_objective(rec[0].interval_widths,
                                              rec[0].control, rec[1].control);
    r.gap        = rec[0].integral_gap;
    r.n_switches = rec[0].n_switches;
    {   // how often the two controls actually differ: the value of independent switching
        int ndiff = 0;
        for (int i = 0; i < rec[0].control.size(); ++i)
            if (rec[0].control(i) != rec[1].control(i)) ++ndiff;
        r.n_split = ndiff;
    }
    r.status     = status;
    return r;
}

//////////////////////////////////////////////////////////////////////////
///////////////////  Main: mesh sweep /////////////////////////////////////
//////////////////////////////////////////////////////////////////////////

int main(void)
{
    const int meshes[] = {40, 80, 160, 320, 640};
    const int nmesh = (int)(sizeof(meshes)/sizeof(meshes[0]));

    printf("\n========== Lotka-Volterra fishing with two independent binary controls ==========\n");
    printf("  The prey and the predator are fished independently, w1 and w2 in {0,1}.\n");
    printf("  PSOPT convexifies over the product of the two admissible sets: four\n");
    printf("  weight-controls and one SOS1 constraint, added internally.\n");
    printf("  The single-control problem of examples/lotka is the case w1 == w2, so its\n");
    printf("  optimum Phi* = %.5f is an upper bound for this one.\n\n", PHI_STAR);
    printf("  %6s  %14s  %14s  %14s  %10s  %10s\n",
           "nodes", "relaxed Phi", "integer Phi", "gap to Phi*", "switches", "w1 != w2");
    printf("  --------------------------------------------------------------------------\n");

    for (int m = 0; m < nmesh; ++m) {
        SweepResult r = solve_and_round(meshes[m]);
        printf("  %6d  %14.6f  %14.6f  %14.3e  %10d  %10d\n",
               meshes[m], r.relaxed, r.integer, r.integer - PHI_STAR,
               r.n_switches, r.n_split);
    }

    printf("  --------------------------------------------------------------------------\n");
    printf("  Dynamics written once in terms of w1 and w2; the convexification over the\n");
    printf("  four product modes was performed internally by declare_integer_control, and\n");
    printf("  reconstruct_integer_controls rounded over those modes and decoded each one\n");
    printf("  back into a value for w1 and a value for w2. The last column counts the\n");
    printf("  intervals on which the two differ: those are the ones a single joint\n");
    printf("  control could not have produced, and they are why the integer objective\n");
    printf("  here can fall below the single-control optimum %.5f.\n\n", PHI_STAR);
    return 0;
}

//////////////////////////////////////////////////////////////////////////
///////////////////////      END OF FILE     //////////////////////////////
//////////////////////////////////////////////////////////////////////////
