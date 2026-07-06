//////////////////////////////////////////////////////////////////////////
////////////////           PSOPT  Example             ////////////////////
//////////////////////////////////////////////////////////////////////////
//////// Title:   IN-CORE mixed-integer optimal control via SCIP       ////
//////// Illustrates: algorithm.nlp_method = "SCIP" with an INTEGER    ////
////////              control declared through                         ////
////////              problem.phase[i].integer_controls. PSOPT builds  ////
////////              the (linear) transcription and SCIP solves it to  ///
////////              GLOBAL optimality by branch-and-bound.           ////
//////// Requires:    build with -DPSOPT_WITH_SCIP=ON.                 ////
//////// Problem:     double integrator, rest-to-rest (pos 0->1, vel   ////
////////              0->0) over T=3, with INTEGER thrust u in {-1,0,1} ////
////////              minimizing total actuation INT|u| dt.            ////
//////// Formulation (kept LINEAR so SCIP's model is exact): an        ////
////////              auxiliary continuous control a >= |u| carries    ////
////////              the objective INT a dt, enforced by two linear   ////
////////              path constraints a-u>=0, a+u>=0.                 ////
//////////////////////////////////////////////////////////////////////////

#include "psopt.h"

adouble endpoint_cost(adouble* i, adouble* f, adouble* p, adouble& t0,
                      adouble& tf, adouble* xad, int iphase, Workspace* w)
{ return 0.0; }

adouble integrand_cost(adouble* states, adouble* controls, adouble* p,
                       adouble& time, adouble* xad, int iphase, Workspace* w)
{ return controls[1]; }   // a  (>= |u|)  -> minimise INT a dt  (LINEAR)

void dae(adouble* d, adouble* path, adouble* states, adouble* controls,
         adouble* p, adouble& time, adouble* xad, int iphase, Workspace* w)
{
    adouble vel = states[1];
    adouble u   = controls[0];   // integer thrust
    adouble a   = controls[1];   // continuous, >= |u|
    d[0] = vel;                  // pos' = vel     (linear)
    d[1] = u - 0.03*vel*vel;      // NONLINEAR quadratic drag -> SLP+SCIP
    path[0] = a - u;             // >= 0  => a >= u
    path[1] = a + u;             // >= 0  => a >= -u   (so a >= |u|)
}

void events(adouble* e, adouble* i, adouble* f, adouble* p, adouble& t0,
            adouble& tf, adouble* xad, int iphase, Workspace* w)
{ e[0]=i[0]; e[1]=i[1]; e[2]=f[0]; e[3]=f[1]; }

void linkages(adouble* linkages, adouble* xad, Workspace* w) {}

int main(void)
{
    Alg algorithm; Sol solution; Prob problem;
    problem.name = "Nonlinear MIOC via SCIP (SLP heuristic)";
    problem.outfilename = "scip_nlmiop.txt";
    problem.nphases = 1; problem.nlinkages = 0;
    psopt_level1_setup(problem);

    problem.phases(1).nstates   = 2;   // pos, vel
    problem.phases(1).ncontrols = 2;   // u (integer), a (continuous)
    problem.phases(1).nevents   = 4;
    problem.phases(1).npath     = 2;   // a-u>=0, a+u>=0
    problem.phases(1).nodes     << 16;
    psopt_level2_setup(problem, algorithm);

    problem.phases(1).bounds.lower.states(0) = -10.0; problem.phases(1).bounds.upper.states(0) = 10.0;
    problem.phases(1).bounds.lower.states(1) = -10.0; problem.phases(1).bounds.upper.states(1) = 10.0;
    problem.phases(1).bounds.lower.controls(0) = -1.0; problem.phases(1).bounds.upper.controls(0) = 1.0; // u in {-1,0,1}
    problem.phases(1).bounds.lower.controls(1) =  0.0; problem.phases(1).bounds.upper.controls(1) = 1.0; // a in [0,1]

    problem.phases(1).bounds.lower.path << 0.0, 0.0;      // a-u >= 0, a+u >= 0
    problem.phases(1).bounds.upper.path << 1.0e19, 1.0e19;

    problem.phases(1).bounds.lower.events << 0.0, 0.0, 0.90, -0.15;
    problem.phases(1).bounds.upper.events << 0.0, 0.0, 1.10,  0.15;

    problem.phases(1).bounds.lower.StartTime = 0.0; problem.phases(1).bounds.upper.StartTime = 0.0;
    problem.phases(1).bounds.lower.EndTime   = 3.0; problem.phases(1).bounds.upper.EndTime   = 3.0;

    // Declare control 0 (u) INTEGER — solved by SCIP over {-1,0,1} at every node.
    problem.phases(1).integer_controls.resize(1,1);
    problem.phases(1).integer_controls(0) = 0;

    problem.integrand_cost = &integrand_cost;
    problem.endpoint_cost  = &endpoint_cost;
    problem.dae            = &dae;
    problem.events         = &events;
    problem.linkages       = &linkages;

    problem.phases(1).guess.states   = zeros(2,16);
    problem.phases(1).guess.states.row(0) = linspace(0.0,1.0,16);
    problem.phases(1).guess.controls = zeros(2,16);
    problem.phases(1).guess.time     = linspace(0.0,3.0,16);

    algorithm.nlp_method         = "SCIP";        // <-- mixed-integer backend
    algorithm.collocation_method = "trapezoidal"; // local, linear defects
    algorithm.mesh_refinement    = "manual";
    algorithm.scaling            = "automatic";

    psopt(solution, problem, algorithm);

    DMatrix x = solution.get_states_in_phase(1);
    DMatrix u = solution.get_controls_in_phase(1);
    DMatrix t = solution.get_time_in_phase(1);
    Save(x,"x.dat"); Save(u,"u.dat"); Save(t,"t.dat");
    plot(t,x,problem.name,"time (s)","pos vel","pos vel");
    plot(t,u,problem.name,"time (s)","u (integer) / a","u a");
}
