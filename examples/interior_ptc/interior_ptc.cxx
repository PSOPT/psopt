//////////////////////////////////////////////////////////////////////////
////////////////           PSOPT  Example             ////////////////////
//////////////////////////////////////////////////////////////////////////
//////// Title:   NATIVE interior-point constraint (single phase)      ////
//////// Illustrates: the new problem.interior_point_constraints hook, ////
////////              which enforces a constraint on the state/control ////
////////              interpolated to a specified interior time WITHOUT ///
////////              splitting the phase. Uses the same Lagrange       ////
////////              endpoint-interpolation machinery as Radau/Gauss.  ////
//////// Problem:     minimum-energy double integrator (pos 0->1,       ////
////////              vel 0->0) with a WAYPOINT: pos = 0.5 at 30% of    ////
////////              the (normalized) phase. The waypoint forces early ////
////////              progress and reshapes the min-energy trajectory.  ////
//////////////////////////////////////////////////////////////////////////

#include "psopt.h"

adouble endpoint_cost(adouble* i, adouble* f, adouble* p, adouble& t0,
                      adouble& tf, adouble* xad, int iphase, Workspace* w)
{ return 0.0; }

adouble integrand_cost(adouble* states, adouble* controls, adouble* p,
                       adouble& time, adouble* xad, int iphase, Workspace* w)
{ adouble u = controls[0]; return u*u; }

void dae(adouble* d, adouble* path, adouble* states, adouble* controls,
         adouble* p, adouble& time, adouble* xad, int iphase, Workspace* w)
{
    adouble pos = states[0], vel = states[1], u = controls[0];
    d[0] = vel;
    d[1] = u;
}

void events(adouble* e, adouble* i, adouble* f, adouble* p, adouble& t0,
            adouble& tf, adouble* xad, int iphase, Workspace* w)
{ e[0]=i[0]; e[1]=i[1]; e[2]=f[0]; e[3]=f[1]; }

void linkages(adouble* linkages, adouble* xad, Workspace* w) {}

// NATIVE interior-point constraint: called with the state (xi) and control (ui)
// interpolated at the interior time; index selects which interior point.
void interior_point_constraints(adouble* g, adouble* xi, adouble* ui,
                                adouble* parameters, adouble& time, int index,
                                int iphase, Workspace* w)
{
    // Interior point 0: position waypoint  pos(t_interior) = 0.5
    g[0] = xi[0] - 0.5;
}

int main(void)
{
    Alg algorithm; Sol solution; Prob problem;
    problem.name = "Interior-point (waypoint) constrained double integrator";
    problem.outfilename = "interior_ptc.txt";
    problem.nphases = 1; problem.nlinkages = 0;
    psopt_level1_setup(problem);

    problem.phases(1).nstates    = 2;
    problem.phases(1).ncontrols  = 1;
    problem.phases(1).nevents    = 4;
    problem.phases(1).npath      = 0;
    problem.phases(1).ninterior  = 1;   // <-- one native interior-point constraint
    problem.phases(1).nodes     << 40;
    psopt_level2_setup(problem, algorithm);

    problem.phases(1).bounds.lower.states(0) = -5.0; problem.phases(1).bounds.upper.states(0) = 5.0;
    problem.phases(1).bounds.lower.states(1) = -5.0; problem.phases(1).bounds.upper.states(1) = 5.0;
    problem.phases(1).bounds.lower.controls(0) = -30.0; problem.phases(1).bounds.upper.controls(0) = 30.0;

    problem.phases(1).bounds.lower.events << 0.0, 0.0, 1.0, 0.0;
    problem.phases(1).bounds.upper.events << 0.0, 0.0, 1.0, 0.0;

    // Interior-point constraint: g = 0 (equality) at 30% of the phase.
    problem.phases(1).interior_time(0)        = 0.30;   // normalized [0,1]
    problem.phases(1).bounds.lower.interior(0) = 0.0;
    problem.phases(1).bounds.upper.interior(0) = 0.0;

    problem.phases(1).bounds.lower.StartTime = 0.0; problem.phases(1).bounds.upper.StartTime = 0.0;
    problem.phases(1).bounds.lower.EndTime   = 1.0; problem.phases(1).bounds.upper.EndTime   = 1.0;

    problem.integrand_cost = &integrand_cost;
    problem.endpoint_cost  = &endpoint_cost;
    problem.dae            = &dae;
    problem.events         = &events;
    problem.linkages       = &linkages;
    problem.interior_point_constraints = &interior_point_constraints;  // <-- register the hook

    problem.phases(1).guess.states   = zeros(2,40);
    problem.phases(1).guess.states.row(0) = linspace(0.0,1.0,40);
    problem.phases(1).guess.controls = zeros(1,40);
    problem.phases(1).guess.time     = linspace(0.0,1.0,40);

    algorithm.nlp_method         = "IPOPT";
    algorithm.scaling            = "automatic";
    algorithm.derivatives        = "automatic";
    algorithm.collocation_method = "Legendre";
    algorithm.nlp_iter_max       = 1000;
    algorithm.nlp_tolerance      = 1.e-6;

    psopt(solution, problem, algorithm);

    DMatrix x = solution.get_states_in_phase(1);
    DMatrix u = solution.get_controls_in_phase(1);
    DMatrix t = solution.get_time_in_phase(1);
    Save(x,"x.dat"); Save(u,"u.dat"); Save(t,"t.dat");
    plot(t,x,problem.name,"time (s)","pos vel","pos vel");
    plot(t,u,problem.name,"time (s)","control","u");
}
