//////////////////////////////////////////////////////////////////////////
////////////////           PSOPT  Example             ////////////////////
//////////////////////////////////////////////////////////////////////////
//////// Title:   Time-WINDOWED path constraint                       ////
//////// Illustrates: how to enforce a path constraint only during a  ////
////////              time window [ta,tb]. PSOPT applies path          ////
////////              constraints at EVERY node uniformly, so we gate  ////
////////              the constraint with a smooth time bump that is   ////
////////              ~1 inside [ta,tb] and ~0 outside, relaxing the   ////
////////              constraint away from the window.                 ////
//////// Problem:     double integrator move (pos 0->1, vel 0->0),     ////
////////              minimum control energy, with pos <= 0.35 ONLY    ////
////////              during t in [0.4, 0.6] (a transient keep-under). ////
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
    adouble pos = states[0];
    adouble vel = states[1];
    adouble u   = controls[0];
    d[0] = vel;
    d[1] = u;

    // Time-windowed path constraint:  vel <= VLIM only for t in [ta,tb].
    // gate(t) ~ 1 inside the window, ~0 outside (product of two sharp sigmoids).
    // Velocity peaks (~1.5) mid-move for the unconstrained solution, so a limit
    // active around the peak binds and reshapes the trajectory.
    const double ta = 0.35, tb = 0.65, VLIM = 1.30, s = 200.0, BIG = 50.0;
    adouble gate = (1.0/(1.0+exp(-s*(time-ta)))) * (1.0/(1.0+exp(s*(time-tb))));
    // Outside the window gate~0 => path ~ vel - VLIM - BIG <= 0 (always ok);
    // inside gate~1 => path ~ vel - VLIM <= 0 (enforced).
    path[0] = (vel - VLIM) - BIG*(1.0 - gate);
}

void events(adouble* e, adouble* i, adouble* f, adouble* p, adouble& t0,
            adouble& tf, adouble* xad, int iphase, Workspace* w)
{ e[0]=i[0]; e[1]=i[1]; e[2]=f[0]; e[3]=f[1]; }

void linkages(adouble* linkages, adouble* xad, Workspace* w) {}

int main(void)
{
    Alg algorithm; Sol solution; Prob problem;
    problem.name = "Time-windowed path constraint";
    problem.outfilename = "path_window.txt";
    problem.nphases = 1; problem.nlinkages = 0;
    psopt_level1_setup(problem);

    problem.phases(1).nstates   = 2;
    problem.phases(1).ncontrols = 1;
    problem.phases(1).nevents   = 4;
    problem.phases(1).npath     = 1;   // the gated path constraint
    problem.phases(1).nodes     << 60;
    psopt_level2_setup(problem, algorithm);

    problem.phases(1).bounds.lower.states(0) = -5.0; problem.phases(1).bounds.upper.states(0) = 5.0;
    problem.phases(1).bounds.lower.states(1) = -5.0; problem.phases(1).bounds.upper.states(1) = 5.0;
    problem.phases(1).bounds.lower.controls(0) = -20.0; problem.phases(1).bounds.upper.controls(0) = 20.0;

    // path <= 0  (lower bound -inf, upper 0)
    problem.phases(1).bounds.lower.path(0) = -1.0e19;
    problem.phases(1).bounds.upper.path(0) = 0.0;

    problem.phases(1).bounds.lower.events << 0.0, 0.0, 1.0, 0.0;
    problem.phases(1).bounds.upper.events << 0.0, 0.0, 1.0, 0.0;

    problem.phases(1).bounds.lower.StartTime = 0.0; problem.phases(1).bounds.upper.StartTime = 0.0;
    problem.phases(1).bounds.lower.EndTime   = 1.0; problem.phases(1).bounds.upper.EndTime   = 1.0;

    problem.integrand_cost = &integrand_cost;
    problem.endpoint_cost  = &endpoint_cost;
    problem.dae            = &dae;
    problem.events         = &events;
    problem.linkages       = &linkages;

    problem.phases(1).guess.states   = zeros(2,60);
    problem.phases(1).guess.states.row(0) = linspace(0.0,1.0,60);
    problem.phases(1).guess.controls = zeros(1,60);
    problem.phases(1).guess.time     = linspace(0.0,1.0,60);

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
