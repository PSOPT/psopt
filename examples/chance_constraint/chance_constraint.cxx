//////////////////////////////////////////////////////////////////////////
////////////////           PSOPT  Example             ////////////////////
//////////////////////////////////////////////////////////////////////////
//////// Title:   CHANCE constraint via deterministic tightening      ////
//////// Illustrates: a probabilistic path constraint                 ////
////////              P{ pos(t) <= POSMAX } >= 1 - eps                 ////
////////              enforced deterministically by tightening the     ////
////////              nominal constraint by k*sigma, where             ////
////////              k = Phi^{-1}(1-eps) and sigma is the (here       ////
////////              constant, illustrative) state std-deviation.     ////
////////              Full version: propagate sigma(t) via a Lyapunov  ////
////////              ODE as extra states (see docs/CONSTRAINTS.md).   ////
//////// Problem:     double integrator move (pos 0->1, vel 0->0),     ////
////////              minimum energy, with the chance constraint on    ////
////////              the position overshoot.                          ////
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

    // Chance constraint  P{ vel <= VMAX } >= 1 - eps  becomes, for a Gaussian
    // velocity estimate,  vel + k*sigma <= VMAX: the nominal velocity limit is
    // TIGHTENED by the safety margin k*sigma. (Velocity peaks mid-move at 1.5
    // for the unconstrained solution, so this constraint binds.)
    const double VMAX  = 1.30;    // nominal velocity limit
    const double sigma = 0.03;    // (illustrative) velocity std-dev
    const double k     = 1.6449;  // Phi^{-1}(0.95) -> eps = 5%
    path[0] = vel + k*sigma - VMAX;   // <= 0  => vel <= 1.251
}

void events(adouble* e, adouble* i, adouble* f, adouble* p, adouble& t0,
            adouble& tf, adouble* xad, int iphase, Workspace* w)
{ e[0]=i[0]; e[1]=i[1]; e[2]=f[0]; e[3]=f[1]; }

void linkages(adouble* linkages, adouble* xad, Workspace* w) {}

int main(void)
{
    Alg algorithm; Sol solution; Prob problem;
    problem.name = "Chance-constrained double integrator";
    problem.outfilename = "chance_constraint.txt";
    problem.nphases = 1; problem.nlinkages = 0;
    psopt_level1_setup(problem);

    problem.phases(1).nstates   = 2;
    problem.phases(1).ncontrols = 1;
    problem.phases(1).nevents   = 4;
    problem.phases(1).npath     = 1;   // the tightened (chance) constraint
    problem.phases(1).nodes     << 40;
    psopt_level2_setup(problem, algorithm);

    problem.phases(1).bounds.lower.states(0) = -5.0; problem.phases(1).bounds.upper.states(0) = 5.0;
    problem.phases(1).bounds.lower.states(1) = -5.0; problem.phases(1).bounds.upper.states(1) = 5.0;
    problem.phases(1).bounds.lower.controls(0) = -20.0; problem.phases(1).bounds.upper.controls(0) = 20.0;

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
