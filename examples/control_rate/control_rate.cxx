//////////////////////////////////////////////////////////////////////////
////////////////           PSOPT  Example             ////////////////////
//////////////////////////////////////////////////////////////////////////
//////// Title:   Control-RATE constraint  |du/dt| <= R                ////
//////// Illustrates: how to impose a rate limit on a control when     ////
////////              PSOPT bounds controls but not their derivatives. ////
//////// Technique:   promote the control u to a STATE and introduce   ////
////////              its rate v = du/dt as the new control; then the  ////
////////              rate limit is a simple bound on the new control. ////
//////// Problem:     minimum-energy rest-to-rest move of a double     ////
////////              integrator (pos 0->1, vel 0->0) with |du/dt|<=R. ////
//////////////////////////////////////////////////////////////////////////

#include "psopt.h"

adouble endpoint_cost(adouble* i, adouble* f, adouble* p, adouble& t0,
                      adouble& tf, adouble* xad, int iphase, Workspace* w)
{ return 0.0; }

adouble integrand_cost(adouble* states, adouble* controls, adouble* p,
                       adouble& time, adouble* xad, int iphase, Workspace* w)
{
    // u is now a state (states[2]); minimise control energy int u^2 dt
    adouble u = states[2];
    return u*u;
}

void dae(adouble* d, adouble* path, adouble* states, adouble* controls,
         adouble* p, adouble& time, adouble* xad, int iphase, Workspace* w)
{
    adouble pos = states[0];
    adouble vel = states[1];
    adouble u   = states[2];   // control promoted to a state
    adouble v   = controls[0]; // v = du/dt is the new (rate) control

    d[0] = vel;                // pos' = vel
    d[1] = u;                  // vel' = u
    d[2] = v;                  // u'   = v   (so v is the control rate)
}

void events(adouble* e, adouble* i, adouble* f, adouble* p, adouble& t0,
            adouble& tf, adouble* xad, int iphase, Workspace* w)
{
    e[0] = i[0];  // pos(0)=0
    e[1] = i[1];  // vel(0)=0
    e[2] = i[2];  // u(0)=0
    e[3] = f[0];  // pos(tf)=1
    e[4] = f[1];  // vel(tf)=0
}

void linkages(adouble* linkages, adouble* xad, Workspace* w) {}

int main(void)
{
    Alg algorithm; Sol solution; Prob problem;
    problem.name = "Control-rate constrained double integrator";
    problem.outfilename = "control_rate.txt";
    problem.nphases = 1; problem.nlinkages = 0;
    psopt_level1_setup(problem);

    problem.phases(1).nstates   = 3;   // pos, vel, u
    problem.phases(1).ncontrols = 1;   // v = du/dt
    problem.phases(1).nevents   = 5;
    problem.phases(1).npath     = 0;
    problem.phases(1).nodes     << 40;
    psopt_level2_setup(problem, algorithm);

    const double R = 4.0;   // control-rate limit: |du/dt| <= R

    problem.phases(1).bounds.lower.states(0) = -5.0;   problem.phases(1).bounds.upper.states(0) = 5.0;
    problem.phases(1).bounds.lower.states(1) = -5.0;   problem.phases(1).bounds.upper.states(1) = 5.0;
    problem.phases(1).bounds.lower.states(2) = -10.0;  problem.phases(1).bounds.upper.states(2) = 10.0;
    problem.phases(1).bounds.lower.controls(0) = -R;   problem.phases(1).bounds.upper.controls(0) = R;  // the rate limit

    problem.phases(1).bounds.lower.events << 0.0, 0.0, 0.0, 1.0, 0.0;
    problem.phases(1).bounds.upper.events << 0.0, 0.0, 0.0, 1.0, 0.0;

    problem.phases(1).bounds.lower.StartTime = 0.0;  problem.phases(1).bounds.upper.StartTime = 0.0;
    problem.phases(1).bounds.lower.EndTime   = 2.0;  problem.phases(1).bounds.upper.EndTime   = 2.0;

    problem.integrand_cost = &integrand_cost;
    problem.endpoint_cost  = &endpoint_cost;
    problem.dae            = &dae;
    problem.events         = &events;
    problem.linkages       = &linkages;

    problem.phases(1).guess.states   = zeros(3,40);
    problem.phases(1).guess.states.row(0) = linspace(0.0,1.0,40);
    problem.phases(1).guess.controls = zeros(1,40);
    problem.phases(1).guess.time     = linspace(0.0,2.0,40);

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
    // Row 2 of x is the physical control u; u.dat is its rate du/dt (bounded by R).
    plot(t,x,problem.name,"time (s)","pos vel u","pos vel u");
    plot(t,u,problem.name,"time (s)","du/dt (rate)","v");
}
