//////////////////////////////////////////////////////////////////////////
////////////////           PSOPT  Example             ////////////////////
//////////////////////////////////////////////////////////////////////////
//////// Title:   NATIVE integral (isoperimetric) constraint          ////
//////// Illustrates: the new problem.integral_constraints hook, which ////
////////              declares an integral constraint  I_L <= INT q dt ////
////////              <= I_U directly, WITHOUT the auxiliary-state      ////
////////              reformulation. PSOPT integrates the user integrand ///
////////              q with the collocation quadrature (same weights   ////
////////              as the running cost) and bounds the result.       ////
//////// Problem:     minimum-energy double integrator (pos 0->1,       ////
////////              vel 0->0) with an area constraint INT_0^1 pos dt  ////
////////              = 0.40. The unconstrained optimum has area 0.5,   ////
////////              so the constraint binds and reshapes the path.    ////
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

// NATIVE integral constraint: return the integrand(s) q at each node; PSOPT
// integrates them over the phase and applies the integral bounds.
void integral_constraints(adouble* q, adouble* states, adouble* controls,
                          adouble* parameters, adouble& time, int iphase,
                          Workspace* w)
{
    // Integral 0: area under the position curve, INT pos dt
    q[0] = states[0];
}

int main(void)
{
    Alg algorithm; Sol solution; Prob problem;
    problem.name = "Integral (isoperimetric) constrained double integrator";
    problem.outfilename = "integral_ctc.txt";
    problem.nphases = 1; problem.nlinkages = 0;
    psopt_level1_setup(problem);

    problem.phases(1).nstates    = 2;
    problem.phases(1).ncontrols  = 1;
    problem.phases(1).nevents    = 4;
    problem.phases(1).npath      = 0;
    problem.phases(1).nintegral  = 1;   // <-- one native integral constraint
    problem.phases(1).nodes     << 40;
    psopt_level2_setup(problem, algorithm);

    problem.phases(1).bounds.lower.states(0) = -5.0; problem.phases(1).bounds.upper.states(0) = 5.0;
    problem.phases(1).bounds.lower.states(1) = -5.0; problem.phases(1).bounds.upper.states(1) = 5.0;
    problem.phases(1).bounds.lower.controls(0) = -30.0; problem.phases(1).bounds.upper.controls(0) = 30.0;

    problem.phases(1).bounds.lower.events << 0.0, 0.0, 1.0, 0.0;
    problem.phases(1).bounds.upper.events << 0.0, 0.0, 1.0, 0.0;

    // Integral constraint:  INT_0^1 pos dt = 0.40  (equality)
    problem.phases(1).bounds.lower.integral(0) = 0.40;
    problem.phases(1).bounds.upper.integral(0) = 0.40;

    problem.phases(1).bounds.lower.StartTime = 0.0; problem.phases(1).bounds.upper.StartTime = 0.0;
    problem.phases(1).bounds.lower.EndTime   = 1.0; problem.phases(1).bounds.upper.EndTime   = 1.0;

    problem.integrand_cost = &integrand_cost;
    problem.endpoint_cost  = &endpoint_cost;
    problem.dae            = &dae;
    problem.events         = &events;
    problem.linkages       = &linkages;
    problem.integral_constraints = &integral_constraints;   // <-- register the hook

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
