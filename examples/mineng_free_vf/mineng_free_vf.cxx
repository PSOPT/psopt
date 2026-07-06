//////////////////////////////////////////////////////////////////////////
////////////////           PSOPT  Example             ////////////////////
//////////////////////////////////////////////////////////////////////////
//////// Title:   Min-energy, FREE final velocity (transversality)      ////
//////// Reference: A. Locatelli, "Optimal Control of a Double          ////
////////            Integrator", Springer 2017, Sect. 3.3 (x(tf) not    ////
////////            fully given).                                       ////
//////// Problem:  xdot1=x2, xdot2=u; x(0)=(0,0), x1(tf)=1, x2(tf) FREE,////
////////           tf=1 fixed; minimise J = INT_0^1 u^2/2 dt.          ////
//////// Maximum Principle:  x2(tf) free => TRANSVERSALITY l2(tf)=0, so ////
////////           the H-optimal control u=-l2 satisfies u(tf)=0.       ////
////////           Closed form:  u(t) = 3 - 3 t,  J* = 3/2,             ////
////////           x1(1)=1,  x2(1) = 3/2  (ends moving),  u(1) = 0.     ////
//////////////////////////////////////////////////////////////////////////

#include "psopt.h"

adouble endpoint_cost(adouble* i, adouble* f, adouble* p, adouble& t0,
                      adouble& tf, adouble* xad, int iphase, Workspace* w)
{ return 0.0; }

adouble integrand_cost(adouble* states, adouble* controls, adouble* p,
                       adouble& time, adouble* xad, int iphase, Workspace* w)
{ adouble u = controls[0]; return 0.5*u*u; }

void dae(adouble* d, adouble* path, adouble* states, adouble* controls,
         adouble* p, adouble& time, adouble* xad, int iphase, Workspace* w)
{ d[0] = states[1];  d[1] = controls[0]; }

void events(adouble* e, adouble* i, adouble* f, adouble* p, adouble& t0,
            adouble& tf, adouble* xad, int iphase, Workspace* w)
{ e[0]=i[0]; e[1]=i[1]; e[2]=f[0]; }                     // x(0)=(0,0), x1(tf)=1 (x2 free)

void linkages(adouble* l, adouble* xad, Workspace* w) {}

int main(void)
{
    Alg algorithm; Sol solution; Prob problem;
    problem.name = "Min-energy, free final velocity (Locatelli 3.3)";
    problem.outfilename = "mineng_free_vf.txt";
    problem.nphases = 1; problem.nlinkages = 0;
    psopt_level1_setup(problem);

    const int N = 40;
    problem.phases(1).nstates   = 2;
    problem.phases(1).ncontrols = 1;
    problem.phases(1).nevents   = 3;                     // one fewer: x2(tf) free
    problem.phases(1).npath     = 0;
    problem.phases(1).nodes     << N;
    psopt_level2_setup(problem, algorithm);

    problem.phases(1).bounds.lower.states << -5.0, -5.0;  problem.phases(1).bounds.upper.states << 5.0, 5.0;
    problem.phases(1).bounds.lower.controls(0) = -50.0;   problem.phases(1).bounds.upper.controls(0) = 50.0;
    problem.phases(1).bounds.lower.events << 0.0, 0.0, 1.0;
    problem.phases(1).bounds.upper.events << 0.0, 0.0, 1.0;
    problem.phases(1).bounds.lower.StartTime = 0.0; problem.phases(1).bounds.upper.StartTime = 0.0;
    problem.phases(1).bounds.lower.EndTime   = 1.0; problem.phases(1).bounds.upper.EndTime   = 1.0;

    problem.integrand_cost = &integrand_cost;
    problem.endpoint_cost  = &endpoint_cost;
    problem.dae            = &dae;
    problem.events         = &events;
    problem.linkages       = &linkages;

    problem.phases(1).guess.states   = zeros(2,N);
    problem.phases(1).guess.states.row(0) = linspace(0.0,1.0,N);
    problem.phases(1).guess.controls = zeros(1,N);
    problem.phases(1).guess.time     = linspace(0.0,1.0,N);

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
    // Transversality signature: u(tf) = 0 (l2(tf)=0), final velocity x2(1)=3/2.
    plot(t,x,problem.name,"time (s)","x1 x2","x1 x2");
    plot(t,u,problem.name,"time (s)","u = 3 - 3 t  (u(tf)=0)","u");
}
