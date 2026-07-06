//////////////////////////////////////////////////////////////////////////
////////////////           PSOPT  Example             ////////////////////
//////////////////////////////////////////////////////////////////////////
//////// Title:   Minimum-energy double integrator (simple constraints) ////
//////// Reference: A. Locatelli, "Optimal Control of a Double          ////
////////            Integrator", Springer 2017, Chaps. 3-5.            ////
//////// Problem:  xdot1=x2, xdot2=u; (x1,x2):(0,0)->(1,0), tf=1 fixed; ////
////////           minimise  J = INT_0^1 u^2/2 dt.                     ////
//////// Maximum Principle:  H = u^2/2 + l1 x2 + l2 u, dH/du=0 => u=-l2,////
////////           l1=const, l2 linear => u is LINEAR in t. Closed form:////
////////           u(t) = 6 - 12 t,  x1 = 3t^2-2t^3,  x2 = 6t-6t^2,     ////
////////           J* = 6,  peak speed x2(1/2) = 3/2.                   ////
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
{ e[0]=i[0]; e[1]=i[1]; e[2]=f[0]; e[3]=f[1]; }          // (0,0) -> (1,0)

void linkages(adouble* l, adouble* xad, Workspace* w) {}

int main(void)
{
    Alg algorithm; Sol solution; Prob problem;
    problem.name = "Minimum-energy double integrator (Locatelli Ch 3)";
    problem.outfilename = "mineng_di.txt";
    problem.nphases = 1; problem.nlinkages = 0;
    psopt_level1_setup(problem);

    const int N = 40;
    problem.phases(1).nstates   = 2;
    problem.phases(1).ncontrols = 1;
    problem.phases(1).nevents   = 4;
    problem.phases(1).npath     = 0;
    problem.phases(1).nodes     << N;
    psopt_level2_setup(problem, algorithm);

    problem.phases(1).bounds.lower.states << -5.0, -5.0;  problem.phases(1).bounds.upper.states << 5.0, 5.0;
    problem.phases(1).bounds.lower.controls(0) = -50.0;   problem.phases(1).bounds.upper.controls(0) = 50.0;
    problem.phases(1).bounds.lower.events << 0.0, 0.0, 1.0, 0.0;
    problem.phases(1).bounds.upper.events << 0.0, 0.0, 1.0, 0.0;
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
    plot(t,x,problem.name,"time (s)","x1 x2","x1 x2");
    plot(t,u,problem.name,"time (s)","u = 6 - 12 t","u");
}
