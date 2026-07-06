//////////////////////////////////////////////////////////////////////////
////////////////           PSOPT  Example             ////////////////////
//////////////////////////////////////////////////////////////////////////
//////// Title:   Min-energy to a TARGET SET (nonstandard final state)  ////
//////// Reference: A. Locatelli, "Optimal Control of a Double          ////
////////            Integrator", Springer 2017, Ch 6 (final state on a  ////
////////            variety, not a point).                              ////
//////// Problem:  xdot1=x2, xdot2=u; x(0)=(0,0), tf=1 fixed; the final ////
////////           state must lie on the LINE  x1(tf) - 2 x2(tf) = 1;   ////
////////           minimise J = INT_0^1 u^2/2 dt.                       ////
//////// Maximum Principle:  the ORTHOGONALITY (transversality) condition////
////////           makes l(tf) parallel to grad(x1-2x2) = (1,-2), i.e.  ////
////////           the costate is orthogonal to the target line. Closed  ///
////////           form:  u(t) = -(3/7)(1+t),  J* = 3/14,               ////
////////           x1(1) = -2/7,  x2(1) = -9/14  (on the line).         ////
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
{ e[0]=i[0]; e[1]=i[1]; e[2] = f[0] - 2.0*f[1]; }        // (0,0) start; x1(tf)-2 x2(tf)=1

void linkages(adouble* l, adouble* xad, Workspace* w) {}

int main(void)
{
    Alg algorithm; Sol solution; Prob problem;
    problem.name = "Min-energy to a target line (Locatelli Ch 6)";
    problem.outfilename = "mineng_target_set.txt";
    problem.nphases = 1; problem.nlinkages = 0;
    psopt_level1_setup(problem);

    const int N = 40;
    problem.phases(1).nstates   = 2;
    problem.phases(1).ncontrols = 1;
    problem.phases(1).nevents   = 3;
    problem.phases(1).npath     = 0;
    problem.phases(1).nodes     << N;
    psopt_level2_setup(problem, algorithm);

    problem.phases(1).bounds.lower.states << -5.0, -5.0;  problem.phases(1).bounds.upper.states << 5.0, 5.0;
    problem.phases(1).bounds.lower.controls(0) = -50.0;   problem.phases(1).bounds.upper.controls(0) = 50.0;
    problem.phases(1).bounds.lower.events << 0.0, 0.0, 1.0;   // x1(0)=0, x2(0)=0, line=1
    problem.phases(1).bounds.upper.events << 0.0, 0.0, 1.0;
    problem.phases(1).bounds.lower.StartTime = 0.0; problem.phases(1).bounds.upper.StartTime = 0.0;
    problem.phases(1).bounds.lower.EndTime   = 1.0; problem.phases(1).bounds.upper.EndTime   = 1.0;

    problem.integrand_cost = &integrand_cost;
    problem.endpoint_cost  = &endpoint_cost;
    problem.dae            = &dae;
    problem.events         = &events;
    problem.linkages       = &linkages;

    problem.phases(1).guess.states   = zeros(2,N);
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
    // Final state lands on the line x1-2 x2 = 1 (costate orthogonal to it).
    plot(t,x,problem.name,"time (s)","x1 x2","x1 x2");
    plot(t,u,problem.name,"time (s)","u = -(3/7)(1+t)","u");
}
