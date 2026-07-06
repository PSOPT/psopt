//////////////////////////////////////////////////////////////////////////
////////////////           PSOPT  Example             ////////////////////
//////////////////////////////////////////////////////////////////////////
//////// Title:   SINGULAR-ARC double integrator                      ////
//////// Reference: A. Locatelli, "Optimal Control of a Double         ////
////////            Integrator", Springer 2017, Problem 11.2-1a.       ////
//////// Problem:  xdot1=x2, xdot2=u, |u|<=1; x(0)=(2,0), x(tf)=(5,0), ////
////////           tf FREE; minimise  J = INT_0^tf (1 + x2^2)/2 dt.    ////
//////// Maximum Principle:  H = (1+x2^2)/2 + l1 x2 + l2 u is LINEAR    ////
////////           in u, so u = -sign(l2) is bang-bang UNLESS the       ///
////////           switching function l2(t) vanishes on an interval --  ///
////////           a SINGULAR arc, where u is set by higher-order       ////
////////           conditions (d^2/dt^2 (dH/du)=0 => u=0, holding x2=1).///
//////// Analytical optimum:  u = +1 on [0,1], u = 0 (singular, x2=1)   ////
////////           on [1,3], u = -1 on [3,4];  tf* = 4,  J* = 10/3.     ////
//////////////////////////////////////////////////////////////////////////

#include "psopt.h"

adouble endpoint_cost(adouble* i, adouble* f, adouble* p, adouble& t0,
                      adouble& tf, adouble* xad, int iphase, Workspace* w)
{ return 0.0; }

adouble integrand_cost(adouble* states, adouble* controls, adouble* p,
                       adouble& time, adouble* xad, int iphase, Workspace* w)
{ adouble x2 = states[1]; return 0.5*(1.0 + x2*x2); }   // (1 + x2^2)/2

void dae(adouble* d, adouble* path, adouble* states, adouble* controls,
         adouble* p, adouble& time, adouble* xad, int iphase, Workspace* w)
{ d[0] = states[1];  d[1] = controls[0]; }               // double integrator

void events(adouble* e, adouble* i, adouble* f, adouble* p, adouble& t0,
            adouble& tf, adouble* xad, int iphase, Workspace* w)
{ e[0]=i[0]; e[1]=i[1]; e[2]=f[0]; e[3]=f[1]; }          // (2,0) -> (5,0)

void linkages(adouble* linkages, adouble* xad, Workspace* w) {}

int main(void)
{
    Alg algorithm; Sol solution; Prob problem;
    problem.name = "Singular-arc double integrator (Locatelli 11.2-1a)";
    problem.outfilename = "singular_di.txt";
    problem.nphases = 1; problem.nlinkages = 0;
    psopt_level1_setup(problem);

    const int N = 60;
    problem.phases(1).nstates   = 2;
    problem.phases(1).ncontrols = 1;
    problem.phases(1).nevents   = 4;
    problem.phases(1).npath     = 0;
    problem.phases(1).nodes     << N;
    psopt_level2_setup(problem, algorithm);

    problem.phases(1).bounds.lower.states << -10.0, -5.0;
    problem.phases(1).bounds.upper.states <<  10.0,  5.0;
    problem.phases(1).bounds.lower.controls(0) = -1.0;   // |u| <= 1
    problem.phases(1).bounds.upper.controls(0) =  1.0;
    problem.phases(1).bounds.lower.events << 2.0, 0.0, 5.0, 0.0;
    problem.phases(1).bounds.upper.events << 2.0, 0.0, 5.0, 0.0;
    problem.phases(1).bounds.lower.StartTime = 0.0; problem.phases(1).bounds.upper.StartTime = 0.0;
    problem.phases(1).bounds.lower.EndTime   = 2.0; problem.phases(1).bounds.upper.EndTime   = 8.0;

    problem.integrand_cost = &integrand_cost;
    problem.endpoint_cost  = &endpoint_cost;
    problem.dae            = &dae;
    problem.events         = &events;
    problem.linkages       = &linkages;

    // Analytical bang-singular-bang warm start (tf=4: +1 on [0,1], 0 on [1,3], -1 on [3,4]).
    MatrixXd xg(2,N), ug(1,N), tg(1,N);
    for (int k = 0; k < N; ++k) {
        double t = 4.0 * k / (N - 1);
        tg(0,k) = t;
        if (t < 1.0)      { ug(0,k)= 1.0; xg(0,k)=2.0+0.5*t*t;         xg(1,k)=t; }
        else if (t < 3.0) { ug(0,k)= 0.0; xg(0,k)=2.5+(t-1.0);          xg(1,k)=1.0; }
        else              { double s=t-3.0; ug(0,k)=-1.0; xg(0,k)=4.5+s-0.5*s*s; xg(1,k)=1.0-s; }
    }
    problem.phases(1).guess.states   = xg;
    problem.phases(1).guess.controls = ug;
    problem.phases(1).guess.time     = tg;

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
    // The middle arc holds x2 == 1 (the singular coast); J* = 10/3, tf* = 4.
    plot(t,x,problem.name,"time (s)","x1 x2","x1 x2");
    plot(t,u,problem.name,"time (s)","u (bang-singular-bang)","u");
}
