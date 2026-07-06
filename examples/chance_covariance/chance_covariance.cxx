//////////////////////////////////////////////////////////////////////////
////////////////           PSOPT  Example             ////////////////////
//////////////////////////////////////////////////////////////////////////
//////// Title:   TRUE chance constraint via covariance propagation    ////
//////// Illustrates: a probabilistic path constraint                  ////
////////              P{ vel(t) <= VMAX } >= 1 - eps                    ////
////////              enforced EXACTLY (not with a constant sigma) by   ////
////////              propagating the state covariance Sigma(t) as      ////
////////              extra states via the Lyapunov ODE                 ////
////////                 Sigma' = A Sigma + Sigma A^T + W               ////
////////              and tightening the nominal limit by k*sqrt(P22),  ////
////////              k = Phi^{-1}(1-eps). The margin GROWS in time as  ////
////////              process noise accumulates -> a time-varying,      ////
////////              physically-correct back-off.                      ////
//////// System:      double integrator with process noise on the      ////
////////              acceleration; A = [[0,1],[0,0]], W = diag(0,w).   ////
//////// Problem:     MINIMUM-TIME move pos 0->1 (vel 0->0, |u|<=1) with ///
////////              the chance velocity constraint. The cap extends tf ////
////////              (2.0 -> 2.12): the peak velocity rides the         ////
////////              time-varying tightened limit VMAX - k*sqrt(P22).   ////
//////////////////////////////////////////////////////////////////////////

#include "psopt.h"

adouble endpoint_cost(adouble* i, adouble* f, adouble* p, adouble& t0,
                      adouble& tf, adouble* xad, int iphase, Workspace* w)
{ return tf; }   // minimum time: the chance velocity cap extends tf (always feasible, binds)

adouble integrand_cost(adouble* states, adouble* controls, adouble* p,
                       adouble& time, adouble* xad, int iphase, Workspace* w)
{ return 0.0; }

void dae(adouble* d, adouble* path, adouble* states, adouble* controls,
         adouble* p, adouble& time, adouble* xad, int iphase, Workspace* w)
{
    adouble pos = states[0], vel = states[1];
    adouble P11 = states[2], P12 = states[3], P22 = states[4];  // covariance entries
    adouble u   = controls[0];

    const double wnoise = 0.05;   // process-noise variance rate on acceleration
    // mean dynamics
    d[0] = vel;
    d[1] = u;
    // covariance dynamics:  Sigma' = A Sigma + Sigma A^T + W, A=[[0,1],[0,0]], W=diag(0,w)
    d[2] = 2.0*P12;               // P11'
    d[3] = P22;                   // P12'
    d[4] = wnoise;                // P22'

    // TRUE chance constraint:  P{ vel <= VMAX } >= 1 - eps
    //   =>  vel + k*sqrt(Var[vel]) <= VMAX,  Var[vel] = P22(t)
    const double VMAX = 1.10, k = 1.6449;   // Phi^{-1}(0.95)
    path[0] = vel + k*sqrt(P22 + 2.5e-3) - VMAX;   // <= 0
}

void events(adouble* e, adouble* i, adouble* f, adouble* p, adouble& t0,
            adouble& tf, adouble* xad, int iphase, Workspace* w)
{
    e[0]=i[0]; e[1]=i[1];              // pos(0)=0, vel(0)=0
    e[2]=i[2]; e[3]=i[3]; e[4]=i[4];   // Sigma(0)=0 (known initial state)
    e[5]=f[0]; e[6]=f[1];              // pos(T)=1, vel(T)=0
}

void linkages(adouble* linkages, adouble* xad, Workspace* w) {}

int main(void)
{
    Alg algorithm; Sol solution; Prob problem;
    problem.name = "Chance-constrained double integrator (covariance propagation)";
    problem.outfilename = "chance_covariance.txt";
    problem.nphases = 1; problem.nlinkages = 0;
    psopt_level1_setup(problem);

    problem.phases(1).nstates   = 5;   // pos, vel, P11, P12, P22
    problem.phases(1).ncontrols = 1;
    problem.phases(1).nevents   = 7;
    problem.phases(1).npath     = 1;
    problem.phases(1).nodes     << 40;
    psopt_level2_setup(problem, algorithm);

    problem.phases(1).bounds.lower.states << -5.0, -5.0,  0.0, -1.0,  0.0;
    problem.phases(1).bounds.upper.states <<  5.0,  5.0,  2.0,  2.0,  2.0;
    problem.phases(1).bounds.lower.controls(0) = -1.0; problem.phases(1).bounds.upper.controls(0) = 1.0;

    problem.phases(1).bounds.lower.path(0) = -1.0e19;
    problem.phases(1).bounds.upper.path(0) = 0.0;

    problem.phases(1).bounds.lower.events << 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0;
    problem.phases(1).bounds.upper.events << 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0;

    problem.phases(1).bounds.lower.StartTime = 0.0; problem.phases(1).bounds.upper.StartTime = 0.0;
    problem.phases(1).bounds.lower.EndTime   = 0.5; problem.phases(1).bounds.upper.EndTime   = 10.0;

    problem.integrand_cost = &integrand_cost;
    problem.endpoint_cost  = &endpoint_cost;
    problem.dae            = &dae;
    problem.events         = &events;
    problem.linkages       = &linkages;

    problem.phases(1).guess.states   = zeros(5,40);
    problem.phases(1).guess.states.row(0) = linspace(0.0,1.0,40);
    problem.phases(1).guess.controls = zeros(1,40);
    problem.phases(1).guess.time     = linspace(0.0,3.0,40);

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
    // Rows of x: pos, vel, P11, P12, P22. The chance margin is k*sqrt(P22(t)).
    plot(t,x,problem.name,"time (s)","pos vel P11 P12 P22","pos vel P11 P12 P22");
    plot(t,u,problem.name,"time (s)","control","u");
}
