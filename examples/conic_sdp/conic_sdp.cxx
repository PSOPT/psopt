//////////////////////////////////////////////////////////////////////////
////////////////           PSOPT  Example             ////////////////////
//////////////////////////////////////////////////////////////////////////
//////// Title:   SEMIDEFINITE (SDP) matrix constraint in optimal ctrl  ////
//////// Illustrates: a linear-matrix-inequality  M(x) >= 0  (positive  ////
////////              semidefinite) enforced through its leading        ////
////////              principal minors (Sylvester's criterion) as       ////
////////              smooth scalar path constraints.                   ////
////////              For M = [[1, b],[b, c]]:  M >= 0  <=>  c - b^2 >=0 ////
//////// Problem:     drive the off-diagonal b from 0 to 0.8 while       ////
////////              MINIMISING the diagonal c(T) at final time, subject///
////////              to M >= 0. The minor constraint c >= b^2 pushes    ////
////////              c(T) down to the PSD boundary c = b^2 = 0.64.      ////
//////////////////////////////////////////////////////////////////////////

#include "psopt.h"

adouble endpoint_cost(adouble* i, adouble* f, adouble* p, adouble& t0,
                      adouble& tf, adouble* xad, int iphase, Workspace* w)
{ return f[1]; }   // minimise c(T): pushes the matrix toward the PSD boundary

adouble integrand_cost(adouble* states, adouble* controls, adouble* p,
                       adouble& time, adouble* xad, int iphase, Workspace* w)
{ adouble ub=controls[0], uc=controls[1]; return 0.01*(ub*ub + uc*uc); }

void dae(adouble* d, adouble* path, adouble* states, adouble* controls,
         adouble* p, adouble& time, adouble* xad, int iphase, Workspace* w)
{
    adouble b = states[0], c = states[1];   // M = [[1, b],[b, c]]
    adouble ub = controls[0], uc = controls[1];
    d[0] = ub;   d[1] = uc;

    // M >= 0  (2x2, leading minors):  M11 = 1 > 0  and  det = c - b^2 >= 0.
    path[0] = b*b - c;                       // <= 0   =>   c >= b^2   =>   M >= 0
}

void events(adouble* e, adouble* i, adouble* f, adouble* p, adouble& t0,
            adouble& tf, adouble* xad, int iphase, Workspace* w)
{
    e[0]=i[0]; e[1]=i[1];   // b(0)=0, c(0)=1  (M0 = I, PSD)
    e[2]=f[0];             // b(T)=0.8
}

void linkages(adouble* linkages, adouble* xad, Workspace* w) {}

int main(void)
{
    Alg algorithm; Sol solution; Prob problem;
    problem.name = "Semidefinite-constrained optimal control (PSD via minors)";
    problem.outfilename = "conic_sdp.txt";
    problem.nphases = 1; problem.nlinkages = 0;
    psopt_level1_setup(problem);

    problem.phases(1).nstates   = 2;   // b, c
    problem.phases(1).ncontrols = 2;   // ub, uc
    problem.phases(1).nevents   = 3;
    problem.phases(1).npath     = 1;   // det minor >= 0
    problem.phases(1).nodes     << 30;
    psopt_level2_setup(problem, algorithm);

    problem.phases(1).bounds.lower.states << -2.0, -0.1;
    problem.phases(1).bounds.upper.states <<  2.0,  3.0;
    problem.phases(1).bounds.lower.controls << -5.0, -5.0;
    problem.phases(1).bounds.upper.controls <<  5.0,  5.0;

    problem.phases(1).bounds.lower.path(0) = -1.0e19;
    problem.phases(1).bounds.upper.path(0) =  0.0;   // b^2 - c <= 0

    problem.phases(1).bounds.lower.events << 0.0, 1.0, 0.8;
    problem.phases(1).bounds.upper.events << 0.0, 1.0, 0.8;

    problem.phases(1).bounds.lower.StartTime = 0.0; problem.phases(1).bounds.upper.StartTime = 0.0;
    problem.phases(1).bounds.lower.EndTime   = 1.0; problem.phases(1).bounds.upper.EndTime   = 1.0;

    problem.integrand_cost = &integrand_cost;
    problem.endpoint_cost  = &endpoint_cost;
    problem.dae            = &dae;
    problem.events         = &events;
    problem.linkages       = &linkages;

    problem.phases(1).guess.states   = zeros(2,30);
    problem.phases(1).guess.states.row(0) = linspace(0.0,0.8,30);
    problem.phases(1).guess.states.row(1) = linspace(1.0,1.0,30);
    problem.phases(1).guess.controls = zeros(2,30);
    problem.phases(1).guess.time     = linspace(0.0,1.0,30);

    algorithm.nlp_method         = "IPOPT";
    algorithm.scaling            = "automatic";
    algorithm.derivatives        = "automatic";
    algorithm.collocation_method = "Legendre";
    algorithm.nlp_iter_max       = 1000;
    algorithm.nlp_tolerance      = 1.e-6;
    algorithm.mesh_refinement    = "automatic";

    psopt(solution, problem, algorithm);

    DMatrix x = solution.get_states_in_phase(1);
    DMatrix u = solution.get_controls_in_phase(1);
    DMatrix t = solution.get_time_in_phase(1);
    Save(x,"x.dat"); Save(u,"u.dat"); Save(t,"t.dat");
    // det(M)(t) = c - b^2 >= 0 throughout; at T it rides the PSD boundary (~0).
    plot(t,x,problem.name,"time (s)","b c","b c");
    plot(t,u,problem.name,"time (s)","controls","ub uc");
}
