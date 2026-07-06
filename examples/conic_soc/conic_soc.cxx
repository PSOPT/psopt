//////////////////////////////////////////////////////////////////////////
////////////////           PSOPT  Example             ////////////////////
//////////////////////////////////////////////////////////////////////////
//////// Title:   SECOND-ORDER CONE (SOC) constraints in optimal control ///
//////// Illustrates: a Lorentz / second-order cone constraint          ////
////////              || (ax,ay) ||_2 <= amax   as a smooth CONVEX path  ////
////////              constraint. The norm form sqrt(ax^2+ay^2) - amax   ////
////////              (tiny apex smoothing) conditions better than the   ////
////////              squared form for IPOPT.                            ////
//////// Problem:     minimum-energy 2D point mass from (0,0) at rest to   ////
////////              (1, 0.5) at rest (T=2.5); the thrust vector rides  ////
////////              the second-order cone (||a|| = amax on the arc).   ////
//////////////////////////////////////////////////////////////////////////

#include "psopt.h"

adouble endpoint_cost(adouble* i, adouble* f, adouble* p, adouble& t0,
                      adouble& tf, adouble* xad, int iphase, Workspace* w)
{ return 0.0; }

adouble integrand_cost(adouble* states, adouble* controls, adouble* p,
                       adouble& time, adouble* xad, int iphase, Workspace* w)
{ adouble ax=controls[0], ay=controls[1]; return ax*ax + ay*ay; }  // minimum control energy

void dae(adouble* d, adouble* path, adouble* states, adouble* controls,
         adouble* p, adouble& time, adouble* xad, int iphase, Workspace* w)
{
    adouble vx = states[2], vy = states[3];
    adouble ax = controls[0], ay = controls[1];   // thrust vector
    d[0] = vx;   d[1] = vy;
    d[2] = ax;   d[3] = ay;

    const double amax = 0.8;
    // thrust second-order cone in squared (smooth, convex) form
    path[0] = sqrt(ax*ax + ay*ay + 1.0e-8) - amax;   // ||(ax,ay)||_2 <= amax (norm form)
}

void events(adouble* e, adouble* i, adouble* f, adouble* p, adouble& t0,
            adouble& tf, adouble* xad, int iphase, Workspace* w)
{
    e[0]=i[0]; e[1]=i[1]; e[2]=i[2]; e[3]=i[3];   // start at (0,0), at rest
    e[4]=f[0]; e[5]=f[1]; e[6]=f[2]; e[7]=f[3];   // reach (1,1), at rest
}

void linkages(adouble* linkages, adouble* xad, Workspace* w) {}

int main(void)
{
    Alg algorithm; Sol solution; Prob problem;
    problem.name = "Second-order cone constrained minimum-energy point mass";
    problem.outfilename = "conic_soc.txt";
    problem.nphases = 1; problem.nlinkages = 0;
    psopt_level1_setup(problem);

    problem.phases(1).nstates   = 4;   // x, y, vx, vy
    problem.phases(1).ncontrols = 2;   // ax, ay
    problem.phases(1).nevents   = 8;
    problem.phases(1).npath     = 1;   // thrust cone
    problem.phases(1).nodes     << 40;
    psopt_level2_setup(problem, algorithm);

    problem.phases(1).bounds.lower.states << -5.0, -5.0, -5.0, -5.0;
    problem.phases(1).bounds.upper.states <<  5.0,  5.0,  5.0,  5.0;
    problem.phases(1).bounds.lower.controls << -2.0, -2.0;
    problem.phases(1).bounds.upper.controls <<  2.0,  2.0;

    problem.phases(1).bounds.lower.path(0) = -1.0e19;
    problem.phases(1).bounds.upper.path(0) =  0.0;

    problem.phases(1).bounds.lower.events << 0.0, 0.0, 0.0, 0.0, 1.0, 0.5, 0.0, 0.0;
    problem.phases(1).bounds.upper.events << 0.0, 0.0, 0.0, 0.0, 1.0, 0.5, 0.0, 0.0;

    problem.phases(1).bounds.lower.StartTime = 0.0; problem.phases(1).bounds.upper.StartTime = 0.0;
    problem.phases(1).bounds.lower.EndTime   = 2.5; problem.phases(1).bounds.upper.EndTime   = 2.5;

    problem.integrand_cost = &integrand_cost;
    problem.endpoint_cost  = &endpoint_cost;
    problem.dae            = &dae;
    problem.events         = &events;
    problem.linkages       = &linkages;

    problem.phases(1).guess.states   = zeros(4,40);
    problem.phases(1).guess.states.row(0) = linspace(0.0,1.0,40);
    problem.phases(1).guess.states.row(1) = linspace(0.0,0.5,40);
    problem.phases(1).guess.controls = zeros(2,40);
    problem.phases(1).guess.time     = linspace(0.0,2.5,40);

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
    plot(t,x,problem.name,"time (s)","x y vx vy","x y vx vy");
    plot(t,u,problem.name,"time (s)","thrust","ax ay");
}
