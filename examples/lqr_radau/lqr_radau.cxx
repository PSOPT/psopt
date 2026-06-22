// lqr_radau.cxx  --  Radau end-to-end validation inside PSOPT (fixed terminal).
//
//   min  1/2 \int_0^T (x^2 + u^2) dt,   x' = u,   x(0)=1,   x(T)=xf
//   analytic:  x = A e^t + B e^{-t},  u = A e^t - B e^{-t},  lambda = -u,
//   with A = (xf - e^{-T})/(e^T - e^{-T}),  B = 1 - A.   (u(T) != 0, stresses the
//   terminal-control interpolation pin.)

#include "psopt.h"
#include <cmath>
using namespace std;

#ifndef COLLOC
#define COLLOC "Radau"
#endif
static const double Tf = 2.0;
static const double XF = 0.5;

adouble endpoint_cost(adouble* x0, adouble* xf, adouble* p,
                      adouble& t0, adouble& tf, adouble* xad, int iphase, Workspace* ws)
{ return 0.0; }

adouble integrand_cost(adouble* states, adouble* controls, adouble* p,
                       adouble& time, adouble* xad, int iphase, Workspace* ws)
{ adouble x = states[0], u = controls[0]; return 0.5*(x*x + u*u); }

void dae(adouble* derivatives, adouble* path, adouble* states, adouble* controls,
         adouble* p, adouble& time, adouble* xad, int iphase, Workspace* ws)
{ derivatives[0] = controls[0]; }

void events(adouble* e, adouble* x0, adouble* xf, adouble* p,
            adouble& t0, adouble& tf, adouble* xad, int iphase, Workspace* ws)
{ e[0] = x0[0]; e[1] = xf[0]; }

void linkages(adouble* linkages, adouble* xad, Workspace* ws) {}

int main(void)
{
    Alg algorithm; Sol solution; Prob problem;
    problem.name = "LQR Radau validation"; problem.outfilename = "lqr.txt";
    problem.nphases = 1; problem.nlinkages = 0;
    psopt_level1_setup(problem);

    problem.phases(1).nstates = 1; problem.phases(1).ncontrols = 1;
    problem.phases(1).nevents = 2; problem.phases(1).npath = 0;
    problem.phases(1).nodes = (RowVectorXi(1) << 30).finished();
    psopt_level2_setup(problem, algorithm);

    problem.phases(1).bounds.lower.states(0)   = -50.0;
    problem.phases(1).bounds.upper.states(0)   =  50.0;
    problem.phases(1).bounds.lower.controls(0) = -50.0;
    problem.phases(1).bounds.upper.controls(0) =  50.0;
    problem.phases(1).bounds.lower.events(0) = 1.0;  problem.phases(1).bounds.upper.events(0) = 1.0;
    problem.phases(1).bounds.lower.events(1) = XF;   problem.phases(1).bounds.upper.events(1) = XF;
    problem.phases(1).bounds.lower.StartTime = 0.0;  problem.phases(1).bounds.upper.StartTime = 0.0;
    problem.phases(1).bounds.lower.EndTime   = Tf;   problem.phases(1).bounds.upper.EndTime   = Tf;

    problem.integrand_cost = &integrand_cost; problem.endpoint_cost = &endpoint_cost;
    problem.dae = &dae; problem.events = &events; problem.linkages = &linkages;

    problem.phases(1).guess.controls = zeros(1,40);
    problem.phases(1).guess.states   = linspace(1.0,XF,40);
    problem.phases(1).guess.time     = linspace(0.0,Tf,40);

    algorithm.nlp_method = "IPOPT"; algorithm.scaling = "automatic";
    algorithm.derivatives = "automatic"; algorithm.collocation_method = COLLOC;
    algorithm.mesh_refinement = "manual"; algorithm.nlp_iter_max = 1000;
    algorithm.nlp_tolerance = 1.e-9;

    psopt(solution, problem, algorithm);

    MatrixXd x = solution.get_states_in_phase(1);
    MatrixXd u = solution.get_controls_in_phase(1);
    MatrixXd t = solution.get_time_in_phase(1);
    MatrixXd lam = solution.get_dual_costates_in_phase(1);

    double A = (XF - exp(-Tf))/(exp(Tf) - exp(-Tf));
    double B = 1.0 - A;
    auto xe = [&](double tt){ return A*exp(tt) + B*exp(-tt); };
    auto ue = [&](double tt){ return A*exp(tt) - B*exp(-tt); };
    auto le = [&](double tt){ return -(A*exp(tt) - B*exp(-tt)); };

    int n = t.cols();
    printf("\nmethod = %s,  nodes = %d\n", COLLOC, n);
    double xm=0, um=0, lm=0;
    for (int k=0;k<n;k++){
        double tk=t(0,k);
        xm = max(xm, fabs(x(0,k)-xe(tk)));
        if (k<n-1) { um = max(um, fabs(u(0,k)-ue(tk))); lm = max(lm, fabs(lam(0,k)-le(tk))); }
    }
    printf("max |x - x_exact|             = %.3e  (all nodes)\n", xm);
    printf("max |u - u_exact|             = %.3e  (collocation nodes)\n", um);
    printf("max |lambda - lambda_exact|   = %.3e  (collocation nodes)\n", lm);
    printf("\nterminal-control pin check (non-collocated node, tau=+1):\n");
    printf("   u_pinned(T) = %.8f   u_exact(T) = %.8f   err = %.2e\n",
           u(0,n-1), ue(Tf), fabs(u(0,n-1)-ue(Tf)));
    printf("terminal-costate recovery check (non-collocated node):\n");
    printf("   lambda(T)   = %.8f   lambda_exact(T) = %.8f   err = %.2e\n",
           lam(0,n-1), le(Tf), fabs(lam(0,n-1)-le(Tf)));
    return 0;
}
