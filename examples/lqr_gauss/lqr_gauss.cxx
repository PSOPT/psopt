// lqr_gauss.cxx -- Legendre-Gauss collocation validation inside PSOPT (increment 1).
//
//   min 1/2 \int_0^T (x^2 + u^2) dt,  x'=u,  x(0)=1,  x(T)=0.5,  T=2
//   analytic: x = A e^t + B e^{-t},  u = A e^t - B e^{-t},
//   A = (xf - e^{-T})/(e^T - e^{-T}), B = 1 - A.
//
// Increment 1 checks state (all stored nodes) and control (Gauss nodes only; the
// initial node is non-collocated and its control is unpinned for now). The fixed
// terminal couples the appended terminal-state variable x_f (via the event) to the
// solution, so matching the analytic trajectory also exercises the Gauss-quadrature
// defining constraint.

#include "psopt.h"
#include <cmath>
using namespace std;

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
    problem.name = "LQR Gauss validation"; problem.outfilename = "lqrg.txt";
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
    algorithm.derivatives = "automatic"; algorithm.collocation_method = "Gauss";
    algorithm.mesh_refinement = "manual"; algorithm.nlp_iter_max = 1000;
    algorithm.nlp_tolerance = 1.e-9;

    psopt(solution, problem, algorithm);

    MatrixXd x = solution.get_states_in_phase(1);
    MatrixXd u = solution.get_controls_in_phase(1);
    MatrixXd t = solution.get_time_in_phase(1);

    double A = (XF - exp(-Tf))/(exp(Tf) - exp(-Tf));
    double B = 1.0 - A;
    auto xe = [&](double tt){ return A*exp(tt) + B*exp(-tt); };
    auto ue = [&](double tt){ return A*exp(tt) - B*exp(-tt); };

    int n = t.cols();
    printf("\nmethod = Gauss,  stored nodes = %d  (node 0 = initial, non-collocated)\n", n);
    double xm=0, um=0;
    for (int k=0;k<n;k++){
        double tk=t(0,k);
        xm = max(xm, fabs(x(0,k)-xe(tk)));
        if (k>0) um = max(um, fabs(u(0,k)-ue(tk)));   // exclude non-collocated initial control
    }
    printf("max |x - x_exact|       = %.3e  (all stored nodes)\n", xm);
    printf("max |u - u_exact|       = %.3e  (Gauss nodes, excl. initial)\n", um);
    printf("state at t=0            = %.8f  (should be 1.0)\n", x(0,0));
    printf("state at last Gauss pt  = %.8f  (t=%.4f, exact %.8f)\n",
           x(0,n-1), t(0,n-1), xe(t(0,n-1)));
    return 0;
}
