//////////////////////////////////////////////////////////////////////////
////////////////           PSOPT  Example             ////////////////////
//////////////////////////////////////////////////////////////////////////
//////// Title:   Minimum-energy double integrator, and a check on   //////
////////          the costates PSOPT returns for it                  //////
//////// Reference: A. Locatelli, "Optimal Control of a Double       //////
////////            Integrator", Springer 2017, Chaps. 3-5.          //////
//////////////////////////////////////////////////////////////////////////
////////
//////// Problem:  xdot1 = x2,  xdot2 = u,   (x1,x2): (0,0) -> (1,0),
////////           tf = 1 fixed,   minimise  J = INT_0^1 u^2/2 dt.
////////
//////// Maximum principle, with the convention H = L + lambda^T f that
//////// PSOPT and its Hamiltonian use:
////////
////////     H          = u^2/2 + l1 x2 + l2 u
////////     dH/du = 0  =>  u = -l2
////////     l1dot      = -dH/dx1 = 0        =>  l1 constant
////////     l2dot      = -dH/dx2 = -l1      =>  l2 linear in t
////////
//////// so the optimal control is linear, and imposing the four boundary
//////// conditions on the resulting cubic x1 gives the whole solution in
//////// closed form:
////////
////////     u*(t)  = 6 - 12 t          x1*(t) = 3t^2 - 2t^3
////////     x2*(t) = 6t - 6t^2         J*     = 6
////////     l1*(t) = -12               l2*(t) = 12t - 6
////////     H*     = -18   (constant: the problem is autonomous)
////////
//////// and the peak speed is x2*(1/2) = 3/2.
////////
//////// WHY THIS EXAMPLE EXISTS. The states are cubic, the control linear
//////// and the cost integrand quadratic, so the Legendre, Radau, Gauss and
//////// Hermite-Simpson transcriptions are all EXACT on this problem: the
//////// discretized problem has the same solution as the continuous one.
//////// That makes it a poor benchmark for accuracy and an excellent one for
//////// the covector mapping, because the discretization contributes nothing
//////// and any error left in the reported costate has nowhere to hide.
////////
//////// The run therefore prints, alongside the usual state and control
//////// errors, the two costates against their closed forms and two checks
//////// that need no closed form at all -- the stationarity residual
//////// dH/du = u + l2, formed from PSOPT's own control and costate, and the
//////// constancy of the Hamiltonian. Those two are what a user has on a
//////// problem whose answer is not known, and this example shows what they
//////// look like when everything is right.
////////
//////// Usage:  ./mineng_di [collocation_method] [nodes]
////////         collocation_method: Legendre (default), Chebyshev, Radau,
////////                             Gauss, trapezoidal or Hermite-Simpson
////////         nodes:              default 40
////////
//////// The trapezoidal scheme is second order and is NOT exact here, so it
//////// is the one run that shows a discretization error rather than a
//////// covector-mapping one; its errors fall like 1/N^2 and reading them is
//////// part of the point.
//////////////////////////////////////////////////////////////////////////

#include "psopt.h"

#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <string>

using namespace std;
using namespace PSOPT;

//////////////////////////////////////////////////////////////////////////
///////////////////  The closed-form solution  ///////////////////////////
//////////////////////////////////////////////////////////////////////////

static double u_exact  (double t) { return 6.0 - 12.0*t;      }
static double x1_exact (double t) { return 3.0*t*t - 2.0*t*t*t; }
static double x2_exact (double t) { return 6.0*t - 6.0*t*t;   }
static double l1_exact (double t) { (void) t; return -12.0;   }
static double l2_exact (double t) { return 12.0*t - 6.0;      }

static const double J_EXACT = 6.0;
static const double H_EXACT = -18.0;

//////////////////////////////////////////////////////////////////////////
///////////////////  Problem functions  //////////////////////////////////
//////////////////////////////////////////////////////////////////////////

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

//////////////////////////////////////////////////////////////////////////

int main(int argc, char** argv)
{
    const string method = (argc > 1) ? argv[1] : "Legendre";
    const int    N      = (argc > 2) ? atoi(argv[2]) : 40;

    Alg algorithm; Sol solution; Prob problem;
    problem.name = "Minimum-energy double integrator (Locatelli Ch 3)";
    problem.outfilename = "mineng_di.txt";
    problem.nphases = 1; problem.nlinkages = 0;
    psopt_level1_setup(problem);

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
    algorithm.collocation_method = method;
    algorithm.nlp_iter_max       = 1000;
    // Tighter than the 1e-6 an ordinary example would ask for. The point here is
    // to measure the covector mapping, and at 1e-6 the NLP's own termination is
    // what the costate error would be reporting.
    algorithm.nlp_tolerance      = 1.e-10;

    if ( psopt(solution, problem, algorithm) != 0 || solution.error_flag ) {
        printf("\nThe problem was not solved: %s\n", solution.error_msg.c_str());
        return 1;
    }

    MatrixXd x   = solution.get_states_in_phase(1);
    MatrixXd u   = solution.get_controls_in_phase(1);
    MatrixXd t   = solution.get_time_in_phase(1);
    MatrixXd lam = solution.get_dual_costates_in_phase(1);
    MatrixXd H   = solution.get_dual_hamiltonian_in_phase(1);

    const long M = t.cols();

    //////////////////////////////////////////////////////////////////////
    ///////  What the solution should be, and what it is  ////////////////
    //////////////////////////////////////////////////////////////////////

    MatrixXd exact(4, M);                       // l1*, l2* and, below, the reported pair
    double ex1 = 0.0, ex2 = 0.0, eu = 0.0;
    double el1 = 0.0, el2 = 0.0, stat = 0.0;
    double Hmin =  1.0e300, Hmax = -1.0e300, eH = 0.0;

    for (long k = 0; k < M; k++) {
        const double tk = t(0,k);

        ex1  = max(ex1,  fabs(x(0,k) - x1_exact(tk)));
        ex2  = max(ex2,  fabs(x(1,k) - x2_exact(tk)));
        eu   = max(eu,   fabs(u(0,k) - u_exact (tk)));
        el1  = max(el1,  fabs(lam(0,k) - l1_exact(tk)));
        el2  = max(el2,  fabs(lam(1,k) - l2_exact(tk)));

        // dH/du = u + l2, formed from PSOPT's own output and needing no closed form
        stat = max(stat, fabs(u(0,k) + lam(1,k)));

        const double Hk = H(0,k);
        Hmin = min(Hmin, Hk);  Hmax = max(Hmax, Hk);
        eH   = max(eH, fabs(Hk - H_EXACT));

        exact(0,k) = lam(0,k);  exact(1,k) = lam(1,k);
        exact(2,k) = l1_exact(tk);  exact(3,k) = l2_exact(tk);
    }

    printf("\n\n");
    printf("=====================================================================\n");
    printf("  Minimum-energy double integrator: the solution against its closed\n");
    printf("  form, and the costates against the maximum principle.\n");
    printf("  Collocation method: %-16s   nodes: %ld\n", method.c_str(), M);
    printf("=====================================================================\n");
    printf("\n  The primal solution\n");
    printf("    J                                %.12f   (exact %.1f)\n", solution.cost, J_EXACT);
    printf("    |J - J*|                         %10.3e\n", fabs(solution.cost - J_EXACT));
    printf("    max |x1 - x1*|                   %10.3e\n", ex1);
    printf("    max |x2 - x2*|                   %10.3e\n", ex2);
    printf("    max |u  - u* |                   %10.3e   (|u*| <= 6)\n", eu);

    printf("\n  The costates, against l1* = -12 and l2* = 12t - 6\n");
    printf("    max |l1 - l1*|                   %10.3e   (|l1*| = 12)\n", el1);
    printf("    max |l2 - l2*|                   %10.3e   (|l2*| <= 6)\n", el2);
    printf("    l1 at the first and last node    %12.6f  %12.6f\n", lam(0,0), lam(0,M-1));
    printf("    l2 at the first and last node    %12.6f  %12.6f   (exact -6, +6)\n",
           lam(1,0), lam(1,M-1));

    printf("\n  Two checks that need no closed form\n");
    printf("    max |dH/du| = max |u + l2|       %10.3e\n", stat);
    printf("    H over the mesh                  [%.8f, %.8f]\n", Hmin, Hmax);
    printf("    H spread (H is constant here)    %10.3e\n", Hmax - Hmin);
    printf("    max |H - H*|                     %10.3e   (H* = -18)\n", eH);

    printf("\n  Reading these numbers\n");
    printf("    The states are cubic, the control linear and the integrand quadratic,\n");
    printf("    so Legendre, Chebyshev, Radau, Gauss and Hermite-Simpson are all EXACT\n");
    printf("    on this problem and should report every error above near the level of\n");
    printf("    the linear algebra, 1e-10 to 1e-6. Anything larger from one of those\n");
    printf("    five is the covector mapping and not the discretization, since the\n");
    printf("    discretization has no error here to be confused with.\n");
    printf("\n");
    printf("    The trapezoidal scheme is the exception. It is second order and not\n");
    printf("    exact here, so its errors are a genuine discretization error and fall\n");
    printf("    like 1/N^2: about 3e-2 on the costates at N = 40, a quarter of that at\n");
    printf("    N = 80. Run it at a few node counts and watch the ratio.\n");
    printf("=====================================================================\n\n");

    //////////////////////////////////////////////////////////////////////
    ///////  Files and figures  //////////////////////////////////////////
    //////////////////////////////////////////////////////////////////////

    Save(x,"x.dat"); Save(u,"u.dat"); Save(t,"t.dat"); Save(lam,"lambda.dat");

    plot(t,x,problem.name + ": states","time (s)","x1 x2","x1 x2");
    plot(t,u,problem.name + ": control","time (s)","u = 6 - 12 t","u");
    plot(t,exact,problem.name + ": costates against the closed form",
         "time (s)","costates","l1 l2 l1* l2*");
    plot(t,H,problem.name + ": Hamiltonian","time (s)","H","H");

    plot(t,x,problem.name + ": states","time (s)","x1 x2","x1 x2",
         "pdf","mineng_di_states.pdf");
    plot(t,u,problem.name + ": control","time (s)","u = 6 - 12 t","u",
         "pdf","mineng_di_control.pdf");
    plot(t,exact,problem.name + ": costates against the closed form",
         "time (s)","costates","l1 l2 l1* l2*",
         "pdf","mineng_di_costates.pdf");

    return 0;
}
