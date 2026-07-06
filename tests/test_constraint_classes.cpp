//////////////////////////////////////////////////////////////////////////////
// test_constraint_classes.cpp
//
// Extends the double-integrator cross-comparison to the CONSTRAINT-CLASS axis:
// each native/advanced constraint type is exercised on a double integrator and
// (a) its constraint is checked satisfied, (b) the result is invariant across
// derivative modes {analytic, finite-diff} -- a mutual cross-check.
//
//   - PATH / state constraint  -> Bryson-Denham, analytical J* = 4.
//   - INTERIOR-POINT (native)  -> waypoint x1(tf/2) = 0.7 hit; cost invariant.
//   - INTEGRAL (native)        -> INT x1 dt = 0.7 satisfied; cost invariant.
//////////////////////////////////////////////////////////////////////////////

#include "gtest/gtest.h"
#include <psopt.h>
#include <cmath>

// ================= (1) PATH / state constraint: Bryson-Denham ================
namespace cc_path {
adouble endpoint_cost(adouble*, adouble*, adouble*, adouble&, adouble&, adouble*, int, Workspace*) { return 0.0; }
adouble integrand_cost(adouble*, adouble* c, adouble*, adouble&, adouble*, int, Workspace*) { adouble u=c[0]; return 0.5*u*u; }
void dae(adouble* d, adouble* path, adouble* s, adouble* c, adouble*, adouble&, adouble*, int, Workspace*)
{ d[0]=s[1]; d[1]=c[0]; path[0]=s[0]; }                          // state x1 (constrained <= 1/9)
void events(adouble* e, adouble* i, adouble* f, adouble*, adouble&, adouble&, adouble*, int, Workspace*)
{ e[0]=i[0]; e[1]=i[1]; e[2]=f[0]; e[3]=f[1]; }                  // x1(0)=0,v(0)=1,x1(1)=0,v(1)=-1
void linkages(adouble*, adouble*, Workspace*) {}
}

TEST(ConstraintClasses, PathState_BrysonDenham_J4)
{
    using namespace cc_path;
    const int N = 40; const double L = 1.0/9.0;
    Alg algorithm; Sol solution; Prob problem;
    problem.name="bryson-denham"; problem.outfilename="cc.txt"; problem.nphases=1; problem.nlinkages=0;
    psopt_level1_setup(problem);
    problem.phases(1).nstates=2; problem.phases(1).ncontrols=1; problem.phases(1).nevents=4;
    problem.phases(1).npath=1; problem.phases(1).nodes<<N;
    psopt_level2_setup(problem, algorithm);
    problem.phases(1).bounds.lower.states << -1.0,-10.0;  problem.phases(1).bounds.upper.states << 1.0,10.0;
    problem.phases(1).bounds.lower.controls(0)=-20.0;     problem.phases(1).bounds.upper.controls(0)=20.0;
    problem.phases(1).bounds.lower.path(0)=-1e19;         problem.phases(1).bounds.upper.path(0)=L;   // x1 <= 1/9
    problem.phases(1).bounds.lower.events << 0.0,1.0,0.0,-1.0;
    problem.phases(1).bounds.upper.events << 0.0,1.0,0.0,-1.0;
    problem.phases(1).bounds.lower.StartTime=0.0; problem.phases(1).bounds.upper.StartTime=0.0;
    problem.phases(1).bounds.lower.EndTime=1.0;   problem.phases(1).bounds.upper.EndTime=1.0;
    problem.integrand_cost=&integrand_cost; problem.endpoint_cost=&endpoint_cost;
    problem.dae=&dae; problem.events=&events; problem.linkages=&linkages;
    problem.phases(1).guess.states=zeros(2,N);
    problem.phases(1).guess.states.row(1)=linspace(1.0,-1.0,N);
    problem.phases(1).guess.controls=zeros(1,N); problem.phases(1).guess.time=linspace(0.0,1.0,N);
    algorithm.nlp_method="IPOPT"; algorithm.scaling="automatic"; algorithm.derivatives="automatic";
    algorithm.collocation_method="Legendre"; algorithm.nlp_iter_max=1000; algorithm.nlp_tolerance=1e-6;
    psopt(solution, problem, algorithm);
    ASSERT_EQ(solution.error_flag, false);
    DMatrix x=solution.get_states_in_phase(1); int n=(int)x.cols();
    EXPECT_NEAR(solution.get_cost(), 4.0, 5.0e-3) << "Bryson-Denham J* = 4/(9L) = 4";
    double x1max=-1e9; for(int k=0;k<n;++k) x1max=std::max(x1max,x(0,k));
    EXPECT_LE(x1max, L + 1.0e-4) << "state path constraint x1 <= 1/9 respected";
}

// ================= min-energy double integrator with a native constraint ====
// Shared: (0,0)->(1,0), tf=1, min INT u^2/2 dt. Returns cost via `deriv` mode.
namespace cc_di {
adouble endpoint_cost(adouble*, adouble*, adouble*, adouble&, adouble&, adouble*, int, Workspace*) { return 0.0; }
adouble integrand_cost(adouble*, adouble* c, adouble*, adouble&, adouble*, int, Workspace*) { adouble u=c[0]; return 0.5*u*u; }
void dae(adouble* d, adouble*, adouble* s, adouble* c, adouble*, adouble&, adouble*, int, Workspace*) { d[0]=s[1]; d[1]=c[0]; }
void events(adouble* e, adouble* i, adouble* f, adouble*, adouble&, adouble&, adouble*, int, Workspace*)
{ e[0]=i[0]; e[1]=i[1]; e[2]=f[0]; e[3]=f[1]; }
void linkages(adouble*, adouble*, Workspace*) {}
// interior-point waypoint: x1 at the interior time must equal 0.7
void interior(adouble* g, adouble* xi, adouble*, adouble*, adouble&, int, int, Workspace*) { g[0]=xi[0]-0.7; }
// integral: integrand is x1 (=> INT x1 dt bounded)
void integral(adouble* q, adouble* s, adouble*, adouble*, adouble&, int, Workspace*) { q[0]=s[0]; }
}

// Build+solve the min-energy DI with an optional native interior/integral constraint.
static double solve_native(const char* deriv, int which /*1=interior,2=integral*/, double* cons_out)
{
    using namespace cc_di;
    const int N=40; Alg algorithm; Sol solution; Prob problem;
    problem.name="native"; problem.outfilename="cc.txt"; problem.nphases=1; problem.nlinkages=0;
    psopt_level1_setup(problem);
    problem.phases(1).nstates=2; problem.phases(1).ncontrols=1; problem.phases(1).nevents=4; problem.phases(1).npath=0;
    if (which==1) problem.phases(1).ninterior=1; else problem.phases(1).nintegral=1;
    problem.phases(1).nodes<<N;
    psopt_level2_setup(problem, algorithm);
    problem.phases(1).bounds.lower.states << -5.0,-5.0;  problem.phases(1).bounds.upper.states << 5.0,5.0;
    problem.phases(1).bounds.lower.controls(0)=-50.0;    problem.phases(1).bounds.upper.controls(0)=50.0;
    problem.phases(1).bounds.lower.events << 0.0,0.0,1.0,0.0;
    problem.phases(1).bounds.upper.events << 0.0,0.0,1.0,0.0;
    problem.phases(1).bounds.lower.StartTime=0.0; problem.phases(1).bounds.upper.StartTime=0.0;
    problem.phases(1).bounds.lower.EndTime=1.0;   problem.phases(1).bounds.upper.EndTime=1.0;
    if (which==1) { problem.phases(1).interior_time(0)=0.5;
                    problem.phases(1).bounds.lower.interior(0)=0.0; problem.phases(1).bounds.upper.interior(0)=0.0;
                    problem.interior_point_constraints=&interior; }
    else          { problem.phases(1).bounds.lower.integral(0)=0.7; problem.phases(1).bounds.upper.integral(0)=0.7;
                    problem.integral_constraints=&integral; }
    problem.integrand_cost=&integrand_cost; problem.endpoint_cost=&endpoint_cost;
    problem.dae=&dae; problem.events=&events; problem.linkages=&linkages;
    problem.phases(1).guess.states=zeros(2,N);
    problem.phases(1).guess.states.row(0)=linspace(0.0,1.0,N);
    problem.phases(1).guess.controls=zeros(1,N); problem.phases(1).guess.time=linspace(0.0,1.0,N);
    algorithm.nlp_method="IPOPT"; algorithm.scaling="automatic"; algorithm.derivatives=deriv;
    algorithm.collocation_method="Legendre"; algorithm.nlp_iter_max=1000; algorithm.nlp_tolerance=1e-6;
    psopt(solution, problem, algorithm);
    EXPECT_EQ(solution.error_flag, false);
    DMatrix x=solution.get_states_in_phase(1), t=solution.get_time_in_phase(1);
    int n=(int)t.cols();
    if (cons_out) {
        if (which==1) { // x1 at t=0.5, linearly interpolated between the bracketing
                        // nodes (LGL nodes rarely land exactly on the interior time)
            *cons_out = x(0,0);
            for(int k=0;k<n-1;++k) if (t(0,k)<=0.5 && t(0,k+1)>=0.5) {
                double a=(0.5-t(0,k))/(t(0,k+1)-t(0,k)); *cons_out=(1.0-a)*x(0,k)+a*x(0,k+1); break; } }
        else { double s=0; for(int k=1;k<n;++k) s+=0.5*(x(0,k)+x(0,k-1))*(t(0,k)-t(0,k-1)); *cons_out=s; } // trapz INT x1 dt
    }
    return solution.get_cost();
}

TEST(ConstraintClasses, InteriorPoint_Waypoint_Invariant)
{
    double c_a, c_fd, wp_a, wp_fd;
    c_a  = solve_native("automatic", 1, &wp_a);
    c_fd = solve_native("numerical", 1, &wp_fd);
    EXPECT_NEAR(wp_a, 0.7, 2.0e-2) << "interior waypoint x1(tf/2) = 0.7 hit";
    EXPECT_GT(c_a, 6.0)            << "waypoint above the natural 0.5 raises the cost";
    EXPECT_NEAR(c_a, c_fd, 1.0e-2) << "cost invariant: analytic == finite-diff derivatives";
}

TEST(ConstraintClasses, Integral_Isoperimetric_Invariant)
{
    double c_a, c_fd, I_a, I_fd;
    c_a  = solve_native("automatic", 2, &I_a);
    c_fd = solve_native("numerical", 2, &I_fd);
    EXPECT_NEAR(I_a, 0.7, 1.0e-2)  << "integral constraint INT x1 dt = 0.7 satisfied";
    EXPECT_GT(c_a, 6.0)            << "INT x1 dt = 0.7 (vs natural 0.5) raises the cost";
    EXPECT_NEAR(c_a, c_fd, 1.0e-2) << "cost invariant: analytic == finite-diff derivatives";
}
