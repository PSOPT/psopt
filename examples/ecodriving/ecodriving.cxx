//////////////////////////////////////////////////////////////////////////
////////////////           PSOPT  Example             ////////////////////
//////////////////////////////////////////////////////////////////////////
//////// Title:   Eco-driving an electric vehicle through a corner,   /////
////////          trading energy against what the driver wants        /////
//////// Reference: R. Lot, J. Fleming, B. Chen and S. Evangelou,     /////
////////            "Eco-driving optimal control for electric         /////
////////            vehicles with driver preferences", Transportation /////
////////            Engineering 19 (2025) 100302,                     /////
////////            doi:10.1016/j.treng.2025.100302 (CC BY 4.0).      /////
////////            Section 4.1, the cornering case.                  /////
//////////////////////////////////////////////////////////////////////////
////////
//////// THE PROBLEM. A front-wheel-drive battery electric vehicle drives
//////// 800 m of level road: two straights joined by a 90 degree curve of
//////// radius 30 m. It starts and finishes at the same speed, which is
//////// itself free, so the run is one period of driving a road with such a
//////// curve every 800 m. The final time is free.
////////
////////     xdot = v
////////     m vdot = Fp - Fr - Ff - (1/2) rho Cd A v^2 - m g c_rr
////////
//////// with three non-negative controls: Fp the propulsive force, Fr the
//////// regenerative braking force and Ff the friction braking force.
//////// Splitting the braking in two is what lets the powertrain losses and
//////// the brake-disc losses be told apart, and it is what makes the
//////// regenerative bias constraint expressible.
////////
//////// WHAT IS BEING MINIMISED. Two things at once, blended by a weight:
////////
////////     J = int_0^T [ (1-w) L'_d  +  (w/L0)(L_e + L_v + L_i) ] dt
////////
//////// L'_d is a model of what a human driver wants -- a preferred speed,
//////// gentle accelerations, a comfortable cornering speed -- calibrated to
//////// behave like the Intelligent Driver Model. L_e, L_v and L_i are the
//////// powertrain, vehicle and idle losses, whose integral is the energy
//////// drawn from the battery. At w = 0 the car is driven as a person would
//////// drive it; at w = 1 as an energy-minimising machine would; in between
//////// it traces a Pareto frontier of energy against journey time.
////////
//////// WHY THIS EXAMPLE EXISTS. Three reasons.
////////
//////// First, it is a small problem with an unusual boundary condition: the
//////// initial speed is free and the final speed is only required to equal
//////// it. That is a periodicity condition rather than a boundary value, and
//////// it is what makes the energy balance exact -- with no change in kinetic
//////// energy over the run, minimising the losses IS minimising the energy.
////////
//////// Second, it carries a modelling trick worth teaching. The regenerative
//////// braking bias is properly an equality with a min() in it, because the
//////// front/rear brake split is fixed until the motor saturates. Written
//////// that way the problem is hard to solve. Written as the inequality
//////// Fr <= beta_r (Fr + Ff) it is easy, and the inequality is active at the
//////// optimum whenever the motor is not saturated, because regenerating is
//////// cheaper than heating a brake disc. The example checks afterwards that
//////// it was active, which is the only thing that makes the relaxation
//////// legitimate.
////////
//////// Third, the energy is integrated twice, along two paths that are equal
//////// only if the dynamics hold:
////////
////////     E1 = int (Fm v + L_e) dt        (from the powertrain side)
////////     E2 = int (L_v + L_e) dt         (from the dissipation side)
////////
//////// They differ by the change in kinetic energy at EVERY instant, and the
//////// largest discrepancy over the mesh is a check on the whole transcription
//////// that needs no reference solution and no closed form. It is what a user
//////// has on a problem nobody has solved before.
////////
//////// Note the word "mesh". At t = T the two integrals agree to machine
//////// precision whatever the mesh, because the periodicity condition makes
//////// both sides zero and the transcription enforces that exactly. The
//////// endpoint is therefore not a check at all; it reports success on any
//////// discretisation whatever. A conservation law evaluated where the
//////// boundary conditions already force agreement measures nothing.
////////
//////// Usage:   ecodriving [w] [L0] [nodes] [collocation_method]
////////          ecodriving pareto [nodes]
////////          ecodriving corner [nodes]
////////          ecodriving eq28   [nodes]
////////
////////     w        energy weight in [0,1] (default 0)
////////     L0       normalising constant in kW (default 6)
////////     nodes    number of collocation nodes (default 60)
////////
//////// pareto traces the frontier from w = 0 upwards by continuation, each
//////// solve started from the one before it, and then traces it again with
//////// every point started from the same cold guess, so that the difference
//////// can be seen rather than asserted. On this problem there is none: all
//////// thirteen points land in the same place either way. That is a result
//////// and not a reason to skip the check -- the failure it guards against
//////// shows up as structure in a plot and in no exit status, so the only way
//////// to know it did not happen is to have looked. (Contrast
//////// examples/lts_costates, whose sweep must NOT be continued, because
//////// there the differences between the points are the measurement.)
//////// corner varies the corner geometry, which the source paper gives only
//////// as "approximately 30 m".
//////// eq28 settles which reading of the source paper's driver cost was
//////// intended, by solving w = 0 under both.
////////
//////// Every file a run writes carries that run's weight in its name, so no
//////// sweep can overwrite the single run's output.
////////
//////// WHAT WE CHOSE, WHERE THE PAPER IS SILENT. Stated here so that anyone
//////// comparing against the paper knows what is theirs and what is ours.
////////
////////  1. v_d = 70 km/h. Table 1 of the paper gives 31.2 m/s, but Section
////////     4.1 says "v_d = 70 km/h in the driver model", and the reported
////////     mean speed of 62.1 km/h can only belong to the latter. Table 1
////////     appears to hold the motorway value used in their Section 4.4.
////////  2. The acceleration penalties of the paper's Eq. (28) are printed
////////     with m*u_a and m*u_b in the denominators. Since Fp/m IS the
////////     acceleration, that makes each term identically one; they are read
////////     here as m*a and m*b, the preferred accelerations of Table 1, as
////////     in the paper's own Eq. (14).
////////  3. The corner is taken as exactly R = 30 m over exactly 90 degrees,
////////     with equal straights, and the curvature is blended over a
////////     transition length CORNER_EPS. A road cannot step in curvature and
////////     neither can a differentiable path constraint. Run "corner" to see
////////     what these two choices are worth.
////////  4. L0 is not given in the paper at all. It does not change the set of
////////     solutions -- dividing by (1-w) shows the problem depends only on
////////     lambda = w/((1-w)L0) -- so it only decides which w labels which
////////     point of the frontier. 6 kW is chosen here because it spreads the
////////     frontier evenly over w in [0,1].
////////  5. The road speed limit is set to 90 km/h. The paper says only that
////////     it is "just above the legal speed limit" and inactive outside the
////////     corner; since the driver's speed penalty is symmetric about v_d,
////////     any value comfortably above 70 km/h gives the same answer.
////////  6. The reported energy excludes the idle loss, although the cost
////////     includes it. See the note by ENERGY_EXCLUDES_IDLE below.
////////
//////////////////////////////////////////////////////////////////////////
////////     Copyright (c) Victor M. Becerra, 2026                   /////
//////////////////////////////////////////////////////////////////////////
//////// This is part of the PSOPT software library, which ///////////////
//////// is distributed under the terms of the GNU Lesser ////////////////
//////// General Public License (LGPL)                    ////////////////
//////////////////////////////////////////////////////////////////////////

#include "psopt.h"

#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <string>

using namespace std;
using namespace PSOPT;

//////////////////////////////////////////////////////////////////////////
///////////////////  Vehicle data (paper, Table 1)  //////////////////////
//////////////////////////////////////////////////////////////////////////

static const double MASS      = 1500.0;      // kg
static const double GRAV      = 9.81;        // m/s^2
static const double C_RR      = 0.01;        // rolling resistance coefficient
static const double RHO_CD_A  = 0.86;        // kg/m,  drag force = 0.5*RHO_CD_A*v^2
static const double R_WHEEL   = 0.330;       // m
static const double N_GEAR    = 3.8;         // reduction gear ratio
static const double P_M_MAX   = 80.0e3;      // W,  motor peak power
static const double T_M_MAX   = 280.0;       // Nm, motor maximum torque
static const double P_R_MAX   = 80.0e3;      // W,  generator peak power
static const double T_R_MAX   = 280.0;       // Nm, generator maximum torque
static const double BETA_R    = 0.5;         // regenerative braking bias
static const double L_IDLE    = 500.0;       // W,  auxiliary/idle power

// Driver model (paper, Table 1), with the desired speed of Section 4.1.
static const double ACC_PREF  = 2.5;         // m/s^2, preferred acceleration
static const double BRK_PREF  = 3.0;         // m/s^2, preferred deceleration
static const double V_DES     = 70.0/3.6;    // m/s,   preferred speed (Sec. 4.1)
static const double IDM_DELTA = 4.0;         // IDM acceleration exponent
static const double LAT_MAX   = 3.0;         // m/s^2, comfortable lateral acc.
static const double CURV_MARG = 0.002;       // 1/m,   curvature safety margin

// Powertrain loss polynomial (paper, Table 2). Fm in kN, v in m/s, L_e in kW.
static const double A01 = 0.0207, A11 = 0.0308, A20 = 0.0206;
static const double A30 = 0.00167, A21 = 0.0279;

// Force limits implied by the torque limits.
static const double FP_MAX = N_GEAR*T_M_MAX/R_WHEEL;   // N
static const double FR_MAX = N_GEAR*T_R_MAX/R_WHEEL;   // N
static const double FF_MAX = 12000.0;                  // N, ~0.8 g, never active
static const double V_LIMIT = 90.0/3.6;                // m/s, road speed limit

//////////////////////////////////////////////////////////////////////////
///////////////////  The road  ///////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////

static const double ROAD_LENGTH = 800.0;     // m

static double CORNER_R   = 30.0;             // m,  corner radius
static double CORNER_EPS = 5.0;              // m,  curvature transition length

// Curvature as a function of distance: zero, then 1/R through a 90 degree arc
// placed centrally, then zero again, with the two edges blended over
// CORNER_EPS. The blend is not cosmetic. A step in curvature is a step in the
// speed limit, which is a path constraint with no derivative, and it is also
// not a road anyone could build: a real corner is entered on a transition
// spiral. CORNER_EPS is the length of that transition.
static adouble road_curvature(adouble x)
{
    const double arc = 0.5*PSOPT::pi*CORNER_R;          // 90 degrees of arc
    const double x1  = 0.5*(ROAD_LENGTH - arc);
    const double x2  = 0.5*(ROAD_LENGTH + arc);

    return (1.0/CORNER_R)*( smooth_heaviside(x - x1, CORNER_EPS)
                          - smooth_heaviside(x - x2, CORNER_EPS) );
}

//////////////////////////////////////////////////////////////////////////
///////////////////  The two objectives  /////////////////////////////////
//////////////////////////////////////////////////////////////////////////

// Weights actually used by the cost callback. Held here rather than passed
// because PSOPT's callbacks are plain functions; run_one() sets them.
static double W_DRIVER = 1.0;                // (1-w)
static double W_ENERGY = 0.0;                // w/L0, with L in kW

// The transcription. Hermite-Simpson by default: see the note in setup().
static string COLLOCATION = "Hermite-Simpson";

// A multiplier on the three acceleration penalties of the driver cost. It
// exists to settle a question about the source paper, whose Eq. (28) divides
// the forces by m*u_a and m*u_b, where u_a and u_b are the accelerations
// themselves -- which makes every one of those terms identically one, and so
// no penalty at all. Set to 1 the terms are read as intended, with the
// PREFERRED accelerations a and b in the denominators; set to 0 they are read
// as printed. See the "eq28" mode.
static double ACCEL_PENALTY = 1.0;

// Powertrain losses, paper Eq. (23). Fm in newtons here; converted inside.
static adouble powertrain_loss_kW(adouble Fm_N, adouble v)
{
    adouble Fm = Fm_N/1000.0;                // kN, the units Table 2 is fitted in
    return A01*v + A11*Fm*v + A20*Fm*Fm + A30*Fm*Fm*Fm + A21*Fm*Fm*v;
}

// Vehicle losses, paper Eq. (22): friction braking, rolling resistance, drag.
static adouble vehicle_loss_W(adouble Ff, adouble v)
{
    return (Ff + MASS*GRAV*C_RR + 0.5*RHO_CD_A*v*v)*v;
}

//////////////////////////////////////////////////////////////////////////
///////////////////  Define the end point (Mayer) cost function //////////
//////////////////////////////////////////////////////////////////////////

adouble endpoint_cost(adouble* initial_states, adouble* final_states,
                      adouble* parameters, adouble& t0, adouble& tf,
                      adouble* xad, int iphase, Workspace* workspace)
{
    return 0.0;                              // everything is in the integrand
}

//////////////////////////////////////////////////////////////////////////
///////////////////  Define the integrand (Lagrange) cost function  //////
//////////////////////////////////////////////////////////////////////////

adouble integrand_cost(adouble* states, adouble* controls, adouble* parameters,
                       adouble& time, adouble* xad, int iphase,
                       Workspace* workspace)
{
    adouble v  = states[1];
    adouble Fp = controls[0];
    adouble Fr = controls[1];
    adouble Ff = controls[2];

    // What the driver wants, paper Eq. (28). The headway term is absent: this
    // scenario has no other traffic. The denominators are the PREFERRED
    // accelerations -- see note 2 in the header.
    adouble r  = v/V_DES - 1.0;
    adouble Ld = IDM_DELTA*IDM_DELTA*r*r
               + ACCEL_PENALTY*( (Fp/(MASS*ACC_PREF))*(Fp/(MASS*ACC_PREF))
                               + (Fr/(MASS*BRK_PREF))*(Fr/(MASS*BRK_PREF))
                               + (Ff/(MASS*BRK_PREF))*(Ff/(MASS*BRK_PREF)) );

    // What it costs, in kW.
    adouble Fm = Fp - Fr;
    adouble L  = powertrain_loss_kW(Fm, v)
               + (vehicle_loss_W(Ff, v) + L_IDLE)/1000.0;

    return W_DRIVER*Ld + W_ENERGY*L;
}

//////////////////////////////////////////////////////////////////////////
///////////////////  Define the DAE's ////////////////////////////////////
//////////////////////////////////////////////////////////////////////////

void dae(adouble* derivatives, adouble* path, adouble* states,
         adouble* controls, adouble* parameters, adouble& time,
         adouble* xad, int iphase, Workspace* workspace)
{
    adouble x  = states[0];
    adouble v  = states[1];
    adouble Fp = controls[0];
    adouble Fr = controls[1];
    adouble Ff = controls[2];

    adouble Fm = Fp - Fr;                              // net motor force
    adouble Le = powertrain_loss_kW(Fm, v);            // kW
    adouble Lv = vehicle_loss_W(Ff, v);                // W

    derivatives[0] = v;
    derivatives[1] = ( Fp - Fr - Ff
                     - 0.5*RHO_CD_A*v*v - MASS*GRAV*C_RR )/MASS;

    // The same energy along two different paths, in kWh. E1 counts what the
    // battery supplies; E2 counts what is dissipated. They differ by the change
    // in kinetic energy, which the boundary conditions force to zero, so their
    // disagreement at t = T measures the transcription and nothing else.
    derivatives[2] = (Fm*v + 1000.0*Le)/3.6e6;         // E1
    derivatives[3] = (Lv   + 1000.0*Le)/3.6e6;         // E2

    // Cornering, paper Eq. (17): the lateral acceleration a driver will accept
    // falls as the speed rises, because the estimate of the curvature gets
    // worse. Written in this form rather than as v <= sqrt(Gamma/(kappa+Delta))
    // it is a polynomial in v, which is better conditioned and has no square
    // root to differentiate near zero.
    path[0] = (road_curvature(x) + CURV_MARG)*v*v;     // <= LAT_MAX

    // The motor cannot deliver or absorb more than its rated power.
    path[1] = Fm*v;                                    // in [-P_R_MAX, P_M_MAX]

    // The regenerative braking bias, paper Eq. (27), relaxed from the equality
    // (26). With beta_r = 1/2 this is simply Fr <= Ff: the front axle may not
    // take more than half the braking, or the car becomes directionally
    // unstable under braking. Whether it is active at the optimum is checked
    // after the solve.
    path[2] = (1.0 - BETA_R)*Fr - BETA_R*Ff;           // <= 0
}

//////////////////////////////////////////////////////////////////////////
///////////////////  Define the events function //////////////////////////
//////////////////////////////////////////////////////////////////////////

void events(adouble* e, adouble* initial_states, adouble* final_states,
            adouble* parameters, adouble& t0, adouble& tf, adouble* xad,
            int iphase, Workspace* workspace)
{
    e[0] = initial_states[0];                          // x(0) = 0
    e[1] = final_states[0];                            // x(T) = l
    e[2] = final_states[1] - initial_states[1];        // v(T) = v(0), both free
    e[3] = initial_states[2];                          // E1(0) = 0
    e[4] = initial_states[3];                          // E2(0) = 0
}

//////////////////////////////////////////////////////////////////////////
///////////////////  Define the phase linkages function //////////////////
//////////////////////////////////////////////////////////////////////////

void linkages(adouble* linkages, adouble* xad, Workspace* workspace)
{
    // Single phase: no linkages
}

//////////////////////////////////////////////////////////////////////////
///////////////////  One run  ////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////

// ENERGY_EXCLUDES_IDLE. The paper's Table 3 reports an energy that leaves the
// idle loss out, although the cost function it minimises puts the idle loss in.
// The paper does not say so; it follows from its own Table 1. Discard the whole
// powertrain model -- take the drivetrain to be lossless -- and keep only
// rolling resistance and drag, and the energy per km of ANY trajectory of mean
// speed vbar over a level road with v(0) = v(T) is at least
//
//     m g c_rr + (1/2) rho Cd A vbar^2 + L_i/vbar,
//
// which at the reported 42.0 km/h is 0.06904 kWh/km. Table 3 reports 0.068. A
// lossless car cannot reach that with the idle loss counted, so it is not
// counted. Reported here the same way, for comparability.
struct Run
{
    bool   ok;
    double w, L0;
    double T;              // s,     journey time
    double vmean;          // km/h,  l/T
    double v0;             // km/h,  the free periodic speed
    double energy;         // kWh/km, idle excluded, as the paper reports it
    double energy_idle;    // kWh/km, idle included, as the cost counts it
    double balance;        // kWh,   |E1(T) - E2(T)|, should be zero
    double bias_gap;       // N,     max slack in the relaxed bias constraint
    double regen_frac;     // fraction of braking taken regeneratively
    double max_lat;        // m/s^2, largest lateral acceleration reached
};

static void setup(Prob& problem, Alg& algorithm, int nnodes)
{
    problem.name        = "Eco-driving through a corner";
    problem.outfilename = "ecodriving.txt";

    problem.nphases     = 1;
    problem.nlinkages   = 0;
    psopt_level1_setup(problem);

    problem.phases(1).nstates   = 4;
    problem.phases(1).ncontrols = 3;
    problem.phases(1).nevents   = 5;
    problem.phases(1).npath     = 3;
    problem.phases(1).nodes     = (RowVectorXi(1) << nnodes).finished();
    psopt_level2_setup(problem, algorithm);

    problem.phases(1).bounds.lower.states << 0.0, 1.0, -5.0, -5.0;
    problem.phases(1).bounds.upper.states << ROAD_LENGTH, V_LIMIT, 5.0, 5.0;

    problem.phases(1).bounds.lower.controls << 0.0, 0.0, 0.0;
    problem.phases(1).bounds.upper.controls << FP_MAX, FR_MAX, FF_MAX;

    problem.phases(1).bounds.lower.events << 0.0, ROAD_LENGTH, 0.0, 0.0, 0.0;
    problem.phases(1).bounds.upper.events << 0.0, ROAD_LENGTH, 0.0, 0.0, 0.0;

    problem.phases(1).bounds.lower.path << 0.0,      -P_R_MAX, -1.0e6;
    problem.phases(1).bounds.upper.path << LAT_MAX,   P_M_MAX,  0.0;

    problem.phases(1).bounds.lower.StartTime = 0.0;
    problem.phases(1).bounds.upper.StartTime = 0.0;
    problem.phases(1).bounds.lower.EndTime   = 20.0;
    problem.phases(1).bounds.upper.EndTime   = 200.0;

    problem.integrand_cost = &integrand_cost;
    problem.endpoint_cost  = &endpoint_cost;
    problem.dae            = &dae;
    problem.events         = &events;
    problem.linkages       = &linkages;

    // A cold guess: cruise the whole road at the preferred speed. It violates
    // the cornering constraint, which is the point -- the solver has to find
    // the slowing-down for itself.
    const int ng = 20;
    const double Tg = ROAD_LENGTH/V_DES;
    MatrixXd xg(4, ng);
    xg << linspace(0.0, ROAD_LENGTH, ng),
          V_DES*ones(1, ng),
          linspace(0.0, 0.1, ng),
          linspace(0.0, 0.1, ng);
    MatrixXd ug(3, ng);
    ug << (MASS*GRAV*C_RR + 0.5*RHO_CD_A*V_DES*V_DES)*ones(1, ng),
          zeros(1, ng),
          zeros(1, ng);
    problem.phases(1).guess.states   = xg;
    problem.phases(1).guess.controls = ug;
    problem.phases(1).guess.time     = linspace(0.0, Tg, ng);

    algorithm.nlp_method         = "IPOPT";
    algorithm.scaling            = "automatic";
    algorithm.derivatives        = "automatic";
    algorithm.collocation_method = COLLOCATION;
    algorithm.nlp_iter_max       = 3000;
    algorithm.nlp_tolerance      = 1.e-8;
    algorithm.mesh_refinement    = "automatic";
    algorithm.ode_tolerance      = 1.e-6;
    algorithm.mr_max_iterations  = 8;

    // WHY LOCAL COLLOCATION. The solution has two corners in it: the speed is
    // flat, falls hard to the cornering limit, holds, then rises hard. A global
    // pseudospectral scheme approximates that with one polynomial over the whole
    // journey, and rings on the straights either side of the bend -- about one
    // per cent in speed at 60 nodes, which is the size of the differences this
    // study is about. Worse, refining it makes matters worse rather than
    // better, because raising the degree of a single polynomial through a
    // corner is the classic way to provoke a Gibbs oscillation: run this
    // example with "Legendre" and watch the discretisation error rise from the
    // fourth mesh iteration onwards. Hermite-Simpson refines by subdividing
    // instead, which is what a solution with corners in it needs.
}

// Read everything the study wants out of a finished solve.
static Run harvest(Sol& solution, double w, double L0)
{
    Run r;
    r.ok = true; r.w = w; r.L0 = L0;

    MatrixXd t = solution.get_time_in_phase(1);
    MatrixXd x = solution.get_states_in_phase(1);
    MatrixXd u = solution.get_controls_in_phase(1);

    const long   nc = t.cols();
    const double T  = t(0, nc-1);
    const double E1 = x(2, nc-1);

    r.T       = T;
    r.vmean   = 3.6*ROAD_LENGTH/T;
    r.v0      = 3.6*x(1, 0);
    r.energy  = E1/(ROAD_LENGTH/1000.0);
    r.energy_idle = r.energy + (L_IDLE*T/3.6e6)/(ROAD_LENGTH/1000.0);

    // The energy identity, checked ALONG the trajectory rather than at its end.
    // E1 - E2 should equal the change in kinetic energy at every instant. At
    // t = T it does so trivially, because the periodicity condition makes both
    // sides zero and the transcription enforces that exactly; the endpoint is
    // therefore no test at all. The largest discrepancy over the mesh is one,
    // and it is a test of the discretisation, not of the model: it needs no
    // reference solution, no closed form and no second code.
    const double v0 = x(1, 0);
    r.balance = 0.0;
    for (long k = 0; k < nc; k++) {
        const double dke = 0.5*MASS*(x(1,k)*x(1,k) - v0*v0)/3.6e6;   // kWh
        const double err = fabs((x(2,k) - x(3,k)) - dke);
        if (err > r.balance) r.balance = err;
    }

    // Was the relaxed bias constraint active whenever the brakes were on? The
    // relaxation is legitimate only where it was -- and the paper's own
    // derivation says it will be slack exactly where the motor saturates in
    // torque or in power, since then the friction brakes must take up the rest.
    // So the slack is measured only where the motor is NOT saturated; anywhere
    // else a gap is the model working as intended, not the relaxation failing.
    r.bias_gap   = 0.0;
    r.max_lat    = 0.0;
    double wreg = 0.0, wtot = 0.0;
    for (long k = 0; k < nc && k < u.cols(); k++) {
        const double Fp = u(0, k), Fr = u(1, k), Ff = u(2, k);
        const double v = x(1, k), xx = x(0, k);
        const double braking = Fr + Ff;
        const bool saturated = (Fr > 0.999*FR_MAX)
                            || (fabs((Fp - Fr)*v) > 0.999*P_R_MAX);
        if (braking > 1.0 && !saturated) {
            const double gap = BETA_R*Ff - (1.0 - BETA_R)*Fr;
            if (gap > r.bias_gap) r.bias_gap = gap;
        }
        wreg += Fr*v; wtot += braking*v;          // braking POWER, not force

        // Recompute the curvature outside the AD tape, for reporting.
        const double arc = 0.5*M_PI*CORNER_R;
        const double x1  = 0.5*(ROAD_LENGTH - arc), x2 = 0.5*(ROAD_LENGTH + arc);
        const double kap = (1.0/CORNER_R)*( 0.5*(1.0 + tanh((xx-x1)/CORNER_EPS))
                                          - 0.5*(1.0 + tanh((xx-x2)/CORNER_EPS)) );
        const double lat = (kap + CURV_MARG)*v*v;
        if (lat > r.max_lat) r.max_lat = lat;
    }
    r.regen_frac = (wtot > 0.0) ? wreg/wtot : 0.0;
    return r;
}

// Solve once. If continue_from is non-null the guess is taken from it, which is
// what makes a scan of w one calculation instead of many.
static Run run_one(Prob& problem, Alg& algorithm, Sol& solution,
                   double w, double L0, const string& tag, bool write_files)
{
    W_DRIVER = 1.0 - w;
    W_ENERGY = w/L0;

    problem.outfilename = "ecodriving" + tag + ".txt";

    Run r; r.ok = false; r.w = w; r.L0 = L0;
    if (psopt(solution, problem, algorithm) != 0 || solution.error_flag) {
        printf("\n ecodriving: the solve did not succeed at w = %g: %s\n",
               w, solution.error_msg.c_str());
        return r;
    }

    r = harvest(solution, w, L0);

    if (write_files) {
        MatrixXd t = solution.get_time_in_phase(1);
        MatrixXd x = solution.get_states_in_phase(1);
        MatrixXd u = solution.get_controls_in_phase(1);
        const string dat = "ecodriving" + tag + ".dat";
        FILE* fp = fopen(dat.c_str(), "w");
        fprintf(fp, "# eco-driving, w = %g, L0 = %g kW, R = %g m, eps = %g m\n",
                w, L0, CORNER_R, CORNER_EPS);
        fprintf(fp, "# t  x  v  E1  E2  Fp  Fr  Ff  kappa  lateral_acc\n");
        for (long k = 0; k < t.cols(); k++) {
            const double xx = x(0,k), v = x(1,k);
            const double arc = 0.5*M_PI*CORNER_R;
            const double x1 = 0.5*(ROAD_LENGTH-arc), x2 = 0.5*(ROAD_LENGTH+arc);
            const double kap = (1.0/CORNER_R)*( 0.5*(1.0+tanh((xx-x1)/CORNER_EPS))
                                              - 0.5*(1.0+tanh((xx-x2)/CORNER_EPS)) );
            fprintf(fp, "%.8e  %.8e  %.8e  %.8e  %.8e", t(0,k), xx, v, x(2,k), x(3,k));
            for (int j = 0; j < 3; j++)
                fprintf(fp, "  %.8e", (k < u.cols()) ? u(j,k) : 0.0);
            fprintf(fp, "  %.8e  %.8e\n", kap, (kap+CURV_MARG)*v*v);
        }
        fclose(fp);
        printf("\n  wrote %s\n", dat.c_str());
    }
    return r;
}

//////////////////////////////////////////////////////////////////////////
///////////////////  Reporting  //////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////

// The floor of the previous note, evaluated here so that the run can be
// compared against something that does not depend on the powertrain model.
static double energy_floor(double vmean_kmh, bool with_idle)
{
    const double v = vmean_kmh/3.6;
    double per_m = MASS*GRAV*C_RR + 0.5*RHO_CD_A*v*v;
    if (with_idle) per_m += L_IDLE/v;
    return per_m*1000.0/3.6e6;                 // J/m -> kWh/km
}

// The energy per km of driving the whole road at one constant speed, from the
// published loss model. This needs no optimal control at all -- it is a
// one-dimensional function that can be minimised by hand -- and on a straight
// road it is what the w = 1 answer must be. It is therefore the sharpest check
// available on the whole implementation: if the optimal control solution does
// not reproduce it, one of them is wrong.
static double constant_speed_energy(double v, bool with_idle)   // v in m/s
{
    const double Fm_N = MASS*GRAV*C_RR + 0.5*RHO_CD_A*v*v;      // steady state
    const double Fm   = Fm_N/1000.0;                            // kN
    const double Le_W = 1000.0*(A01*v + A11*Fm*v + A20*Fm*Fm
                              + A30*Fm*Fm*Fm + A21*Fm*Fm*v);
    double per_m = (Fm_N*v + Le_W + (with_idle ? L_IDLE : 0.0))/v;
    return per_m*1000.0/3.6e6;                 // J/m -> kWh/km
}

static double best_constant_speed(void)        // m/s, by golden-section search
{
    double lo = 2.0, hi = 25.0;
    for (int i = 0; i < 300; i++) {
        const double m1 = lo + (hi-lo)/3.0, m2 = hi - (hi-lo)/3.0;
        if (constant_speed_energy(m1, true) < constant_speed_energy(m2, true))
            hi = m2; else lo = m1;
    }
    return 0.5*(lo+hi);
}

static void report_one(const Run& r)
{
    printf("\n=====================================================================\n");
    printf("  Eco-driving through a corner: w = %g, L0 = %g kW\n", r.w, r.L0);
    printf("=====================================================================\n");
    printf("  Journey\n");
    printf("    time             %8.3f s\n", r.T);
    printf("    mean speed       %8.2f km/h\n", r.vmean);
    printf("    periodic speed   %8.2f km/h   (free: v(0) = v(T))\n", r.v0);
    printf("    peak lateral acc %8.3f m/s^2  (limit %.1f)\n", r.max_lat, LAT_MAX);

    printf("\n  Energy\n");
    printf("    idle excluded    %8.4f kWh/km   (as the paper reports it)\n", r.energy);
    printf("    idle included    %8.4f kWh/km   (as the cost counts it)\n", r.energy_idle);
    printf("    lossless floor   %8.4f kWh/km   at this mean speed, idle excluded\n",
           energy_floor(r.vmean, false));
    printf("    lossless floor   %8.4f kWh/km   at this mean speed, idle included\n",
           energy_floor(r.vmean, true));
    printf("    The floors use nothing but the mass, the rolling and drag\n");
    printf("    coefficients and the idle power. No trajectory can beat them.\n");

    const double vb = best_constant_speed();
    printf("\n  The one-dimensional problem underneath this one\n");
    printf("    best constant speed  %8.2f km/h at %.4f kWh/km including idle\n",
           3.6*vb, constant_speed_energy(vb, true));
    printf("    corner speed limit   %8.2f km/h\n",
           3.6*sqrt(LAT_MAX/(1.0/CORNER_R + CURV_MARG)));
    printf("    Where the corner limit is above the best constant speed, the\n");
    printf("    cornering constraint never binds and the minimum-energy answer\n");
    printf("    IS that constant speed. Minimising a function of one variable by\n");
    printf("    hand and solving the optimal control problem must then agree,\n");
    printf("    and it is worth checking that they do before believing either.\n");

    printf("\n  Checks that need no reference solution\n");
    printf("    energy identity  %8.2e kWh    max over the mesh of\n", r.balance);
    printf("                                    |(E1 - E2) - dKE|, the same\n");
    printf("                                    energy counted two ways\n");
    printf("    bias slack       %8.2e N      largest slack in Fr <= Ff where\n", r.bias_gap);
    printf("                                    the motor is NOT saturated; the\n");
    printf("                                    relaxation of the equality is\n");
    printf("                                    legitimate only if this is small\n");
    printf("    regen fraction   %8.3f        of the braking ENERGY recovered\n", r.regen_frac);
    printf("=====================================================================\n\n");
}

//////////////////////////////////////////////////////////////////////////
///////////////////  The Pareto frontier  ////////////////////////////////
//////////////////////////////////////////////////////////////////////////

// A tag that identifies a run by its weight, so that no run of the scan can
// write over the file of another, or over the single run's.
//
// The weight is written as three digits rather than as a decimal -- _w050 and
// not _w0.5 -- because PSOPT builds the name of its mesh statistics file by
// cutting the output file name at the first dot, so "ecodriving_w0.5.txt" and
// "ecodriving_w0.txt" produce the same mesh statistics file and the later run
// silently replaces the earlier one. That is the same hazard the naming scheme
// exists to prevent, arriving through a different door: it is not enough for a
// name to be distinct, it must still be distinct after everything downstream
// has finished trimming it.
static string weight_tag(double w)
{
    char buf[32];
    snprintf(buf, sizeof buf, "_w%03d", (int) lround(w*100.0));
    return string(buf);
}

static int run_pareto(int nnodes, double L0)
{
    static const double ws[] = { 0.0, 0.05, 0.1, 0.2, 0.3, 0.4, 0.5,
                                 0.6, 0.7, 0.8, 0.9, 0.95, 1.0 };
    const int nw = (int)(sizeof ws / sizeof ws[0]);

    Alg algorithm; Prob problem; Sol solution;
    setup(problem, algorithm, nnodes);
    algorithm.print_level = 0;

    Run runs[32], cold[32];
    bool complete = true;

    // Traced by continuation, from the naturalistic end upwards. These points
    // are meant to lie on one curve, so each solve starts from the one before
    // it; started cold, a few of them would settle on a different branch and
    // appear on the plot as structure that is not there.
    for (int k = 0; k < nw; k++) {
        runs[k] = run_one(problem, algorithm, solution, ws[k], L0,
                          weight_tag(ws[k]), true);
        if (!runs[k].ok) { complete = false; continue; }
        set_guess_from_solution(problem, solution);
    }

    // The same frontier again, every point from the same cold guess, so that
    // the difference the continuation makes can be seen rather than asserted.
    // This second sweep is what Section 8.6.10 of the companion book is about;
    // it is not how the frontier should be computed.
    for (int k = 0; k < nw; k++) {
        Alg alg2; Prob prob2; Sol sol2;
        setup(prob2, alg2, nnodes);
        alg2.print_level = 0;
        cold[k] = run_one(prob2, alg2, sol2, ws[k], L0, "_cold", false);
    }

    // The frontier itself, in a form a plotting script can read.
    FILE* fd = fopen("ecodriving_pareto.dat", "w");
    if (fd) {
        fprintf(fd, "# w  lambda  mean_speed_kmh  energy_kWh_per_km  "
                    "floor_kWh_per_km  time_s  cold_mean_kmh  cold_energy\n");
        for (int k = 0; k < nw; k++) {
            if (!runs[k].ok) continue;
            const double lam = (ws[k] < 1.0) ? ws[k]/((1.0-ws[k])*L0) : -1.0;
            fprintf(fd, "%g  %.6e  %.6e  %.6e  %.6e  %.6e  %.6e  %.6e\n",
                    ws[k], lam, runs[k].vmean, runs[k].energy,
                    energy_floor(runs[k].vmean, false), runs[k].T,
                    cold[k].ok ? cold[k].vmean  : 0.0,
                    cold[k].ok ? cold[k].energy : 0.0);
        }
        fclose(fd);
    }

    FILE* fp = fopen("ecodriving_pareto.txt", "w");
    for (int pass = 0; pass < 2; pass++) {
        FILE* out = (pass == 0) ? stdout : fp;
        if (out == NULL) continue;
        fprintf(out, "=====================================================================\n");
        fprintf(out, "  Pareto frontier: energy against mean speed, L0 = %g kW\n", L0);
        fprintf(out, "  %d nodes. Traced by continuation from w = 0 upwards.\n", nnodes);
        fprintf(out, "=====================================================================\n");
        fprintf(out, "  %6s %8s %10s %10s %10s %10s | %10s %10s\n",
                "w", "lambda", "mean km/h", "kWh/km", "floor", "|E1-E2|",
                "cold km/h", "cold kWh/km");
        int ndiff = 0;
        for (int k = 0; k < nw; k++) {
            if (!runs[k].ok) { fprintf(out, "  %6g   did not solve\n", ws[k]); continue; }
            const double lam = (ws[k] < 1.0) ? ws[k]/((1.0-ws[k])*L0) : INFINITY;
            const bool moved = cold[k].ok
                             && fabs(cold[k].energy - runs[k].energy) > 0.005*runs[k].energy;
            if (moved) ndiff++;
            fprintf(out, "  %6g %8.4f %10.2f %10.4f %10.4f %10.2e | %10.2f %10.4f%s\n",
                    ws[k], lam, runs[k].vmean, runs[k].energy,
                    energy_floor(runs[k].vmean, false), runs[k].balance,
                    cold[k].ok ? cold[k].vmean : 0.0,
                    cold[k].ok ? cold[k].energy : 0.0,
                    cold[k].ok ? (moved ? "  <-- differs" : "") : "  <-- failed");
        }
        fprintf(out, "\n  Points where the cold start landed somewhere else: %d of %d.\n",
                ndiff, nw);
        if (ndiff == 0) {
            fprintf(out, "  None, on this problem: the guess is a good one and the frontier\n");
            fprintf(out, "  is a single well-behaved branch, so continuation buys nothing\n");
            fprintf(out, "  here. That is worth knowing and is not an argument for skipping\n");
            fprintf(out, "  the check. A scan that has drifted onto another branch reports\n");
            fprintf(out, "  success at every point, and the damage appears only as a bump in\n");
            fprintf(out, "  a plot that was the object of the exercise. The second column\n");
            fprintf(out, "  costs one extra sweep and is the only thing that distinguishes\n");
            fprintf(out, "  'the scan is clean' from 'nobody looked'.\n");
        } else {
            fprintf(out, "  A cold-started scan reports success at every point either way.\n");
            fprintf(out, "  What distinguishes the two columns is not an exit status but a\n");
            fprintf(out, "  plot, and on a plot a point that has moved looks like structure.\n");
        }
        fprintf(out, "\n  The published endpoints, Lot et al. (2025) Table 3, cornering:\n");
        fprintf(out, "    w = 0   0.153 kWh/km at 62.1 km/h\n");
        fprintf(out, "    w = 1   0.068 kWh/km at 42.0 km/h\n");
        fprintf(out, "  L0 is not given in that paper and does not change this frontier;\n");
        fprintf(out, "  it only decides which w labels which point of it. Compare the\n");
        fprintf(out, "  curve and its endpoints, not the rows.\n");
        fprintf(out, "=====================================================================\n");
    }
    if (fp) { fclose(fp); printf("\n  wrote ecodriving_pareto.txt\n\n"); }

    return complete ? 0 : 1;
}

//////////////////////////////////////////////////////////////////////////
///////////////////  What the corner assumption is worth  ////////////////
//////////////////////////////////////////////////////////////////////////

// The paper gives the corner as "approximately 30 m" and says nothing about how
// the curvature is entered. Both are ours. This measures what they are worth,
// so that the comparison against the published numbers can be read with the
// right number of significant figures.
static int run_corner(int nnodes, double L0)
{
    static const double radii[] = { 28.0, 30.0, 32.0 };
    static const double epss[]  = { 2.0, 5.0, 10.0, 20.0 };
    const int nr = 3, ne = 4;
    const double R_keep = CORNER_R, E_keep = CORNER_EPS;

    printf("\n=====================================================================\n");
    printf("  What the corner assumptions are worth (w = 0 and w = 1)\n");
    printf("  The paper gives 'approximately 30 m' and no transition at all.\n");
    printf("=====================================================================\n");
    printf("  %6s %6s | %10s %10s | %10s %10s\n",
           "R [m]", "eps", "w=0 km/h", "w=0 kWh/km", "w=1 km/h", "w=1 kWh/km");

    for (int i = 0; i < nr; i++) {
        for (int j = 0; j < ne; j++) {
            CORNER_R = radii[i]; CORNER_EPS = epss[j];

            Alg algorithm; Prob problem; Sol solution;
            setup(problem, algorithm, nnodes);
            algorithm.print_level = 0;

            char tag[64];
            snprintf(tag, sizeof tag, "_R%g_e%g", CORNER_R, CORNER_EPS);
            Run r0 = run_one(problem, algorithm, solution, 0.0, L0, tag, false);
            if (r0.ok) set_guess_from_solution(problem, solution);
            Run r1 = run_one(problem, algorithm, solution, 1.0, L0, tag, false);

            printf("  %6.1f %6.1f | %10.2f %10.4f | %10.2f %10.4f\n",
                   CORNER_R, CORNER_EPS,
                   r0.ok ? r0.vmean : 0.0, r0.ok ? r0.energy : 0.0,
                   r1.ok ? r1.vmean : 0.0, r1.ok ? r1.energy : 0.0);
            fflush(stdout);
        }
    }
    CORNER_R = R_keep; CORNER_EPS = E_keep;
    printf("=====================================================================\n\n");
    return 0;
}

//////////////////////////////////////////////////////////////////////////
///////////////////  Which reading of Eq. (28)?  /////////////////////////
//////////////////////////////////////////////////////////////////////////

// The naturalistic end of the frontier is set entirely by the driver cost, so
// it is the place where a misreading of that cost shows up. Solve w = 0 both
// ways and compare against the published point.
static int run_eq28(int nnodes)
{
    printf("\n=====================================================================\n");
    printf("  Two readings of the source paper's Eq. (28), at w = 0\n");
    printf("=====================================================================\n");
    printf("  %-34s %10s %10s\n", "reading", "km/h", "kWh/km");

    const double keep = ACCEL_PENALTY;
    for (int k = 0; k < 2; k++) {
        ACCEL_PENALTY = (k == 0) ? 1.0 : 0.0;

        Alg algorithm; Prob problem; Sol solution;
        setup(problem, algorithm, nnodes);
        algorithm.print_level = 0;

        const string tag = (k == 0) ? "_eq28_intended" : "_eq28_asprinted";
        Run r = run_one(problem, algorithm, solution, 0.0, 6.0, tag, true);
        printf("  %-34s %10.2f %10.4f  regen %.3f\n",
               (k == 0) ? "denominators m*a, m*b (intended)"
                        : "denominators m*u_a, m*u_b (printed)",
               r.ok ? r.vmean : 0.0, r.ok ? r.energy : 0.0,
               r.ok ? r.regen_frac : 0.0);
        fflush(stdout);
    }
    ACCEL_PENALTY = keep;

    printf("  %-34s %10.2f %10.4f\n", "Lot et al. (2025), Table 3", 62.1, 0.153);
    printf("=====================================================================\n\n");
    return 0;
}

//////////////////////////////////////////////////////////////////////////
///////////////////  Define the main routine /////////////////////////////
//////////////////////////////////////////////////////////////////////////

int main(int argc, char* argv[])
{
    string mode = (argc > 1) ? argv[1] : "0";

    if (mode == "pareto") {
        int nnodes = (argc > 2) ? atoi(argv[2]) : 60;
        return run_pareto(nnodes, 6.0);
    }
    if (mode == "corner") {
        int nnodes = (argc > 2) ? atoi(argv[2]) : 60;
        return run_corner(nnodes, 6.0);
    }
    if (mode == "eq28") {
        int nnodes = (argc > 2) ? atoi(argv[2]) : 60;
        return run_eq28(nnodes);
    }

    double w      = atof(mode.c_str());
    double L0     = (argc > 2) ? atof(argv[2]) : 6.0;
    int    nnodes = (argc > 3) ? atoi(argv[3]) : 60;
    if (argc > 4) COLLOCATION = argv[4];

    if (w < 0.0 || w > 1.0 || L0 <= 0.0 || nnodes < 8) {
        printf("usage: %s [w in 0..1] [L0 kW > 0] [nodes >= 8]\n"
               "       %s pareto [nodes]\n"
               "       %s corner [nodes]\n", argv[0], argv[0], argv[0]);
        return 1;
    }

    Alg algorithm; Prob problem; Sol solution;
    setup(problem, algorithm, nnodes);

    const Run r = run_one(problem, algorithm, solution, w, L0, weight_tag(w), true);
    if (!r.ok) return 1;
    report_one(r);

    MatrixXd t = solution.get_time_in_phase(1);
    MatrixXd x = solution.get_states_in_phase(1);
    MatrixXd u = solution.get_controls_in_phase(1);
    MatrixXd v = x.row(1), pos = x.row(0);

    plot(t, v, problem.name + ": speed", "time (s)", "v (m/s)", "v");
    plot(t, u, problem.name + ": forces", "time (s)", "force (N)", "Fp Fr Ff");
    plot(pos, v, problem.name + ": speed against distance",
         "distance (m)", "v (m/s)", "v");
    plot(pos, v, problem.name + ": speed against distance",
         "distance (m)", "v (m/s)", "v", "pdf",
         ("ecodriving" + weight_tag(w) + "_speed.pdf").c_str());

    return 0;
}

//////////////////////////////////////////////////////////////////////////
///////////////////////      END OF FILE     /////////////////////////////
//////////////////////////////////////////////////////////////////////////
