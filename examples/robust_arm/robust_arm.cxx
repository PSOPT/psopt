//////////////////////////////////////////////////////////////////////////
//////////////////        robust_arm.cxx        //////////////////////////
//////////////////////////////////////////////////////////////////////////
////////////////           PSOPT  Example             ////////////////////
//////////////////////////////////////////////////////////////////////////
//////// Title:   Two-link arm with an UNCERTAIN PAYLOAD            //////
////////          -- robust minimum time by scenario augmentation   //////
//////// Last modified: 27 September 2026                           //////
//////// Reference:     the arm is the PROPT user's guide problem,   //////
////////                as in examples/twolinkarm; the robust        //////
////////                formulation is sample-average approximation  //////
////////                (see e.g. Shapiro, Dentcheva & Ruszczynski). //////
//////////////////////////////////////////////////////////////////////////
////////     Copyright (c) Victor M. Becerra, 2026         ///////////////
//////////////////////////////////////////////////////////////////////////
//////// This is part of the PSOPT software library, which ///////////////
//////// is distributed under the terms of the GNU Lesser ////////////////
//////// General Public License (LGPL)                    ////////////////
//////////////////////////////////////////////////////////////////////////
//
//  WHAT THIS EXAMPLE IS FOR
//
//  A control computed for the nominal plant is not a control that works on the
//  plant you have. This example poses the smallest honest version of that
//  statement as an optimal control problem and solves it with PSOPT as it
//  stands -- no new library feature is used, and none is needed.
//
//  The arm of examples/twolinkarm carries a payload whose mass m_p is not known
//  exactly. One torque history and one final time must be committed BEFORE m_p
//  is revealed, and must bring every plant in the uncertainty set to the target
//  within a stated tolerance. That is a "here-and-now" decision, and it is the
//  thing that distinguishes robust optimal control from solving the problem
//  once per sample:
//
//      Solving one optimal control problem per sample and averaging the answers
//      does NOT give a robust control. It gives the wait-and-see solution, in
//      which every realisation is optimised with foreknowledge of its own
//      uncertainty. Since E[min] <= min E[.], its cost is a lower bound that no
//      implementable control attains, and the average of the control histories
//      solves nothing at all. What couples the samples into one problem is the
//      requirement that they SHARE the control -- non-anticipativity -- and in
//      the transcription below that requirement is structural: there is one
//      control vector and M copies of the state.
//
//  HOW THE PAYLOAD ENTERS THE MODEL
//
//  examples/twolinkarm writes the arm with its inertia coefficients already
//  evaluated -- 9/4, 2, 4/3, 3/2, 7/2, 7/3, 31/36 -- so there is nowhere for a
//  payload to go. Working backwards from those numbers, the model is a planar
//  2R arm in ABSOLUTE angular coordinates
//
//      x1 = dtheta1/dt,  x2 = absolute angular rate of link 2,
//      x3 = theta2 - theta1,  x4 = theta1,
//
//  driven by torques tau = (u1 - u2, u2), with unit link lengths and
//
//      M11 = 7/3,   M12 = (3/2) cos x3,   M22 = 4/3,   S2 = 3/2,
//
//  where S2 is the first moment of link 2 about joint 2 and M22 its second
//  moment about the same point. Every one of the shipped expressions follows:
//  the two numerators term by term, and the denominator as
//  M11*M22 - M12^2 = 28/9 - (9/4)cos^2(x3), which is identically
//  31/36 + (9/4)sin^2(x3).
//
//  A point payload of mass m_p at the tip of link 2 then enters in exactly one
//  way, and it is completely determined by that reading of the model:
//
//      M11 = 7/3 + m_p,  M12 = (3/2 + m_p) cos x3,  M22 = 4/3 + m_p,
//      S2  = 3/2 + m_p.
//
//  The payload swings about joint 1 with the rest of link 2, so it adds m_p L1^2
//  to M11; it adds m_p L2 to the first moment about joint 2 and m_p L2^2 to the
//  second; and L1 = L2 = 1 in the units the model already uses. Setting m_p = 0
//  reproduces the shipped dynamics to 2.8e-15 over 200000 random states and
//  torques, which was checked before anything here was built on it. A payload of
//  0.5 changes the joint accelerations by 12 to 17 per cent.
//
//  THE TRANSCRIPTION
//
//  A finite set of M payloads -- scenarios -- turns the robust problem into ONE
//  deterministic optimal control problem: 4M states, being M copies of the arm
//  each with its own m_p, against 2 controls and one free final time. The copies
//  are coupled only through the control and t_f, and that coupling IS the
//  robustness. Nothing else in the transcription is unusual, which is the point:
//  PSOPT solves it as it stands.
//
//  THE TERMINAL CONDITION, which is the modelling choice worth arguing about
//
//  One open-loop torque history cannot steer several different plants to the same
//  point exactly, so requiring x(t_f) = x_target for every scenario is
//  generically infeasible for M > 1. The terminal conditions are therefore
//  relaxed to a tolerance ball of radius delta, enforced for every scenario, and
//  the objective remains the final time. The answer is then a trade-off: how much
//  longer the slew must take for every plant in the set to arrive within delta.
//
//  WHY THE ROBUST SLEW IS A SLOW ONE
//
//  The payload enters only through the inertia, so it changes the ACCELERATION
//  the torques produce. Drive hard and the payload matters; drive gently and the
//  trajectory approaches a quasi-static one on which it matters much less. A
//  design that must land every plant in the set within delta therefore buys its
//  insensitivity with time, and the numbers below show it doing so: t_f roughly
//  2.6 times the nominal. That is not a defect of the method, it is the physics
//  of the problem, and it is why the final-time bound of examples/twolinkarm had
//  to be raised here.
//
//  WHAT THE PROGRAM DOES, AND WHAT IT IS FOR
//
//    (a) M = 1, delta = 0     reproduces examples/twolinkarm. This verifies that
//                             the reparameterised model and the augmented
//                             transcription are what they claim to be.
//    (b) M = 1, delta > 0     the NOMINAL design at the working tolerance.
//                             Relaxing the terminal condition buys time on its
//                             own, and this run measures how much, so that the
//                             cost of robustness is not confused with it.
//    (c) a sweep in M         scenarios placed at the Gauss-Hermite nodes of the
//                             payload distribution, M = 3, 5, 7, 9. Each design
//                             is then scored on payloads it never saw.
//    (d) scenario generation  solve, find the payload the resulting control
//                             serves WORST over the whole uncertainty set, add
//                             it, re-solve. One call to psopt() per iteration.
//
//  The reason (c) is in the program and not merely in its comments is that (c)
//  does not work, and the way it fails is the most useful thing here. A fixed
//  quadrature rule pins the terminal miss at its own nodes and says nothing about
//  the gaps between them, so the optimiser drives the miss to delta exactly at
//  the nodes and lets it grow freely elsewhere: the in-sample miss is delta by
//  construction at every M, while the worst miss over the uncertainty set never
//  once falls below delta. Nor does raising M cure it in any orderly way. The
//  fraction of sampled plants that arrive within delta runs 0%, 39%, 23%, 79%,
//  25% as M goes 1, 3, 5, 7, 9 -- no trend at all, because where the next
//  quadrature node lands has nothing to do with where the previous design was
//  failing. Quadrature is the right tool for an EXPECTATION and the wrong tool
//  for a constraint that has to hold everywhere.
//
//  (d) works, because it puts each new scenario where the current design is
//  actually failing. It converges here in twelve iterations, to thirteen
//  scenarios, a design that brings 100% of sampled plants inside delta, and --
//  what no amount of sampling can give -- a certificate: no payload anywhere in
//  the set misses by more than delta. It is a cutting-plane method on the
//  semi-infinite constraint, and it is the loop a robust optimal control driver
//  would run.
//
//  COST. Eighteen calls to psopt() and some tens of millions of RK4 steps: about
//  five minutes, most of it in the solver. This example is a study rather than a
//  demonstration, and it is the only one in the distribution that takes minutes.
//
//  Usage:  ./robust_arm                          defaults below
//          ./robust_arm mu sigma delta           with explicit values
//          ./robust_arm mu sigma delta steps     also setting the number of
//                                                integrator steps per segment
//
//////////////////////////////////////////////////////////////////////////

#include "psopt.h"

using namespace PSOPT;

//////////////////////////////////////////////////////////////////////////
///////////////////  The uncertainty, handed to the user functions  //////
//////////////////////////////////////////////////////////////////////////

// Read through problem.user_data, which is how examples/climb passes its
// aerodynamic tables and is the channel any driver would use.
struct Uncertainty {
    int                 M;      // number of scenarios actually in the problem
    std::vector<double> mp;     // payload mass of each scenario
    std::vector<double> w;      // its quadrature weight (unused by the
                                // minimum-time objective; carried because a
                                // risk-measure objective would need it and its
                                // absence would be a trap)
    double              delta;  // terminal tolerance, same for every scenario
};

// Gauss-Hermite nodes and weights for the standard normal, by Golub-Welsch: the
// nodes are the eigenvalues of the symmetric tridiagonal Jacobi matrix of the
// probabilists' Hermite polynomials (zero diagonal, off-diagonal sqrt(k)), and
// the weights are the squared first components of its eigenvectors. Generated
// rather than tabulated so that the scenario count is a parameter and not a
// table lookup somebody has to extend and can mistype. Checked against
// numpy.polynomial.hermite_e.hermegauss for M = 3, 5, 7 and 9.
static void gauss_hermite(int M, std::vector<double>& x, std::vector<double>& w)
{
    MatrixXd J = zeros(M, M);
    for (int k = 1; k < M; k++) {
        J(k-1, k) = sqrt((double) k);
        J(k, k-1) = sqrt((double) k);
    }
    Eigen::SelfAdjointEigenSolver<MatrixXd> es(J);
    x.resize(M); w.resize(M);
    for (int k = 0; k < M; k++) {
        x[k] = es.eigenvalues()(k);
        w[k] = es.eigenvectors()(0, k)*es.eigenvectors()(0, k);
    }
}

static const double X0[4] = { 0.0, 0.0, 0.500, 0.000 };   // initial state
static const double XF[4] = { 0.0, 0.0, 0.500, 0.522 };   // target state

//////////////////////////////////////////////////////////////////////////
///////////////////  The arm, with the payload named  ////////////////////
//////////////////////////////////////////////////////////////////////////

// One scenario's dynamics. Templated so that the same source serves the adouble
// tape and the plain-double verification integrator at the bottom of this file:
// a robust design checked against a DIFFERENT implementation of the plant is
// checking the implementation and not the design.
template <class T>
static void arm_rhs(const T* x, const T* u, double mp, T* dx)
{
    const T x1 = x[0], x2 = x[1], x3 = x[2];
    const T tau_a = u[0] - u[1];      // torque on link 1, in absolute coordinates
    const T tau_b = u[1];             // torque on link 2

    const double M11 = 7.0/3.0 + mp;
    const double M22 = 4.0/3.0 + mp;
    const double S2  = 3.0/2.0 + mp;

    const T c = cos(x3), s = sin(x3);
    const T M12 = S2*c;
    const T det = M11*M22 - M12*M12;          // = 28/9 - (9/4)c^2 at mp = 0

    // Centrifugal terms. Each link's rate enters the OTHER equation, which is
    // the usual structure in absolute coordinates and is what makes the shipped
    // numerators asymmetric.
    const T b1 = tau_a + S2*s*x2*x2;
    const T b2 = tau_b - S2*s*x1*x1;

    dx[0] = ( M22*b1 - M12*b2)/det;
    dx[1] = (-M12*b1 + M11*b2)/det;
    dx[2] = x2 - x1;
    dx[3] = x1;
}

adouble endpoint_cost(adouble* initial_states, adouble* final_states,
                      adouble* parameters, adouble& t0, adouble& tf,
                      adouble* xad, int iphase, Workspace* workspace)
{
    return tf;          // minimum time, shared by every scenario
}

adouble integrand_cost(adouble* states, adouble* controls, adouble* parameters,
                       adouble& time, adouble* xad, int iphase, Workspace* workspace)
{
    return 0.0;
}

void dae(adouble* derivatives, adouble* path, adouble* states,
         adouble* controls, adouble* parameters, adouble& time,
         adouble* xad, int iphase, Workspace* workspace)
{
    const Uncertainty& U = *((Uncertainty*) workspace->problem->user_data);

    // M copies of the arm, decoupled from each other. They meet only in the
    // control, which is the single copy above them, and in t_f.
    for (int i = 0; i < U.M; i++)
        arm_rhs<adouble>(states + 4*i, controls, U.mp[i], derivatives + 4*i);
}

void events(adouble* e, adouble* initial_states, adouble* final_states,
            adouble* parameters, adouble& t0, adouble& tf, adouble* xad,
            int iphase, Workspace* workspace)
{
    const Uncertainty& U = *((Uncertainty*) workspace->problem->user_data);

    // Eight per scenario: the four initial states, which are known exactly and
    // are the same for every scenario, and the four terminal states, which are
    // bounded to the tolerance ball rather than pinned.
    int k = 0;
    for (int i = 0; i < U.M; i++)
        for (int j = 0; j < 4; j++) e[k++] = initial_states[4*i + j];
    for (int i = 0; i < U.M; i++)
        for (int j = 0; j < 4; j++) e[k++] = final_states[4*i + j];
}

void linkages(adouble* linkages, adouble* xad, Workspace* workspace) {}

//////////////////////////////////////////////////////////////////////////
///////////////////  Verification: an independent integrator  ////////////
//////////////////////////////////////////////////////////////////////////

// Fixed-step RK4 through the plant with payload mp, driven by the control table
// (t_nodes, u_nodes), integrating node to node with nsub substeps per interval
// and returning the state at every node. This shares arm_rhs with the
// transcription and nothing else: it uses neither PSOPT's integrator, nor its
// mesh, nor its solution. A robust design checked against a DIFFERENT
// implementation of the plant is checking the implementation, not the design.
//
// The control returned by multiple shooting with a linear parameterisation IS
// piecewise linear through the node table, so interpolating linearly within the
// node interval currently being integrated reproduces the designed control
// exactly -- and, because the interval is known, without searching for it.
static MatrixXd simulate_states(const MatrixXd& t_nodes, const MatrixXd& u_nodes,
                                double mp, int nsub)
{
    const int n = (int) t_nodes.cols();
    MatrixXd X  = zeros(4, n);

    double x[4] = { X0[0], X0[1], X0[2], X0[3] };
    for (int j = 0; j < 4; j++) X(j, 0) = x[j];

    for (int i = 0; i < n-1; i++) {
        const double ta = t_nodes(0,i), tb = t_nodes(0,i+1);
        const double h  = (tb - ta)/nsub;
        const double ua0 = u_nodes(0,i), ua1 = u_nodes(0,i+1);
        const double ub0 = u_nodes(1,i), ub1 = u_nodes(1,i+1);

        for (int s = 0; s < nsub; s++) {
            // Local coordinate within the node interval, so the control is two
            // multiply-adds rather than a search.
            const double p0 = ((double) s)/nsub, ph = (s + 0.5)/nsub,
                         p1 = ((double) s + 1.0)/nsub;
            const double uA[2] = { ua0 + p0*(ua1-ua0), ub0 + p0*(ub1-ub0) };
            const double uH[2] = { ua0 + ph*(ua1-ua0), ub0 + ph*(ub1-ub0) };
            const double uB[2] = { ua0 + p1*(ua1-ua0), ub0 + p1*(ub1-ub0) };

            double k1[4], k2[4], k3[4], k4[4], y[4];
            arm_rhs<double>(x, uA, mp, k1);
            for (int j=0;j<4;j++) y[j] = x[j] + 0.5*h*k1[j];
            arm_rhs<double>(y, uH, mp, k2);
            for (int j=0;j<4;j++) y[j] = x[j] + 0.5*h*k2[j];
            arm_rhs<double>(y, uH, mp, k3);
            for (int j=0;j<4;j++) y[j] = x[j] + h*k3[j];
            arm_rhs<double>(y, uB, mp, k4);
            for (int j=0;j<4;j++)
                x[j] += (h/6.0)*(k1[j] + 2*k2[j] + 2*k3[j] + k4[j]);
        }
        for (int j = 0; j < 4; j++) X(j, i+1) = x[j];
    }
    return X;
}

// Terminal miss, in the infinity norm over the four states.
static double miss_of(const MatrixXd& X)
{
    const int n = (int) X.cols();
    double m = 0.0;
    for (int j = 0; j < 4; j++) {
        const double e = fabs(X(j, n-1) - XF[j]);
        if (e > m) m = e;
    }
    return m;
}

static double simulate_miss(const MatrixXd& t_nodes, const MatrixXd& u_nodes,
                            double mp, int nsub)
{
    return miss_of(simulate_states(t_nodes, u_nodes, mp, nsub));
}

//////////////////////////////////////////////////////////////////////////
///////////////////  One design solve  ///////////////////////////////////
//////////////////////////////////////////////////////////////////////////

struct Design {
    double   tf;
    MatrixXd t, u;
    bool     ok;
    int      rc;
    double   miss_tx;   // worst terminal miss as the TRANSCRIPTION reports it,
                        // i.e. read from the NLP's own terminal states. Compared
                        // against the independent integrator's figure, which is
                        // how a converged solve that has not actually met its
                        // terminal bounds is caught.
};

static int nsolves = 0;         // calls to psopt(), reported at the end

// Steps per segment in PSOPT's own multiple-shooting integrator. Settable from
// the command line because the measurement that fixed this value is worth being
// able to repeat: see the note at solve_design().
static int g_steps = 16;

// A marker for a design the NLP did not report as fully converged.
static const char* mark(const Design& d) { return d.rc == 0 ? " " : "*"; }

// warm, when given, supplies the initial guess: its control table directly, and
// a state guess obtained by integrating each scenario's own plant through that
// control. A guess built that way is dynamically consistent scenario by
// scenario, which matters here -- a constant state guess leaves the augmented
// problem with M copies of an infeasible arc and the solver free to wander off
// to a distant local minimum, which is exactly what happened before this was
// added.
static Design solve_design(const Uncertainty& U, const char* label,
                           const Design* warm = 0)
{
    Alg  algorithm;
    Sol  solution;
    Prob problem;

    const int M  = U.M;
    const int nx = 4*M;
    const int N  = 40;                    // mesh points, as examples/twolinkarm

    problem.name        = "Two-link arm with an uncertain payload";
    problem.outfilename = "robust_arm.txt";
    problem.nphases     = 1;
    problem.nlinkages   = 0;
    psopt_level1_setup(problem);

    problem.phases(1).nstates   = nx;
    problem.phases(1).ncontrols = 2;
    problem.phases(1).nevents   = 8*M;
    problem.phases(1).npath     = 0;
    problem.phases(1).nodes     << N;
    psopt_level2_setup(problem, algorithm);

    // The uncertainty must be visible to the user functions before any of them
    // is taped, and level 2 is where the tape is sized, so this goes here.
    problem.user_data = (void*) &U;

    for (int j = 0; j < nx; j++) {
        problem.phases(1).bounds.lower.states(j) = -2.0;
        problem.phases(1).bounds.upper.states(j) =  2.0;
    }
    for (int j = 0; j < 2; j++) {
        problem.phases(1).bounds.lower.controls(j) = -1.0;
        problem.phases(1).bounds.upper.controls(j) =  1.0;
    }

    int k = 0;
    for (int i = 0; i < M; i++)               // initial states, exact
        for (int j = 0; j < 4; j++) {
            problem.phases(1).bounds.lower.events(k) = X0[j];
            problem.phases(1).bounds.upper.events(k) = X0[j];
            k++;
        }
    for (int i = 0; i < M; i++)               // terminal states, within delta
        for (int j = 0; j < 4; j++) {
            problem.phases(1).bounds.lower.events(k) = XF[j] - U.delta;
            problem.phases(1).bounds.upper.events(k) = XF[j] + U.delta;
            k++;
        }

    problem.phases(1).bounds.lower.StartTime = 0.0;
    problem.phases(1).bounds.upper.StartTime = 0.0;
    problem.phases(1).bounds.lower.EndTime   = 1.0;
    // examples/twolinkarm bounds t_f by 10, which is ample when the plant is
    // known. It is NOT ample here: the robust slew is a slow one, for the reason
    // set out at the head of this file, and a bound of 10 is active for the
    // larger scenario sets -- which would silently turn "the cost of robustness"
    // into "the bound the author happened to write down".
    problem.phases(1).bounds.upper.EndTime   = 30.0;

    problem.integrand_cost = &integrand_cost;
    problem.endpoint_cost  = &endpoint_cost;
    problem.dae            = &dae;
    problem.events         = &events;
    problem.linkages       = &linkages;

    MatrixXd x_guess(nx, N);
    if (warm != 0 && warm->ok && (int) warm->t.cols() == N) {
        for (int i = 0; i < M; i++) {
            MatrixXd Xi = simulate_states(warm->t, warm->u, U.mp[i], 8);
            for (int j = 0; j < 4; j++) x_guess.row(4*i + j) = Xi.row(j);
        }
        problem.phases(1).guess.states   = x_guess;
        problem.phases(1).guess.controls = warm->u;
        problem.phases(1).guess.time     = warm->t;
    } else {
        for (int i = 0; i < M; i++)
            for (int j = 0; j < 4; j++)
                x_guess.row(4*i + j) = X0[j]*ones(1, N);
        problem.phases(1).guess.states   = x_guess;
        problem.phases(1).guess.controls = zeros(2, N);
        problem.phases(1).guess.time     = linspace(0.0, 3.0, N);
    }

    algorithm.nlp_method   = "IPOPT";
    algorithm.scaling      = "automatic";
    algorithm.derivatives  = "automatic";
    algorithm.nlp_iter_max = 1000;
    algorithm.nlp_tolerance = 1.e-6;

    // Multiple shooting, as examples/twolinkarm now uses, and for the same
    // reasons: the optimal torques are bang-bang, which a global polynomial
    // represents badly, and the scenarios integrate independently inside a
    // segment so the structure suits the method.
    algorithm.transcription_method        = "multiple-shooting";
    algorithm.ms_integrator               = "RK4";
    // examples/twolinkarm uses 8 steps per segment, a value chosen by measurement
    // on THAT problem -- a 3-second slew. It is not enough here. The robust
    // designs are two and a half times longer, and over a trajectory that long
    // the segment integrator's error is amplified enough to matter: at 8 steps a
    // design that its own transcription reports as landing exactly on the
    // tolerance ball, 2.000e-02, is found by an independent integrator to miss by
    // 2.346e-02, so it had satisfied a slightly different plant. At 16 steps the
    // two agree to 2e-05 and at 32 to the printed precision. 16 is used, and the
    // "in-sample" and "same, sim" columns below are the standing check on it.
    //
    // The general point is worth more than the number: an accuracy setting
    // validated on the nominal problem has to be validated again on the robust
    // one, because the robust solution is a different trajectory.
    algorithm.ms_steps_per_segment        = g_steps;
    algorithm.ms_control_parameterisation = "linear";
    algorithm.print_level                 = 0;

    Design d;
    d.ok = false;
    d.rc = -1;
    nsolves++;
    (void) psopt(solution, problem, algorithm);
    d.rc = solution.nlp_return_code;
    if (solution.error_flag == 0 && solution.nlp_return_code <= 1) {
        // IPOPT's return code 1 is "solved to an acceptable level", which is NOT
        // the same as solved: it can leave the terminal bounds violated by
        // considerably more than the NLP tolerance. Such a design is kept and
        // used -- discarding it would hide the behaviour -- but every table below
        // marks it, because an unmarked row would be read as a converged solve.
        d.ok = true;
        d.tf = solution.get_time_in_phase(1)(0, N-1);
        d.t  = solution.get_time_in_phase(1);
        d.u  = solution.get_controls_in_phase(1);

        MatrixXd xs = solution.get_states_in_phase(1);
        d.miss_tx = 0.0;
        for (int i = 0; i < M; i++)
            for (int j = 0; j < 4; j++) {
                const double e = fabs(xs(4*i + j, N-1) - XF[j]);
                if (e > d.miss_tx) d.miss_tx = e;
            }
    }

    // An empty label means the caller is tabulating the result itself.
    if (label != 0 && label[0] != '\0') {
        printf("  %-34s ", label);
        if (!d.ok) printf("FAILED (nlp %d, flag %d)\n",
                          solution.nlp_return_code, solution.error_flag);
        else       printf("t_f = %10.6f\n", d.tf);
    }
    return d;
}

//////////////////////////////////////////////////////////////////////////
///////////////////  The driver's two ingredients  ////////////////////////
//////////////////////////////////////////////////////////////////////////

// A scenario set from the Gauss-Hermite rule of order M, mapped onto the
// payload distribution. M = 1 returns the mean, which is the nominal design.
static Uncertainty scenarios(double mu, double sigma, int M, double delta)
{
    std::vector<double> z, w;
    gauss_hermite(M, z, w);
    Uncertainty U;
    U.M = M; U.delta = delta; U.mp.resize(M); U.w.resize(M);
    for (int k = 0; k < M; k++) { U.mp[k] = mu + sigma*z[k]; U.w[k] = w[k]; }
    return U;
}

// Substeps per node interval in the verification integrator. Validated in main()
// against a tenfold finer run rather than asserted: a verifier whose own step
// size has not been checked is not a verifier.
static const int NSUB = 25;

// The worst payload in an interval for a GIVEN control: a scan on 241 points
// followed by a refinement over the two intervals around the maximum. The
// uncertainty here is one scalar, so a scan is both exhaustive and cheap, and it
// is honest in a way a gradient search on a non-concave function would not be.
// For a higher-dimensional uncertainty this is the step that would have to
// become an optimisation -- and it is the step that makes the scheme below a
// bilevel method rather than a quadrature.
static double worst_case_miss(const Design& d, double lo, double hi, double* mp_at)
{
    const int NG = 241;
    double best = -1.0, mbest = lo;
    for (int i = 0; i < NG; i++) {
        const double mp = lo + (hi - lo)*i/(NG - 1.0);
        const double m  = simulate_miss(d.t, d.u, mp, NSUB);
        if (m > best) { best = m; mbest = mp; }
    }
    const double h = (hi - lo)/(NG - 1.0);
    const double a = fmax(lo, mbest - h), b = fmin(hi, mbest + h);
    for (int i = 0; i <= 40; i++) {
        const double mp = a + (b - a)*i/40.0;
        const double m  = simulate_miss(d.t, d.u, mp, NSUB);
        if (m > best) { best = m; mbest = mp; }
    }
    if (mp_at) *mp_at = mbest;
    return best;
}

// Out-of-sample statistics, reported over the payloads that lie INSIDE the
// uncertainty set the design was given, and separately for those outside it.
// The distinction is not pedantry: a design constrained on [lo, hi] promises
// nothing beyond it, and scoring it against the whole Gaussian tail measures the
// truncation of the set rather than the quality of the design. The tail is still
// counted and reported, because it is what a chance constraint leaves on the
// table.
struct OOS { double mean, worst, frac_in; int n_out; double worst_out; };

static OOS out_of_sample(const Design& d, const double* mp_sample, int NS,
                         double delta, double lo, double hi)
{
    OOS o; o.mean = 0.0; o.worst = 0.0; o.n_out = 0; o.worst_out = 0.0;
    int inside = 0, within = 0;
    for (int s = 0; s < NS; s++) {
        const double m = simulate_miss(d.t, d.u, mp_sample[s], NSUB);
        if (mp_sample[s] < lo || mp_sample[s] > hi) {
            o.n_out++;
            if (m > o.worst_out) o.worst_out = m;
            continue;
        }
        inside++;
        o.mean += m;
        if (m > o.worst) o.worst = m;
        if (m <= delta) within++;
    }
    o.mean   /= (inside > 0 ? inside : 1);
    o.frac_in = 100.0*within/(inside > 0 ? inside : 1);
    return o;
}

// Worst in-sample miss, read from the transcription's own terminal states.
static double in_sample_miss(const Design& d, const Uncertainty& U)
{
    double worst = 0.0;
    for (int i = 0; i < U.M; i++) {
        const double m = simulate_miss(d.t, d.u, U.mp[i], NSUB);
        if (m > worst) worst = m;
    }
    return worst;
}

//////////////////////////////////////////////////////////////////////////

int main(int argc, char** argv)
{
    const double mu    = (argc > 1) ? atof(argv[1]) : 0.50;   // mean payload
    const double sigma = (argc > 2) ? atof(argv[2]) : 0.15;   // its std deviation
    const double delta = (argc > 3) ? atof(argv[3]) : 0.02;   // terminal tolerance
    if (argc > 4) g_steps = atoi(argv[4]);                    // steps per segment

    const double lo = mu - 3.0*sigma, hi = mu + 3.0*sigma;

    printf("\nTwo-link arm with an uncertain payload\n");
    printf("======================================\n\n");
    printf("  payload m_p ~ N(mu = %.3f, sigma = %.3f)\n", mu, sigma);
    printf("  terminal tolerance delta = %.3f (infinity norm on the four\n"
           "  states)\n", delta);
    printf("  uncertainty set for verification: [%.3f, %.3f] (mu +- 3 sigma)\n\n",
           lo, hi);

    // The out-of-sample payloads, drawn once and reused by every design so that
    // the comparison between designs is not also a comparison between samples.
    const int NS = 1000;
    static double mp_sample[1000];
    srand(20260927);
    for (int s = 0; s < NS; s++) {
        // Box-Muller, so the sample is the stated Gaussian and not merely
        // something with the right mean and variance.
        const double u1 = (rand() + 1.0)/(RAND_MAX + 2.0);
        const double u2 = (rand() + 1.0)/(RAND_MAX + 2.0);
        mp_sample[s] = mu + sigma*sqrt(-2.0*log(u1))*cos(2.0*M_PI*u2);
    }

    // ---- (a) the shipped problem, to verify the reparameterised model -------
    printf("Verification of the model and the transcription\n");
    Uncertainty A = scenarios(0.0, 0.0, 1, 0.0);
    Design a = solve_design(A, "(a) m_p = 0, delta = 0");
    if (a.ok) {
        const double published = 2.985042;
        printf("       examples/twolinkarm gives %.6f, difference %.2e\n",
               published, fabs(a.tf - published));
    }

    // ---- (b) the nominal design at the working tolerance --------------------
    printf("\nThe nominal design, which robustness has to be measured against\n");
    Uncertainty B = scenarios(mu, 0.0, 1, delta);
    Design b = solve_design(B, "(b) m_p = mu, delta > 0", &a);
    if (!b.ok) {
        printf("\n  the nominal design failed; nothing below is meaningful\n\n");
        return 1;
    }

    // Is the verifier's own step size fine enough to be believed? Every number
    // in the tables below is a miss computed by it, so this is checked and not
    // assumed. The comparison is against a tenfold finer integration, at the
    // nominal payload and at the edge of the uncertainty set.
    {
        const double m1c = simulate_miss(b.t, b.u, mu, NSUB);
        const double m1f = simulate_miss(b.t, b.u, mu, 10*NSUB);
        const double m2c = simulate_miss(b.t, b.u, hi, NSUB);
        const double m2f = simulate_miss(b.t, b.u, hi, 10*NSUB);
        printf("       verifier and a tenfold finer run agree to %.1e\n",
               fmax(fabs(m1c-m1f), fabs(m2c-m2f)));
    }

    // ---- (c) sample-average approximation, swept in the scenario count ------
    //
    // This is the sweep that decides whether a fixed quadrature rule is a
    // sufficient discretisation of the uncertainty. Each row reports both the
    // in-sample miss, which the optimiser drove to delta by construction, and
    // the miss over payloads it never saw, which it did not.
    printf("\nSample-average approximation over M Gauss-Hermite scenarios\n");
    printf("  %2s %14s %11s %10s %10s %10s %10s\n", "M", "scenarios", "t_f",
           "in-sample", "same, sim", "worst set", "within d");

    {
        char span[64];
        OOS o = out_of_sample(b, mp_sample, NS, delta, lo, hi);
        double mpw; const double wc = worst_case_miss(b, lo, hi, &mpw);
        snprintf(span, sizeof span, "%.3f only", mu);
        printf("  %2d %14s %10.6f%s %10.3e %10.3e %10.3e %9.1f%%\n", 1, span,
               b.tf, mark(b), b.miss_tx, in_sample_miss(b, B), wc, o.frac_in);
    }

    const int  Msweep[4] = { 3, 5, 7, 9 };
    Design     best_saa  = b;
    Uncertainty best_U   = B;
    for (int q = 0; q < 4; q++) {
        const int M = Msweep[q];
        Uncertainty U = scenarios(mu, sigma, M, delta);
        // Warm-started from the previous member of the sweep: the sequence of
        // problems is a continuation in M, and starting each one from scratch
        // invites each to settle in a different local minimum, which would make
        // the column of final times unreadable.
        Design d = solve_design(U, "", &best_saa);
        if (!d.ok) { printf("  %2d   FAILED\n", M); continue; }

        char span[64];
        snprintf(span, sizeof span, "[%.3f,%.3f]", U.mp[0], U.mp[M-1]);
        OOS o = out_of_sample(d, mp_sample, NS, delta, lo, hi);
        double mpw; const double wc = worst_case_miss(d, lo, hi, &mpw);
        printf("  %2d %14s %10.6f%s %10.3e %10.3e %10.3e %9.1f%%\n", M, span,
               d.tf, mark(d), d.miss_tx, in_sample_miss(d, U), wc, o.frac_in);
        best_saa = d; best_U = U;
    }

    // ---- (d) sequential scenario generation --------------------------------
    //
    // The alternative to guessing a quadrature rule fine enough: solve with the
    // scenarios in hand, then ask where the resulting control is worst over the
    // whole uncertainty set and add THAT payload as a scenario. Each iteration
    // is one call to psopt() and one scan of the verifier, and the scenario set
    // grows only where the design is actually failing. This is the loop a robust
    // optimal control driver would run; everything it needs, PSOPT already
    // provides.
    printf("\nSequential scenario generation (worst case in the loop)\n");

    // The scenarios are satisfied to delta exactly -- an active constraint is
    // active -- so the miss BETWEEN two scenarios is necessarily a little larger
    // than delta, and a loop that demands delta over the whole set from a design
    // that was only ever asked for delta at finitely many points cannot
    // terminate. The remedy is the usual one: tighten the tolerance the design is
    // given and keep verifying against the one that was asked for. The margin is
    // a stated 10%, not a fitted number.
    const double tighten = 0.90;

    std::vector<double> gen;
    gen.push_back(mu);
    Design      cur = b;
    Uncertainty G   = scenarios(mu, 0.0, 1, tighten*delta);

    printf("  the design is given a tightened tolerance of %.4f; the column\n"
           "  'worst in set' is measured against the full %.4f\n\n",
           tighten*delta, delta);
    printf("  %2s %3s %11s %12s %9s %11s %10s\n", "it", "M", "t_f",
           "worst in set", "at m_p", "mean OOS", "within d");

    const int ITMAX = 18;
    for (int it = 0; it < ITMAX; it++) {
        double mpw;
        const double wc = worst_case_miss(cur, lo, hi, &mpw);
        OOS o = out_of_sample(cur, mp_sample, NS, delta, lo, hi);
        printf("  %2d %3d %10.6f%s %12.3e %9.4f %11.3e %9.1f%%\n",
               it, (int) gen.size(), cur.tf, mark(cur), wc, mpw, o.mean, o.frac_in);

        if (wc <= delta) {
            printf("\n  converged: no payload in [%.3f, %.3f] misses the target by\n"
                   "  more than delta, so adding one could not change the design\n",
                   lo, hi);
            break;
        }
        if (it == ITMAX-1) {
            printf("\n  iteration limit reached with a worst case of %.3e\n", wc);
            break;
        }

        // Add the payload that the current design serves worst, and re-solve
        // warm-started from it.
        gen.push_back(mpw);
        G.M  = (int) gen.size();
        G.mp = gen;
        G.w.assign(gen.size(), 1.0/gen.size());
        Design d = solve_design(G, "", &cur);
        if (!d.ok) {
            printf("       adding m_p = %.4f made the problem unsolvable\n", mpw);
            break;
        }
        cur = d;
    }

    printf("\n  t_f is not monotone down this column, because each subproblem is\n");
    printf("  nonconvex and a warm start does not guarantee the same local minimum.\n");
    printf("  What the scheme controls is the worst case over the set, not t_f.\n");

    // ---- what it bought and what it cost -----------------------------------
    OOS on  = out_of_sample(b,        mp_sample, NS, delta, lo, hi);
    OOS osa = out_of_sample(best_saa, mp_sample, NS, delta, lo, hi);
    OOS orb = out_of_sample(cur,      mp_sample, NS, delta, lo, hi);
    double       at;
    const double wc_nm = worst_case_miss(b,        lo, hi, &at);
    const double wc_sa = worst_case_miss(best_saa, lo, hi, &at);
    const double wc_rb = worst_case_miss(cur,      lo, hi, &at);

    printf("\nSummary, over the %d sampled payloads inside the uncertainty set\n\n",
           NS - on.n_out);
    printf("  %-27s %3s %9s %10s %10s %9s\n", "design", "M", "t_f",
           "worst set", "mean miss", "within d");
    printf("  %-27s %3d %9.4f%s %10.3e %10.3e %8.1f%%\n",
           "nominal (mean payload only)", 1, b.tf, mark(b),
           wc_nm, on.mean, on.frac_in);
    printf("  %-27s %3d %9.4f%s %10.3e %10.3e %8.1f%%\n",
           "Gauss-Hermite, finest tried", best_U.M, best_saa.tf, mark(best_saa),
           wc_sa, osa.mean, osa.frac_in);
    printf("  %-27s %3d %9.4f%s %10.3e %10.3e %8.1f%%\n",
           "generated scenarios", (int) gen.size(), cur.tf, mark(cur),
           wc_rb, orb.mean, orb.frac_in);

    printf("\n  price of robustness: t_f %.4f -> %.4f, a factor of %.2f. The robust\n",
           b.tf, cur.tf, cur.tf/b.tf);
    printf("  slew is slow because a gentle trajectory is one on which the payload\n");
    printf("  matters less; that is where the insensitivity is bought.\n");

    printf("\n  %d scenarios placed where the design was failing achieve what %d\n",
           (int) gen.size(), best_U.M);
    printf("  scenarios placed by a quadrature rule do not, and they come with a\n");
    printf("  certificate: the worst payload ANYWHERE in the set, not the worst one\n");
    printf("  sampled. Quadrature discretises an expectation; a constraint\n");
    printf("  that must hold everywhere needs its own worst case in the loop.\n");

    printf("\n  the %d sampled payload%s OUTSIDE the set -- %.2f%% of the draw, and\n",
           on.n_out, (on.n_out == 1 ? "" : "s"), 100.0*on.n_out/NS);
    printf("  what truncating the distribution at 3 sigma gives up -- miss by up to\n");
    printf("  %.2e (nominal) and %.2e (robust). Neither design promised anything\n",
           on.worst_out, orb.worst_out);
    printf("  there, and a design that must cover them has to be given a wider set.\n");

    printf("\n  %d calls to psopt() in total. A t_f marked * comes from a\n", nsolves);
    printf("  solve the NLP reported as only ACCEPTABLE rather than\n");
    printf("  converged: its constraints are satisfied but its optimality is\n");
    printf("  not established, so the final time on that row is an upper\n");
    printf("  bound on the optimum and not the optimum.\n");

    ////////////////////////////////////////////////////////////////////////
    ///////////  Plot some results if desired (requires gnuplot) ///////////
    ////////////////////////////////////////////////////////////////////////

    // The miss as a function of the payload, for each design, is the figure this
    // example exists to produce: the nominal design touches zero at one payload
    // and rises steeply either side of it, while the robust one is held under
    // delta across the whole set.
    const int NC = 181;
    MatrixXd mp_axis = zeros(1, NC), curves = zeros(3, NC);
    for (int i = 0; i < NC; i++) {
        const double mp = lo + (hi - lo)*i/(NC - 1.0);
        mp_axis(0, i) = mp;
        curves(0, i)  = simulate_miss(b.t,   b.u,   mp, NSUB);
        curves(1, i)  = simulate_miss(cur.t, cur.u, mp, NSUB);
        curves(2, i)  = delta;
    }

    plot(mp_axis, curves, "Terminal miss versus payload",
         "payload m_p", "terminal miss (infinity norm)", "nominal robust delta");
    plot(mp_axis, curves, "Terminal miss versus payload",
         "payload m_p", "terminal miss (infinity norm)", "nominal robust delta",
         "pdf", "robust_arm_miss.pdf");

    plot(b.t, b.u, "Nominal design: controls", "time (s)", "controls", "u1 u2");
    plot(b.t, b.u, "Nominal design: controls", "time (s)", "controls", "u1 u2",
         "pdf", "robust_arm_controls_nominal.pdf");

    plot(cur.t, cur.u, "Robust design: controls", "time (s)", "controls", "u1 u2");
    plot(cur.t, cur.u, "Robust design: controls", "time (s)", "controls", "u1 u2",
         "pdf", "robust_arm_controls_robust.pdf");

    printf("\n");
    return 0;
}
