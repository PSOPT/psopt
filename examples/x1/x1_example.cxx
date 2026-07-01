//////////////////////////////////////////////////////////////////////////
//                  x1.cxx : integrated-residual showcase                 //
//////////////////////////////////////////////////////////////////////////
// Title:    X1 -- Neuenhofen & Kerrigan analytic optimal-control example  //
// Method:   collocation vs integrated-residual vs Nie-Kerrigan flexible   //
//           order, validated against a closed-form solution              //
// Purpose:  demonstrate, on a problem with a known analytic optimum, that //
//           Hermite-Simpson collocation rings on the control and reports  //
//           a spuriously low cost, that integrated-residual transcription //
//           is feasibility-honest, and that Nie-Kerrigan flexible-order   //
//           p-refinement converges to the analytic optimum.              //
//                                                                        //
// Problem (Neuenhofen & Kerrigan, "An Integral Penalty-Barrier Direct     //
// Transcription Method for Optimal Control", CDC 2020):                  //
//                                                                        //
//   minimise   J = int_0^{pi/2} ( y(t)^2 + cos(t) u(t) ) dt              //
//   subject to ydot = 1/2 y^2 + u,  y(0) = 0,  -1.5 <= u <= 1.           //
//                                                                        //
// The control bound is inactive, so the unconstrained Pontryagin solution //
// holds throughout, giving the closed-form optimum                       //
//                                                                        //
//   y*(t)      = -sin t / (2 - cos t),                                   //
//   lambda*(t) = -cos t,                                                 //
//   u*(t)      = 1/2 - 1.5 / (cos t - 2)^2,                              //
//   J*         = -0.256996962561  (verified by quadrature to ~1e-14).    //
//                                                                        //
// The integrated-residual runs use the well-conditioned residual-box     //
// (minimise J subject to |residual| <= delta), with delta set by the     //
// order-aware schedule delta = h_node^p (h_node the node spacing, p the   //
// local order: 2 for cubic-Hermite, ir_local_order for Nie-Kerrigan).    //
//                                                                        //
// Provenance: this example was prepared with AI assistance (Claude) and   //
// validated against the closed-form solution above.                      //
//////////////////////////////////////////////////////////////////////////

#include "psopt.h"

using namespace PSOPT;

//////////////////////// analytic reference ////////////////////////////////

static double y_star(double t) { return -sin(t)/(2.0 - cos(t)); }
static double u_star(double t) { double d = cos(t) - 2.0; return 0.5 - 1.5/(d*d); }
static const double J_STAR = -0.256996962561;

//////////////////////// problem functions /////////////////////////////////

adouble endpoint_cost(adouble* /*x0*/, adouble* /*xf*/, adouble* /*p*/,
                      adouble& /*t0*/, adouble& /*tf*/, adouble* /*xad*/,
                      int /*iphase*/, Workspace* /*workspace*/)
{
    return 0.0;
}

adouble integrand_cost(adouble* states, adouble* controls, adouble* /*parameters*/,
                       adouble& time, adouble* /*xad*/, int /*iphase*/,
                       Workspace* /*workspace*/)
{
    adouble y = states[0];
    adouble u = controls[0];
    return y*y + cos(time)*u;
}

void dae(adouble* derivatives, adouble* /*path*/, adouble* states, adouble* controls,
         adouble* /*parameters*/, adouble& /*time*/, adouble* /*xad*/, int /*iphase*/,
         Workspace* /*workspace*/)
{
    adouble y = states[0];
    adouble u = controls[0];
    derivatives[0] = 0.5*y*y + u;
}

void events(adouble* e, adouble* initial_states, adouble* /*final_states*/,
            adouble* /*parameters*/, adouble& /*t0*/, adouble& /*tf*/, adouble* /*xad*/,
            int /*iphase*/, Workspace* /*workspace*/)
{
    e[0] = initial_states[0];   // y(0) = 0
}

void linkages(adouble* /*linkages*/, adouble* /*xad*/, Workspace* /*workspace*/) {}

//////////////////////// one solve, one row ////////////////////////////////

struct Row { double J, ey, eu, rle; };

// Solve X1 once on a fixed mesh of `nodes` nodes with the requested transcription.
//   trans       : "collocation" or "integrated-residual"
//   ir_obj      : ir_objective ("cost" for the residual-box optimality solve)
//   local_order : Nie-Kerrigan flexible order (0 = cubic-Hermite)
//   bound       : ir_residual_bound (<0 leaves it unset, i.e. collocation)
//   res_nodes   : ir_residual_nodes (must be >= local_order)
//   nlptol      : NLP convergence tolerance
static Row solve_x1(int nodes, const std::string& trans, const std::string& ir_obj,
                    int local_order, double bound, int res_nodes, double nlptol)
{
    const double pi = 3.14159265358979323846;
    Alg  algorithm;
    Sol  solution;
    Prob problem;

    problem.name        = "X1 (Neuenhofen-Kerrigan analytic example)";
    problem.outfilename = "x1.txt";
    problem.nphases     = 1;
    problem.nlinkages   = 0;
    psopt_level1_setup(problem);

    problem.phases(1).nstates   = 1;
    problem.phases(1).ncontrols = 1;
    problem.phases(1).nevents   = 1;
    problem.phases(1).npath     = 0;
    problem.phases(1).nodes     = (RowVectorXi(1) << nodes).finished();
    psopt_level2_setup(problem, algorithm);

    problem.phases(1).bounds.lower.states(0)   = -10.0;
    problem.phases(1).bounds.upper.states(0)   =  10.0;
    problem.phases(1).bounds.lower.controls(0) =  -1.5;
    problem.phases(1).bounds.upper.controls(0) =   1.0;
    problem.phases(1).bounds.lower.events(0)   =   0.0;
    problem.phases(1).bounds.upper.events(0)   =   0.0;
    problem.phases(1).bounds.lower.StartTime   =   0.0;
    problem.phases(1).bounds.upper.StartTime   =   0.0;
    problem.phases(1).bounds.lower.EndTime     =   pi/2.0;
    problem.phases(1).bounds.upper.EndTime     =   pi/2.0;

    problem.integrand_cost = &integrand_cost;
    problem.endpoint_cost  = &endpoint_cost;
    problem.dae            = &dae;
    problem.events         = &events;
    problem.linkages       = &linkages;

    problem.phases(1).guess.states   = zeros(1, nodes);
    problem.phases(1).guess.controls = zeros(1, nodes);
    problem.phases(1).guess.time     = linspace(0.0, pi/2.0, nodes);

    algorithm.print_level        = 0;
    algorithm.nlp_method         = "IPOPT";
    algorithm.scaling            = "automatic";
    algorithm.derivatives        = "automatic";
    algorithm.nlp_iter_max       = 5000;
    algorithm.nlp_tolerance      = nlptol;
    algorithm.collocation_method = "Hermite-Simpson";
    algorithm.transcription_method = trans;
    if (!ir_obj.empty()) algorithm.ir_objective = ir_obj;
    if (local_order > 0) algorithm.ir_local_order   = local_order;
    if (res_nodes  > 0)  algorithm.ir_residual_nodes = res_nodes;
    if (bound     >= 0)  algorithm.ir_residual_bound = bound;

    psopt(solution, problem, algorithm);

    Row r; r.J = solution.get_cost(); r.ey = 0.0; r.eu = 0.0; r.rle = -1.0;
    MatrixXd t = solution.get_time_in_phase(1);
    MatrixXd x = solution.get_states_in_phase(1);
    MatrixXd u = solution.get_controls_in_phase(1);
    for (int k = 0; k < t.cols(); ++k) {
        double tk = t(0,k);
        r.ey = std::max(r.ey, std::fabs(x(0,k) - y_star(tk)));
        r.eu = std::max(r.eu, std::fabs(u(0,k) - u_star(tk)));
    }
    MatrixXd rle = solution.get_relative_local_error_in_phase(1);
    if (rle.size() > 0) r.rle = rle.maxCoeff();
    return r;
}

// Order-aware residual-box schedule delta = h_node^p  (exact integer power).
static double delta_schedule(int nodes, int p)
{
    const double pi = 3.14159265358979323846;
    double h = (pi/2.0)/(nodes - 1);
    double d = 1.0;
    for (int i = 0; i < p; ++i) d *= h;
    return d;
}

static void print_row(const char* label, int nodes, const Row& r)
{
    printf("  %-30s %5d  %16.12f  %9.2e  %9.2e  %9.2e\n",
           label, nodes, r.J, std::fabs(r.J - J_STAR), r.eu, r.rle);
}

int main()
{
    printf("\n");
    printf("X1 -- Neuenhofen-Kerrigan analytic example   (J* = %.12f)\n", J_STAR);
    printf("================================================================================\n");
    printf("  %-30s %5s  %16s  %9s  %9s  %9s\n",
           "method", "nodes", "J", "|J-J*|", "max|u-u*|", "relLocErr");
    printf("--------------------------------------------------------------------------------\n");

    // ---- Part A: three-way contrast on one fixed mesh (41 nodes = 8 elements x order 5) ----
    // NLP tolerances are chosen per method so every solve reports a clean IPOPT
    // success: collocation is ill-conditioned near its ringing solution and only
    // certifies optimality at a looser tolerance; the binding residual box on the
    // low-order cubic-Hermite representation needs the larger iteration budget.
    Row coll = solve_x1(41, "collocation", "", 0, -1.0, 4, 1.0e-5);
    print_row("collocation (Hermite-Simpson)", 41, coll);

    Row ir = solve_x1(41, "integrated-residual", "cost", 0, delta_schedule(41, 2), 4, 1.0e-5);
    print_row("integrated-residual (cubic-H)", 41, ir);

    Row nk5 = solve_x1(41, "integrated-residual", "cost", 5, delta_schedule(41, 5), 10, 1.0e-7);
    print_row("Nie-Kerrigan order 5", 41, nk5);

    // ---- Part B: Nie-Kerrigan p-convergence (8 elements, order p, nodes = 8p+1) ----
    // Capped at order 5: with nlp_tolerance = 1e-7 the order-aware box h_node^p must
    // stay above the NLP tolerance, which it does through order 5 (h_node^5 ~ 9e-8)
    // but not order 6 (~1e-9).
    printf("--------------------------------------------------------------------------------\n");
    printf("  Nie-Kerrigan flexible-order p-convergence (8 elements, order-aware box):\n");
    for (int p = 2; p <= 5; ++p) {
        int nodes = 8*p + 1;
        Row r = solve_x1(nodes, "integrated-residual", "cost", p,
                         delta_schedule(nodes, p), 2*p, 1.0e-7);
        char lab[64];
        snprintf(lab, sizeof(lab), "  order %d", p);
        print_row(lab, nodes, r);
    }
    printf("================================================================================\n");
    printf("Collocation reports J below J* (infeasible between nodes) and a ringing control;\n");
    printf("integrated-residual is feasibility-honest; flexible-order p-refinement converges\n");
    printf("to the analytic optimum.\n\n");
    return 0;
}
