// _psopt_driver.cpp  -- Layer A (B-1): minimal single-phase driver.
// Exposes solve_single_phase(spec) to Python: builds the PSOPT Prob/Alg in the
// correct level1->dims->level2->bounds->register->guess->algorithm order, dlopens a
// JIT-compiled math .so (Layer B) and registers its extern "C" user functions as
// PSOPT function pointers, solves, and returns the solution as NumPy arrays.
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/eigen.h>
#include <dlfcn.h>
#include <stdexcept>
#include <string>
#include "psopt.h"
#include "integer_controls.h"
#include "integer_parameters.h"
#include <vector>

namespace py = pybind11;
using namespace PSOPT;

// Exact PSOPT user-function pointer signatures (from psopt.h).
typedef void   (*dae_t)      (adouble*,adouble*,adouble*,adouble*,adouble*,adouble&,adouble*,int,Workspace*);
typedef adouble(*endpoint_t) (adouble*,adouble*,adouble*,adouble&,adouble&,adouble*,int,Workspace*);
typedef adouble(*integrand_t)(adouble*,adouble*,adouble*,adouble&,adouble*,int,Workspace*);
typedef void   (*events_t)   (adouble*,adouble*,adouble*,adouble*,adouble&,adouble&,adouble*,int,Workspace*);
typedef void   (*linkages_t) (adouble*,adouble*,Workspace*);
typedef void   (*observation_t)(adouble*,adouble*,adouble*,adouble*,adouble&,int,adouble*,int,Workspace*);

static void* must_sym(void* h, const char* n) {
    dlerror();
    void* s = dlsym(h, n);
    if (!s) throw std::runtime_error(std::string("symbol not found in math library: ") + n);
    return s;
}

static void apply_algorithm(Alg& a, py::dict o) {
    auto S = [&](const char* k, std::string& f){ if (o.contains(k)) f = py::cast<std::string>(o[k]); };
    auto I = [&](const char* k, int& f){ if (o.contains(k)) f = py::cast<int>(o[k]); };
    auto D = [&](const char* k, double& f){ if (o.contains(k)) f = py::cast<double>(o[k]); };
    auto B = [&](const char* k, bool& f){ if (o.contains(k)) f = py::cast<bool>(o[k]); };
    // core
    S("nlp_method", a.nlp_method); S("scaling", a.scaling); S("derivatives", a.derivatives);
    S("collocation_method", a.collocation_method); I("nlp_iter_max", a.nlp_iter_max);
    D("nlp_tolerance", a.nlp_tolerance);
    // general extras (quiet mode, linear solver, hessian, diagnostics)
    S("constraint_scaling", a.constraint_scaling); S("hessian", a.hessian);
    S("ipopt_linear_solver", a.ipopt_linear_solver); I("print_level", a.print_level);
    D("ipopt_max_cpu_time", a.ipopt_max_cpu_time); I("diagnostic_level", a.diagnostic_level);
    // hp-adaptive mesh refinement
    S("mesh_refinement", a.mesh_refinement); I("mr_max_iterations", a.mr_max_iterations);
    D("ode_tolerance", a.ode_tolerance); D("mr_max_growth_factor", a.mr_max_growth_factor);
    I("mr_min_order", a.mr_min_order); I("mr_max_order", a.mr_max_order);
    // integrated-residual transcription / Nie-Kerrigan flexible-order
    S("transcription_method", a.transcription_method); I("ir_residual_nodes", a.ir_residual_nodes);
    D("ir_regularization", a.ir_regularization); S("ir_objective", a.ir_objective);
    D("ir_residual_bound", a.ir_residual_bound); B("ir_dair", a.ir_dair);
    D("ir_dair_delta_factor", a.ir_dair_delta_factor); I("ir_local_order", a.ir_local_order);
    S("ir_include_path", a.ir_include_path); D("ir_path_weight", a.ir_path_weight);
    S("ir_residual_scaling", a.ir_residual_scaling);
    B("ir_element_local_controls", a.ir_element_local_controls);
    // The flexible mesh. Reachable from Python for the same reason every other
    // integrated-residual option is: an option the C++ interface has and the Python one
    // does not is an option a Python user cannot discover exists.
    B("ir_flexible_mesh", a.ir_flexible_mesh);
    D("ir_min_element_fraction", a.ir_min_element_fraction);
    // Multiple shooting. transcription_method above is what selects it, and without these
    // seven the selection reaches Python and nothing that configures it does -- so a Python
    // user would get the default scheme and step count, a piecewise-constant control, path
    // constraints at the boundaries only and a fixed partition, with no way to change any of
    // them and no sign that there was anything to change.
    I("ms_steps_per_segment", a.ms_steps_per_segment);
    S("ms_integrator", a.ms_integrator);
    S("ms_control_parameterisation", a.ms_control_parameterisation);
    I("ms_path_samples", a.ms_path_samples);
    B("ms_flexible_segments", a.ms_flexible_segments);
    D("ms_min_segment_fraction", a.ms_min_segment_fraction);
    D("ms_refine_tolerance", a.ms_refine_tolerance);
    I("ms_algebraic_iterations", a.ms_algebraic_iterations);
    B("ms_adaptive_steps", a.ms_adaptive_steps);
    I("ms_max_steps_per_segment", a.ms_max_steps_per_segment);
    I("ms_implicit_iterations", a.ms_implicit_iterations);
    // PSOPT's own SQP solver. nlp_method = "SQP" was already reachable from Python and
    // not one of the seven settings that configure it, which is the same defect the
    // multiple-shooting block above describes: a Python user could select the solver and
    // then had no way to choose its QP backend, its strategy, or its trust region, and
    // no sign from the interface that there was anything to choose.
    S("qp_solver", a.qp_solver);
    S("qp_restoration", a.qp_restoration);
    S("sqp_strategy", a.sqp_strategy);
    I("qp_iter_max", a.qp_iter_max);
    S("trust_region", a.trust_region);
    D("trust_region_radius", a.trust_region_radius);
    S("elastic_penalty", a.elastic_penalty);
    // Parameter estimation. The Python interface has carried nobserved, nsamples and the
    // observation function since it was written, so a Python user could pose an
    // estimation problem and get the parameters -- and could not ask for the covariance,
    // the confidence intervals or the residual norm that say whether those parameters
    // mean anything. Both halves are here now: this option asks PSOPT to compute them,
    // and the statistics are returned in the solution.
    S("parameter_statistics", a.parameter_statistics);
    I("parameter_estimation_norm", a.parameter_estimation_norm);
    // The rest of Alg, so that the two interfaces offer the same problem to the solver.
    S("objective_form", a.objective_form);
    S("defect_scaling", a.defect_scaling);
    S("diff_matrix", a.diff_matrix);
    D("jac_sparsity_ratio", a.jac_sparsity_ratio);
    D("hess_sparsity_ratio", a.hess_sparsity_ratio);
    I("save_sparsity_pattern", a.save_sparsity_pattern);
    I("nsteps_error_integration", a.nsteps_error_integration);
    D("mr_kappa", a.mr_kappa);
    I("mr_M1", a.mr_M1);
    I("mr_switch_detection", a.mr_switch_detection);
    I("switch_order", a.switch_order);
    B("hessian_verify", a.hessian_verify);
    S("on_error", a.on_error);
    I("max_integer_combinations", a.max_integer_combinations);
}

// ---- what the solve produced, beyond the trajectory ------------------------
//
// Everything below was reachable from C++ and from nothing in Python. The costates are
// the most conspicuous: PSOPT computes the discrete adjoint of every phase, the examples
// and the book use it to check a solution against the maximum principle, and a Python
// user could not see it. The same is true of the constraint multipliers, of the local
// discretisation error, and -- most importantly -- of whether the solve worked at all.
//
// That last point deserves its own sentence. solution.get_cost() returns a number
// whatever happened, so a Python script that read only the objective could not tell a
// converged solve from one that hit its iteration limit, and had no way to find out.
// nlp_return_code, error_flag and error_msg are now returned on every solve, and the
// Python Solution object turns them into a success flag.
static void pack_status(py::dict& out, Sol& solution)
{
    out["nlp_return_code"] = solution.nlp_return_code;
    out["error_flag"]      = solution.error_flag;
    out["error_msg"]       = solution.error_msg;
    out["cpu_time"]        = solution.cpu_time;
    out["mesh_refinement_iterations"] = solution.mesh_refinement_iterations;

    // One entry per mesh-refinement iteration: the mesh it used, its size, how much work
    // it cost and the discretisation error it achieved. mesh_stats is allocated by the
    // solve and indexed from zero over the iterations actually performed.
    py::list ms;
    if (solution.mesh_stats != NULL && solution.mesh_refinement_iterations > 0) {
        const int n = solution.mesh_refinement_iterations;
        for (int k = 0; k < n; ++k) {
            py::dict d;
            d["method"]           = solution.mesh_stats[k].method;
            d["nnodes"]           = solution.mesh_stats[k].nnodes;
            d["nvars"]            = solution.mesh_stats[k].nvars;
            d["ncons"]            = solution.mesh_stats[k].ncons;
            d["n_obj_evals"]      = solution.mesh_stats[k].n_obj_evals;
            d["n_con_evals"]      = solution.mesh_stats[k].n_con_evals;
            d["n_jacobian_evals"] = solution.mesh_stats[k].n_jacobian_evals;
            d["n_hessian_evals"]  = solution.mesh_stats[k].n_hessian_evals;
            d["n_ode_rhs_evals"]  = solution.mesh_stats[k].n_ode_rhs_evals;
            d["epsilon_max"]      = solution.mesh_stats[k].epsilon_max;
            d["cpu_time"]         = solution.mesh_stats[k].CPU_time;
            ms.append(d);
        }
    }
    out["mesh_stats"] = ms;
}

// The duals and diagnostics of one phase. Each getter returns a reference into the
// solution, so each is copied into a NumPy array on the way out.
static void pack_phase_duals(py::dict& d, Sol& solution, int iphase, int npath, int nevents)
{
    d["costates"]          = Eigen::MatrixXd(solution.get_dual_costates_in_phase(iphase));
    d["hamiltonian"]       = Eigen::MatrixXd(solution.get_dual_hamiltonian_in_phase(iphase));
    d["terminal_state"]    = Eigen::MatrixXd(solution.get_terminal_state_in_phase(iphase));
    d["terminal_costate"]  = Eigen::MatrixXd(solution.get_dual_terminal_costate_in_phase(iphase));
    d["relative_local_error"] =
        Eigen::MatrixXd(solution.get_relative_local_error_in_phase(iphase));
    if (npath   > 0) d["dual_path"]   = Eigen::MatrixXd(solution.get_dual_path_in_phase(iphase));
    if (nevents > 0) d["dual_events"] = Eigen::MatrixXd(solution.get_dual_events_in_phase(iphase));
}

// The statistics of a parameter-estimation solve, when one was asked for and succeeded.
// parameter_statistics_ok is the solver's own verdict and is passed through rather than
// inferred from whether the arrays are empty: a covariance that could not be formed and
// one that happens to be zero are different things.
static void pack_parameter_statistics(py::dict& out, Sol& solution)
{
    if (!solution.parameter_statistics_ok) return;
    py::dict d;
    d["covariance"]       = Eigen::MatrixXd(solution.parameter_covariance);
    d["confidence_low"]   = Eigen::MatrixXd(solution.parameter_confidence_low);
    d["confidence_high"]  = Eigen::MatrixXd(solution.parameter_confidence_high);
    d["residuals"]        = Eigen::MatrixXd(solution.observation_residuals);
    d["sigma_hat"]        = solution.sigma_hat;
    out["parameter_statistics"] = d;
}


// ---- discrete-valued declarations -----------------------------------------
// Reads the integer_controls / integer_parameters entries of a phase spec and
// records them on the phase. Must run after psopt_level2_setup, since both
// declarations modify problem state that level 2 has already sized, and before
// psopt()/psopt_solve_integer, which consume them.
static void declare_discrete(Prob& problem, int iphase, py::dict p)
{
    if (p.contains("integer_controls")) {
        for (auto item : p["integer_controls"].cast<py::list>()) {
            py::dict d = item.cast<py::dict>();
            std::vector<double> v = py::cast<std::vector<double>>(d["values"]);
            RowVectorXd vals((int) v.size());
            for (size_t j = 0; j < v.size(); ++j) vals((int) j) = v[j];
            declare_integer_control(problem, iphase, py::cast<int>(d["index"]), vals);
        }
    }
    if (p.contains("integer_parameters")) {
        for (auto item : p["integer_parameters"].cast<py::list>()) {
            py::dict d = item.cast<py::dict>();
            std::vector<double> v = py::cast<std::vector<double>>(d["values"]);
            RowVectorXd vals((int) v.size());
            for (size_t j = 0; j < v.size(); ++j) vals((int) j) = v[j];
            declare_integer_parameter(problem, iphase, py::cast<int>(d["index"]), vals);
        }
    }
}

static bool any_integer_parameters(Prob& problem)
{
    for (int i = 1; i <= problem.nphases; ++i)
        if (!problem.phases(i).integer_parameters.empty()) return true;
    return false;
}

// Reconstructions, packed for Python. Integer controls are rounded from the relaxed
// weights; integer parameters are read back from the solution that psopt_solve_integer
// selected.
static py::list pack_integer_controls(Sol& solution, Prob& problem, int iphase)
{
    py::list out;
    std::vector<IntegerControlReconstruction> rec =
        reconstruct_integer_controls(solution, problem, iphase);
    for (size_t k = 0; k < rec.size(); ++k) {
        py::dict d;
        d["control"]         = Eigen::MatrixXd(rec[k].control);
        d["mode_index"]      = Eigen::MatrixXi(rec[k].mode_index);
        d["interval_widths"] = Eigen::MatrixXd(rec[k].interval_widths);
        d["integral_gap"]    = rec[k].integral_gap;
        d["n_switches"]      = rec[k].n_switches;
        out.append(d);
    }
    return out;
}

static py::list pack_integer_parameters(Sol& solution, Prob& problem, int iphase)
{
    py::list out;
    const std::vector<IntegerParameter>& ip = problem.phases(iphase).integer_parameters;
    for (size_t k = 0; k < ip.size(); ++k) {
        IntegerParameterReconstruction r =
            reconstruct_integer_parameter(solution, problem, iphase, ip[k].parameter_index);
        py::dict d;
        d["index"]      = ip[k].parameter_index;
        d["value"]      = r.value;
        d["mode_index"] = r.index;
        out.append(d);
    }
    return out;
}

// Weighting of a parameter-estimation problem. Both were reachable from C++ and from
// nothing in Python, so a Python user could only ever pose the unweighted, unregularised
// least-squares problem. An empty residual_weights is left empty rather than filled with
// ones here, because PSOPT does that itself in validate.cxx and there is no reason for
// two places to know the default.
static void set_estimation_weights(phases_str& ph, py::dict p)
{
    if (p.contains("residual_weights")) {
        Eigen::MatrixXd w = py::cast<Eigen::MatrixXd>(p["residual_weights"]);
        if (w.size() > 0) ph.residual_weights = w;
    }
    if (p.contains("regularization_factor"))
        ph.regularization_factor = py::cast<double>(p["regularization_factor"]);
}

static py::dict solve_single_phase(py::dict spec) {
    Alg algorithm; Sol solution; Prob problem;

    problem.outfilename = py::cast<std::string>(spec["outfilename"]);
    problem.nphases     = 1;
    problem.nlinkages   = 0;
    psopt_level1_setup(problem);

    auto& ph = problem.phases(1);
    ph.nstates     = py::cast<int>(spec["nstates"]);
    ph.ncontrols   = py::cast<int>(spec["ncontrols"]);
    ph.nparameters = py::cast<int>(spec["nparameters"]);
    ph.nevents     = py::cast<int>(spec["nevents"]);
    ph.npath       = py::cast<int>(spec["npath"]);
    // The algebraic declaration of a semi-explicit index-1 DAE. A count, like every other
    // size here; it is read only by the multiple-shooting segment integrator.
    if (spec.contains("nalgebraic")) ph.nalgebraic = py::cast<int>(spec["nalgebraic"]);
    ph.nobserved   = py::cast<int>(spec["nobserved"]);
    ph.nsamples    = py::cast<int>(spec["nsamples"]);
    {
        auto nodes = py::cast<std::vector<int>>(spec["nodes"]);
        ph.nodes.resize(1, (int)nodes.size());
        for (size_t i = 0; i < nodes.size(); ++i) ph.nodes(i) = nodes[i];
    }

    psopt_level2_setup(problem, algorithm);

    // ---- discrete-valued declarations (after level 2, before the solve) ----
    declare_discrete(problem, 1, spec);

    // ---- bounds ----
    auto set_vec = [](MatrixXd& dst, const std::vector<double>& v) {
        for (size_t i = 0; i < v.size(); ++i) dst(i) = v[i];
    };
    set_vec(ph.bounds.lower.states,   py::cast<std::vector<double>>(spec["states_lower"]));
    set_vec(ph.bounds.upper.states,   py::cast<std::vector<double>>(spec["states_upper"]));
    if (ph.ncontrols > 0) {
        set_vec(ph.bounds.lower.controls, py::cast<std::vector<double>>(spec["controls_lower"]));
        set_vec(ph.bounds.upper.controls, py::cast<std::vector<double>>(spec["controls_upper"]));
    }
    if (ph.nevents > 0) {
        set_vec(ph.bounds.lower.events, py::cast<std::vector<double>>(spec["events_lower"]));
        set_vec(ph.bounds.upper.events, py::cast<std::vector<double>>(spec["events_upper"]));
    }
    if (ph.npath > 0) {
        set_vec(ph.bounds.lower.path, py::cast<std::vector<double>>(spec["path_lower"]));
        set_vec(ph.bounds.upper.path, py::cast<std::vector<double>>(spec["path_upper"]));
    }
    if (ph.nparameters > 0) {
        set_vec(ph.bounds.lower.parameters, py::cast<std::vector<double>>(spec["parameters_lower"]));
        set_vec(ph.bounds.upper.parameters, py::cast<std::vector<double>>(spec["parameters_upper"]));
    }
    ph.bounds.lower.StartTime = py::cast<double>(spec["t0_lower"]);
    ph.bounds.upper.StartTime = py::cast<double>(spec["t0_upper"]);
    ph.bounds.lower.EndTime   = py::cast<double>(spec["tf_lower"]);
    ph.bounds.upper.EndTime   = py::cast<double>(spec["tf_upper"]);

    // ---- register JIT-compiled user functions ----
    std::string so = py::cast<std::string>(spec["so_path"]);
    void* h = dlopen(so.c_str(), RTLD_NOW | RTLD_GLOBAL);
    if (!h) throw std::runtime_error(std::string("dlopen failed: ") + dlerror());
    problem.dae            = (dae_t)      must_sym(h, "psopt_dae");
    problem.endpoint_cost  = (endpoint_t) must_sym(h, "psopt_endpoint_cost");
    problem.integrand_cost = (integrand_t)must_sym(h, "psopt_integrand_cost");
    problem.events         = (events_t)   must_sym(h, "psopt_events");
    problem.linkages       = (linkages_t) must_sym(h, "psopt_linkages");
    if (ph.nobserved > 0) {
        problem.observation_function = (observation_t) must_sym(h, "psopt_observation");
        ph.observation_nodes = py::cast<Eigen::MatrixXd>(spec["observation_nodes"]);
        ph.observations      = py::cast<Eigen::MatrixXd>(spec["observations"]);
        set_estimation_weights(ph, spec);
    }

    // ---- guess ----
    ph.guess.states   = py::cast<Eigen::MatrixXd>(spec["guess_states"]);
    if (ph.ncontrols > 0) ph.guess.controls = py::cast<Eigen::MatrixXd>(spec["guess_controls"]);
    if (ph.nparameters > 0) ph.guess.parameters = py::cast<Eigen::MatrixXd>(spec["guess_parameters"]);
    ph.guess.time     = py::cast<Eigen::MatrixXd>(spec["guess_time"]);

    // ---- algorithm ----
    apply_algorithm(algorithm, spec["algorithm"].cast<py::dict>());

    // A declared integer parameter cannot be relaxed, so the problem is solved by
    // enumeration over the admissible combinations; otherwise the ordinary driver runs
    // and any integer controls are handled by the outer convexification inside it.
    if (any_integer_parameters(problem))
        (void) psopt_solve_integer(solution, problem, algorithm);
    else
        psopt(solution, problem, algorithm);

    py::dict out;
    out["objective"] = solution.get_cost();
    out["integer_controls"]   = pack_integer_controls(solution, problem, 1);
    out["integer_parameters"] = pack_integer_parameters(solution, problem, 1);
    out["states"]    = Eigen::MatrixXd(solution.get_states_in_phase(1));
    if (ph.ncontrols > 0)
        out["controls"] = Eigen::MatrixXd(solution.get_controls_in_phase(1));
    out["time"]      = Eigen::MatrixXd(solution.get_time_in_phase(1));
    if (ph.nparameters > 0)
        out["parameters"] = Eigen::MatrixXd(solution.get_parameters_in_phase(1));
    pack_status(out, solution);
    pack_phase_duals(out, solution, 1, ph.npath, ph.nevents);
    pack_parameter_statistics(out, solution);
    return out;
}

static void set_phase_bounds(phases_str& ph, py::dict p) {
    auto set_vec = [](MatrixXd& dst, const std::vector<double>& v) {
        for (size_t i = 0; i < v.size(); ++i) dst(i) = v[i];
    };
    set_vec(ph.bounds.lower.states, py::cast<std::vector<double>>(p["states_lower"]));
    set_vec(ph.bounds.upper.states, py::cast<std::vector<double>>(p["states_upper"]));
    if (ph.ncontrols > 0) {
        set_vec(ph.bounds.lower.controls, py::cast<std::vector<double>>(p["controls_lower"]));
        set_vec(ph.bounds.upper.controls, py::cast<std::vector<double>>(p["controls_upper"]));
    }
    if (ph.nevents > 0) {
        set_vec(ph.bounds.lower.events, py::cast<std::vector<double>>(p["events_lower"]));
        set_vec(ph.bounds.upper.events, py::cast<std::vector<double>>(p["events_upper"]));
    }
    if (ph.npath > 0) {
        set_vec(ph.bounds.lower.path, py::cast<std::vector<double>>(p["path_lower"]));
        set_vec(ph.bounds.upper.path, py::cast<std::vector<double>>(p["path_upper"]));
    }
    // Static parameters. These were handled in the single-phase driver and not here, so a
    // multi-phase problem with parameters was set up with whatever bounds level 2 had left
    // in place -- the sizes were right, the numbers were not, and nothing said so. The
    // Python front end allowed the problem to be posed, which is what made it a defect
    // rather than an absent feature.
    if (ph.nparameters > 0) {
        set_vec(ph.bounds.lower.parameters, py::cast<std::vector<double>>(p["parameters_lower"]));
        set_vec(ph.bounds.upper.parameters, py::cast<std::vector<double>>(p["parameters_upper"]));
    }
    ph.bounds.lower.StartTime = py::cast<double>(p["t0_lower"]);
    ph.bounds.upper.StartTime = py::cast<double>(p["t0_upper"]);
    ph.bounds.lower.EndTime   = py::cast<double>(p["tf_lower"]);
    ph.bounds.upper.EndTime   = py::cast<double>(p["tf_upper"]);
}

static py::dict solve_multiphase(py::dict spec) {
    Alg algorithm; Sol solution; Prob problem;
    problem.outfilename = py::cast<std::string>(spec["outfilename"]);
    int N = py::cast<int>(spec["nphases"]);
    problem.nphases   = N;
    problem.nlinkages = py::cast<int>(spec["nlinkages"]);
    psopt_level1_setup(problem);

    py::list phases = spec["phases"];
    for (int k = 1; k <= N; ++k) {
        py::dict p = phases[k - 1].cast<py::dict>();
        auto& ph = problem.phases(k);
        ph.nstates     = py::cast<int>(p["nstates"]);
        ph.ncontrols   = py::cast<int>(p["ncontrols"]);
        ph.nparameters = py::cast<int>(p["nparameters"]);
        ph.nevents     = py::cast<int>(p["nevents"]);
        ph.npath       = py::cast<int>(p["npath"]);
        // Observations, for a multi-phase parameter-estimation problem. Sized here with
        // everything else, because level 2 reads them.
        ph.nobserved   = py::cast<int>(p["nobserved"]);
        ph.nsamples    = py::cast<int>(p["nsamples"]);
        if (p.contains("nalgebraic")) ph.nalgebraic = py::cast<int>(p["nalgebraic"]);
        auto nodes = py::cast<std::vector<int>>(p["nodes"]);
        ph.nodes.resize(1, (int)nodes.size());
        for (size_t i = 0; i < nodes.size(); ++i) ph.nodes(i) = nodes[i];
    }

    psopt_level2_setup(problem, algorithm);

    for (int k = 1; k <= N; ++k)
        declare_discrete(problem, k, phases[k - 1].cast<py::dict>());

    for (int k = 1; k <= N; ++k) {
        py::dict p = phases[k - 1].cast<py::dict>();
        auto& ph = problem.phases(k);
        set_phase_bounds(ph, p);
        ph.guess.states = py::cast<Eigen::MatrixXd>(p["guess_states"]);
        if (ph.ncontrols > 0) ph.guess.controls = py::cast<Eigen::MatrixXd>(p["guess_controls"]);
        if (ph.nparameters > 0)
            ph.guess.parameters = py::cast<Eigen::MatrixXd>(p["guess_parameters"]);
        ph.guess.time = py::cast<Eigen::MatrixXd>(p["guess_time"]);
    }

    std::string so = py::cast<std::string>(spec["so_path"]);
    void* h = dlopen(so.c_str(), RTLD_NOW | RTLD_GLOBAL);
    if (!h) throw std::runtime_error(std::string("dlopen failed: ") + dlerror());
    problem.dae            = (dae_t)      must_sym(h, "psopt_dae");
    problem.endpoint_cost  = (endpoint_t) must_sym(h, "psopt_endpoint_cost");
    problem.integrand_cost = (integrand_t)must_sym(h, "psopt_integrand_cost");
    problem.events         = (events_t)   must_sym(h, "psopt_events");
    problem.linkages       = (linkages_t) must_sym(h, "psopt_linkages");
    {
        bool any_observed = false;
        for (int k = 1; k <= N; ++k) if (problem.phases(k).nobserved > 0) any_observed = true;
        if (any_observed) {
            problem.observation_function = (observation_t) must_sym(h, "psopt_observation");
            for (int k = 1; k <= N; ++k) {
                py::dict p = phases[k - 1].cast<py::dict>();
                auto& ph = problem.phases(k);
                if (ph.nobserved > 0) {
                    ph.observation_nodes = py::cast<Eigen::MatrixXd>(p["observation_nodes"]);
                    ph.observations      = py::cast<Eigen::MatrixXd>(p["observations"]);
                    set_estimation_weights(ph, p);
                }
            }
        }
    }

    apply_algorithm(algorithm, spec["algorithm"].cast<py::dict>());

    if (any_integer_parameters(problem))
        (void) psopt_solve_integer(solution, problem, algorithm);
    else
        psopt(solution, problem, algorithm);

    py::dict out;
    out["objective"] = solution.get_cost();
    py::list st, ct, tm, ic, ip, pr, du;
    for (int k = 1; k <= N; ++k) {
        auto& ph = problem.phases(k);
        st.append(Eigen::MatrixXd(solution.get_states_in_phase(k)));
        ct.append(Eigen::MatrixXd(solution.get_controls_in_phase(k)));
        tm.append(Eigen::MatrixXd(solution.get_time_in_phase(k)));
        ic.append(pack_integer_controls(solution, problem, k));
        ip.append(pack_integer_parameters(solution, problem, k));
        pr.append(ph.nparameters > 0
                  ? py::object(py::cast(Eigen::MatrixXd(solution.get_parameters_in_phase(k))))
                  : py::none());
        py::dict d;
        pack_phase_duals(d, solution, k, ph.npath, ph.nevents);
        du.append(d);
    }
    out["states"] = st; out["controls"] = ct; out["time"] = tm;
    out["integer_controls"] = ic; out["integer_parameters"] = ip;
    out["parameters"] = pr;
    out["duals"] = du;
    if (problem.nlinkages > 0)
        out["dual_linkages"] = Eigen::MatrixXd(solution.get_dual_linkages());
    pack_status(out, solution);
    pack_parameter_statistics(out, solution);
    return out;
}

PYBIND11_MODULE(_psopt, m) {
    m.doc() = "PSOPT Python binding (B-1 single-phase, B-2 multi-phase)";
    m.def("solve_single_phase", &solve_single_phase,
          "Set up and solve a single-phase PSOPT problem from a spec dict.");
    m.def("solve_multiphase", &solve_multiphase,
          "Set up and solve a multi-phase PSOPT problem (with linkages) from a spec dict.");
}
