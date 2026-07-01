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
    ph.nobserved   = py::cast<int>(spec["nobserved"]);
    ph.nsamples    = py::cast<int>(spec["nsamples"]);
    {
        auto nodes = py::cast<std::vector<int>>(spec["nodes"]);
        ph.nodes.resize(1, (int)nodes.size());
        for (size_t i = 0; i < nodes.size(); ++i) ph.nodes(i) = nodes[i];
    }

    psopt_level2_setup(problem, algorithm);

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
    }

    // ---- guess ----
    ph.guess.states   = py::cast<Eigen::MatrixXd>(spec["guess_states"]);
    if (ph.ncontrols > 0) ph.guess.controls = py::cast<Eigen::MatrixXd>(spec["guess_controls"]);
    if (ph.nparameters > 0) ph.guess.parameters = py::cast<Eigen::MatrixXd>(spec["guess_parameters"]);
    ph.guess.time     = py::cast<Eigen::MatrixXd>(spec["guess_time"]);

    // ---- algorithm ----
    apply_algorithm(algorithm, spec["algorithm"].cast<py::dict>());

    psopt(solution, problem, algorithm);

    py::dict out;
    out["objective"] = solution.get_cost();
    out["states"]    = Eigen::MatrixXd(solution.get_states_in_phase(1));
    if (ph.ncontrols > 0)
        out["controls"] = Eigen::MatrixXd(solution.get_controls_in_phase(1));
    out["time"]      = Eigen::MatrixXd(solution.get_time_in_phase(1));
    if (ph.nparameters > 0)
        out["parameters"] = Eigen::MatrixXd(solution.get_parameters_in_phase(1));
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
        auto nodes = py::cast<std::vector<int>>(p["nodes"]);
        ph.nodes.resize(1, (int)nodes.size());
        for (size_t i = 0; i < nodes.size(); ++i) ph.nodes(i) = nodes[i];
    }

    psopt_level2_setup(problem, algorithm);

    for (int k = 1; k <= N; ++k) {
        py::dict p = phases[k - 1].cast<py::dict>();
        auto& ph = problem.phases(k);
        set_phase_bounds(ph, p);
        ph.guess.states = py::cast<Eigen::MatrixXd>(p["guess_states"]);
        if (ph.ncontrols > 0) ph.guess.controls = py::cast<Eigen::MatrixXd>(p["guess_controls"]);
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

    apply_algorithm(algorithm, spec["algorithm"].cast<py::dict>());

    psopt(solution, problem, algorithm);

    py::dict out;
    out["objective"] = solution.get_cost();
    py::list st, ct, tm;
    for (int k = 1; k <= N; ++k) {
        st.append(Eigen::MatrixXd(solution.get_states_in_phase(k)));
        ct.append(Eigen::MatrixXd(solution.get_controls_in_phase(k)));
        tm.append(Eigen::MatrixXd(solution.get_time_in_phase(k)));
    }
    out["states"] = st; out["controls"] = ct; out["time"] = tm;
    return out;
}

PYBIND11_MODULE(_psopt, m) {
    m.doc() = "PSOPT Python binding (B-1 single-phase, B-2 multi-phase)";
    m.def("solve_single_phase", &solve_single_phase,
          "Set up and solve a single-phase PSOPT problem from a spec dict.");
    m.def("solve_multiphase", &solve_multiphase,
          "Set up and solve a multi-phase PSOPT problem (with linkages) from a spec dict.");
}
