// CASADI_interface.cxx - CasADi NLP backend for PSOPT (guarded by USE_CASADI).
//
// Wraps PSOPT's objective/constraint callbacks as a CasADi Callback and solves
// with casadi::nlpsol("ipopt", ...). The IPOPT linear solver is taken from
// algorithm.ipopt_linear_solver, so a CasADi/IPOPT built with HSL uses ma57/97.
//
// See casadi_psopt.h for scope and limitations.

#include "casadi_psopt.h"

#ifdef USE_CASADI

#include <casadi/casadi.hpp>
#include <vector>

using namespace casadi;

// Forward declaration of PSOPT's constraint-bounds helper (defined in the
// PSOPT sources and reused by the IPOPT interface).
void get_constraint_bounds(double* g_l, double* g_u, Workspace* workspace);

namespace {

// A CasADi Callback exposing PSOPT's NLP: x -> (f, g). Derivatives are obtained
// by CasADi finite differences (enable_fd), so no analytic Jacobian is needed.
class PsoptNlp : public Callback {
public:
  PsoptNlp(int nx, int ng,
           double (*f)(MatrixXd&, Workspace*),
           void   (*g)(MatrixXd&, MatrixXd*, Workspace*),
           Workspace* ws)
    : nx_(nx), ng_(ng), f_(f), g_(g), ws_(ws) {
    Dict opts;
    opts["enable_fd"] = true;          // finite-difference the callback
    construct("psopt_nlp", opts);
  }

  // inputs: x (nx), p (0)  ; outputs: f (1), g (ng)
  casadi_int get_n_in()  override { return 2; }
  casadi_int get_n_out() override { return 2; }
  Sparsity get_sparsity_in(casadi_int i) override {
    return (i == 0) ? Sparsity::dense(nx_, 1) : Sparsity::dense(0, 1);
  }
  Sparsity get_sparsity_out(casadi_int i) override {
    return (i == 0) ? Sparsity::dense(1, 1) : Sparsity::dense(ng_, 1);
  }
  std::string get_name_in(casadi_int i)  override { return (i == 0) ? "x" : "p"; }
  std::string get_name_out(casadi_int i) override { return (i == 0) ? "f" : "g"; }

  std::vector<DM> eval(const std::vector<DM>& arg) const override {
    const std::vector<double> xv = arg[0].nonzeros();
    MatrixXd x(nx_, 1);
    for (int j = 0; j < nx_; ++j) x(j, 0) = xv[j];

    double fval = f_(x, ws_);

    MatrixXd gvec(ng_, 1);
    if (ng_ > 0) g_(x, &gvec, ws_);
    std::vector<double> gv(ng_);
    for (int j = 0; j < ng_; ++j) gv[j] = gvec(j, 0);

    return { DM(fval), DM(gv) };
  }

private:
  int nx_, ng_;
  double (*f_)(MatrixXd&, Workspace*);
  void   (*g_)(MatrixXd&, MatrixXd*, Workspace*);
  Workspace* ws_;
};

} // namespace

int psopt_casadi_solve(
    Alg&       algorithm,
    MatrixXd*  x0,
    double   (*f)(MatrixXd&, Workspace*),
    void     (*g)(MatrixXd&, MatrixXd*, Workspace*),
    int        nlp_ncons,
    int        /*nlp_neq*/,
    MatrixXd*  xlb,
    MatrixXd*  xub,
    MatrixXd*  lambda,
    Workspace* workspace,
    void*      /*user_data*/)
{
  const int nx = workspace->nvars;
  const int ng = nlp_ncons;

  // Build the CasADi NLP function from PSOPT's callbacks.
  PsoptNlp nlp(nx, ng, f, g, workspace);

  // IPOPT options mirrored from the algorithm structure.
  Dict opts;
  opts["ipopt.linear_solver"]  = algorithm.ipopt_linear_solver;   // mumps / ma57 / ...
  opts["ipopt.tol"]            = algorithm.nlp_tolerance;
  opts["ipopt.max_iter"]       = algorithm.nlp_iter_max;
  opts["ipopt.mu_strategy"]    = std::string("adaptive");
  opts["ipopt.hessian_approximation"] = std::string("limited-memory");
  if (algorithm.print_level == 0) opts["ipopt.print_level"] = 0;

  Function solver = nlpsol("psopt_casadi_solver", "ipopt", nlp, opts);

  // Bounds and initial guess.
  std::vector<double> x0v(nx), lbx(nx), ubx(nx), lbg(ng), ubg(ng);
  for (int j = 0; j < nx; ++j) {
    x0v[j] = (*x0)(j, 0);
    lbx[j] = (*xlb)(j, 0);
    ubx[j] = (*xub)(j, 0);
  }
  if (ng > 0) get_constraint_bounds(lbg.data(), ubg.data(), workspace);

  DMDict arg;
  arg["x0"]  = DM(x0v);
  arg["lbx"] = DM(lbx);
  arg["ubx"] = DM(ubx);
  if (ng > 0) { arg["lbg"] = DM(lbg); arg["ubg"] = DM(ubg); }

  DMDict res = solver(arg);

  // Write the solution back.
  const std::vector<double> xopt = res.at("x").nonzeros();
  for (int j = 0; j < nx; ++j) (*x0)(j, 0) = xopt[j];
  if (ng > 0 && lambda != nullptr) {
    const std::vector<double> lg = res.at("lam_g").nonzeros();
    for (int j = 0; j < ng && j < (int)lg.size(); ++j) (*lambda)(j, 0) = lg[j];
  }
  return 0;
}

#endif // USE_CASADI
