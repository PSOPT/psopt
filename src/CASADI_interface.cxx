// CASADI_interface.cxx - CasADi NLP backend for PSOPT (guarded by USE_CASADI).
//
// Solves PSOPT's transcribed NLP through casadi::nlpsol(<plugin>, ...), so the
// solve can use the IPOPT plugin (with MUMPS or a parallel HSL linear solver)
// or alternative plugins CasADi exposes (sqpmethod, fatrop, ...) that PSOPT has
// no native path to. Selected by algorithm.nlp_method=="CASADI" +
// algorithm.casadi_solver (default "ipopt").
//
// Derivatives: CasADi finite-differences the PSOPT objective/constraint
// callbacks (enable_fd). PSOPT's own f/g still honour algorithm.derivatives
// (analytical ADOL-C tapes vs numerical) internally, so this reuses PSOPT's
// transcription; the CasADi layer only adds the outer FD for the NLP Jacobian.
// TODO (follow-up): supply exact first derivatives by wrapping PSOPT's
// IPOPT_PSOPT::eval_grad_f/eval_jac_g via Callback::get_jacobian (needs the
// TNLP driven through its normal Ipopt lifecycle, not called standalone).
//
// See casadi_psopt.h for scope. Requires a CasADi build with the requested
// plugin (the PSOPT superbuild provides IPOPT).

#include "casadi_psopt.h"

#ifdef USE_CASADI

#include <casadi/casadi.hpp>
#include <vector>

using namespace casadi;

void get_constraint_bounds(double* g_l, double* g_u, Workspace* workspace);

namespace {

// PSOPT NLP exposed to CasADi: x -> (f, g). Derivatives via CasADi FD.
class PsoptNlp : public Callback {
public:
  PsoptNlp(int nx, int ng,
           double (*f)(MatrixXd&, Workspace*),
           void   (*g)(MatrixXd&, MatrixXd*, Workspace*),
           Workspace* ws)
    : nx_(nx), ng_(ng), f_(f), g_(g), ws_(ws) {
    Dict opts; opts["enable_fd"] = true;
    construct("psopt_nlp", opts);
  }
  casadi_int get_n_in()  override { return 2; }               // x, p
  casadi_int get_n_out() override { return 2; }               // f, g
  Sparsity get_sparsity_in(casadi_int i)  override { return (i == 0) ? Sparsity::dense(nx_, 1) : Sparsity::dense(0, 1); }
  Sparsity get_sparsity_out(casadi_int i) override { return (i == 0) ? Sparsity::dense(1, 1)  : Sparsity::dense(ng_, 1); }
  std::string get_name_in(casadi_int i)  override { return (i == 0) ? "x" : "p"; }
  std::string get_name_out(casadi_int i) override { return (i == 0) ? "f" : "g"; }
  std::vector<DM> eval(const std::vector<DM>& arg) const override {
    const std::vector<double> xv = arg[0].nonzeros();
    MatrixXd x(nx_, 1);
    for (int j = 0; j < nx_; ++j) x(j, 0) = xv[j];
    double fval = f_(x, ws_);
    std::vector<double> gv(ng_);
    if (ng_ > 0) { MatrixXd gvec(ng_, 1); g_(x, &gvec, ws_); for (int j = 0; j < ng_; ++j) gv[j] = gvec(j, 0); }
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

  PsoptNlp nlp(nx, ng, f, g, workspace);

  std::string plugin = algorithm.casadi_solver.empty() ? std::string("ipopt")
                                                        : algorithm.casadi_solver;
  Dict opts;
  if (plugin == "ipopt") {
    opts["ipopt.linear_solver"]  = algorithm.ipopt_linear_solver;   // mumps / ma97 / ...
    opts["ipopt.tol"]            = algorithm.nlp_tolerance;
    opts["ipopt.max_iter"]       = algorithm.nlp_iter_max;
    opts["ipopt.mu_strategy"]    = std::string("adaptive");
    opts["ipopt.hessian_approximation"] = std::string("limited-memory");
    if (algorithm.print_level == 0) opts["ipopt.print_level"] = 0;
  } else if (plugin == "sqpmethod") {
    // sqpmethod defaults to an EXACT Hessian of the Lagrangian, which is
    // indefinite for most optimal-control NLPs; qpOASES then fails on the
    // non-convex QP subproblem and (with error_on_fail=true) aborts. Use a
    // limited-memory (BFGS, positive-definite) Hessian so the QP stays convex,
    // and tolerate an occasional QP failure instead of crashing.
    opts["hessian_approximation"] = std::string("limited-memory");
    opts["max_iter"]              = algorithm.nlp_iter_max;
    opts["error_on_fail"]         = false;
    Dict qpopts; qpopts["error_on_fail"] = false;
    if (algorithm.print_level == 0) { qpopts["printLevel"] = std::string("none"); opts["print_time"] = false; opts["print_iteration"] = false; }
    opts["qpsol"]         = std::string("qpoases");
    opts["qpsol_options"] = qpopts;
  } else {
    opts["error_on_fail"] = false;
    if (algorithm.print_level == 0) opts["print_time"] = false;
  }

  Function solver = nlpsol("psopt_casadi_solver", plugin, nlp, opts);

  std::vector<double> x0v(nx), lbx(nx), ubx(nx), lbg(ng), ubg(ng);
  for (int j = 0; j < nx; ++j) { x0v[j] = (*x0)(j, 0); lbx[j] = (*xlb)(j, 0); ubx[j] = (*xub)(j, 0); }
  if (ng > 0) get_constraint_bounds(lbg.data(), ubg.data(), workspace);

  DMDict arg;
  arg["x0"]  = DM(x0v);
  arg["lbx"] = DM(lbx);
  arg["ubx"] = DM(ubx);
  if (ng > 0) { arg["lbg"] = DM(lbg); arg["ubg"] = DM(ubg); }

  DMDict res = solver(arg);

  const std::vector<double> xopt = res.at("x").nonzeros();
  for (int j = 0; j < nx; ++j) (*x0)(j, 0) = xopt[j];
  if (ng > 0 && lambda != nullptr) {
    const std::vector<double> lg = res.at("lam_g").nonzeros();
    for (int j = 0; j < ng && j < (int)lg.size(); ++j) (*lambda)(j, 0) = lg[j];
  }
  return 0;
}

#endif // USE_CASADI
