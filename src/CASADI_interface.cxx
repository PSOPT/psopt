// CASADI_interface.cxx - CasADi NLP backend for PSOPT (guarded by USE_CASADI).
//
// Solves PSOPT's transcribed NLP through casadi::nlpsol(<plugin>, ...), so the
// solve can use the IPOPT plugin (with MUMPS or a parallel HSL linear solver)
// or alternative plugins CasADi exposes (sqpmethod, fatrop, ...) that PSOPT has
// no native path to. Selected by algorithm.nlp_method=="CASADI" +
// algorithm.casadi_solver (default "ipopt").
//
// Derivatives: by default (algorithm.derivatives=="automatic") the objective
// gradient and constraint Jacobian are PSOPT's EXACT ADOL-C derivatives, driven
// directly with the Workspace (ff_ad/gg_ad taped on dedicated tags) and handed
// to CasADi via Callback::get_jacobian. This removes the finite-difference
// penalty that made the backend 5-25x slower and time out on large problems
// (bryson_denham: 24s FD -> 3s exact; obstacle 55s -> 11s; both match IPOPT).
// With algorithm.derivatives=="numerical" the FD PsoptNlp path is used instead.
// The Hessian is limited-memory BFGS built on these exact gradients.
//
// See casadi_psopt.h for scope. Requires a CasADi build with the requested
// plugin (the PSOPT superbuild provides IPOPT).

#include "casadi_psopt.h"

#ifdef USE_CASADI

#include <casadi/casadi.hpp>
#include <adolc/adolc.h>
#include <vector>

using namespace casadi;

void get_constraint_bounds(double* g_l, double* g_u, Workspace* workspace);
// PSOPT's objective/constraints in adoubles (from NLP_objective/NLP_constraints).
adouble ff_ad(adouble* xad, Workspace* workspace);
void    gg_ad(adouble* xad, adouble* gad, Workspace* workspace);

namespace {

// Dedicated ADOL-C tape tags for the CasADi backend, kept clear of PSOPT's
// tags (tag_f=1..tag_gc=5) so re-used sparse patterns are never clobbered.
static const int TAG_F = 21, TAG_G = 22;

// ---------------------------------------------------------------------------
// Exact derivatives via PSOPT's ADOL-C tapes, driven DIRECTLY with the
// Workspace (not through a standalone IPOPT_PSOPT, whose get_nlp_info crashed
// out of the Ipopt lifecycle). Removes the finite-difference penalty that made
// the CasADi backend 5-25x slower and time out on large problems.
// ---------------------------------------------------------------------------

// Exact objective gradient d f / d x (1 x n). Tapes ff_ad once on tag_f.
class GradF : public Callback {
public:
  GradF(int n, Workspace* ws) : n_(n), ws_(ws) {
    // Tape the objective once; gradient() re-evaluates it at new x.
    trace_on(TAG_F);
    for (int i = 0; i < n_; ++i) ws_->xad[i] <<= (*ws_->x0)(i, 0);
    adouble fad = ff_ad(ws_->xad, ws_);
    double fval; fad >>= fval;
    trace_off();
    construct("psopt_grad_f");
  }
  casadi_int get_n_in()  override { return 2; }               // x, nominal f
  casadi_int get_n_out() override { return 1; }               // grad (1 x n)
  Sparsity get_sparsity_in(casadi_int i) override {
    return (i == 0) ? Sparsity::dense(n_, 1) : Sparsity::dense(1, 1); }
  Sparsity get_sparsity_out(casadi_int) override { return Sparsity::dense(1, n_); }
  std::vector<DM> eval(const std::vector<DM>& arg) const override {
    std::vector<double> x = arg[0].nonzeros();
    std::vector<double> gf(n_);
    gradient(TAG_F, n_, x.data(), gf.data());
    return { DM(std::vector<double>(gf.begin(), gf.end())).T() };  // row vector
  }
private:
  int n_; Workspace* ws_;
};

// Exact constraint Jacobian d g / d x (m x n, sparse). Tapes gg_ad on tag_g and
// detects the sparsity once; sparse_jac re-evaluates values at new x.
class JacG : public Callback {
public:
  JacG(int n, int m, Workspace* ws) : n_(n), m_(m), ws_(ws) {
    std::vector<double> x0(n_);
    for (int i = 0; i < n_; ++i) x0[i] = (*ws_->x0)(i, 0);
    trace_on(TAG_G);
    for (int i = 0; i < n_; ++i) ws_->xad[i] <<= x0[i];
    gg_ad(ws_->xad, ws_->gad, ws_);
    std::vector<double> g(m_);
    for (int i = 0; i < m_; ++i) ws_->gad[i] >>= g[i];
    trace_off();
    // sparsity pattern (repeat=0) + first values
    int opt[4] = {0, 0, 0, 0};
    unsigned int* rind = nullptr; unsigned int* cind = nullptr; double* vals = nullptr;
    int nnz = 0;
    sparse_jac(TAG_G, m_, n_, 0, x0.data(), &nnz, &rind, &cind, &vals, opt);
    nnz_ = nnz;
    rows_.assign(rind, rind + nnz); cols_.assign(cind, cind + nnz);
    sp_ = Sparsity::triplet(m_, n_, std::vector<casadi_int>(rows_.begin(), rows_.end()),
                                    std::vector<casadi_int>(cols_.begin(), cols_.end()));
    free(rind); free(cind); free(vals);
    construct("psopt_jac_g");
  }
  casadi_int get_n_in()  override { return 2; }               // x, nominal g
  casadi_int get_n_out() override { return 1; }               // jac (m x n)
  Sparsity get_sparsity_in(casadi_int i) override {
    return (i == 0) ? Sparsity::dense(n_, 1) : Sparsity::dense(m_, 1); }
  Sparsity get_sparsity_out(casadi_int) override { return sp_; }
  std::vector<DM> eval(const std::vector<DM>& arg) const override {
    std::vector<double> x = arg[0].nonzeros();
    int opt[4] = {0, 0, 0, 0};
    unsigned int* rind = nullptr; unsigned int* cind = nullptr; double* vals = nullptr;
    int nnz = nnz_;   // repeat=1 requires the previously detected nonzero count
    sparse_jac(TAG_G, m_, n_, 1, x.data(), &nnz, &rind, &cind, &vals, opt);
    // Build via triplet so CasADi places values into its column-major (CCS)
    // nonzero order; repeat=1 keeps (rind,cind) == (rows_,cols_).
    std::vector<casadi_int> R(rows_.begin(), rows_.end()), C(cols_.begin(), cols_.end());
    std::vector<double> V(vals, vals + nnz);
    free(rind); free(cind); free(vals);
    return { DM::triplet(R, C, V, m_, n_) };
  }
private:
  int n_, m_; Workspace* ws_;
  int nnz_ = 0;
  std::vector<int> rows_, cols_;
  Sparsity sp_;
};

// Objective f(x) with an exact Jacobian (GradF). Held as a member so it outlives
// the solve; get_jacobian hands it to CasADi.
class ObjFun : public Callback {
public:
  ObjFun(int n, Workspace* ws,
         double (*f)(MatrixXd&, Workspace*)) : n_(n), ws_(ws), f_(f), gradf_(n, ws) {
    construct("psopt_obj");
  }
  casadi_int get_n_in()  override { return 1; }
  casadi_int get_n_out() override { return 1; }
  Sparsity get_sparsity_in(casadi_int)  override { return Sparsity::dense(n_, 1); }
  Sparsity get_sparsity_out(casadi_int) override { return Sparsity::dense(1, 1); }
  std::vector<DM> eval(const std::vector<DM>& arg) const override {
    std::vector<double> xv = arg[0].nonzeros();
    MatrixXd x(n_, 1); for (int i = 0; i < n_; ++i) x(i, 0) = xv[i];
    return { DM(f_(x, ws_)) };
  }
  bool has_jacobian() const override { return true; }
  Function get_jacobian(const std::string& name, const std::vector<std::string>&,
                        const std::vector<std::string>&, const Dict&) const override {
    return gradf_;
  }
private:
  int n_; Workspace* ws_;
  double (*f_)(MatrixXd&, Workspace*);
  GradF gradf_;
};

// Constraints g(x) with an exact Jacobian (JacG).
class ConFun : public Callback {
public:
  ConFun(int n, int m, Workspace* ws,
         void (*g)(MatrixXd&, MatrixXd*, Workspace*)) : n_(n), m_(m), ws_(ws), g_(g), jacg_(n, m, ws) {
    construct("psopt_con");
  }
  casadi_int get_n_in()  override { return 1; }
  casadi_int get_n_out() override { return 1; }
  Sparsity get_sparsity_in(casadi_int)  override { return Sparsity::dense(n_, 1); }
  Sparsity get_sparsity_out(casadi_int) override { return Sparsity::dense(m_, 1); }
  std::vector<DM> eval(const std::vector<DM>& arg) const override {
    std::vector<double> xv = arg[0].nonzeros();
    MatrixXd x(n_, 1); for (int i = 0; i < n_; ++i) x(i, 0) = xv[i];
    MatrixXd gv(m_, 1); g_(x, &gv, ws_);
    std::vector<double> out(m_); for (int i = 0; i < m_; ++i) out[i] = gv(i, 0);
    return { DM(out) };
  }
  bool has_jacobian() const override { return true; }
  Function get_jacobian(const std::string& name, const std::vector<std::string>&,
                        const std::vector<std::string>&, const Dict&) const override {
    return jacg_;
  }
private:
  int n_, m_; Workspace* ws_;
  void (*g_)(MatrixXd&, MatrixXd*, Workspace*);
  JacG jacg_;
};

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

  // Exact ADOL-C derivatives (default) vs finite differences ("numerical").
  const bool exact = (algorithm.derivatives != "numerical");

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

  // Build the NLP. Exact path: wrap PSOPT's objective/constraints as CasADi
  // Functions with exact Jacobians (ObjFun/ConFun) and assemble an MX NLP so
  // CasADi/IPOPT get analytic first derivatives. FD path: the PsoptNlp callback.
  Function solver;
  ObjFun*   objp = nullptr;
  ConFun*   conp = nullptr;
  PsoptNlp* fdp  = nullptr;
  if (exact) {
    objp = new ObjFun(nx, workspace, f);
    if (ng > 0) conp = new ConFun(nx, ng, workspace, g);
    MX x = MX::sym("x", nx);
    MX fx = objp->operator()(std::vector<MX>{x}).at(0);
    MXDict nlp; nlp["x"] = x; nlp["f"] = fx;
    if (ng > 0) nlp["g"] = conp->operator()(std::vector<MX>{x}).at(0);
    solver = nlpsol("psopt_casadi_solver", plugin, nlp, opts);
  } else {
    fdp = new PsoptNlp(nx, ng, f, g, workspace);
    solver = nlpsol("psopt_casadi_solver", plugin, *fdp, opts);
  }

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
  delete objp; delete conp; delete fdp;
  return 0;
}

#endif // USE_CASADI
