// CASADI_interface.cxx - CasADi NLP backend for PSOPT (guarded by USE_CASADI).
//
// Solves PSOPT's transcribed NLP through casadi::nlpsol("ipopt", ...), so the
// IPOPT linear solver (algorithm.ipopt_linear_solver) can be MUMPS or a
// parallel HSL code (ma86/ma97). EXACT first derivatives are supplied to CasADi
// by reusing PSOPT's own IPOPT_PSOPT TNLP (its eval_grad_f / eval_jac_g compute
// the objective gradient and the sparse constraint Jacobian via ADOL-C), so no
// finite differences are used. The Hessian is limited-memory (as in the native
// IPOPT path), so second derivatives are not needed.
//
// See casadi_psopt.h for scope. Requires a CasADi build with the IPOPT plugin
// (provided by -DPSOPT_SUPERBUILD=ON -DPSOPT_WITH_CASADI=ON).

#include "casadi_psopt.h"

#ifdef USE_CASADI

#include <casadi/casadi.hpp>
#include <vector>

using namespace casadi;
using Ipopt::Index;
using Ipopt::Number;

void get_constraint_bounds(double* g_l, double* g_u, Workspace* workspace);

namespace {

// Exact objective gradient d f / d x (shape 1 x n) via IPOPT_PSOPT::eval_grad_f.
class GradF : public Callback {
public:
  GradF(int n, IPOPT_PSOPT* nlp) : n_(n), nlp_(nlp) { construct("psopt_grad_f"); }
  casadi_int get_n_in()  override { return 2; }               // x, nominal f
  casadi_int get_n_out() override { return 1; }               // jac (1 x n)
  Sparsity get_sparsity_in(casadi_int i) override {
    return (i == 0) ? Sparsity::dense(n_, 1) : Sparsity::dense(1, 1);
  }
  Sparsity get_sparsity_out(casadi_int) override { return Sparsity::dense(1, n_); }
  std::vector<DM> eval(const std::vector<DM>& arg) const override {
    std::vector<double> x = arg[0].nonzeros();
    std::vector<double> gf(n_, 0.0);
    nlp_->eval_grad_f((Index)n_, x.data(), true, gf.data());
    return { DM(reshape(DM(gf), 1, n_)) };
  }
private:
  int n_; IPOPT_PSOPT* nlp_;
};

// Exact constraint Jacobian d g / d x (shape m x n, sparse) via eval_jac_g.
class JacG : public Callback {
public:
  JacG(int n, int m, IPOPT_PSOPT* nlp,
       std::vector<casadi_int> rows, std::vector<casadi_int> cols)
    : n_(n), m_(m), nlp_(nlp), rows_(std::move(rows)), cols_(std::move(cols)) {
    sp_ = Sparsity::triplet(m_, n_, rows_, cols_);
    construct("psopt_jac_g");
  }
  casadi_int get_n_in()  override { return 2; }               // x, nominal g
  casadi_int get_n_out() override { return 1; }               // jac (m x n)
  Sparsity get_sparsity_in(casadi_int i) override {
    return (i == 0) ? Sparsity::dense(n_, 1) : Sparsity::dense(m_, 1);
  }
  Sparsity get_sparsity_out(casadi_int) override { return sp_; }
  std::vector<DM> eval(const std::vector<DM>& arg) const override {
    std::vector<double> x = arg[0].nonzeros();
    const int nnz = (int)rows_.size();
    std::vector<double> vals(nnz, 0.0);
    nlp_->eval_jac_g((Index)n_, x.data(), true, (Index)m_, (Index)nnz,
                     nullptr, nullptr, vals.data());
    // assemble in the sparsity's own nonzero order via triplet construction
    return { DM::triplet(rows_, cols_, std::vector<double>(vals), m_, n_) };
  }
private:
  int n_, m_; IPOPT_PSOPT* nlp_;
  std::vector<casadi_int> rows_, cols_;
  Sparsity sp_;
};

// Objective f(x) with an exact Jacobian (GradF held as a member so it outlives
// the solve; get_jacobian returns it as a Function sharing the same node).
class ObjFun : public Callback {
public:
  ObjFun(int n, IPOPT_PSOPT* nlp,
         double (*f)(MatrixXd&, Workspace*), Workspace* ws)
    : n_(n), nlp_(nlp), f_(f), ws_(ws), gradf_(n, nlp) { construct("psopt_f"); }
  casadi_int get_n_in()  override { return 1; }
  casadi_int get_n_out() override { return 1; }
  Sparsity get_sparsity_in(casadi_int)  override { return Sparsity::dense(n_, 1); }
  Sparsity get_sparsity_out(casadi_int) override { return Sparsity::dense(1, 1); }
  std::string get_name_in(casadi_int)  override { return "x"; }
  std::string get_name_out(casadi_int) override { return "f"; }
  std::vector<DM> eval(const std::vector<DM>& arg) const override {
    std::vector<double> xv = arg[0].nonzeros();
    MatrixXd x(n_, 1);
    for (int j = 0; j < n_; ++j) x(j, 0) = xv[j];
    return { DM(f_(x, ws_)) };
  }
  bool has_jacobian() const override { return true; }
  Function get_jacobian(const std::string&, const std::vector<std::string>&,
                        const std::vector<std::string>&, const Dict&) const override {
    return gradf_;
  }
private:
  int n_; IPOPT_PSOPT* nlp_;
  double (*f_)(MatrixXd&, Workspace*); Workspace* ws_;
  GradF gradf_;
};

// Constraints g(x) with an exact Jacobian (JacG held as a member).
class ConFun : public Callback {
public:
  ConFun(int n, int m, IPOPT_PSOPT* nlp,
         void (*g)(MatrixXd&, MatrixXd*, Workspace*), Workspace* ws,
         std::vector<casadi_int> rows, std::vector<casadi_int> cols)
    : n_(n), m_(m), nlp_(nlp), g_(g), ws_(ws),
      jacg_(n, m, nlp, rows, cols) { construct("psopt_g"); }
  casadi_int get_n_in()  override { return 1; }
  casadi_int get_n_out() override { return 1; }
  Sparsity get_sparsity_in(casadi_int)  override { return Sparsity::dense(n_, 1); }
  Sparsity get_sparsity_out(casadi_int) override { return Sparsity::dense(m_, 1); }
  std::string get_name_in(casadi_int)  override { return "x"; }
  std::string get_name_out(casadi_int) override { return "g"; }
  std::vector<DM> eval(const std::vector<DM>& arg) const override {
    std::vector<double> xv = arg[0].nonzeros();
    MatrixXd x(n_, 1);
    for (int j = 0; j < n_; ++j) x(j, 0) = xv[j];
    MatrixXd gvec(m_, 1);
    if (m_ > 0) g_(x, &gvec, ws_);
    std::vector<double> gv(m_);
    for (int j = 0; j < m_; ++j) gv[j] = gvec(j, 0);
    return { DM(gv) };
  }
  bool has_jacobian() const override { return true; }
  Function get_jacobian(const std::string&, const std::vector<std::string>&,
                        const std::vector<std::string>&, const Dict&) const override {
    return jacg_;
  }
private:
  int n_, m_; IPOPT_PSOPT* nlp_;
  void (*g_)(MatrixXd&, MatrixXd*, Workspace*); Workspace* ws_;
  JacG jacg_;
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
    void*      user_data)
{
  const int nx = workspace->nvars;
  const int ng = nlp_ncons;

  // Reuse PSOPT's exact-derivative TNLP for gradients/Jacobian (ADOL-C inside).
  IPOPT_PSOPT nlpobj(workspace, user_data);
  Index n = nx, m = ng, nnz_jac = 0, nnz_h = 0;
  Ipopt::TNLP::IndexStyleEnum idx;
  nlpobj.get_nlp_info(n, m, nnz_jac, nnz_h, idx);

  // Constraint-Jacobian sparsity pattern (rows/cols) from eval_jac_g(values=NULL).
  std::vector<Index> iRow(nnz_jac), jCol(nnz_jac);
  if (nnz_jac > 0)
    nlpobj.eval_jac_g(n, nullptr, false, m, nnz_jac, iRow.data(), jCol.data(), nullptr);
  std::vector<casadi_int> rows(nnz_jac), cols(nnz_jac);
  for (int k = 0; k < (int)nnz_jac; ++k) { rows[k] = iRow[k]; cols[k] = jCol[k]; }

  // Build the NLP with exact-derivative CasADi functions.
  ObjFun objfun(nx, &nlpobj, f, workspace);
  ConFun confun(nx, ng, &nlpobj, g, workspace, rows, cols);

  MX x = MX::sym("x", nx);
  MX fx = objfun(std::vector<MX>{x}).at(0);
  MX gx = (ng > 0) ? confun(std::vector<MX>{x}).at(0) : MX::zeros(0, 1);
  MXDict nlp = {{"x", x}, {"f", fx}, {"g", gx}};

  // Choose the CasADi nlpsol plugin (net-new: solvers PSOPT has no native path
  // to, e.g. sqpmethod / fatrop). Default "ipopt".
  std::string plugin = algorithm.casadi_solver.empty() ? std::string("ipopt")
                                                        : algorithm.casadi_solver;
  Dict opts;
  if (plugin == "ipopt") {
    opts["ipopt.linear_solver"]  = algorithm.ipopt_linear_solver;   // mumps / ma97 / ...
    opts["ipopt.tol"]            = algorithm.nlp_tolerance;
    opts["ipopt.max_iter"]       = algorithm.nlp_iter_max;
    opts["ipopt.mu_strategy"]    = std::string("adaptive");
    // NOTE: exact Hessian is intentionally NOT requested here. The CasADi path
    // wraps PSOPT's numeric objective/constraint callbacks opaquely, so CasADi
    // cannot build the exact Hessian of the Lagrangian from them; that requires
    // a symbolic (MX/SX) transcription. Exact ADOL-C Hessian is available only
    // via the native nlp_method=="IPOPT" path (algorithm.hessian=="exact").
    opts["ipopt.hessian_approximation"] = std::string("limited-memory");
    if (algorithm.print_level == 0) opts["ipopt.print_level"] = 0;
  } else {
    // Generic options for non-IPOPT plugins (sqpmethod/fatrop/...); these carry
    // their own defaults and option namespaces.
    if (algorithm.print_level == 0) opts["print_time"] = false;
  }

  Function solver = nlpsol("psopt_casadi_solver", plugin, nlp, opts);

  std::vector<double> x0v(nx), lbx(nx), ubx(nx), lbg(ng), ubg(ng);
  for (int j = 0; j < nx; ++j) {
    x0v[j] = (*x0)(j, 0); lbx[j] = (*xlb)(j, 0); ubx[j] = (*xub)(j, 0);
  }
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
