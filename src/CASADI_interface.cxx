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
// The Hessian is limited-memory BFGS by default, but can be switched to an
// exact Hessian when algorithm.hessian=="exact" and the selected CasADi
// plugin is IPOPT.
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
// TAG_HC: weighted-constraint Hessian tape (re-taped per Hessian eval).
static const int TAG_F = 21, TAG_G = 22, TAG_HC = 24;

// Build a symmetric CasADi DM from an ADOL-C lower-triangular sparse_hess result.
static DM symm_from_lower(int n, int nnz, unsigned int* ir, unsigned int* jc, double* v) {
  std::vector<casadi_int> R, C; std::vector<double> V;
  for (int k = 0; k < nnz; ++k) {
    R.push_back(ir[k]); C.push_back(jc[k]); V.push_back(v[k]);
    if (ir[k] != jc[k]) { R.push_back(jc[k]); C.push_back(ir[k]); V.push_back(v[k]); }
  }
  return DM::triplet(R, C, V, n, n);
}

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

// ---------------------------------------------------------------------------
// Reverse-mode (adjoint) overrides so CasADi can assemble the EXACT Hessian of
// the Lagrangian. The adjoint of the objective is f_bar*grad_f, of the
// constraints is jac_g'*g_bar; their Jacobians w.r.t. x are f_bar*Hess_f and
// Hess(g_bar'*g). With f_bar=obj_factor, g_bar=lambda, CasADi's sum of these is
// exactly PSOPT's Lagrangian Hessian (Lagrangian_ad -> sparse_hess).
// ---------------------------------------------------------------------------

// Jacobian of the objective adjoint x_bar = f_bar*grad_f w.r.t. (x, f_nom, f_bar).
// Output (n x (n+2)) = [ f_bar*Hess_f | 0 | grad_f ].
class ObjRevJac : public Callback {
public:
  ObjRevJac(int n) : n_(n) { construct("psopt_obj_revjac"); }
  casadi_int get_n_in()  override { return 4; }  // x, f_nom, f_bar, xbar_nom
  casadi_int get_n_out() override { return 3; }  // d xbar/dx, /d f_nom, /d f_bar
  Sparsity get_sparsity_in(casadi_int i) override {
    if (i == 0 || i == 3) return Sparsity::dense(n_, 1);
    return Sparsity::dense(1, 1);              // f_nom, f_bar are scalars
  }
  Sparsity get_sparsity_out(casadi_int i) override {
    return (i == 0) ? Sparsity::dense(n_, n_) : Sparsity::dense(n_, 1); }
  std::vector<DM> eval(const std::vector<DM>& arg) const override {
    std::vector<double> x = arg[0].nonzeros();
    double f_bar = arg[2].nonzeros()[0];
    std::vector<double> gf(n_); gradient(TAG_F, n_, x.data(), gf.data());
    int opt[2] = {0, 0};
    unsigned int* ir = nullptr; unsigned int* jc = nullptr; double* v = nullptr; int nnz = 0;
    sparse_hess(TAG_F, n_, 0, x.data(), &nnz, &ir, &jc, &v, opt);
    DM Hx = DM::zeros(n_, n_);
    for (int k = 0; k < nnz; ++k) {            // f_bar * Hess_f (symmetric)
      Hx(ir[k], jc[k]) = f_bar * v[k];
      if (ir[k] != jc[k]) Hx(jc[k], ir[k]) = f_bar * v[k];
    }
    free(ir); free(jc); free(v);
    return { Hx, DM::zeros(n_, 1), DM(gf) };   // d/dx, d/d f_nom, d/d f_bar
  }
private:
  int n_;
};

// Objective adjoint: (x, f_nom, f_bar) -> x_bar = f_bar*grad_f.
class ObjRev : public Callback {
public:
  ObjRev(int n) : n_(n), jac_(new ObjRevJac(n)) { construct("psopt_obj_rev"); }
  casadi_int get_n_in()  override { return 3; }  // x, f_nom, f_bar
  casadi_int get_n_out() override { return 1; }  // x_bar
  Sparsity get_sparsity_in(casadi_int i) override {
    return (i == 0) ? Sparsity::dense(n_, 1) : Sparsity::dense(1, 1); }
  Sparsity get_sparsity_out(casadi_int) override { return Sparsity::dense(n_, 1); }
  std::vector<DM> eval(const std::vector<DM>& arg) const override {
    std::vector<double> x = arg[0].nonzeros();
    double f_bar = arg[2].nonzeros()[0];
    std::vector<double> gf(n_); gradient(TAG_F, n_, x.data(), gf.data());
    std::vector<double> xb(n_); for (int i = 0; i < n_; ++i) xb[i] = f_bar * gf[i];
    return { DM(xb) };
  }
  bool has_jacobian() const override { return true; }
  Function get_jacobian(const std::string&, const std::vector<std::string>&,
                        const std::vector<std::string>&, const Dict&) const override { return *jac_; }
private:
  int n_; std::shared_ptr<ObjRevJac> jac_;
};

// Jacobian of the constraint adjoint x_bar = jac_g'*g_bar w.r.t. (x, g_nom, g_bar).
// Output (n x (n+2m)) = [ Hess(g_bar'*g) | 0 | jac_g' ].
class ConRevJac : public Callback {
public:
  ConRevJac(int n, int m, Workspace* ws) : n_(n), m_(m), ws_(ws) { construct("psopt_con_revjac"); }
  casadi_int get_n_in()  override { return 4; }  // x, g_nom, g_bar, xbar_nom
  casadi_int get_n_out() override { return 3; }  // d xbar/dx, /d g_nom, /d g_bar
  Sparsity get_sparsity_in(casadi_int i) override {
    if (i == 0 || i == 3) return Sparsity::dense(n_, 1);
    return Sparsity::dense(m_, 1);             // g_nom, g_bar
  }
  Sparsity get_sparsity_out(casadi_int i) override {
    if (i == 0) return Sparsity::dense(n_, n_);      // d/dx = Hess(g_bar'g)
    if (i == 1) return Sparsity::dense(n_, m_);      // d/d g_nom = 0
    return Sparsity::dense(n_, m_);                  // d/d g_bar = jac_g'
  }
  std::vector<DM> eval(const std::vector<DM>& arg) const override {
    std::vector<double> x = arg[0].nonzeros();
    std::vector<double> gbar = arg[2].nonzeros();
    DM Hx = DM::zeros(n_, n_);
    DM JT = DM::zeros(n_, m_);
    // constraint Jacobian transpose block via TAG_G
    { int opt[4] = {0,0,0,0}; unsigned int* ri=nullptr; unsigned int* ci=nullptr; double* vv=nullptr; int nnz=0;
      sparse_jac(TAG_G, m_, n_, 1, x.data(), &nnz, &ri, &ci, &vv, opt);
      for (int k = 0; k < nnz; ++k) JT(ci[k], ri[k]) = vv[k];   // jac_g'[j,i]
      free(ri); free(ci); free(vv); }
    // Hess(g_bar'*g): re-tape L_c = sum gbar_i g_i, then sparse_hess (TAG_HC)
    { trace_on(TAG_HC);
      for (int i = 0; i < n_; ++i) ws_->xad[i] <<= x[i];
      gg_ad(ws_->xad, ws_->gad, ws_);
      adouble Lc = 0.0; for (int i = 0; i < m_; ++i) Lc += gbar[i] * ws_->gad[i];
      double Lcv; Lc >>= Lcv; trace_off();
      int opt[2] = {0,0}; unsigned int* ir=nullptr; unsigned int* jc=nullptr; double* v=nullptr; int nnz=0;
      sparse_hess(TAG_HC, n_, 0, x.data(), &nnz, &ir, &jc, &v, opt);
      for (int k = 0; k < nnz; ++k) { Hx(ir[k], jc[k]) = v[k]; if (ir[k]!=jc[k]) Hx(jc[k], ir[k]) = v[k]; }
      free(ir); free(jc); free(v); }
    return { Hx, DM::zeros(n_, m_), JT };
  }
private:
  int n_, m_; Workspace* ws_;
};

// Constraint adjoint: (x, g_nom, g_bar) -> x_bar = jac_g'*g_bar.
class ConRev : public Callback {
public:
  ConRev(int n, int m, Workspace* ws) : n_(n), m_(m), ws_(ws), jac_(new ConRevJac(n, m, ws)) {
    construct("psopt_con_rev"); }
  casadi_int get_n_in()  override { return 3; }  // x, g_nom, g_bar
  casadi_int get_n_out() override { return 1; }
  Sparsity get_sparsity_in(casadi_int i) override {
    return (i == 0) ? Sparsity::dense(n_, 1) : Sparsity::dense(m_, 1); }
  Sparsity get_sparsity_out(casadi_int) override { return Sparsity::dense(n_, 1); }
  std::vector<DM> eval(const std::vector<DM>& arg) const override {
    std::vector<double> x = arg[0].nonzeros();
    std::vector<double> gbar = arg[2].nonzeros();
    std::vector<double> xb(n_, 0.0);
    int opt[4] = {0,0,0,0}; unsigned int* ri=nullptr; unsigned int* ci=nullptr; double* vv=nullptr; int nnz=0;
    sparse_jac(TAG_G, m_, n_, 1, x.data(), &nnz, &ri, &ci, &vv, opt);
    for (int k = 0; k < nnz; ++k) xb[ci[k]] += vv[k] * gbar[ri[k]];   // (jac_g' g_bar)
    free(ri); free(ci); free(vv);
    return { DM(xb) };
  }
  bool has_jacobian() const override { return true; }
  Function get_jacobian(const std::string&, const std::vector<std::string>&,
                        const std::vector<std::string>&, const Dict&) const override { return *jac_; }
private:
  int n_, m_; Workspace* ws_; std::shared_ptr<ConRevJac> jac_;
};

// Objective f(x) with an exact Jacobian (GradF). Held as a member so it outlives
// the solve; get_jacobian hands it to CasADi.
class ObjFun : public Callback {
public:
  ObjFun(int n, Workspace* ws,
         double (*f)(MatrixXd&, Workspace*)) : n_(n), ws_(ws), f_(f), gradf_(n, ws), rev_(new ObjRev(n)) {
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
  bool has_reverse(casadi_int nadj) const override { return nadj == 1; }
  Function get_reverse(casadi_int, const std::string&, const std::vector<std::string>&,
                       const std::vector<std::string>&, const Dict&) const override { return *rev_; }
private:
  int n_; Workspace* ws_;
  double (*f_)(MatrixXd&, Workspace*);
  GradF gradf_;
  std::shared_ptr<ObjRev> rev_;
};

// Constraints g(x) with an exact Jacobian (JacG).
class ConFun : public Callback {
public:
  ConFun(int n, int m, Workspace* ws,
         void (*g)(MatrixXd&, MatrixXd*, Workspace*)) : n_(n), m_(m), ws_(ws), g_(g), jacg_(n, m, ws), rev_(new ConRev(n, m, ws)) {
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
  bool has_reverse(casadi_int nadj) const override { return nadj == 1; }
  Function get_reverse(casadi_int, const std::string&, const std::vector<std::string>&,
                       const std::vector<std::string>&, const Dict&) const override { return *rev_; }
private:
  int n_, m_; Workspace* ws_;
  void (*g_)(MatrixXd&, MatrixXd*, Workspace*);
  JacG jacg_;
  std::shared_ptr<ConRev> rev_;
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
  const bool exact_hessian = exact && algorithm.hessian == "exact";

  if (plugin == "ipopt") {
    opts["ipopt.linear_solver"]  = algorithm.ipopt_linear_solver;   // mumps / ma97 / ...
    opts["ipopt.tol"]            = algorithm.nlp_tolerance;
    opts["ipopt.max_iter"]       = algorithm.nlp_iter_max;
    opts["ipopt.mu_strategy"]    = std::string("adaptive");
    opts["ipopt.hessian_approximation"] = exact_hessian ? std::string("exact")
                                                        : std::string("limited-memory");
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
