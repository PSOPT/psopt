// ============================================================================
//  psopt_ad_driver.hpp  —  Level-A AD backend interface (record + sparse Jac/Hess)
//  Results are returned as RAII SparseTriplet (owned vectors); no manual
//  allocation survives at the call sites.
//
//  CppAD (taped) is the production path. AUTODIFF / XAD are forward-mode and
//  build the Jacobian densely via directional sweeps (diagnostic only; no exact
//  Hessian -- the caller falls back to limited-memory BFGS).
// ============================================================================
#ifndef PSOPT_AD_DRIVER_HPP
#define PSOPT_AD_DRIVER_HPP
#include "psopt_ad_backend.hpp"
#include <functional>
#include <vector>
#include <cstdlib>

namespace psopt_ad {

using ADVecFunc = std::function<void(const adouble*, adouble*)>;

struct SparseTriplet {                 // row-compressed entries, RAII
    std::vector<unsigned int> row, col;
    std::vector<double>       val;
    int nnz() const { return (int)val.size(); }
};

struct ADHandle {
    int tag = -1, n = 0, m = 0;
    ADVecFunc f;
#if PSOPT_AD_BACKEND == PSOPT_AD_CPPAD
    CppAD::ADFun<double> fun;
    bool jac_struct=false, hess_struct=false;
    CppAD::sparse_rc<std::vector<size_t> >                       jac_pat, hes_pat;
    CppAD::sparse_rcv<std::vector<size_t>, std::vector<double> > jac_sub, hes_sub;
    CppAD::sparse_jac_work jac_work;
    CppAD::sparse_hes_work hes_work;
#endif
    ADHandle() = default;
    ADHandle(const ADHandle&) = delete;
    ADHandle& operator=(const ADHandle&) = delete;
};

inline void ad_record(ADHandle& h, int n, int m, const double* x0, ADVecFunc f){
    h.n=n; h.m=m; h.f=f;
#if PSOPT_AD_BACKEND == PSOPT_AD_CPPAD
    // Tape with native AD<double> vectors (CppAD requires the exact element type), then
    // hand the user functor an adouble* view of the same storage (layout-identical subclass).
    std::vector<CppAD::AD<double> > x(n), y(m);
    for(int i=0;i<n;i++) x[i]=x0[i];
    CppAD::Independent(x);
    f(reinterpret_cast<const adouble*>(x.data()), reinterpret_cast<adouble*>(y.data()));
    h.fun.Dependent(x,y);
    // Keep cached sparsity patterns across re-records: IPOPT fixes the Hessian/Jacobian
    // structure (from get_nlp_info) and eval_h re-records the lambda-weighted tape every
    // call. Only the work objects are tied to the specific tape and must be rebuilt.
    h.jac_work.clear(); h.hes_work.clear();
#else
    // Forward-mode backends: nothing to record now -- the functor is stored in h.f and
    // re-invoked per direction in ad_sparse_jacobian / ad_gradient.
    (void)x0;
#endif
}

inline SparseTriplet ad_sparse_jacobian(ADHandle& h, const double* x, bool reuse){
    SparseTriplet T;
#if PSOPT_AD_BACKEND == PSOPT_AD_CPPAD
    std::vector<double> xv(x, x+h.n);
    if(!(reuse && h.jac_struct)){
        size_t n=(size_t)h.n;
        CppAD::sparse_rc<std::vector<size_t> > eye(n,n,n);
        for(size_t k=0;k<n;k++) eye.set(k,k,k);
        h.fun.for_jac_sparsity(eye,false,false,false,h.jac_pat);
        h.jac_sub = CppAD::sparse_rcv<std::vector<size_t>,std::vector<double> >(h.jac_pat);
        h.jac_work.clear();
        h.jac_struct=true;
    }
    h.fun.sparse_jac_for((size_t)1, xv, h.jac_sub, h.jac_pat, "cppad", h.jac_work);
    size_t nz=h.jac_sub.nnz();
    T.row.resize(nz); T.col.resize(nz); T.val.resize(nz);
    for(size_t k=0;k<nz;k++){ T.row[k]=(unsigned int)h.jac_sub.row()[k]; T.col[k]=(unsigned int)h.jac_sub.col()[k]; T.val[k]=h.jac_sub.val()[k]; }
#else
    // ---- Forward-mode (autodiff / XAD): dense Jacobian via n directional sweeps ----
    // No tape: re-evaluate the stored functor once per input direction. The (i,j) emission
    // order is deterministic, so the structure pass (reuse=false) and value passes agree.
    (void)reuse;
    int n=h.n, m=h.m;
    std::vector<adouble> X(n), Y(m);
    for(int j=0;j<n;j++) X[j]=adouble(x[j]);
    T.row.reserve((size_t)m*n); T.col.reserve((size_t)m*n); T.val.reserve((size_t)m*n);
    for(int j=0;j<n;j++){
        seed_deriv(X[j], 1.0);
        h.f(X.data(), Y.data());
        for(int i=0;i<m;i++){
            T.row.push_back((unsigned int)i);
            T.col.push_back((unsigned int)j);
            T.val.push_back(get_deriv(Y[i]));
        }
        seed_deriv(X[j], 0.0);
    }
#endif
    return T;
}

inline SparseTriplet ad_sparse_hessian(ADHandle& h, const double* x, bool reuse){
    SparseTriplet T;
#if PSOPT_AD_BACKEND == PSOPT_AD_CPPAD
    std::vector<double> xv(x, x+h.n);
    // IPOPT fixes the Hessian structure from get_nlp_info; detect once per problem size
    // (lambda-independent superset pattern) and reuse it for every eval_h, regardless of the
    // caller's reuse hint. sparse_hes returns exact zeros for entries inactive at this point.
    (void)reuse;
    if(!h.hess_struct || (int)h.hes_pat.nr()!=h.n){
        std::vector<bool> sd((size_t)h.n,true), sr((size_t)h.m,true);
        h.fun.for_hes_sparsity(sd, sr, false, h.hes_pat);
        h.hes_sub = CppAD::sparse_rcv<std::vector<size_t>,std::vector<double> >(h.hes_pat);
        h.hes_work.clear();
        h.hess_struct=true;
    }
    std::vector<double> w((size_t)h.m, 1.0);
    h.fun.sparse_hes(xv, w, h.hes_sub, h.hes_pat, "cppad.symmetric", h.hes_work);
    for(size_t k=0;k<h.hes_sub.nnz();k++){           // keep lower triangle (IPOPT convention)
        size_t r=h.hes_sub.row()[k], c=h.hes_sub.col()[k];
        if(r>=c){ T.row.push_back((unsigned int)r); T.col.push_back((unsigned int)c); T.val.push_back(h.hes_sub.val()[k]); }
    }
#else
    // Forward 1st-order modes provide no exact Hessian; return empty so the caller falls
    // back to limited-memory BFGS (enforced in NLP_interface for non-taped backends).
    (void)x;(void)reuse;
#endif
    return T;
}

inline std::vector<double> ad_gradient(ADHandle& h, const double* x){
    std::vector<double> g(h.n, 0.0);
#if PSOPT_AD_BACKEND == PSOPT_AD_CPPAD
    std::vector<double> xv(x, x+h.n), J = h.fun.Jacobian(xv);
    for(int i=0;i<h.n;i++) g[i]=J[i];
#else
    // Forward-mode gradient of the scalar objective (output 0): n directional sweeps.
    int n=h.n, m=h.m;
    std::vector<adouble> X(n), Y(m);
    for(int j=0;j<n;j++) X[j]=adouble(x[j]);
    for(int j=0;j<n;j++){
        seed_deriv(X[j], 1.0);
        h.f(X.data(), Y.data());
        g[j]=get_deriv(Y[0]);
        seed_deriv(X[j], 0.0);
    }
#endif
    return g;
}

} // namespace psopt_ad
#endif
