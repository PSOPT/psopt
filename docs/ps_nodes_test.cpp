// Standalone validation of Legendre-Gauss (LG) and Legendre-Gauss-Radau (LGR)
// nodes/weights/differentiation-matrix + endpoint interpolation, before wiring
// into PSOPT. Checks against known values.
#include <Eigen/Dense>
#include <Eigen/Eigenvalues>
#include <cstdio>
#include <cmath>
using Eigen::MatrixXd; using Eigen::VectorXd;

static int len(const MatrixXd& x){ return (int)x.size(); }

static void bary_weights(const MatrixXd& x, MatrixXd& b){
  int n=len(x); b.resize(n,1);
  for(int j=0;j<n;j++){ double p=1.0; for(int k=0;k<n;k++) if(k!=j) p*=(x(j)-x(k)); b(j)=1.0/p; }
}
static void bary_diffmat(const MatrixXd& x, MatrixXd& D){
  int n=len(x); D.resize(n,n); MatrixXd b; bary_weights(x,b);
  for(int i=0;i<n;i++){ double diag=0.0;
    for(int j=0;j<n;j++) if(j!=i){ D(i,j)=(b(j)/b(i))/(x(i)-x(j)); diag-=D(i,j);} D(i,i)=diag; }
}
static void lagrange_at(const MatrixXd& x, double t, MatrixXd& L){
  int n=len(x); L.resize(1,n); MatrixXd b; bary_weights(x,b);
  double s=0.0; for(int k=0;k<n;k++){ L(0,k)=b(k)/(t-x(k)); s+=L(0,k);} L/=s;
}

// M = N+1 Legendre-Gauss nodes (Golub-Welsch), interior to (-1,1).
static void lgnodes(int N, MatrixXd& x, MatrixXd& w, MatrixXd& D){
  int M=N+1; MatrixXd J=MatrixXd::Zero(M,M);
  for(int k=1;k<M;k++){ double b=k/std::sqrt(4.0*k*k-1.0); J(k-1,k)=b; J(k,k-1)=b; }
  Eigen::SelfAdjointEigenSolver<MatrixXd> es(J);
  VectorXd ev=es.eigenvalues(); MatrixXd V=es.eigenvectors();
  x.resize(M,1); w.resize(M,1);
  for(int i=0;i<M;i++){ x(i)=ev(i); w(i)=2.0*V(0,i)*V(0,i); }
  bary_diffmat(x,D);
}
// M = N+1 Legendre-Gauss-Radau nodes (fixed at -1), Golub-Welsch-Radau.
static void lgrnodes(int N, MatrixXd& x, MatrixXd& w, MatrixXd& D){
  int M=N+1; double a=-1.0;
  VectorXd beta(M); beta.setZero();
  for(int k=1;k<M;k++) beta(k)=k/std::sqrt(4.0*k*k-1.0);
  // leading (M-1) block T; solve (T - a I) delta = beta_{M-1}^2 e_{M-1}
  int m1=M-1; MatrixXd T=MatrixXd::Zero(m1,m1);
  for(int k=1;k<m1;k++){ T(k-1,k)=beta(k); T(k,k-1)=beta(k); }
  VectorXd rhs=VectorXd::Zero(m1); rhs(m1-1)=beta(M-1)*beta(M-1);
  VectorXd delta=(T - a*MatrixXd::Identity(m1,m1)).colPivHouseholderQr().solve(rhs);
  double alphaM = a + delta(m1-1);
  MatrixXd J=MatrixXd::Zero(M,M);
  for(int k=1;k<M;k++){ J(k-1,k)=beta(k); J(k,k-1)=beta(k); }
  J(M-1,M-1)=alphaM;
  Eigen::SelfAdjointEigenSolver<MatrixXd> es(J);
  VectorXd ev=es.eigenvalues(); MatrixXd V=es.eigenvectors();
  x.resize(M,1); w.resize(M,1);
  for(int i=0;i<M;i++){ x(i)=ev(i); w(i)=2.0*V(0,i)*V(0,i); }
  bary_diffmat(x,D);
}

int main(){
  // 3-point Gauss: nodes +-sqrt(3/5),0 ; weights 5/9,8/9,5/9
  MatrixXd x,w,D; lgnodes(2,x,w,D);
  std::printf("LG(3) nodes:  %.6f %.6f %.6f  (exp -0.774597 0 0.774597)\n", x(0),x(1),x(2));
  std::printf("LG(3) wts:    %.6f %.6f %.6f  (exp 0.555556 0.888889 0.555556)\n", w(0),w(1),w(2));
  // differentiation exactness: d/dt of f=t^2 is 2t at nodes
  VectorXd f(3); for(int i=0;i<3;i++) f(i)=x(i)*x(i);
  VectorXd df=D*f; std::printf("LG(3) D on t^2: %.4f %.4f %.4f (exp 2x= %.4f %.4f %.4f)\n",
      df(0),df(1),df(2), 2*x(0),2*x(1),2*x(2));
  // 3-point Radau (fixed -1): nodes -1, 0.289898, 0.689898 ; weights 2/9, ...
  MatrixXd xr,wr,Dr; lgrnodes(2,xr,wr,Dr);
  std::printf("LGR(3) nodes: %.6f %.6f %.6f  (exp -1 0.289898 0.689898)\n", xr(0),xr(1),xr(2));
  std::printf("LGR(3) wts:   %.6f %.6f %.6f  (sum exp 2; 0.222222 1.024972 0.752806)\n", wr(0),wr(1),wr(2));
  std::printf("LGR wt sum=%.6f (exp 2)\n", wr.sum());
  // endpoint interpolation: interpolate f=t^2 to t=+1 from LG nodes -> should be 1
  MatrixXd Ltf; lagrange_at(x, 1.0, Ltf);
  double fend=0; for(int k=0;k<3;k++) fend+=Ltf(0,k)*f(k);
  std::printf("LG interp t^2 to t=+1: %.6f (exp 1.0)\n", fend);
  return 0;
}
