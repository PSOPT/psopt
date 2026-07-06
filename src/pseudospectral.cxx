/*********************************************************************************************

This file is part of the PSOPT library, a software tool for computational optimal control

Copyright (C) 2009-2020 Victor M. Becerra

This library is free software; you can redistribute it and/or
modify it under the terms of the GNU Lesser General Public
License as published by the Free Software Foundation; either
version 2.1 of the License, or (at your option) any later version.

This library is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public
License along with this library; if not, write to the Free Software
Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA,
or visit http://www.gnu.org/licenses/

Author:    Professor Victor M. Becerra
Address:   University of Portsmouth
           School of Energy and Electronic Engineering
           Portsmouth PO1 3DJ
           United Kingdom
e-mail:    v.m.becerra@ieee.org

**********************************************************************************************/


#include "psopt.h"

using namespace Eigen;

using namespace PSOPT;

double delta(long l, long N)
{
      double delta_l;
      if (l == 0 || l==N)
            delta_l = 0.5;
      else
            delta_l = 1.0;
      return delta_l;
}

void diffmat_central_differences(MatrixXd& D, MatrixXd & x )
{

int i;

double h1, h2, h;

int N1 = length(x);

int N = N1-1;

D.resize(N1,N1);
D.setZero();

h = x(1)-x(0);

D(0,0)  = -1/h;

D(0,1)  =  1/h;

for (i=1;i<N;i++) { // EIGEN_UPDATE
    h2 = x(i+1)-x(i);
    h1 = x(i)-x(i-1);
    D(i,i+1) = 1/(h1+h2);
    D(i,i-1) = -1/(h1+h2);
}

h = x(N-1)-x(N);

D(N,N-1)  =1/h;

D(N,N)=-1/h;

}

void diffmat_lagrange3pt(MatrixXd& D, MatrixXd & x )
{

int i;

double  h;

int N1 = length(x);

int N = N1-1;


D.resize(N1,N1);
D.setZero();


  h = x(1)-x(0);
//D(1,1)  = -1/h;
  D(0,0)  = -1/h;
//D(1,2)  =  1/h;
  D(0,1)  =  1/h;

for (i=1;i<N;i++) { // EIGEN_UPDATE
    // 3 Point differentiation based on lagrange polynomial interpolation
    D(i,i-1)   =  (x(i)-x(i+1))/((x(i-1)-x(i))*(x(i-1)-x(i+1)));
    D(i,i)     =  (2*x(i)-x(i-1)-x(i+1))/((x(i)-x(i-1))*(x(i)-x(i+1)));
    D(i,i+1)   =  (x(i)-x(i-1))/((x(i+1)-x(i-1))*(x(i+1)-x(i)));

}


h = x(N-1)-x(N);

D(N,N-1)  =1/h;

D(N,N)=-1/h;

}


void diffmat_lagrange5pt(MatrixXd& D, MatrixXd & x )
{

int i;

double h;

int N1 = length(x);

int N = N1-1;


D.resize(N1,N1);
D.setZero();


  h = x(1)-x(0);

  D(0,0)  = -1/h;

  D(0,1)  = 1/h;

    // 3 Point differentiation based on lagrange polynomial interpolation
i=1; // EIGEN_UPDATE
    D(i,i-1)   =  (x(i)-x(i+1))/((x(i-1)-x(i))*(x(i-1)-x(i+1)));
    D(i,i)     =  (2*x(i)-x(i-1)-x(i+1))/((x(i)-x(i-1))*(x(i)-x(i+1)));
    D(i,i+1)   =  (x(i)-x(i-1))/((x(i+1)-x(i-1))*(x(i+1)-x(i)));

i=N-1; // EIGEN_UPDATE
    D(i,i-1)   =  (x(i)-x(i+1))/((x(i-1)-x(i))*(x(i-1)-x(i+1)));
    D(i,i)     =  (2*x(i)-x(i-1)-x(i+1))/((x(i)-x(i-1))*(x(i)-x(i+1)));
    D(i,i+1)   =  (x(i)-x(i-1))/((x(i+1)-x(i-1))*(x(i+1)-x(i)));

for (i=2;i<(N-1);i++) {  // EIGEN_UPDATE
    // 5 Point central differentiation based on lagrange polynomial interpolation
    double x1 = x(i-2);
    double x2 = x(i-1);
    double x3 = x(i);
    double x4 = x(i+1);
    double x5 = x(i+2);
    D(i,i-2)   =  -((x2-x3)*(x3-x4)*(x3-x5))/((x1-x2)*(x1-x3)*(x1-x4)*(x1-x5));
    D(i,i-1)   =  ((x1-x3)*(x3-x4)*(x3-x5))/((x1-x2)*(x2-x3)*(x2-x4)*(x2-x5));
    D(i,i)     =  (x3*(4*x3*x3+2*x4*x5-3*x3*(x4+x5))+x2*(-3*x3*x3-x4*x5+2*x3*(x4+x5))+x1*(-3*x3*x3+2*x3*x4+2*x3*x5-x4*x5-x2*(-2*x3+x4+x5)))/((x1-x3)*(x2-x3)*(x3-x4)*(x3-x5));
    D(i,i+1)   =  ((x1-x3)*(-x2+x3)*(x3-x5))/((x1-x4)*(-x2+x4)*(-x3+x4)*(x4-x5));
    D(i,i+2)   =  ((-x1+x3)*(-x2+x3)*(x3-x4))/((-x1+x5)*(-x2+x5)*(-x3+x5)*(-x4+x5));

}

h = x(N-1)-x(N);

D(N,N-1)  =1/h;

D(N,N)=-1/h;

}

void legendre_points(int N, MatrixXd& x, MatrixXd& w)
{
// Finds the roots of the Legendre polynomials in (-1,1), also known as the Legendre points,
// together with the corresponding weights for quadrature.
//
      MatrixXd beta;

      int N1 = N+1;

      int i;

      MatrixXd T1(N1,N1), T2(N1,N1), T(N1,N1);

      MatrixXd V(N1,N1);

      RowVectorXi indx(N1,1);

      T1.setZero();
      T2.setZero();


      MatrixXd C(N1-1,1);
      
      for(i=0; i< N1-1; i++) C(i) = (double) i;
      
      C = 2.0*C;
      
      C = C.array().pow(-2.0);
      
      C = ones(N,1)-C;
      
      C = C.array().sqrt();

      beta = 0.5*elemDivision( ones(N,1) , C );

      for(i=0;i<N1-2;i++) { // EIGEN_UPDATE
  	    T1(i,i+1) = beta(i);
            T2(i+1,i) = beta(i);
      }
      T  = T1 + T2;

      EigenSolver<MatrixXd> es(T);
      
      x = es.eigenvalues().real();
      V = es.eigenvectors().real();

      sort_vector(x,indx);

     MatrixXd V1 = V.row(0);
     for (i=0; i<length(x); i++) {
        w(i) = 2.0*pow(V1(indx(i)), 2.0);     
     }

}

void lglnodes(int N, MatrixXd& x, MatrixXd& w, MatrixXd& P, MatrixXd& D, Workspace* workspace)
{
// Computes the Legendre-Gauss-Lobatto nodes, weights and the LGL Vandermonde
// matrix. The LGL nodes are the zeros of (1-x^2)*P'_N(x).
//
// Reference on LGL nodes and weights:
//   C. Canuto, M. Y. Hussaini, A. Quarteroni, T. A. Tang, "Spectral Methods
//   in Fluid Dynamics," Section 2.3. Springer-Verlag 1987
//


  long N1 = N+1;
  long k;
  long l,i,j;

  MatrixXd xold, X, Xdiff, L;
  // x.resize(1,N1);
  x.resize(N1,1);
  for (i=0; i<N+1; i++) {
    x(i) = cos((pi*i)/N);   
  }
//  Use the Chebyshev-Gauss-Lobatto nodes as the first guess

  
//  The Legendre Vandermonde Matrix

  P.resize(N1,N1);

  w.resize(N1,1);
  P.setZero();
  w.setZero();
// Compute P_(N) using the recursion relation
// Compute its first and second derivatives and
// update x using the Newton-Raphson method.
  xold = 2*ones(N1,1);
  while ( Max( Abs( x-xold ) ) > PSOPT_extras::GetEPS() )
  {

       xold = x;

       P.col(0)=ones(N1,1); P.col(1)=x;

       for ( k = 1; k<N; k++ ) // EIGEN_UPDATE
       {
	          double kd = (double) k;

             P.col(k+1) = ( elemProduct((2*kd+1)*x,P.col(k)) - (kd)*P.col(k-1) )/(kd+1);

       }

       double N1d = (double) N1;

       x = xold-elemDivision(  elemProduct(x,P.col(N1-1))-P.col(N-1) , N1d*P.col(N1-1) );


  }

  for( k=0; k<N1; k++) // EIGEN_UPDATE
  {
     w(k) = 2/((N*N1)*P(k,N1-1)*P(k,N1-1) );
  }

  if ( workspace->differential_defects == "standard") {

	// Compute now the differentiation matrix D using the standard formula

        X.resize(N1,N1);
        X.setZero();
	for(k=0;k<N1;k++) // EIGEN_UPDATE
	{

      X.col(k) = x;
	}

        Xdiff.resize(N1,N1);

        for(i=0;i<N1;i++) {  // EIGEN_UPDATE
                for(j=0;j<N1;j++) {
                         double identity;
                         if (i==j) identity = 1.0; else identity=0.0;
                         Xdiff(i,j) = X(i,j)-X(j,i) + identity;
                }
        }


        L.resize(N1,N1);
        L.setZero();
	for(k=0;k<N1;k++) // EIGEN_UPDATE
	{

      L.col(k) =  P.col(N1-1);
	}

	for (k=0;k<N1*N1;k=k+N1+1) // EIGEN_UPDATE
	{
	     L(k)=1;
	}


   D = elemDivision(L, elemProduct(Xdiff, (L).transpose() ) );


	for (k=0; k<N1*N1; k=k+N1+1)  // EIGEN_UPDATE
	{
	    D(k) = 0.0;
	}


   D(0) = ((double) N1*N)/4.0;


   D(N1*N1-1)=-((double) N1*N)/4.0;

        for(i=0;i<N1;i++) { // EIGEN_UPDATE
                for(j=0;j<N1;j++) {
                     D(i,j) = -D(i,j);
                }
        }

  }

  else if (workspace->differential_defects == "reduced-roundoff") {

  // Use equation 2_4_35 from the book Canuto et al (2006)
	for(j=0;j<=N;j++) {
		for(l=0;l<=N;l++)  {
		if(j!=l) {

         D(j,l) = delta(l,N)*pow(-1.0,(double) (j+l))/(delta(j,N)*(x(j)-x(l)));
		}
		else {
			double sum=0.0;
			long i;
			for(i=0;i<=N;i++) {
				if (i != j) {

               sum += delta(i,N)*pow(-1.0,(double) (i+j))/(delta(j,N)*(x(j)-x(i)));
				}
			}

         D(j,l) = -sum;
		}
		}
	}

        for(i=0;i<N1;i++) { // EIGEN_UNDATE
                for(j=0;j<N1;j++) {
                     D(i,j) = -D(i,j);
                }
        }

  }

  else if (workspace->differential_defects == "central-differences") {
      diffmat_central_differences( D, x );
        for(i=0;i<N1;i++) { // EIGEN_UPDATE
                for(j=0;j<N1;j++) {
                     D(i,j) = -D(i,j);
                }
        }
  }

  else if (workspace->differential_defects == "Lagrange-3pt") {
      diffmat_lagrange3pt( D, x );
        for(i=0;i<N1;i++) {  // EIGEN_UPDATE
                for(j=0;j<N1;j++) {
                     D(i,j) = -D(i,j);
                }
        }
  }

  return;



}

double cbar(long j,long N)
{
    double retval;

    if (j==0 || j==N)
        retval = 2.0;
    else
        retval = 1.0;

    return retval;
}


void cglnodes(int N, MatrixXd& x, MatrixXd& w,  MatrixXd& D, Workspace* workspace)
{
// Computes the Chebyshev-Gauss-Lobatto nodes, weights and the first order
// differentiation matrix
//
// Reference
//   C. Canuto, M. Y. Hussaini, A. Quarteroni, T. A. Tang, "Spectral Methods
//   in Fluid Dynamics," Springer-Verlag 1987
//


  long N1 = N+1;
  long k;
  long i,j,l;

//  compute the Chebyshev-Gauss-Lobatto nodes

  x.resize(N1,1);
  for (i=0; i<N+1; i++) {
    x(i) = cos((pi*i)/N); 
  }

  w = zeros(N1,1);



  w(0)=pi/(2*N);

  for( k=1; k<N1-1; k++)  // EIGEN_UPDATE
  {
     w(k) = pi/N;
  }



 w(N1-1)=pi/(2*N);

  // Compute now the differentiation matrix D

  if ( workspace->differential_defects == "standard" ) {

	for(j=0;j<=N;j++) {
		for(l=0;l<=N;l++)  {
		if(j!=l)

			D(j,l)= -cbar(j,N)/(2*cbar(l,N))*pow(-1.0,(int) (j+l))/( sin((j+l)*pi/(2*N))*sin((j-l)*pi/(2*N)));
			else if ( j==l && 1<=j && l <= N-1)

         D(j,l) =  -x(j)/(2*pow( sin(j*pi/N), 2));
			else if ( j==0 && l==0 )

         D(j,l)=(2*pow(N,2.0)+1.0)/6.0;
			else if ( j==N && l==N)

         D(j,l)=-(2*pow(N,2.0)+1.0)/6.0;
		}
	}


        for(i=0;i<N1;i++) {  // EIGEN_UPDATE
                for(j=0;j<N1;j++) {
                     D(i,j) = -D(i,j);
                }
        }

  }

  else if ( workspace->differential_defects == "reduced-roundoff") {

  // Use equation No. 2_4_35 from Canuto et al (2006).
	for(j=0;j<=N;j++) {
		for(l=0;l<=N;l++)  {
		if(j!=l) {

         D(j,l) = delta(l,N)*pow(-1.0,(int) (j+l))/(delta(j,N)*(x(j)-x(l)));
		}
		else {
			double sum=0.0;
			long i;
			for(i=0;i<=N;i++) {
				if (i != j) {
               sum += delta(i,N)*pow(-1.0,(int) (i+j))/(delta(j,N)*(x(j)-x(i)));
				}
			}

         D(j,l) = -sum;
		}
		}
	}

        for(i=0;i<N1;i++) {  // EIGEN_UPDATE
                for(j=0;j<N1;j++) {
                     D(i,j) = -D(i,j);
                }
        }

  }

  else if (workspace->differential_defects == "central-differences") {
      diffmat_central_differences( D, x );

        for(i=0;i<N1;i++) {   // EIGEN_UPDATE
                for(j=0;j<N1;j++) {
                     D(i,j) = -D(i,j);
                }
        }
  }


  return;

}


// ---------------------------------------------------------------------------
// Legendre-Gauss (LG) and Legendre-Gauss-Radau (LGR) pseudospectral node sets.
// Both build a square barycentric differentiation matrix over M=N+1 nodes and
// return Gauss/Radau quadrature weights. Unlike LGL these do not include both
// (LG: neither; LGR: only t=-1) endpoints, so boundary conditions must be
// interpolated to the endpoints (see lagrange_endpoint_weights + NLP_constraints).
// The node/weight math is validated against exact values in ps_nodes_test.cpp.
// ---------------------------------------------------------------------------

// Barycentric differentiation matrix D (M x M) for a distinct node set x.
static void barycentric_diffmat(const MatrixXd& x, MatrixXd& D)
{
    int n = (int) x.size();
    MatrixXd b(n,1);
    for (int j=0;j<n;j++){ double p=1.0; for(int k=0;k<n;k++) if(k!=j) p*=(x(j)-x(k)); b(j)=1.0/p; }
    D.resize(n,n);
    for (int i=0;i<n;i++){
        double diag=0.0;
        for (int j=0;j<n;j++) if(j!=i){ D(i,j)=(b(j)/b(i))/(x(i)-x(j)); diag-=D(i,j); }
        D(i,i)=diag;
    }
}

// Lagrange interpolation weights L (1 x M) so that f(t) ~= sum_k L(k) f(x_k).
// If t coincides with a node, returns the corresponding unit vector (so LGL/CGL,
// whose endpoints ARE nodes, reduce to selecting that node's value).
void lagrange_endpoint_weights(const MatrixXd& x, double t, MatrixXd& L)
{
    int n = (int) x.size();
    L.resize(1,n); L.setZero();
    for (int k=0;k<n;k++) if (fabs(t-x(k)) < 1e-12) { L(0,k)=1.0; return; }
    MatrixXd b(n,1);
    for (int j=0;j<n;j++){ double p=1.0; for(int kk=0;kk<n;kk++) if(kk!=j) p*=(x(j)-x(kk)); b(j)=1.0/p; }
    double s=0.0;
    for (int k=0;k<n;k++){ L(0,k)=b(k)/(t-x(k)); s+=L(0,k); }
    for (int k=0;k<n;k++) L(0,k)/=s;
}

// M = N+1 Legendre-Gauss nodes/weights (Golub-Welsch) + barycentric diff matrix.
void lgnodes(int N, MatrixXd& x, MatrixXd& w, MatrixXd& D, Workspace* /*workspace*/)
{
    int M = N+1;
    MatrixXd J = MatrixXd::Zero(M,M);
    for (int k=1;k<M;k++){ double b=k/sqrt(4.0*k*k-1.0); J(k-1,k)=b; J(k,k-1)=b; }
    Eigen::SelfAdjointEigenSolver<MatrixXd> es(J);
    Eigen::VectorXd ev=es.eigenvalues(); MatrixXd V=es.eigenvectors();
    x.resize(M,1); w.resize(M,1);
    for (int i=0;i<M;i++){ x(i)=ev(i); w(i)=2.0*V(0,i)*V(0,i); }   // ascending
    barycentric_diffmat(x,D);
}

// M = N+1 Legendre-Gauss-Radau nodes/weights (fixed at t=-1) + diff matrix,
// via the Golub-Welsch-Radau modification of the last Jacobi diagonal.
void lgrnodes(int N, MatrixXd& x, MatrixXd& w, MatrixXd& D, Workspace* /*workspace*/)
{
    int M = N+1; double a=-1.0;
    Eigen::VectorXd beta(M); beta.setZero();
    for (int k=1;k<M;k++) beta(k)=k/sqrt(4.0*k*k-1.0);
    int m1=M-1;
    MatrixXd T = MatrixXd::Zero(m1,m1);
    for (int k=1;k<m1;k++){ T(k-1,k)=beta(k); T(k,k-1)=beta(k); }
    Eigen::VectorXd rhs=Eigen::VectorXd::Zero(m1); rhs(m1-1)=beta(M-1)*beta(M-1);
    Eigen::VectorXd delta=(T - a*MatrixXd::Identity(m1,m1)).colPivHouseholderQr().solve(rhs);
    double alphaM = a + delta(m1-1);
    MatrixXd J = MatrixXd::Zero(M,M);
    for (int k=1;k<M;k++){ J(k-1,k)=beta(k); J(k,k-1)=beta(k); }
    J(M-1,M-1)=alphaM;
    Eigen::SelfAdjointEigenSolver<MatrixXd> es(J);
    Eigen::VectorXd ev=es.eigenvalues(); MatrixXd V=es.eigenvectors();
    x.resize(M,1); w.resize(M,1);
    for (int i=0;i<M;i++){ x(i)=ev(i); w(i)=2.0*V(0,i)*V(0,i); }
    barycentric_diffmat(x,D);
}
