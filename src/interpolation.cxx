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


void lagrange_interpolation(MatrixXd& y, MatrixXd& x, MatrixXd& pointx, MatrixXd& pointy)
{
//
//       This function approximates a point-defined function using Lagrange polynomial interpolation
//
//       lagrange_interpolation(Y, X,POINTX,POINTY) approximates the function defined by the points:
//       P1=(POINTX(1),POINTY(1)), P2=(POINTX(2),POINTY(2)), ..., PN(POINTX(N),POINTY(N))
//       and calculate it in each element of X


   int i,j;

   long n = length(pointx);

   MatrixXd L(n, x.cols() );

   for(i=0;i<n;i++) {  // EIGEN_UPDATE: indices i and j shifted by -1
	  for(j=0;j<x.cols();j++) {
            L(i,j) = 1.0;
	  }
   }

   for (i=0;i<n;i++) {
      for (j=0;j<n;j++) { // EIGEN_UPDATE: indices i and j shifted by -1
           if (i != j) {
                  L.row(i) = elemProduct( L.row(i), (  x-pointx(j)*ones(1,length(x)) )  )/(pointx(i)-pointx(j));
           }
      }
   }

   y.resize(1,length(x));
   y = zeros(1,length(x));
   for (i=0;i<n;i++) {  // EIGEN_UPDATE index i shifted by -1
   	y = y+pointy(i)*L.row(i);
   }
}


void lagrange_interpolation_ad(adouble* y, adouble& x, adouble* pointx, adouble* pointy, int npoints, Workspace* workspace)
{
//
//       This function approximates a point-defined function using Lagrange polynomial interpolation
//       (version for automatic differentiation)
//       lagrange_interpolation(Y, X,POINTX,POINTY) approximates the function defined by the points:
//       P1=(POINTX(1),POINTY(1)), P2=(POINTX(2),POINTY(2)), ..., PN(POINTX(N),POINTY(N))
//       and calculate it at point x


   int i,j;

   int n = npoints;

   adouble* L = workspace->L_ad_tmp;

   for(i=0;i<npoints;i++) L[i] = 1.0;

   for (i=0;i<n;i++) {
      for (j=0;j<n;j++) {
           if (i != j) {
                L[i] = ( L[i]*(  x-pointx[j] )  )/(pointx[i]-pointx[j]);
           }
      }
   }

   *y = 0;

   for (i=0;i<n;i++) {
   	*y = (*y)+pointy[i]*L[i];
   }

}


void linear_interpolation(adouble* y, adouble& x, adouble* pointx, adouble* pointy, int npoints)
{
//    Linear interpolation from point values (version for automatic differentiation)

   int i,j;
   bool jdone = false;

   if( x < pointx[0] )
   {
      j = 0;
      jdone = true;
   }



  if ( x > pointx[npoints-1] )
   {
      j=npoints-2;
      jdone = true;
   }

  if (!jdone) {
     for(i=0;i<npoints-1;i++) {
       if( x >= pointx[i] && x <=pointx[i+1] )
       {
         j=i;
         break;
       }
     }
  }

   *y = pointy[j] + (x-pointx[j])*(pointy[j+1]-pointy[j])/(pointx[j+1]-pointx[j]);

}

void linear_interpolation(MatrixXd& y, double x, MatrixXd& pointx, MatrixXd& pointy, int npoints)
{
//    Linear interpolation from point values

   int i,j;
   bool jdone = false;

   if( x < pointx(0) )  // EIGEN_UPDATE
   {
      j = 0;
      jdone = true;
   }



  if ( x > pointx(npoints-1) ) // EIGEN_UPDATE
   {
      j=npoints-2;
      jdone = true;
   }

  if (!jdone) {
     for(i=0;i<npoints-1;i++) {  // EIGEN_UPDATE
       if( x >= pointx(i) && x <=pointx(i+1) )
       {
         j=i;
         break;
       }
     }
  }

 y = pointy.col(j) + (x-pointx(j))*(pointy.col(j+1)-pointy.col(j))/(pointx(j+1)-pointx(j));

}

void linear_interpolation(adouble* y, adouble x, MatrixXd& pointx, MatrixXd& pointy, int npoints)
{
//    Linear interpolation from point values

   int i,j;
   bool jdone = false;

   if( x.value() < pointx(0) ) // EIGEN_UPDATE
   {
      j = 0;
      jdone = true;
   }



  if ( x.value() > pointx(npoints-1 ) )  // EIGEN_UPDATE
   {
      j=npoints-2;
      jdone = true;
   }

  if (!jdone) {
     for(i=0;i<npoints-1;i++) {   // EIGEN_UPDATE
       if( x.value() >= pointx(i) && x <=pointx(i+1) )
       {
         j=i;
         break;
       }
     }
  }

  *y = pointy(j) + (x-pointx(j))*(pointy(j+1)-pointy(j))/(pointx(j+1)-pointx(j));

}


void linear_interpolation(MatrixXd& y, MatrixXd& x, MatrixXd& pointx, MatrixXd& pointy, int npoints)
{
//    Linear interpolation from point values

   int i,j=0, k;
   MatrixXd pyj(  pointy.rows(),1 );
   MatrixXd pyj1( pointy.rows(),1 );

   for(k=0;k<length(x);k++)  {  // EIGEN_UPDATE

		bool jdone = false;

		if( x(k) < pointx(0) ) // EIGEN_UPDATE
		{
	   	j = 0;
	   	jdone = true;
		}



		if ( x(k) > pointx(npoints-1) ) // EIGEN_UPDATE
		{
	    	j=npoints-2;
	    	jdone = true;
		}

		if (!jdone) {
	   	for(i=0;i<npoints-1;i++) {
	     		if( x(k) >= pointx(i) && x(k) <=pointx(i+1) )
	     		{
		      	j=i;
	       		break;
	     		}
	   	}
		}


   	y.resize( pointy.rows(), length(x) );


   	pyj  = pointy.col(j);
 
   	pyj1 = pointy.col(j+1);


    	y.col(k) = pyj + (x(k)-pointx(j))*(pyj1-pyj)/(pointx(j+1)-pointx(j));

  }

}


void smooth_linear_interpolation(adouble* y, adouble& x, MatrixXd& Xdata, MatrixXd& Ydata, int n)
{
  // Given MatrixXd objects Xdata and Ydata that contain a tabulated function, Ydata(i) = f(Xdata(i)), with
  // i=1,...n, and  x(j+1)>x(j),   given a value of x, this function returns the interpolated value y
  // using a smoothed linear interpolation
  int i;
  adouble Si, Li;
  double a, m;
  double epsilon = 0.1;

  *y=0.0;

  if ( length(Xdata) != length(Ydata) ) {
     error_message("The legnth of vectors Xdata and Ydata must be the same in function smooth_linear_interpolation()");
  }

  if( x.value() < Xdata(0) || x.value() > Xdata(length(Xdata)-1))
  {
      error_message("Extrapolation not allowed in function smooth_linear_interpolation()");
  }


  for (i=0; i<=n-2; i++) {  // EIGEN_UPDATE: i index shifted by -1

    if ( Xdata(i+1)<= Xdata(i) ) error_message("Xdata vector must be strictly monotonic increasing in smooth_linear_interpolation()");
    if (i==0) {
       a  = epsilon* fabs( Xdata(i)-Xdata(i+1));
    }
    else {
       a = epsilon * min( fabs( Xdata(i-1)-Xdata(i)), fabs( Xdata(i)-Xdata(i+1)) );
    }

    Si = smooth_heaviside( (x-Xdata(i)), a ) - smooth_heaviside( (x-Xdata(i+1)), a );

    m  = (Ydata(i+1)-Ydata(i))/(Xdata(i+1)-Xdata(i));

    Li = Ydata(i) + m*(x-Xdata(i));

    *y += Si*Li;

  }


}



void smooth_linear_interpolation(adouble* y, adouble& x, adouble* Xdata, adouble* Ydata, int n)
{
  // Given MatrixXd objects Xdata and Ydata that contain a tabulated function, Ydata(i) = f(Xdata(i)), with
  // i=1,...n, and  x(j+1)>x(j),   given a value of x, this function returns the interpolated value y
  // using a smoothed linear interpolation
  int i;
  adouble Si, Li;
  double a;
  adouble m;
  double epsilon = 0.1;

  *y=0.0;

  if( x.value() < Xdata[0].value() || x.value() > Xdata[n-1].value())
  {
//      error_message("Extrapolation not allowed in function smooth_linear_interpolation()");
  }


  for (i=1; i<=n-1; i++) {

    if ( Xdata[i]<= Xdata[i-1] ) error_message("Xdata vector must be strictly monotonic increasing in smooth_linear_interpolation()");
    if (i==1) {
       a  = epsilon* fabs( Xdata[i-1].value()-Xdata[i].value());
    }
    else {
      a = epsilon * min( fabs( Xdata[i-2].value()-Xdata[i-1].value()), fabs( Xdata[i-1].value()-Xdata[i].value()) );
    }

    Si = smooth_heaviside( (x-Xdata[i-1]), a ) - smooth_heaviside( (x-Xdata[i]), a );

    m  = (Ydata[i]-Ydata[i-1])/(Xdata[i]-Xdata[i-1]);

    Li = Ydata[i-1] + m*(x-Xdata[i-1]);

    *y += Si*Li;

  }

}


void zoh_interpolation(adouble* y, adouble x, MatrixXd& pointx, MatrixXd& pointy, int npoints)
{
//    Zero Order Hold interpolation from point values

   int i,j;
   bool jdone = false;

   if( x.value() < pointx(0) ) // EIGEN_UPDATE
   {
      j = 0;
      jdone = true;
   }



  if ( x.value() > pointx(npoints-1) )  // EIGEN_UPDATE
   {
      j=npoints-2;
      jdone = true;
   }

  if (!jdone) {
     for(i=0;i<npoints-1;i++) { // EIGEN_UPDATE: index i shifted by -1
       if( x.value() >= pointx(i) && x <=pointx(i+1) )
       {
         j=i;
         break;
       }
     }
  }

  *y = pointy(j);

}


void spline_2d_interpolation(adouble* z, adouble& x, adouble& y, MatrixXd& X, MatrixXd& Y, MatrixXd& Z, Workspace* workspace)
{
//    Spline 2d interpolation
//    Output: z: a pointer to a single adouble variable. z must point to a valid address on entry to the function
//    Inputs:
//    The adouble variables (x,y) are the values at which the interpolated value of the function is returned.
//    X is a vector of dimension nxpoints x 1
//    Y is a vector of dimension nypoints x 1
//    Z is a matrix of dimensions nxpoints x nypoints
//    Each element Z(i,j) corresponds to the pair ( X(i), Y(j) )
//    Method: Cubic spline 2D interpolation
//    The function does not deal with sparse data.
   int i,j;
   long nxpoints = length(X);
   long nypoints = length(Y);
   long nrowsZ   = Z.rows();
   long ncolsZ   = Z.cols();

   if ( nrowsZ != nxpoints ) {
       error_message("Number of rows of matrix Z must be equal to the length of vector X in call to bilinear_interpolation()");
   }
   if ( ncolsZ != nypoints )  {
         error_message("Number of columns of matrix Z must be equal to the length of vector Y in call to bilinear_interpolation()");
   }

   AutoDiffMatrix Arow(1,nypoints);
   AutoDiffMatrix yval(1,nypoints);
   AutoDiffMatrix xval(1,nxpoints);
   AutoDiffMatrix Zval(1,nxpoints);

   adouble *zval = Zval.GetPr();

   for(j=0;j<nypoints;j++) { // EIGEN_UPDATE: index j shifted by -1
          yval(0,j) = Y(j);
   }

   for(j=0;j<nxpoints;j++) { // EIGEN_UPDATE: index j shifted by -1
          xval(0,j) = X(j);
   }

   for (i=0;i<nxpoints;i++) {   // EIGEN_UPDATE: indices i,j shifted by -1

       for(j=0;j<nypoints;j++) {
          Arow(0,j) = Z(i,j);
       }

       spline_interpolation( &zval[i], y, yval.GetPr(), Arow.GetPr(), nypoints, workspace);

   }

   spline_interpolation( z, x, xval.GetPr(), zval, nxpoints, workspace);

}

void smooth_bilinear_interpolation(adouble* z, adouble& x, adouble& y, MatrixXd& X, MatrixXd& Y, MatrixXd& Z)
{
//    Spline 2d interpolation
//    Output: z: a pointer to a single adouble variable. z must point to a valid address on entry to the function
//    Inputs:
//    The adouble variables (x,y) are the values at which the interpolated value of the function is returned.
//    X is a vector of dimension nxpoints x 1
//    Y is a vector of dimension nypoints x 1
//    Z is a matrix of dimensions nxpoints x nypoints
//    Each element Z(i,j) corresponds to the pair ( X(i), Y(j) )
//    Method: Cubic spline 2D interpolation
//    The function does not deal with sparse data.
   int i,j;
   long nxpoints = length(X);
   long nypoints = length(Y);
   long nrowsZ   = Z.rows();
   long ncolsZ   = Z.cols();

   if ( nrowsZ != nxpoints ) {
       error_message("Number of rows of matrix Z must be equal to the length of vector X in call to bilinear_interpolation()");
   }
   if ( ncolsZ != nypoints )  {
         error_message("Number of columns of matrix Z must be equal to the length of vector Y in call to bilinear_interpolation()");
   }

   AutoDiffMatrix Arow(1,nypoints);
   AutoDiffMatrix yval(1,nypoints);
   AutoDiffMatrix xval(1,nxpoints);
   AutoDiffMatrix Zval(1,nxpoints);

   adouble *zval = Zval.GetPr();

   for(j=0;j<nypoints;j++) {  // EIGEN_UPDATE: index j shifted by -1
          yval(0,j) = Y(j);
   }

   for(j=0;j<nxpoints;j++) {  // EIGEN_UPDATE: index j shifted by -1
          xval(0,j) = X(j);
   }

   for (i=0;i<nxpoints;i++) {  // EIGEN_UPDATE: indices i,j shifted by -1

       for(j=0;j<nypoints;j++) {
          Arow(0,j) = Z(i,j);
       }

       smooth_linear_interpolation( &zval[i], y, yval.GetPr(), Arow.GetPr(), nypoints);

   }

   smooth_linear_interpolation( z, x, xval.GetPr(), zval, nxpoints);

}


void bilinear_interpolation(adouble* z, adouble& x, adouble& y, MatrixXd& X, MatrixXd& Y, MatrixXd& Z)
{
//    Linear 2d interpolation
//    Output: z: a pointer to a single adouble variable. z must point to a valid address on entry to the function
//    Inputs:
//    The adouble variables (x,y) are the values at which the interpolated value of the function is returned.
//    X is a vector of dimension nxpoints x 1
//    Y is a vector of dimension nypoints x 1
//    Z is a matrix of dimensions nxpoints x nypoints
//    Each element Z(i,j) corresponds to the pair ( X(i), Y(j) )
//    Method: Classical bilinear interpolation.
//    The function does not allow extrapolation. An error message is printed if the input pair x,y is out of range.
//    The function does not deal with sparse data.


   long i,jx, jy;
   bool jxdone = false;
   bool jydone = false;
   long nxpoints = length(X);
   long nypoints = length(Y);
   long nrowsZ   = Z.rows();
   long ncolsZ   = Z.cols();
   double x1,x2,y1,y2;
   double z11,z12,z21,z22;

   if ( nrowsZ != nxpoints ) {
       error_message("Number of rows of matrix Z must be equal to the length of vector X in call to bilinear_interpolation()");
   }
   if ( ncolsZ != nypoints )  {
         error_message("Number of columns of matrix Z must be equal to the length of vector Y in call to bilinear_interpolation()");
   }

   if( x.value() < X(0) )  // EIGEN_UPDATE
   {
      jx = 0;
      jxdone = true;

//      error_message("Extrapolation not allowed in function bilinear_interpolation()");
   }


  if ( x.value() > X(nxpoints-1) ) // EIGEN_UPDATE
   {
      jx=nxpoints-2;
      jxdone = true;

//      error_message("Extrapolation not allowed in function bilinear_interpolation()");

   }

  if (!jxdone) {
     for(i=0;i<nxpoints-1;i++) {  // EIGEN_UPDATE: index i shifted by -1
       if( x.value() >= X(i) && x <=X(i+1) )
       {
         jx=i;
         break;
       }
     }
  }

  if( y.value() < Y(0) ) // EIGEN_UPDATE
   {
      jy = 0;
      jydone = true;
//      error_message("Extrapolation not allowed in function bilinear_interpolation()");
   }



  if ( x.value() > Y(nypoints-1) )
   {
      jy=nypoints-2;
      jydone = true;
//      error_message("Extrapolation not allowed in function bilinear_interpolation()");
   }

  if (!jydone) {
     for(i=0;i<nypoints-1;i++) {
       if( y.value() >= Y(i) && y <=Y(i+1) )
       {
         jy=i;
         break;
       }
     }
  }

  x1  = X(jx);        x2 = X(jx+1);
  y1  = Y(jy);        y2 = Y(jy+1);
  z11 = Z(jx,jy);    z12 = Z(jx,jy+1);
  z21 = Z(jx+1,jy);  z22 = Z(jx+1,jy+1);

  *z = z11*(x2-x)*(y2-y)/((x2-x1)*(y2-y1));
  *z+= z21*(x-x1)*(y2-y)/((x2-x1)*(y2-y1));
  *z+= z12*(x2-x)*(y-y1)/((x2-x1)*(y2-y1));
  *z+= z22*(x-x1)*(y-y1)/((x2-x1)*(y2-y1));

}


void spline_second_derivative(adouble* x, adouble* y, int n,  adouble* d2y, Workspace* workspace)
// Given arrays x[] and y[] that contain a tabulated function, y[i] = f(x[i]), with
// i=0,...n-1, and  x[j+1]>x[j],  this function returns an array d2y[] that contains
// the second derivatives of the interpolating cubic polynomial at the points x[i].
//
// Reference: Burden and Faires (2005) "Numerical Analysis". Thompson.
{
      int i,j;

      adouble *c = d2y;
      adouble hi;
      adouble him1;
      adouble alphai;
      adouble li = 0;

      adouble* mu= workspace->u_spline;
      adouble*  z= workspace->z_spline;

      mu[0] = 0.0;
      z[0]  = 0.0;

      for(i=1; i<n-1;i++) {
	 		him1 = x[i]-x[i-1];
	 		hi   = x[i+1]-x[i];
	 		alphai = 3.0/hi*(y[i+1]-y[i])-3.0/him1*(y[i]-y[i-1]);
	 		li = 2*(x[i+1]-x[i-1])-him1*mu[i-1];
	 		mu[i] = hi/li;
	 		z[i] = (alphai-him1*z[i-1])/li;
      }

      c[n-1] = 0.0;

      for(j=n-2;j>=0;j--) {
	  		c[j] = z[j]-mu[j]*c[j+1];
      }
      for(j=1;j<n-1;j++) {
         c[j] = 2*c[j];
      }


}




void spline_second_derivative(MatrixXd& xdata, MatrixXd& ydata, int n,  MatrixXd& d2y)
// Given arrays x[] and y[] that contain a tabulated function, y[i] = f(x[i]), with
// i=0,...n-1, and  x[j+1]>x[j],  this function returns an array d2y[] that contains
// the second derivatives of the interpolating cubic polynomial at the points x[i].
//
// Reference: Burden and Faires (2005) "Numerical Analysis". Thompson.
{
      int i,j;


      double *c = &d2y(0);
      double hi;
      double him1;
      double alphai;
      double li = 0;

      double *x=&xdata(0);


      double *y=&ydata(0);

      MatrixXd MU(1,n);
      MatrixXd Z(1,n);


      double* mu = &MU(0);

      double* z  = &Z(0);

      mu[0] = 0.0;
      z[0]  = 0.0;

      for(i=1; i<n-1;i++) {
	      him1 = x[i]-x[i-1];
	      hi   = x[i+1]-x[i];
	      alphai = 3.0/hi*(y[i+1]-y[i])-3.0/him1*(y[i]-y[i-1]);
	      li = 2*(x[i+1]-x[i-1])-him1*mu[i-1];
	      mu[i] = hi/li;
	      z[i] = (alphai-him1*z[i-1])/li;
      }

      c[n-1] = 0.0;

      for(j=n-2;j>=0;j--) {
	       c[j] = z[j]-mu[j]*c[j+1];
      }
      for(j=1;j<n-1;j++) {
          c[j] = 2*c[j];
      }


}


void spline_interpolation(adouble* y, adouble& x, adouble* xdata, adouble* ydata, int n, Workspace* workspace)
//   Given the arrays xdata[i] and ydata[i], i = 0,...n-1, which tabulate a function, with xdata[i] < xdata[i+1],
//   and given the array d2y[i], which is the output from function spline_second_derivative(),
//   and given a value of x, this function returns the interpolated value y using (natural) cubic-spline interpolation
//
//   Reference: Burden and Faires (2005) "Numerical Analysis". Thompson.
{

   int kleft,kright,k;
   adouble h,A,B,C,D;
   adouble *d2y = workspace->y2a_spline;
   kleft=1;

   spline_second_derivative(xdata, ydata, n, d2y, workspace );

   kright=n;
   while (kright-kleft > 1) {
      k=(int) ( (kright+kleft)/2 );
      if (xdata[k-1] > x) kright=k;
      else kleft=k;
   }
   h=xdata[kright-1]-xdata[kleft-1];
   if (h == 0.0) error_message("Bad xdata input to routine spline_interpolation()");
   A=(xdata[kright-1]-x)/h;
   B=(x-xdata[kleft-1])/h;
   C=(pow(A,3)-A)*(h*h)/6.0;
   D=(pow(B,3)-B)*(h*h)/6.0;
   // Evaluate the cubic spline polynomial
   *y=A*ydata[kleft-1]+B*ydata[kright-1]+C*d2y[kleft-1]+D*d2y[kright-1];
}

void spline_interpolation(adouble* y, adouble& x, MatrixXd& Xdata, MatrixXd& Ydata, int n)
//   Given the arrays xdata[i] and ydata[i], i = 0,...n-1, which tabulate a function, with xdata[i] < xdata[i+1],
//   and given the array d2y[i], which is the output from function spline_second_derivative(),
//   and given a value of x, this function returns the interpolated value y using (natural) cubic-spline interpolation
//
//   Reference: Burden and Faires (2005) "Numerical Analysis". Thompson.
{

   int kleft,kright,k;
   adouble h,A,B,C,D;
   MatrixXd D2Y(1,n);

   double* d2y = &D2Y(0);

   double *xdata = &Xdata(0);

   double *ydata = &Ydata(0);
   kleft=1;

   spline_second_derivative(Xdata, Ydata, n, D2Y );

   kright=n;
   while (kright-kleft > 1) {
      k=(int) ( (kright+kleft)/2 );
      if (xdata[k-1] > x) kright=k;
      else kleft=k;
   }
   h=xdata[kright-1]-xdata[kleft-1];
   if (h == 0.0) error_message("Bad xdata input to routine spline_interpolation()");
   A=(xdata[kright-1]-x)/h;
   B=(x-xdata[kleft-1])/h;
   C=(pow(A,3)-A)*(h*h)/6.0;
   D=(pow(B,3)-B)*(h*h)/6.0;
   // Evaluate the cubic spline polynomial
   *y=A*ydata[kleft-1]+B*ydata[kright-1]+C*d2y[kleft-1]+D*d2y[kright-1];
}






void spline_interpolation(MatrixXd& Y, MatrixXd& X, MatrixXd& Xdata, MatrixXd& Ydata, int n)
//   Given the arrays xdata[i] and ydata[i], i = 0,...n-1, which tabulate a function, with xdata[i] < xdata[i+1],
//   and given the array d2y[i], which is the output from function spline_second_derivative(),
//   and given a vector of independent varaibles X, this function returns a vector of interpolated values Y
//   using (natural) cubic-spline interpolation
//
//   Reference: Burden and Faires (2005) "Numerical Analysis". Thompson.
{



     double *xdata = &Xdata(0);

     double *ydata = &Ydata(0);
   int kleft,kright,k;
   double h,A,B,C,D;
   MatrixXd D2Y(1,n);

   double* d2y = &D2Y(0);
   int i;

   spline_second_derivative(Xdata, Ydata, n, D2Y );


   for(i=0;i< length(X); i++) { // EIGEN_UPDATE
      kleft=1;
      kright=n;
      while (kright-kleft > 1) {
	  k=(int) ( (kright+kleft)/2 );
	  if (xdata[k-1] > X(i)) kright=k;
	  else kleft=k;
      }
      h=xdata[kright-1]-xdata[kleft-1];
      if (h == 0.0) error_message("Bad xdata input to routine spline_interpolation()");
      A=(xdata[kright-1]-X(i))/h;
      B=(X(i)-xdata[kleft-1])/h;
      C=(pow(A,3)-A)*(h*h)/6.0;
      D=(pow(B,3)-B)*(h*h)/6.0;
      // Evaluate the cubic spline polynomial
      Y(i)=A*ydata[kleft-1]+B*ydata[kright-1]+C*d2y[kleft-1]+D*d2y[kright-1];


   } // End for loop
}





void spline_interpolation_with_second_derivative(adouble* y, adouble& x, adouble* xdata, adouble* ydata, adouble* y2d, int n, Workspace* workspace)
//   Given the arrays xdata[i] and ydata[i], i = 0,...n-1, which tabulate a function, with xdata[i] < xdata[i+1],
//   and given the array d2y[i], which is the output from function spline_second_derivative(),
//   and given a value of x, this function returns the interpolated value y using (natural) cubic-spline interpolation
//
//   Reference: Burden and Faires (2005) "Numerical Analysis". Thompson.
{

   int kleft,kright,k;
   adouble h,A,B,C,D;
   adouble *d2y = workspace->y2a_spline;
   kleft=1;

   spline_second_derivative(xdata, ydata, n, d2y, workspace );

   kright=n;
   while (kright-kleft > 1) {
      k=(int) ( (kright+kleft)/2 );
      if (xdata[k-1] > x) kright=k;
      else kleft=k;
   }
   h=xdata[kright-1]-xdata[kleft-1];
   if (h == 0.0) error_message("Bad xdata input to routine spline_interpolation()");
   A=(xdata[kright-1]-x)/h;
   B=(x-xdata[kleft-1])/h;
   C=(pow(A,3)-A)*(h*h)/6.0;
   D=(pow(B,3)-B)*(h*h)/6.0;
   // Evaluate the cubic spline polynomial
   *y=A*ydata[kleft-1]+B*ydata[kright-1]+C*d2y[kleft-1]+D*d2y[kright-1];
}


void resample_trajectory(MatrixXd& Y,  MatrixXd& X, MatrixXd& Ydata, MatrixXd& Xdata )
{
//  This function resamples a trajectory (Xdata,Ydata) given new values of the time vector Xdata using
//   natural cubic spline interpolation. The interpolated values are returned in Y.
//   Xdata has dimensions 1 x N
//   Ydata has dimensions ny x N
//   X has dimension 1 x M
//   On output, Y has dimensions ny x M
//   Xdata and X should be monotonically increasing vectors.
    int i;

    long ny = Ydata.rows();

    long lx = length(X);

    long n = length(Xdata);

    MatrixXd Yi, Yidata;

    if ( (X(0) < Xdata(0)) || X(length(X)-1) > Xdata(length(Xdata)-1) ) {                   // EIGEN_UPDATE
         error_message("No extrapolation is allowed in function resample_trajectory()");
    }

    Y.resize(ny,lx);
    Yi.resize(1,lx);
    Yidata.resize(1,n);

    for(i=0;i<ny;i++) {

         Yidata = Ydata.row(i);
      	spline_interpolation( Yi, X, Xdata, Yidata, n);

         Y.row(i) = Yi;
    }
}


