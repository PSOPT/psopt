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



void lagrange_interpolation(DMatrix& y, DMatrix& x, DMatrix& pointx, DMatrix& pointy)
{
//
//       This function approximates a point-defined function using Lagrange polynomial interpolation
//
//       lagrange_interpolation(Y, X,POINTX,POINTY) approximates the function defined by the points:
//       P1=(POINTX(1),POINTY(1)), P2=(POINTX(2),POINTY(2)), ..., PN(POINTX(N),POINTY(N))
//       and calculate it in each element of X


   int i,j;

   long n = length(pointx);

   DMatrix L(n, x.GetNoCols() );

   for(i=1;i<=n;i++) {
	for(j=1;j<=x.GetNoCols();j++) {
            L(i,j) = 1.0;
	}
    }

   for (i=1;i<=n;i++) {
      for (j=1;j<=n;j++) {
           if (i != j) {
                L(i,colon()) = ( L(i,colon())&(  x-pointx(j)*ones(1,length(x)) )  )/(pointx(i)-pointx(j));
           }
      }
   }

   y.Resize(1,length(x));
   y.FillWithZeros();
   for (i=1;i<=n;i++) {
   	y = y+pointy(i)*L(i,colon());
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

void linear_interpolation(DMatrix& y, double x, DMatrix& pointx, DMatrix& pointy, int npoints)
{
//    Linear interpolation from point values

   int i,j;
   bool jdone = false;

   if( x < pointx(1) )
   {
      j = 1;
      jdone = true;
   }



  if ( x > pointx(npoints) )
   {
      j=npoints-1;
      jdone = true;
   }

  if (!jdone) {
     for(i=1;i<npoints;i++) {
       if( x >= pointx(i) && x <=pointx(i+1) )
       {
         j=i;
         break;
       }
     }
  }

  y = pointy(colon(),j) + (x-pointx(j))*(pointy(colon(),j+1)-pointy(colon(),j))/(pointx(j+1)-pointx(j));

}

void linear_interpolation(adouble* y, adouble x, DMatrix& pointx, DMatrix& pointy, int npoints)
{
//    Linear interpolation from point values

   int i,j;
   bool jdone = false;

   if( x.value() < pointx(1) )
   {
      j = 1;
      jdone = true;
   }



  if ( x.value() > pointx(npoints) )
   {
      j=npoints-1;
      jdone = true;
   }

  if (!jdone) {
     for(i=1;i<npoints;i++) {
       if( x.value() >= pointx(i) && x <=pointx(i+1) )
       {
         j=i;
         break;
       }
     }
  }

  *y = pointy(j) + (x-pointx(j))*(pointy(j+1)-pointy(j))/(pointx(j+1)-pointx(j));

}


void linear_interpolation(DMatrix& y, DMatrix& x, DMatrix& pointx, DMatrix& pointy, int npoints)
{
//    Linear interpolation from point values

   int i,j=0, k;
   DMatrix pyj( pointy.GetNoRows(),1);
   DMatrix pyj1(pointy.GetNoRows(),1);

   for(k=1;k<=length(x);k++)  {

	bool jdone = false;

	if( x(k) < pointx(1) )
	{
	   j = 1;
	   jdone = true;
	}



	if ( x(k) > pointx(npoints) )
	{
	    j=npoints-1;
	    jdone = true;
	}

	if (!jdone) {
	   for(i=1;i<npoints;i++) {
	     if( x(k) >= pointx(i) && x(k) <=pointx(i+1) )
	     {
		j=i;
		break;
	     }
	   }
	}

//        DMatrix temp= pointy(colon(),j) + (x(k)-pointx(j))*(pointy(colon(), j+1)-pointy(colon(),j))/(pointx(j+1)-pointx(j));

        y.Resize( pointy.GetNoRows(), length(x) );

	pyj  = pointy(colon(),j);
	pyj1 = pointy(colon(),j+1);

	y(colon(),k) = pyj + (x(k)-pointx(j))*(pyj1-pyj)/(pointx(j+1)-pointx(j));

  }

}


void smooth_linear_interpolation(adouble* y, adouble& x, DMatrix& Xdata, DMatrix& Ydata, int n)
{
  // Given DMatrix objects Xdata and Ydata that contain a tabulated function, Ydata(i) = f(Xdata(i)), with
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

  if( x.value() < Xdata(1) || x.value() > Xdata("end"))
  {
//      error_message("Extrapolation not allowed in function smooth_linear_interpolation()");
  }


  for (i=1; i<=n-1; i++) {

    if ( Xdata(i+1)<= Xdata(i) ) error_message("Xdata vector must be strictly monotonic increasing in smooth_linear_interpolation()");
    if (i==1) {
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

  Si = smooth_heaviside( (x-(Xdata(1)-100*a)), a ) - smooth_heaviside( (x-Xdata(1)), a );
  *y+=Si*Ydata(1);
  Si = smooth_heaviside( (x-Xdata(n)), a ) - smooth_heaviside( (x-(Xdata(n)+100*a)), a );
  *y+=Si*Ydata(n);

}



void smooth_linear_interpolation(adouble* y, adouble& x, adouble* Xdata, adouble* Ydata, int n)
{
  // Given DMatrix objects Xdata and Ydata that contain a tabulated function, Ydata(i) = f(Xdata(i)), with
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


//  Si = smooth_heaviside( (x-(Xdata[0]-100*a)), a ) - smooth_heaviside( (x-Xdata[0]), a );
//  *y+=Si*Ydata[0];
//  Si = smooth_heaviside( (x-Xdata[n-1]), a ) - smooth_heaviside( (x-(Xdata[n-1]+100*a)), a );
//  *y+=Si*Ydata[n-1];

}


void zoh_interpolation(adouble* y, adouble x, DMatrix& pointx, DMatrix& pointy, int npoints)
{
//    Zero Order Hold interpolation from point values

   int i,j;
   bool jdone = false;

   if( x.value() < pointx(1) )
   {
      j = 1;
      jdone = true;
   }



  if ( x.value() > pointx(npoints) )
   {
      j=npoints-1;
      jdone = true;
   }

  if (!jdone) {
     for(i=1;i<npoints;i++) {
       if( x.value() >= pointx(i) && x <=pointx(i+1) )
       {
         j=i;
         break;
       }
     }
  }

  *y = pointy(j);

}


void spline_2d_interpolation(adouble* z, adouble& x, adouble& y, DMatrix& X, DMatrix& Y, DMatrix& Z, Workspace* workspace)
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
   int i,j, jx, jy;
   bool jxdone = false;
   bool jydone = false;
   int nxpoints = length(X);
   int nypoints = length(Y);
   int nrowsZ   = Z.GetNoRows();
   int ncolsZ   = Z.GetNoCols();

   if ( nrowsZ != nxpoints ) {
       error_message("Number of rows of matrix Z must be equal to the length of vector X in call to bilinear_interpolation()");
   }
   if ( ncolsZ != nypoints )  {
         error_message("Number of columns of matrix Z must be equal to the length of vector Y in call to bilinear_interpolation()");
   }

   ADMatrix Arow(1,nypoints);
   ADMatrix yval(1,nypoints);
   ADMatrix xval(1,nxpoints);
   ADMatrix Zval(1,nxpoints);

   adouble *zval = Zval.GetPr();

   for(j=1;j<=nypoints;j++) {
          yval(1,j) = Y(j);
   }

   for(j=1;j<=nxpoints;j++) {
          xval(1,j) = X(j);
   }

   for (i=1;i<=nxpoints;i++) {

       for(j=1;j<=nypoints;j++) {
          Arow(1,j) = Z(i,j);
       }

       spline_interpolation( &zval[i-1], y, yval.GetPr(), Arow.GetPr(), nypoints, workspace);

   }

   spline_interpolation( z, x, xval.GetPr(), zval, nxpoints, workspace);

}

void smooth_bilinear_interpolation(adouble* z, adouble& x, adouble& y, DMatrix& X, DMatrix& Y, DMatrix& Z)
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
   int i,j, jx, jy;
   bool jxdone = false;
   bool jydone = false;
   int nxpoints = length(X);
   int nypoints = length(Y);
   int nrowsZ   = Z.GetNoRows();
   int ncolsZ   = Z.GetNoCols();

   if ( nrowsZ != nxpoints ) {
       error_message("Number of rows of matrix Z must be equal to the length of vector X in call to bilinear_interpolation()");
   }
   if ( ncolsZ != nypoints )  {
         error_message("Number of columns of matrix Z must be equal to the length of vector Y in call to bilinear_interpolation()");
   }

   ADMatrix Arow(1,nypoints);
   ADMatrix yval(1,nypoints);
   ADMatrix xval(1,nxpoints);
   ADMatrix Zval(1,nxpoints);

   adouble *zval = Zval.GetPr();

   for(j=1;j<=nypoints;j++) {
          yval(1,j) = Y(j);
   }

   for(j=1;j<=nxpoints;j++) {
          xval(1,j) = X(j);
   }

   for (i=1;i<=nxpoints;i++) {

       for(j=1;j<=nypoints;j++) {
          Arow(1,j) = Z(i,j);
       }

       smooth_linear_interpolation( &zval[i-1], y, yval.GetPr(), Arow.GetPr(), nypoints);

   }

   smooth_linear_interpolation( z, x, xval.GetPr(), zval, nxpoints);

}


void bilinear_interpolation(adouble* z, adouble& x, adouble& y, DMatrix& X, DMatrix& Y, DMatrix& Z)
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


   int i,jx, jy;
   bool jxdone = false;
   bool jydone = false;
   int nxpoints = length(X);
   int nypoints = length(Y);
   int nrowsZ   = Z.GetNoRows();
   int ncolsZ   = Z.GetNoCols();
   double x1,x2,y1,y2;
   double z11,z12,z21,z22;

   if ( nrowsZ != nxpoints ) {
       error_message("Number of rows of matrix Z must be equal to the length of vector X in call to bilinear_interpolation()");
   }
   if ( ncolsZ != nypoints )  {
         error_message("Number of columns of matrix Z must be equal to the length of vector Y in call to bilinear_interpolation()");
   }

   if( x.value() < X(1) )
   {
      jx = 1;
      jxdone = true;

//      error_message("Extrapolation not allowed in function bilinear_interpolation()");
   }


  if ( x.value() > X(nxpoints) )
   {
      jx=nxpoints-1;
      jxdone = true;

//      error_message("Extrapolation not allowed in function bilinear_interpolation()");

   }

  if (!jxdone) {
     for(i=1;i<nxpoints;i++) {
       if( x.value() >= X(i) && x <=X(i+1) )
       {
         jx=i;
         break;
       }
     }
  }

  if( y.value() < Y(1) )
   {
      jy = 1;
      jydone = true;
//      error_message("Extrapolation not allowed in function bilinear_interpolation()");
   }



  if ( x.value() > Y(nypoints) )
   {
      jy=nypoints-1;
      jydone = true;
//      error_message("Extrapolation not allowed in function bilinear_interpolation()");
   }

  if (!jydone) {
     for(i=1;i<nypoints;i++) {
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




void spline_second_derivative(DMatrix& xdata, DMatrix& ydata, int n,  DMatrix& d2y)
// Given arrays x[] and y[] that contain a tabulated function, y[i] = f(x[i]), with
// i=0,...n-1, and  x[j+1]>x[j],  this function returns an array d2y[] that contains
// the second derivatives of the interpolating cubic polynomial at the points x[i].
//
// Reference: Burden and Faires (2005) "Numerical Analysis". Thompson.
{
      int i,j;

      double *c = d2y.GetPr();
      double hi;
      double him1;
      double alphai;
      double li = 0;
      double *x=xdata.GetPr();
      double *y=ydata.GetPr();

      DMatrix MU(1,n);
      DMatrix Z(1,n);

      double* mu= MU.GetPr();
      double*  z= Z.GetPr();

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

void spline_interpolation(adouble* y, adouble& x, DMatrix& Xdata, DMatrix& Ydata, int n)
//   Given the arrays xdata[i] and ydata[i], i = 0,...n-1, which tabulate a function, with xdata[i] < xdata[i+1],
//   and given the array d2y[i], which is the output from function spline_second_derivative(),
//   and given a value of x, this function returns the interpolated value y using (natural) cubic-spline interpolation
//
//   Reference: Burden and Faires (2005) "Numerical Analysis". Thompson.
{

   int kleft,kright,k;
   adouble h,A,B,C,D;
   DMatrix D2Y(1,n);
   double* d2y = D2Y.GetPr();
   double *xdata = Xdata.GetPr();
   double *ydata = Ydata.GetPr();

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






void spline_interpolation(DMatrix& Y, DMatrix& X, DMatrix& Xdata, DMatrix& Ydata, int n)
//   Given the arrays xdata[i] and ydata[i], i = 0,...n-1, which tabulate a function, with xdata[i] < xdata[i+1],
//   and given the array d2y[i], which is the output from function spline_second_derivative(),
//   and given a vector of independent varaibles X, this function returns a vector of interpolated values Y
//   using (natural) cubic-spline interpolation
//
//   Reference: Burden and Faires (2005) "Numerical Analysis". Thompson.
{


   double *xdata = Xdata.GetPr();
   double *ydata = Ydata.GetPr();
   int kleft,kright,k;
   double h,A,B,C,D;
   DMatrix D2Y(1,n);
   double* d2y = D2Y.GetPr();
   int i;

   spline_second_derivative(Xdata, Ydata, n, D2Y );


   for(i=1;i<= length(X); i++) {
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


void resample_trajectory(DMatrix& Y,  DMatrix& X, DMatrix& Ydata, DMatrix& Xdata )
{
//  This function resamples a trajectory (Xdata,Ydata) given new values of the time vector Xdata using
//   natural cubic spline interpolation. The interpolated values are returned in Y.
//   Xdata has dimensions 1 x N
//   Ydata has dimensions ny x N
//   X has dimension 1 x M
//   On output, Y has dimensions ny x M
//   Xdata and X should be monotonically increasing vectors.
    int i;

    int ny = Ydata.GetNoRows();

    int lx = length(X);

    int n = length(Xdata);

    DMatrix Yi, Yidata;

    if ( (X(1) < Xdata(1)) || X("end") > Xdata("end") ) {
         error_message("No extrapolation is allowed in function resample_trajectory()");
    }

    Y.Resize(ny,lx);
    Yi.Resize(1,lx);
    Yidata.Resize(1,n);

    for(i=1;i<=ny;i++) {
        Yidata = Ydata(i,colon());
      	spline_interpolation( Yi, X, Xdata, Yidata, n);
        Y(i,colon()) = Yi;
    }
}


