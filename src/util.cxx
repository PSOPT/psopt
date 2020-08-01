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

using namespace std;
using namespace Eigen;

adouble dot(adouble* x, adouble* y, int n)
{
   // Dot product of two vectors of adoubles with n elements
   adouble retval = 0.0;

   for(int j=0; j<n; j++)
        retval += x[j]*y[j];

   return (retval);
}

void cross(adouble* x, adouble* y, adouble* z)
{
   // returns in z the cross product of two vectors of adoubles x and y.

   z[0] = x[1]*y[2]-x[2]*y[1];
   z[1] = x[2]*y[0]-x[0]*y[2];
   z[2] = x[0]*y[1]-x[1]*y[0];
}


adouble smooth_heaviside(adouble x, double a)
{
   // Smooth Heaviside function
   return (0.5*(1.0+tanh( x/a ) ));

//   if (x>=0) return 1.0;
//   else return 0.0;
}

adouble smooth_sign(adouble x, double a)
{
   // Smooth sign function

   return ( tanh( x/a ) );
}

double smooth_fabs(double x, double eps)
{
// This is a differentiable approximation to the absolute value function
//      return (1.0-exp( -(1.0/eps)*fabs(x)))*fabs(x);

   return (sqrt( x*x + eps*eps ) );

}

adouble smooth_fabs(adouble x, double eps)
{
// This is a differentiable approximation to the absolute value function
//      return (1.0-exp( -(1.0/eps)*fabs(x)))*fabs(x);
     return (sqrt( x*x + eps*eps ) );

}

void chronometer_tic(Workspace* workspace)
{
   workspace->start_ticks = clock();
}

double chronometer_toc(Workspace* workspace)
{
// Chronometer stops

   clock_t  stop_ticks     = clock();
   clock_t  elapsed_ticks;
   double elapsed_time;

   elapsed_ticks= stop_ticks-workspace->start_ticks;

   elapsed_time = ((double) elapsed_ticks/CLOCKS_PER_SEC);

   return (elapsed_time);

}


void get_local_time( string& date_and_time)
{
     char*  dt;
     time_t timer;
     time( &timer );
     dt = asctime( localtime( &timer ) );
     date_and_time = dt;
}

double inverse_twotailed_t_cdf(double A, int ndf)
{


MatrixXd TwoTailed_T_table(56,6);

TwoTailed_T_table <<
1.,	 	3.078,	6.314,	12.706,	31.820,	63.657,
2.,	 	1.886,	2.920,	4.303,	6.965,	9.925,
3.,	 	1.638,	2.353,	3.182,	4.541,	5.841,
4.,	 	1.533,	2.132,	2.776,	3.747,	4.604,
5.,	 	1.476,	2.015,	2.571,	3.365,	4.032,
6.,	 	1.440,	1.943,	2.447,	3.143,	3.707,
7.,	 	1.415,	1.895,	2.365,	2.998,	3.499,
8.,	 	1.397,	1.860,	2.306,	2.897,	3.355,
9.,	 	1.383,	1.833,	2.262,	2.821,	3.250,
10.,	 	1.372,	1.812,	2.228,	2.764,	3.169,
11.,	 	1.363,	1.796,	2.201,	2.718,	3.106,
12.,	 	1.356,	1.782,	2.179,	2.681,	3.055,
13.,	 	1.350,	1.771,	2.160,	2.650,	3.012,
14.,	 	1.345,	1.761,	2.145,	2.625,	2.977,
15.,	 	1.341,	1.753,	2.131,	2.602,	2.947,
16.,	 	1.337,	1.746,	2.120,	2.584,	2.921,
17.,	 	1.333,	1.740,	2.110,	2.567,	2.898,
18.,	 	1.330,	1.734,	2.101,	2.552,	2.878,
19.,	 	1.328,	1.729,	2.093,	2.539,	2.861,
20.,	 	1.325,	1.725,	2.086,	2.528,	2.845,
21.,	 	1.323,	1.721,	2.080,	2.518,	2.831,
22.,	 	1.321,	1.717,	2.074,	2.508,	2.819,
23.,	 	1.319,	1.714,	2.069,	2.500,	2.807,
24.,	 	1.318,	1.711,	2.064,	2.492,	2.797,
25.,	 	1.316,	1.708,	2.060,	2.485,	2.787,
26.,	 	1.315,	1.706,	2.056,	2.479,	2.779,
27.,	 	1.314,	1.703,	2.052,	2.473,	2.771,
28.,	 	1.313,	1.701,	2.048,	2.467,	2.763,
29.,	 	1.311,	1.699,	2.045,	2.462,	2.756,
30.,	 	1.310,	1.697,	2.042,	2.457,	2.750,
31.,	 	1.309,	1.695,	2.040,	2.453,	2.744,
32.,	 	1.309,	1.694,	2.037,	2.449,	2.738,
33.,	 	1.308,	1.692,	2.035,	2.445,	2.733,
34.,	 	1.307,	1.691,	2.032,	2.441,	2.728,
35.,	 	1.306,	1.690,	2.030,	2.438,	2.724,
36.,	 	1.306,	1.688,	2.028,	2.434,	2.719,
37.,	 	1.305,	1.687,	2.026,	2.431,	2.715,
38.,	 	1.304,	1.686,	2.024,	2.429,	2.712,
39.,	 	1.304,	1.685,	2.023,	2.426,	2.708,
40.,	 	1.303,	1.684,	2.021,	2.423,	2.704,
42.,	 	1.302,	1.682,	2.018,	2.418,	2.698,
44.,	 	1.301,	1.680,	2.015,	2.414,	2.692,
46.,	 	1.300,	1.679,	2.013,	2.410,	2.687,
48.,	 	1.299,	1.677,	2.011,	2.407,	2.682,
50.,	 	1.299,	1.676,	2.009,	2.403,	2.678,
60.,	 	1.296,	1.671,	2.000,	2.390,	2.660,
70.,	 	1.294,	1.667,	1.994,	2.381,	2.648,
80.,	 	1.292,	1.664,	1.990,	2.374,	2.639,
90.,	 	1.291,	1.662,	1.987,	2.369,	2.632,
100.,	 	1.290,	1.660,	1.984,	2.364,	2.626,
120.,	 	1.289,	1.658,	1.980,	2.358,	2.617,
150.,	 	1.287,	1.655,	1.976,	2.351,	2.609,
200.,	 	1.286,	1.652,	1.972,	2.345,	2.601,
300.,	 	1.284,	1.650,	1.968,	2.339,	2.592,
500.,	 	1.283,	1.648,	1.965,	2.334,	2.586,
Infinity,	  	1.282,	1.645,	1.960, 	2.326,	2.576;


    MatrixXd pointx;
    MatrixXd pointy;
    MatrixXd y(1,1);

    pointx = TwoTailed_T_table.col(0); // EIGEN_UPDATE:   (colon(),1);

    if (A==0.80) {
      pointy = TwoTailed_T_table.col(1); // EIGEN_UPDATE: (colon(),2);
    }
    else if( A==0.90 ) {
    pointy = TwoTailed_T_table.col(2); 
    }
    else if (A==0.95) {
    pointy = TwoTailed_T_table.col(3);
    }
    else if (A == 0.98) {
    pointy = TwoTailed_T_table.col(4);
    }
    else if (A== 0.99) {
        pointy = TwoTailed_T_table.col(5);
    }
    else {
        error_message("\n A value not supported by inverse_twotailed_t_cdf()");
    }

    long npoints = pointy.rows();

    pointx = pointx.transpose().eval();
    pointy = pointy.transpose().eval();

    linear_interpolation(y, (double) ndf, pointx, pointy, (int) npoints);

    return y(1,1);

}


double nint(double x)
{
 return ceil(x + 0.5) - (fmod(x*0.5 + 0.25, 1.0) != 0);
}



void product_ad(const MatrixXd& A, const adouble* x, int nx, adouble* y)
{
// Multiplies a constant matrix stored in MatrixXd object A (nxn) by adouble vector
// stored in array x, and returns the result in adouble array y.

  long  i,k,na;

  long  rowA;

  adouble sum;

  na=A.rows();

  const double* apr = &A(0); // EIGEN_UPDATE .GetConstPr();

  if (A.cols() != nx )

  { error_message("Incompatible matrix dimensions in product_ad()"); }


  for (i=0; i< na; i++)
  {

      sum=0.0;

      rowA = i;

      for (k=0; k<nx; k++)
      {
          if (apr[ rowA ]!=0)
            sum += apr[ rowA ]*x[ k ];
          rowA += A.rows();

      }

      y[i] = sum;

  }

  return;

}

void sum_ad(const adouble* a, const adouble*b, int n, adouble* c)
{
  // Matrix or vector sum C=A+B optimised for columnwise storage in a single array,
  // suitable for automatic differentiation.
  // Inputs: a: array with the elements of matrix A stored column by column
  //         b: array with the elements of matrix B stored column by column
  //         n: total number of elments of matrix A
  // Output: c: array with the result of the operation, stored column by column.
  int i;
  for (i=0;i<n;i++) {
     c[i] = a[i]+b[i];
  }
  return;
}


void subtract_ad(const adouble* a, const adouble*b, int n, adouble* c)
{
  // Matrix or vector subtraction C=A-B optimised for columnwise storage in a single array,
  // suitable for automatic differentiation.
  // Inputs: a: array with the elements of matrix A stored column by column
  //         b: array with the elements of matrix B stored column by column
  //         n: total number of elements of matrix A
  // Output: c: array with the result of the operation, stored column by column.
  int i;
  for (i=0;i<n;i++) {
     c[i] = a[i]-b[i];
  }
  return;
}


void inverse_ad(adouble* minput, int n, adouble* minv)
{
// Inverts the nxn matrix minput (which is stored column by column in a single array)
// and stores the result in array minv (also with colunnwise storage)
// Array minv must be pre-allocated with length(n*n).
// Inputs:  minput: n x n matrix stored columnwise
//          n: row dimension if input matrix
// Output:
//          minv: n x n matrix stored columnwise with the inverse of minput.
//
  int i,p,k;

  adouble pivot;

  for (i=0; i< n*n; i++)
	minv[i] = minput[i];


  adouble *a = minv;

  for (p=0;p<n;p++)
  {

    pivot = a[p*n+p];

    for ( i=0;i<n;i++)
    {

      if ( i!=p )
      {

        for ( k=0; k<n; k++)
        {

          if ( k!=p )
          {

            if (  pivot==0.0) error_message("singular matrix in inverse_ad()");


              a[k*n+i] +=  -a[p*n+i]*a[k*n+p]/pivot;

          }

        }

      }

    }

    for (i=0;i<n;i++)
    {

      if (i!=p)
      {


          a[p*n+i] = -a[p*n+i]/pivot;


          a[i*n+p]=a[i*n+p]/pivot;


      }

    }


    if (a[p*n+p]==0.0) error_message("singular matrix");


    a[p*n+p]=1/pivot;


  }
}

#define MTX( a, i, j, n )  (a)[((j))*(n) + (i) ]    // EIGEN_UPDATE

void transpose_ad(adouble* Apr, int na, int ma,  adouble* Atpr)
{
    int i,j;

    for(i=0;i<na;i++) {   // EIGEN_UPDATE
        for(j=0;j<ma;j++) {
	     MTX( Atpr, j, i, ma ) = MTX( Apr, i, j, na );
	}
    }
}

#undef MTX

void product_ad(adouble* Apr,adouble* Bpr, int na, int ma, int nb, int mb, adouble* ABpr)
{

  // Matrix product AB=A*B optimised for columnwise storage in a single array,
  // suitable for automatic differentiation.
  // Inputs: A: array with the elements of matrix A stored column by column
  //         B: array with the elements of matrix B stored column by column
  //         na: number of rows of matrix A
  //         ma: number of columns of matrix A
  //         nb: number of rows of matrix B
  //         mb: number of columns of matrix B
  // Output: AB: array with the result of the operation, stored column by column.

  int  i,j,k;

  int  rowA, colB;

  adouble sum;

  if (ma != nb ) error_message("Dimension error in product_ad()");



  for (i=0; i< na; i++)
  {


    for (j=0; j<mb; j++)
    {

      colB = j*nb;

      sum=0.0;

      rowA = i;

      for (k=0; k<nb; k++)
      {

          sum += Apr[ rowA ]*Bpr[ colB + k ];
          rowA += na;

      }

      ABpr[j*na+i] = sum;

    }

  }


}



AutoDiffMatrix::AutoDiffMatrix(int ninput, int minput)
{
   n = ninput;
   m = minput;

   a = new adouble[n*m];
   if (a==NULL) error_message("Allocation error in AutoDiffMatrix::AutoDiffMatrix()");


}

AutoDiffMatrix::AutoDiffMatrix(int ninput)
{
   n = ninput;
   m = 1;

   a = new adouble[n*m];

}

AutoDiffMatrix::AutoDiffMatrix()
{
   n = 0;
   m = 0;

   a = NULL;

}

adouble& AutoDiffMatrix::operator() (int i, int j)
{
  if ( i<0 || j<0 || i>n-1 || j>m-1 )  // EIGEN_UPDATE
     error_message("\nIndex error in AutoDiffMatrix::operator()");

  return a[((j)*n+i) ];
}

adouble& AutoDiffMatrix::elem(int i, int j)
{
  if ( i<0 || j<0 || i>n-1 || j>m-1 )  // EIGEN_UPDATE
     error_message("\nIndex error in AutoDiffMatrix::operator()");

  return a[((j)*n+i) ];
}

adouble& AutoDiffMatrix::operator() (int i)
{
  if ( i<0 || i>n*m-1 )  // EIGEN_UPDATE
     error_message("\nIndex error in AutoDiffMatrix::operator()");

  return a[ i ];
}

adouble AutoDiffMatrix::operator() (int i, int j) const
{
  if ( i<0 || j<0 || i>n-1 || j>m-1 )   // EIGEN_UPDATE
     error_message("\nIndex error in AutoDiffMatrix::operator()");

  return a[((j)*n+i)];
}

adouble AutoDiffMatrix::elem(int i, int j) const
{
  if ( i<0 || j<0 || i>n-1 || j>m-1 )  // EIGEN_UPDATE
     error_message("\nIndex error in AutoDiffMatrix::operator()");

  return a[((j)*n+i)  ];
}

adouble AutoDiffMatrix::operator() (int i) const
{
  if ( i<0 || i>n*m-1 )  // EIGEN_UPDATE
     error_message("\nIndex error in AutoDiffMatrix::operator()");

  return a[ i ];
}

AutoDiffMatrix::~AutoDiffMatrix()
{
   delete[] a;
}


void inverse_ad(const AutoDiffMatrix& minput, AutoDiffMatrix* minv)
{
  int i,p,k,n,m;

  adouble pivot;

  n = minput.rows();
  m = minput.cols();

  if (n!=m) { error_message("\n Matrix must be square in inverse_ad()"); }

  for( i=0;i<n;i++) {
    for(k=0;k<n;k++) {
         (*minv)(i,k) = minput(i,k);
    }
  }

  adouble *a = minv->GetPr();

  for (p=0;p<n;p++)
  {

    pivot = a[p*n+p];

    for ( i=0;i<n;i++)
    {

      if ( i!=p )
      {

        for ( k=0; k<n; k++)
        {

          if ( k!=p )
          {

            if (  pivot==0.0) error_message("singular matrix in inverse_ad()");
//            Temp->elem(i,k)=Temp->elem(i,k)-Temp->elem(i,p)*Temp->elem(p,k)/Temp->elem(p,p);

              a[k*n+i] +=  -a[p*n+i]*a[k*n+p]/pivot;

          }

        }

      }

    }

    for (i=0;i<n;i++)
    {

      if (i!=p)
      {

//        Temp->elem(i,p)=-Temp->elem(i,p)/Temp->elem(p,p);
          a[p*n+i] = -a[p*n+i]/pivot;

//        Temp->elem(p,i)=Temp->elem(p,i)/Temp->elem(p,p);
          a[i*n+p]=a[i*n+p]/pivot;


      }

    }


    if (a[p*n+p]==0.0) error_message("singular matrix");

//    Temp->elem(p,p)=1/Temp->elem(p,p);
    a[p*n+p]=1/pivot;


  }
}


void product_ad(const AutoDiffMatrix& A,const AutoDiffMatrix& B, AutoDiffMatrix* AB)
{

  // Matrix product optimised for vector columnwise storage

  int  i,j,k,na,nb,mb;


  int  rowA, colB;

  adouble sum;

  na=A.rows();

  mb=B.cols();

  nb=B.rows();

  const adouble* apr = A.GetConstPr();

  const adouble* bpr = B.GetConstPr();

  if (A.cols() != B.rows() )

  { error_message("Incoherent matrix dimensions in product_ad()"); }

  for (i=0; i< na; i++)
  {


    for (j=0; j<mb; j++)
    {

      colB = j*B.rows();

      sum=0.0;

      rowA = i;

      for (k=0; k<nb; k++)
      {

          sum += apr[ rowA ]*bpr[ colB + k ];
          rowA += A.rows();

      }

      (*AB)(i,j) = sum;

    }

  }


}


void subtract_ad(const AutoDiffMatrix& A,const AutoDiffMatrix& B, AutoDiffMatrix* AmB)
{

  int i;

  int n = A.rows();

  int m = B.cols();

  int nelem = n*m;

  const adouble *b = B.GetConstPr();

  const adouble *a = A.GetConstPr();


  if (m != B.cols() || n != B.rows())

  { error_message("Incoherent matrix dimensions in subtract_ad()"); }



  for (i=0; i< nelem; i++)
  {

     AmB->GetPr()[i] = a[i] - b[i];

  }

}

void sum_ad(const AutoDiffMatrix& A,const AutoDiffMatrix& B, AutoDiffMatrix* ApB)
{
  int i;

  int n = A.rows();

  int m = B.cols();

  int nelem = n*m;

  const adouble *b = B.GetConstPr();

  const adouble *a = A.GetConstPr();


  if (m != B.cols() || n != B.rows())

  { error_message("Incoherent matrix dimensions in add_ad()"); }



  for (i=0; i< nelem; i++)
  {

     ApB->GetPr()[i] = a[i] + b[i];

  }
}

void AutoDiffMatrix::Print(const char* text) const
{
// Prints an AutoDiffMatrix
// Eg. A.Print(" A = ");

  int i,j;

  fprintf(stderr,"\n%s\n",text);
  for (i=0;i<n;i++)  // EIGEN_UPDATE
  {
    for (j=0;j<m;j++)
    {
      fprintf(stderr," %f ", (this->elem(i,j).value()) );
      fprintf(stderr,"\t");
    }
    fprintf(stderr,"\n");
  }



}

void AutoDiffMatrix::resize(int nn, int mm)
{
    if (a!=NULL) delete[] a;

    n = nn;
    m = mm;

    a = new adouble[n*m];
}

void AutoDiffMatrix::setZero()
{
    for(int i=0;i< n*m; i++) {
       a[i] = 0.0;
    }
}

void AutoDiffMatrix::setOne()
{
    for(int i=0;i< n*m; i++) {
       a[i] = 1.0;
    }
}

double convert_to_original_time(double tbar,double t0,double tf)
{
    return (  (tf+t0)/2.0 + (tf-t0)*tbar/2.0 );
}

adouble convert_to_original_time_ad(double tbar,adouble& t0,adouble& tf)
{

    adouble retval;
    
    retval = (tf+t0)/2.0 + (tf-t0)*(tbar/2.0);

    return retval;
}


bool use_local_collocation(Alg & algorithm)
{

     if (algorithm.collocation_method == "trapezoidal" || algorithm.collocation_method=="Hermite-Simpson")
         return true;
     else return false;

}




bool use_global_collocation(Alg & algorithm)
{
     return (!use_local_collocation(algorithm));
}




bool need_midpoint_controls(Alg& algorithm, Workspace* workspace)
{
    bool retval;
    if ( workspace->differential_defects == "Hermite-Simpson" )
        retval = true;
    else
	retval = false;
    return (retval);
}



#define MTR(a,i,j,n)  (a)[ (j)*(n) + (i) ]

void mtrx_mul_trans(adouble* a,double* b,adouble* ab,int na,int ma,int nb,int mb)
{
  int i,j,k;
  adouble sum;

  /* multiplies matrix a[na,ma] times the transpose of matrix b[nb,mb] and stores
     the result in matrix a_times_b[na,nb]
     example of call
		     mtrx_mult(m1,m2,mr,2,2,2,1)
  */
  if (ma!=mb) {error_message("matrix dimension error in mtrx_mul()"); return;}
  for (i=0; i<na; i++)
  {
    for (j=0; j<mb; j++)
    {
      sum=0.0;
      for (k=0; k<nb; k++)
      {
        if (MTR(b,j,k,nb)!=0)
	sum=sum+MTR(a,i,k,na)*MTR(b,j,k,nb);
      }
      MTR(ab,i,j,na)=sum;
    }
  }

}


bool useAutomaticDifferentiation(Alg& algorithm)
{
   if ( algorithm.derivatives=="automatic" )
     return true;
   else {return false;}
}


void clip_vector_given_bounds(MatrixXd& xp, MatrixXd& xlb, MatrixXd& xub)
{
    long nvars = xp.rows();
    for(int i=0; i<nvars;i++) { // EIGEN_UPDATE
         if( xp(i) > xub(i) )
            xp(i)=xub(i);
         else if ( xp(i) < xlb(i) )
            xp(i)=xlb(i);
    }
}




void copy_decision_variables(Sol& solution, MatrixXd& x, Prob& problem, Alg& algorithm, Workspace* workspace)
{
   // Assign the contents of the decision vector x into the
   // relevant variables of the solution structure.

   int i;

   int iphase_offset = 0;

   for (i=0; i< problem.nphases; i++) {

	MatrixXd& control_scaling  = problem.phase[i].scale.controls;
	MatrixXd& state_scaling    = problem.phase[i].scale.states;
   MatrixXd& param_scaling    = problem.phase[i].scale.parameters;
   double   time_scaling      = problem.phase[i].scale.time;

	int k;

	int norder        = problem.phase[i].current_number_of_intervals;
	int ncontrols     = problem.phase[i].ncontrols;
	int nstates       = problem.phase[i].nstates;
   int nparam        = problem.phase[i].nparameters;
	int offset1       = ncontrols*(norder+1);
   int offset2       = (ncontrols+nstates)*(norder+1);

	int nvars_phase_i = get_nvars_phase_i(problem,i,workspace);

//	double    t0  = x(iphase_offset + nvars_phase_i-1)/time_scaling;
             double    t0  = x(iphase_offset + nvars_phase_i-2)/time_scaling;
//	double    tf  = x(iphase_offset + nvars_phase_i)/time_scaling;
             double    tf  = x(iphase_offset + nvars_phase_i-1)/time_scaling;
	for (k=0; k<norder+1; k++) {  // EIGEN_UPDATE
                if (ncontrols>0) {
//		    (solution.controls[i])(colon(),k) = elemDivision(x(colon(iphase_offset+(k-1)*ncontrols+1,iphase_offset+k*ncontrols) ) , control_scaling);
          (solution.controls[i]).col(k) = elemDivision(x.block(iphase_offset+(k)*ncontrols, 0, ncontrols,1) , control_scaling);
                }
//		(solution.states[i])(colon(),k)   = elemDivision(x(colon(iphase_offset+(k-1)*nstates+1+offset1, iphase_offset + k*nstates+offset1)), state_scaling);
      (solution.states[i]).col(k)   = elemDivision(x.block(iphase_offset+(k)*nstates+offset1, 0, nstates, 1), state_scaling);
  	   (solution.nodes[i])(0,k)          =  convert_to_original_time( (workspace->snodes[i])(k), t0, tf );
	}

//        solution.parameters[i] = elemDivision( x(colon(iphase_offset+offset2+1, iphase_offset+offset2+nparam)), param_scaling);
        solution.parameters[i] = elemDivision( x.block(iphase_offset+offset2,0,nparam,1), param_scaling);

        iphase_offset += nvars_phase_i;

  }

}



void sort_vector(MatrixXd& A, RowVectorXi& sindex)
{
	/* Sorts vector A in ascending order              */
   /* Returns the array of sorted indices  */

   double *a;
   long   n = MAX(A.rows(), A.cols());

   if ( A.rows()!=1 && A.cols() != 1 ) {

      error_message("Argument must be a column or row vector in sort_vector() ");

   }


   long i,j,inc;
   long vi;
   double v;
   inc=1;

   // a = A.a - 1;
   a = &A(0)-1;

   sindex.resize(n);
   
   int* pindx = &sindex(0); pindx--;

   if (1) { for( i=1; i<= n; i++ )  pindx[i] = (int) i; }


        do {
                inc *= 3;
                inc++;
        } while (inc <= n);
        do {
                inc /= 3;
                for (i=inc+1;i<=n;i++) {
                        v=a[i];
                        vi=(int) pindx[i];
                        j=i;
                        while (a[j-inc] > v) {
                                a[j]=a[j-inc];
                                pindx[j]=pindx[j-inc];
                                j -= inc;
                                if (j <= inc) break;
                        }
                        a[j]=v;
                        pindx[j]= (int) vi;

                }
        } while (inc > 1);


   for( i=0; i< n; i++ )  sindex(i) = sindex(i)-1; 

}

void sort(MatrixXd& A)
{
   double *a;
   RowVectorXi sindex;
   long   n = MAX(A.rows(), A.cols());

   if ( A.rows()!=1 && A.cols() != 1 ) {

      error_message("Argument must be a column or row vector in sort_vector() ");

   }


   long i,j,inc;
   int vi;
   double v;
   inc=1;

   // a = A.a - 1;
   a = &A(0)-1;

   sindex.resize(A.rows(),A.cols());
   
   int* pindx = &sindex(0); pindx--;

   if (1) { for( i=1; i<= n; i++ )  pindx[i] = (int) i; }


        do {
                inc *= 3;
                inc++;
        } while (inc <= n);
        do {
                inc /= 3;
                for (i=inc+1;i<=n;i++) {
                        v=a[i];
                        vi=(int) pindx[i];
                        j=i;
                        while (a[j-inc] > v) {
                                a[j]=a[j-inc];
                                pindx[j]=pindx[j-inc];
                                j -= inc;
                                if (j <= inc) break;
                        }
                        a[j]=v;
                        pindx[j]= (int) vi;

                }
        } while (inc > 1);


}

void rearrange_vector(MatrixXd& m, RowVectorXi& sindex)
{
	
	bool rowVector = false;
	 
	if ( m.rows()!=1 && m.cols() != 1 ) {

      error_message("Argument must be a column or row vector in sort() ");

   }
   
   if (m.rows()==1) rowVector=true;
	
   MatrixXd A = stack_columns(m);
   
   m=A;
	
	for (long i =0; i< sindex.size(); i++) {
	
	  m(i) = A(sindex(i));
	  
	}
	
	if (rowVector==true) m=m.transpose().eval();

}

void Save(const MatrixXd& m, const char* filename)
{
  long i,j;

  FILE *fp;

  if ( (fp = fopen(filename,"w")) == NULL )

  {  error_message( "Error opening file in Save()"); }


    for (i=0;i<m.rows();i++) { for (j=0;j<m.cols();j++) {

    fprintf(fp,"%g\t", m(i,j)); } fprintf(fp,"\n"); }

  fclose(fp);
  
}

MatrixXd linspace(double X1, double X2, long N)
{
     MatrixXd  m;
     VectorXd  v1;
     v1 = VectorXd::LinSpaced(N,X1,X2);
     m = v1;
     long lm = length(m);
     m.resize(1,lm);
     return m;
}

double sum(MatrixXd& A)
{
   return A.sum();
}


MatrixXd reshape(MatrixXd& M1, int nrow, int ncol)
{
    Map<MatrixXd> M2(M1.data(), nrow, ncol);
    
    return M2;
}

MatrixXd zeros(long nrows, long ncols)
{
    MatrixXd m(nrows, ncols);
    
    m.setZero();
    
    return m;
}

MatrixXd eye(long n)
{
	 MatrixXd m(n, n);
    
    m.setIdentity(n,n);
    
    return m;
}

MatrixXd ones(long nrows, long ncols)
{
    MatrixXd m(nrows, ncols);
    
    m.setOnes();
    
    return m;

}

MatrixXd elemProduct(const MatrixXd& m1,const MatrixXd& m2)
{
	 MatrixXd m;
	 
	 m = m1.cwiseProduct(m2);
	 
	 return m;

}

MatrixXd elemDivision(const MatrixXd& m1,const MatrixXd& m2)
{
	 MatrixXd m;
	 
	 m = m1.cwiseQuotient(m2);
	 
	 return m;
}

bool any(const MatrixXd& m)
{
   return m.any();
}

MatrixXd load_data(const char* filename, long nrows, long ncols)
{
  long i,j;

  FILE *fp;
  
  MatrixXd m;

  if ( (fp = fopen(filename,"r")) == NULL )

  {  error_message( "Error opening file in load_data()"); }


  for (i=0;i<nrows;i++) { 
      for (j=0;j<ncols;j++) {

           fscanf(fp,"%lf", &m(i,j)); 
      } 
  }

  fclose(fp);
  
  return m;
}

MatrixXd Abs(const MatrixXd& m)
{
   return m.array().abs();
}

MatrixXd tra(const MatrixXd& m)
{
	return m.transpose();
}

MatrixXd stack_columns(const MatrixXd& m)
{
     MatrixXd m1;
     
     m1.resize(m.rows()*m.cols(),1);
     
     long count=0;
     
     for (long i=0;i<m.rows();i++) { 
        for (long j=0;j<m.cols();j++) {
           m1(count,0) = m(i,j);
           count++; 
        } 
     }
     
     return m1;
}

long length(const MatrixXd& m)
{
    return (m.rows()*m.cols());
}

double Max(const MatrixXd& m)
{
    return m.array().maxCoeff();
}

double Max(const MatrixXd& m, long* i)
{
	  long irow, icol;
	 
     if ( m.rows()!=1 && m.cols() != 1 ) {

          error_message("Argument must be a column or row vector in sort() ");

     }	 
	 
    double maxval = m.maxCoeff(&irow,&icol);
    
    
    if (irow == 0) {
       *i = icol;
    }
    else if(icol==0) {
       *i = irow;    
    }
    
    return maxval;
}

double Min(const MatrixXd& m)
{
   return m.array().minCoeff();
}

double Min(const MatrixXd& m, long* i)
{
	 long irow, icol;
    double minval = m.minCoeff(&irow,&icol);
    
    *i = irow*icol-1;
    
    return minval;
}

double mean(const MatrixXd& m)
{
  return m.mean();
}

double MaxAbs(const MatrixXd& m)
{
    return m.lpNorm<Infinity>();
}



MatrixXd sum_columns(const MatrixXd& A)  // matrix sum function
{
// Returns a row vector containing the sum of values of each column

  long i, j;

  double sum;

  MatrixXd Temp;

  sum=0.0;

  Temp.resize( 1, A.cols() );


  for ( j = 0; j<A.cols(); j++ ) {


    for( i=0 ; i< A.rows(); i++ )
    {

        sum += A(i,j);

    }

    Temp(0,j) = sum;

    sum=0;

  }

  return Temp;

}




TripletSparseMatrix::TripletSparseMatrix(void)
{
// Default constructor
   asize = 0;
   a     = 0;
   n     = 0;
   m     = 0;
   nz    = 0;

   a = NULL;
   RowIndx = NULL;
   ColIndx = NULL;
}


TripletSparseMatrix::TripletSparseMatrix( int nn, int mm, int nzz)
{
    n = nn;
    m = mm;
    nz= nzz;
    asize = nz;
    a       = new double[nz];
    RowIndx = new int[nz];
    ColIndx = new int[nz];

    for(int k=0;k<nz;k++)  {
       RowIndx[k]=0;
       ColIndx[k]=0;
       a[k]      = 0.0;
    }

}

TripletSparseMatrix::TripletSparseMatrix(double* aArg, int nArg, int mArg, int nzArg, int* RowIndxArg, int* ColIndxArg)
// Constructor using arrays
{
   int i;
   n     = nArg;
   m     = mArg;
   nz    = nzArg;
   asize = nz;

   a       = new double[nz];
   RowIndx = new int[nz];
   ColIndx = new int[nz];


   for(i=0;i<nz;i++)
   {
         a[i]       = aArg[i];
         ColIndx[i] = ColIndxArg[i];
         RowIndx[i] = RowIndxArg[i];
   }


}



void TripletSparseMatrix::resize(int nnew, int mnew, int nznew)
{

    double *anew;
    int    *RowIndxNew;
    int    *ColIndxNew;
    int    dmax;

    anew = new double[nznew];
    RowIndxNew = new int[nznew];
    ColIndxNew = new int[nznew];

    dmax = MIN( nznew, nz );

    if (a!= NULL)
    {
        memcpy(anew, a, dmax*sizeof(double) );
        delete a;
    }

    if (RowIndx!=NULL) {
    	memcpy(RowIndxNew, RowIndx, dmax*sizeof(int) );
	delete RowIndx;
    }

    if (ColIndx!=NULL) {
    	memcpy(ColIndxNew, ColIndx, dmax*sizeof(int) );
        delete ColIndx;
    }

    a       = anew;
    RowIndx = RowIndxNew;
    ColIndx = ColIndxNew;

    n  = nnew;
    m  = mnew;
    nz = nznew;

}



TripletSparseMatrix::TripletSparseMatrix( const TripletSparseMatrix& A) // copy constructor
{
   asize = A.nz;
   a     = A.a;
   n     = A.n;
   m     = A.m;
   nz    = A.nz;

   a = new double[nz];
   RowIndx = new int[nz];
   ColIndx = new int[nz];

   memcpy(a, A.a, nz*sizeof(double) );
   memcpy(ColIndx, A.ColIndx, nz*sizeof(int) );
   memcpy(RowIndx, A.RowIndx, nz*sizeof(int) );
}

// Destructor
TripletSparseMatrix::~TripletSparseMatrix()
{
    if (a!=NULL)
 	delete a;
    if (RowIndx!= NULL)
        delete RowIndx;
    if (ColIndx!= NULL)
        delete ColIndx;

}




void TripletSparseMatrix::InsertNonZero(int i, int j, double val)
{
    int k;

    double* anew;
    int* RowNew;
    int* ColNew;
    int nznew=nz+1;
    int eflag = 0;

    anew =  new double[nznew];
    RowNew= new int[nznew];
    ColNew= new int[nznew];

    for (k=0; k< nz; k++)
    {
        if (RowIndx[k]==i && ColIndx[k]==j) {
               a[k]= val;
               eflag = 1;
               break;
        }
    }

    if (!eflag) {
    	memcpy(anew,   a,       nz*sizeof(double) );
    	memcpy(RowNew, RowIndx, nz*sizeof(int) );
    	memcpy(ColNew, ColIndx, nz*sizeof(int) );

    	delete a;
    	delete RowIndx;
    	delete ColIndx;

    	a = anew;
    	RowIndx = RowNew;
    	ColIndx = ColNew;
    	nz = nznew;
    	asize = nznew;

    	  RowIndx[nznew-1] = i;
        ColIndx[nznew-1] = j;
        a[nznew-1]       = val;
    }

    // Now resize the matrix if necessary
    for (k=0;k<nz;k++) {
        if (RowIndx[k]>n)
           n = RowIndx[k];
        if (ColIndx[k]>m)
           m = ColIndx[k];
    }

}


// Operators

TripletSparseMatrix& TripletSparseMatrix::operator += (const TripletSparseMatrix &rval)
{
    int i;

    if ( (this->n != rval.n) || (this->m != rval.m) )
        fprintf(stderr,"\nTripletSparseMatrix::operator += error: matrix dimensions do not agree");

    for (i=0; i<rval.nz; i++)
    {
        if ( (*this)(rval.RowIndx[i],rval.ColIndx[i])!=0.0 )
               this->a[i] += rval.a[i];
        else
               this->InsertNonZero( rval.RowIndx[i], rval.ColIndx[i], rval.a[i] );
    }
    this->Compress();
    return (*this);
}


TripletSparseMatrix& TripletSparseMatrix::operator+ (const TripletSparseMatrix& Other_matrix) const
{
   TripletSparseMatrix* Result;

   Result = new TripletSparseMatrix;

   (*Result)=(*this);

   return ((*Result)+=Other_matrix);

}

TripletSparseMatrix& TripletSparseMatrix::operator -= (const TripletSparseMatrix &rval)
{
    int i;

    if ( (this->n != rval.n) || (this->m != rval.m) )
        fprintf(stderr,"\nTripletSparseMatrix::operator -= error: matrix dimensions do not agree");

    for (i=0; i<rval.nz; i++)
    {
        if ( (*this)(rval.RowIndx[i],rval.ColIndx[i])!=0.0 )
               this->a[i] -= rval.a[i];
        else
               this->InsertNonZero( rval.RowIndx[i], rval.ColIndx[i], -rval.a[i] );
    }
    this->Compress();
    return (*this);
}


TripletSparseMatrix& TripletSparseMatrix::operator- (const TripletSparseMatrix& Other_matrix) const
{
   TripletSparseMatrix* Result;

   Result = new TripletSparseMatrix;

   (*Result)=(*this);

   return ((*Result)-=Other_matrix);
}




TripletSparseMatrix& TripletSparseMatrix::operator*= (double Arg)
{
   int i;
   for (i=0; i<nz; i++)
   {
      a[i]*= Arg;
   }
   return (*this);
}


TripletSparseMatrix& TripletSparseMatrix::operator* (double Arg) const
{
   TripletSparseMatrix* Result;

   Result = new TripletSparseMatrix;

   (*Result)=(*this);

   return ( (*Result)*=Arg );
}

TripletSparseMatrix& operator *(double Arg, const TripletSparseMatrix& A)
{
   return (A*Arg);
}



TripletSparseMatrix& elemProduct(const TripletSparseMatrix A, const TripletSparseMatrix& B)
{
          int i, j, k;

          double sum;

          double *a = A.GetPr();

          int *RowIndx = A.GetRowIndx_C_Array();
          int *ColIndx = A.GetColIndx_C_Array();

          TripletSparseMatrix *sp = new TripletSparseMatrix( A.rows(), B.cols(), 0 );

          if (A.rows() != B.rows() ||  A.cols() != B.cols())
                sp_error_message("Inconsistent matrix dimensions in TripletSparseMatrix elemProduct()");



          for (k = 0; k < A.GetNonZero(); k++) {
                i = RowIndx[k];
                j = ColIndx[k];
                sum = a[k]*B(i,j);
                if (sum!=0.0) sp->InsertNonZero(i, j, sum);
          }

          return (*sp);
}


TripletSparseMatrix& TripletSparseMatrix::operator/= (double Arg)
{
   int i;

   if (Arg==0.0) sp_error_message("Division by zero in TripletSparseMatrix::operator/");

   for (i=0; i<nz; i++)
   {
      a[i]/= Arg;
   }
   return (*this);
}

TripletSparseMatrix& TripletSparseMatrix::operator/ (double Arg) const
{
   TripletSparseMatrix* Result;

   if (Arg==0.0) sp_error_message("Division by zero in TripletSparseMatrix::operator/");

   Result = new TripletSparseMatrix;

   (*Result)=(*this);

   return ( (*Result)/=Arg );
}

TripletSparseMatrix& TripletSparseMatrix::operator= (const TripletSparseMatrix& Other_matrix)
{
     this->resize(Other_matrix.n, Other_matrix.m, Other_matrix.nz);
     memcpy(a, Other_matrix.a, nz*sizeof(double) );
     memcpy(RowIndx, Other_matrix.RowIndx, nz*sizeof(int) );
     memcpy(ColIndx, Other_matrix.ColIndx, nz*sizeof(int) );
     return (*this);
}

double& TripletSparseMatrix::operator() (int i,int j)
{
   int k;

    int ii;

    int eflag = 0;

    for (k=0; k< nz; k++)
    {
       if (RowIndx[k]==i && ColIndx[k]==j) {
              ii = k;
              eflag = 1;
              break;
       }

    }

    if( !eflag ) {

    	for (k=0; k<nz; k++) {

        	  if (RowIndx[k]==0 || ColIndx[k]==0) {
           		RowIndx[k]=i;
           		ColIndx[k]=j;
                        ii=k;
                        eflag=1;
                        break;
       		  }
        }

    }

    if (!eflag) {
       this->InsertNonZero(i,j,0.0);
       ii = nz-1;
    }

    return a[ii];

}

double TripletSparseMatrix::operator() (int row,int col) const
{
      int i;
      double retval = 0.0;

      if ( (row>n) || (row<0) || (col<0) || (col>m) ) {
          error_message("Out of range index in TripletSparseMatrix::operator()");
      }

      for (i=0; i< nz; i++) {
          if (RowIndx[i]==row && ColIndx[i]==col)
          {
               retval = a[i];
               break;
          }
      }

      return retval;
}



void TripletSparseMatrix::Print(const char* text)
{
   int i;

   fprintf(stderr,"\nSparse matrix %s",text);
   fprintf(stderr,"\nNumber of rows: %li", n);
   fprintf(stderr,"\nNumber of columns: %li", m);
   fprintf(stderr,"\nNumber of non-zero elements: %li", nz);
   if (n*m!=0) fprintf(stderr,"\nDensity: %lf%%", (  ((double) nz*100)/(double) (n*m)  ) );
   if (nz>0) fprintf(stderr,"\n(Row,Col)\tValue");

   for (i=0; i<nz; i++)
   {
       fprintf(stderr,"\n(%li,%li)\t\t%e", RowIndx[i], ColIndx[i], a[i]);

   }
   fprintf(stderr,"\n");

}

void TripletSparseMatrix::SaveSparsityPattern(const char* text) const
{
   int i, j;
   FILE *outfile;

   outfile = fopen(text,"w");

   for (i=0; i<n; i++)  {
       for(j=0; j<m; j++)  {
           if ( (*this)(i,j)!=0.0 )
               fprintf(outfile,"%c ",'*');
           else
               fprintf(outfile,"%c ",'0');
       }
       fprintf(outfile,"\n");
   }
   fclose(outfile);
}


void TripletSparseMatrix::Compress()
{

    int i, j;
    int deleted_count=0;

    double* anew;
    int* RowNew;
    int* ColNew;
    int nznew;
    int nztemp=nz;

    for (i=0; i< nztemp; i++)
    {
        if (fabs(a[i])<=(PSOPT_extras::GetEPS()) ) {
             for(j=i; j<nztemp-1; j++)
	     {
                    a[j] = a[j+1];
                    RowIndx[j]=RowIndx[j+1];
                    ColIndx[j]=ColIndx[j+1];
	     }
 	     deleted_count++;
             nztemp--;
             i--;
        }
    }

    nznew = nz-deleted_count;
    anew = new double[nznew];
    RowNew= new int[nznew];
    ColNew= new int[nznew];

    memcpy(anew,a, nznew*sizeof(double) );
    memcpy(RowNew, RowIndx, nznew*sizeof(int) );
    memcpy(ColNew, ColIndx, nznew*sizeof(int) );

    delete a;
    delete RowIndx;
    delete ColIndx;

    a = anew;
    RowIndx = RowNew;
    ColIndx = ColNew;

    nz = nznew;
    asize = nznew;

}



TripletSparseMatrix& speye(int n)
{
// Returns a sparse identity matrix
     TripletSparseMatrix *sp = new TripletSparseMatrix;
     int i;
     sp->resize(n, n, n);

     for (i=0;i<sp->GetNonZero();i++) {
         sp->GetRowIndx_C_Array()[i]=i;
         sp->GetColIndx_C_Array()[i]=i;
         sp->GetPr()[i]=1.0;
     }

     return (*sp);
}


void TripletSparseMatrix::Load(const char* fname)
{
  int nrow, ncol, nnz;

  int k;

  FILE *fp;

  if ( (fp = fopen(fname,"r")) == NULL )

  {  sp_error_message( "Error opening file in TripletSparseMatrix::Load()"); }


   fscanf(fp,"%li", &nrow);
   fscanf(fp,"%li", &ncol);
   fscanf(fp,"%li", &nnz);

   this->resize(nrow, ncol, nnz);

   for (k=0;k<nnz;k++) {
            fscanf(fp,"%li",  &RowIndx[k] );
            fscanf(fp,"%li",  &ColIndx[k] );
            fscanf(fp,"%lf", &a[k]);
   }


  fclose(fp);

}

void TripletSparseMatrix::Save(const char* fname) const
{

  int k;

  FILE *fp;

  if ( (fp = fopen(fname,"w")) == NULL )

  {  sp_error_message( "Error opening file in TripletSparseMatrix::Save()"); }


   fprintf(fp,"%li\t%li\t%li\n", n, m, nz);

   for (k=0;k<nz;k++) {
            fprintf(fp,"%li\t%li\t%e\n",  RowIndx[k], ColIndx[k], a[k] );
   }

   fclose(fp);

}



TripletSparseMatrix& tra(const TripletSparseMatrix& A)
{
// Returns a Sparse matrix object with the transposition of input sparse matrix A.

   TripletSparseMatrix *sp = new TripletSparseMatrix(A.cols(), A.rows(), A.GetNonZero() );

   memcpy(sp->RowIndx, A.ColIndx, A.nz*sizeof(int) );
   memcpy(sp->ColIndx, A.RowIndx, A.nz*sizeof(int) );
   memcpy(sp->a      , A.a,       A.nz*sizeof(double) );

   return (*sp);

}





MatrixXd TripletSparseMatrix::col(int j) const
{
// Returns a dense column vector (MatrixXd object) with the j-th column of the calling
// sparse matrix.
  MatrixXd Temp;

  Temp.resize( this->rows(), 1 );
  Temp.setZero();
  for(int i=0;i<nz;i++) {
       if (ColIndx[i]==j) {
       		Temp( RowIndx[i], 0 ) = a[i];
       }
  }
  return Temp;

}

MatrixXd TripletSparseMatrix::row(int i) const
{
// Returns a dense row vector (MatrixXd object) with the i-th row of the calling
// sparse matrix.
  MatrixXd Temp;

  Temp.resize( 1, this->cols());
  Temp.setZero();
  for(int j=0;j<nz;j++) {
       if (RowIndx[j]==i) {
       		Temp(0, ColIndx[j] ) = a[i];
       }
  }
  return Temp;
}




void TripletSparseMatrix::Transpose()
{
// Transposes a sparse matrix object
   int s;
   int *temp = new int[nz];


   memcpy(temp, RowIndx,    nz*sizeof(int) );
   memcpy(RowIndx, ColIndx, nz*sizeof(int) );
   memcpy(ColIndx, temp,    nz*sizeof(int) );

   s = n;
   n = m;
   m = s;

   return;

}


void sp_error_message(const char *error_text)
{
     error_message( error_text );
}


bool isEmpty(const MatrixXd& m)     {  
	return !((bool) m.rows()&&m.cols()); 
}


bool isSymmetric(const MatrixXd& m) {  
	return m.isApprox(m.transpose()); 
}


DEC_THREAD int     PSOPT_extras::errorFlag   = 0;	   

DEC_THREAD int     PSOPT_extras::print_level = 1;        

DEC_THREAD clock_t PSOPT_extras::start_clock = 0;



void PSOPT_extras::SetPrintLevel( int plevel )
{
	 print_level=plevel;
	 return;
}


int PSOPT_extras::PrintLevel() {
	return print_level;
}

#ifndef OUTPUT_STREAM
#define OUTPUT_STREAM stderr
#endif

void error_message(const char *error_text)
{
// Error handling routine

  string m1;
  string m2;

  m1 =  "\n**** Run time error";
  m1 += "\n**** To trace this error, set up your debugger to break at";
  m1 += "\n**** function 'error_message', then do a backtrace.";
  m1 += "\n**** A diagnostic message is given below:";
  m2 = error_text;
  m2 =  "\n**** ====> " + m2 + " <====\n\n";

  if ( PSOPT_extras::PrintLevel() ) {

  fprintf(OUTPUT_STREAM,"%s", m1.c_str());
  fprintf(OUTPUT_STREAM,"%s", m2.c_str());
  FILE* err_file = fopen("error_message.txt","w");
  fprintf(err_file,"%s", m1.c_str() );
  fprintf(err_file,"%s", m2.c_str() );
  fclose(err_file);


 }

  PSOPT_extras::RiseErrorFlag();

  throw ErrorHandler(m1+m2);


}

ErrorHandler::ErrorHandler(const string m)
{
  error_message = m;
}


void PSOPT_extras::tic(void)
{
// Chronometer starts
   PSOPT_extras::SetStartTicks(clock());

}

double PSOPT_extras::toc()
{
// Chronometer stops

   clock_t  stop_ticks     = clock();
   clock_t  elapsed_ticks;
   double elapsed_time;

   elapsed_ticks= stop_ticks-PSOPT_extras::GetStartTicks();

   elapsed_time = ((double) elapsed_ticks/CLOCKS_PER_SEC);

   if (1) {

       if (PSOPT_extras::PrintLevel() ) {
       	fprintf(stderr,"\n Elapsed time is %e seconds\n", elapsed_time );
       }

   }

   return (elapsed_time);

}

