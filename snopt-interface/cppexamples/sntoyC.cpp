#include <stdio.h>
#include <string.h>
#include <iostream>

#include "sntoyC.hpp"
#include "snoptProblem.hpp"

using namespace std;

void toyusrfun ( int *mode,  int *nnObj, int *nnCon,
		 int *nnJac, int *nnL,   int *negCon, double x[],
		 double *fObj,  double gObj[],
		 double fCon[], double gCon[], int *Status,
		 char    *cu, int *lencu,
		 int    iu[], int *leniu,
		 double ru[], int *lenru )
{
  //==================================================================
  // Computes the nonlinear objective and constraint terms for the toy
  // problem featured in the SnoptA users guide.
  // m = 3, n = 2.
  //
  //   Minimize     x(2)
  //
  //   subject to   x(1)**2      + 4 x(2)**2  <= 4,
  //               (x(1) - 2)**2 +   x(2)**2  <= 5,
  //                x(1) >= 0.
  //
  //==================================================================

  if ( *mode == 0 || *mode == 2 ) {
    *fObj   =  0;
    fCon[0] =  x[0]*x[0] + 4*x[1]*x[1];
    fCon[1] = (x[0] - 2)*(x[0] - 2) + x[1]*x[1];
  }

  if ( *mode == 1 || *mode == 2 ) {
    // gObj not set; no nonlinear variables in objective
    gCon[0] = 2*x[0];
    gCon[1] = 2*(x[0] - 2);
    gCon[2] = 8*x[1];
    gCon[3] = 2*x[1];
  }

}

int main( int argc, char **argv)
{
  snoptProblemC ToyProb("ToyC");

  int n     =  2;
  int m     =  3;
  int ne    =  5;
  int nnCon =  2;
  int nnObj =  0;
  int nnJac =  2;

  int    *indJ = new int[ne];
  int    *locJ = new int[n+1];
  double *valJ = new double[ne];

  double *x  = new double[n+m];
  double *bl = new double[n+m];
  double *bu = new double[n+m];
  double *pi = new double[m];
  double *rc = new double[n+m];
  int    *hs = new    int[n+m];

  int    iObj    = 2;
  double ObjAdd  = 0;

  int Cold = 0, Basis = 1, Warm = 2;


  // Set the upper and lower bounds.
  bl[0] =     0;  bu[0] = 1e20;
  bl[1] = -1e20;  bu[1] = 1e20;
  bl[2] = -1e20;  bu[2] =    4;
  bl[3] = -1e20;  bu[3] =    5;
  bl[4] = -1e20;  bu[4] = 1e20;

  // Initialize states, x and multipliers
  for ( int i = 0; i < n+m; i++ ) {
    hs[i] = 0;
     x[i] = 0;
    rc[i] = 0;
  }

  for ( int i = 0; i < m; i++ ) {
    pi[i] = 0;
  }

  x[0]    = 1.0;
  x[1]    = 1.0;

  // Set up the Jacobian matrix
  // Column 1
  locJ[0] = 0;

  indJ[0] = 0;
  valJ[0] = 0;

  indJ[1] = 1;
  valJ[1] = 0;

  // Column 2
  locJ[1] = 2;

  indJ[2] = 0;
  valJ[2] = 0;

  indJ[3] = 1;
  valJ[3] = 0;

  indJ[4] = 2;
  valJ[4] = 1;

  locJ[2] = 5;

  ToyProb.setProblemSize ( m, n, nnCon, nnJac, nnObj );
  ToyProb.setObjective   ( iObj, ObjAdd );
  ToyProb.setJ           ( ne, valJ, indJ, locJ );
  ToyProb.setX           ( bl, bu, x, pi, rc, hs );
  ToyProb.setUserFun     ( toyusrfun );

  ToyProb.setSpecsFile   ( "sntoy.spc" );
  ToyProb.setIntParameter( "Verify level", 3 );
  ToyProb.setIntParameter( "Derivative option", 3 );

  ToyProb.solve          ( Cold );

  delete []indJ;  delete []locJ; delete []valJ;

  delete []x;     delete []bl;   delete []bu;
  delete []pi;    delete []rc;   delete []hs;

}
