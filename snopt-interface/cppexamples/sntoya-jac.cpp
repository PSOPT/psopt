#include <stdio.h>
#include <string.h>
#include <iostream>

#include "sntoyA.hpp"
#include "snoptProblem.hpp"

using namespace std;

void toyusrf_(int    *Status, int *n,    double x[],
	      int    *needF,  int *neF,  double F[],
	      int    *needG,  int *neG,  double G[],
	      char      *cu,  int *lencu,
	      int    iu[],    int *leniu,
	      double ru[],    int *lenru )
{
  //==================================================================
  // Computes the nonlinear objective and constraint terms for the toy
  // problem featured in the SnoptA users guide.
  // neF = 3, n = 2.
  //
  //   Minimize     x(2)
  //
  //   subject to   x(1)**2      + 4 x(2)**2  <= 4,
  //               (x(1) - 2)**2 +   x(2)**2  <= 5,
  //                x(1) >= 0.
  //
  //==================================================================
  F[0] =  x[1];
  F[1] =  x[0]*x[0] + 4*x[1]*x[1];
  F[2] = (x[0] - 2)*(x[0] - 2) + x[1]*x[1];
}

void toyusrfg_( int    *Status, int *n,    double x[],
		int    *needF,  int *neF,  double F[],
		int    *needG,  int *neG,  double G[],
		char      *cu,  int *lencu,
		int    iu[],    int *leniu,
		double ru[],    int *lenru )
{
  //==================================================================
  // Computes the nonlinear objective and constraint terms for the toy
  // problem featured in the SnoptA users guide.
  // neF = 3, n = 2.
  //
  //   Minimize     x(2)
  //
  //   subject to   x(1)**2      + 4 x(2)**2  <= 4,
  //               (x(1) - 2)**2 +   x(2)**2  <= 5,
  //                x(1) >= 0.
  //
  // The triples (g(k),iGfun(k),jGvar(k)), k = 1:neG, define
  // the sparsity pattern and values of the nonlinear elements
  // of the Jacobian.
  //==================================================================

  if ( *needF > 0 ) {
    F[0] =  x[1]; //  Objective row
    F[1] =  x[0]*x[0] + 4*x[1]*x[1];
    F[2] = (x[0] - 2)*(x[0] - 2) + x[1]*x[1];
  }


  if ( *needG > 0 ) {
    // iGfun[0] = 1
    // jGvar[0] = 0
    G[0] = 2*x[0];

    // iGfun[1] = 1
    // jGvar[1] = 1
    G[1] = 8*x[1];

    // iGfun[2] = 2
    // jGvar[2] = 0
    G[2] = 2*(x[0] - 2);

    // iGfun[3] = 2
    // jGvar[3] = 1
    G[3] = 2*x[1];
  }
}


int main( int argc, char **argv)
{
  snoptProblemA ToyProb;

  // Allocate and initialize;
  int n     =  2;
  int neF   =  3;

  int neA, neG;

  double *x      = new double[n];
  double *xlow   = new double[n];
  double *xupp   = new double[n];
  double *xmul   = new double[n];
  int    *xstate = new    int[n];

  double *F      = new double[neF];
  double *Flow   = new double[neF];
  double *Fupp   = new double[neF];
  double *Fmul   = new double[neF];
  int    *Fstate = new int[neF];

  int    ObjRow  = 0;
  double ObjAdd  = 0;

  int Cold = 0, Basis = 1, Warm = 2;


  // Set the upper and lower bounds.
  xlow[0]   =  0.0;  xlow[1]   = -1e20;
  xupp[0]   = 1e20;  xupp[1]   =  1e20;
  xstate[0] =    0;  xstate[1] =  0;

  Flow[0] = -1e20; Flow[1] = -1e20; Flow[2] = -1e20;
  Fupp[0] =  1e20; Fupp[1] =   4.0; Fupp[2] =  5.0;
  x[0]    = 1.0;
  x[1]    = 1.0;


  // Load the data for ToyProb ...
  ToyProb.setProbName   ("Toy0");
  ToyProb.setPrintFile  ( "Toy0.out" );

  ToyProb.setProblemSize( n, neF );
  ToyProb.setObjective  ( ObjRow, ObjAdd );
  ToyProb.setX          ( x, xlow, xupp, xmul, xstate );
  ToyProb.setF          ( F, Flow, Fupp, Fmul, Fstate );

  ToyProb.setUserFun    ( toyusrf_ );

  ToyProb.setIntParameter( "Derivative option", 0 );
  ToyProb.setIntParameter( "Verify level ", 3 );

  // Make an explcit call to compute the Jacobian.
  // neA and neG are returned.
  ToyProb.computeJac(neA,neG);
  printf("neA: %d; neG: %d", neA, neG);

  ToyProb.solve( Cold );


  for (int i = 0; i < n; i++ ){
    cout << "x = " << x[i] << " xstate = " << xstate[i] << endl;
  }
  for (int i = 0; i < neF; i++ ){
    cout << "F = " << F[i] << " Fstate = " << Fstate[i] << endl;
  }


  delete []x;      delete []xlow;   delete []xupp;
  delete []xmul;   delete []xstate;

  delete []F;      delete []Flow;   delete []Fupp;
  delete []Fmul;   delete []Fstate;

}
