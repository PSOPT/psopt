#include <stdlib.h>
#include <stdio.h>
#include "snopt_cwrap.h"

/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
/*  toyusrfun    main                                                          */
/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */

void toyusrfun
( int *mode, int *nnObj, int *nnCon, int *nnJac, int *nnL,
  int  *negCon,  double x[],
  double *fObj,  double gObj[],
  double  fCon[],  double gCon[],  int *nState,
  char   cu[], int   *lencu,
  int    iu[], int   *leniu,
  double ru[], int   *lenru )
{
  if ( *mode == 0 || *mode == 2 ) {
    *fObj = 0;
    fCon[0] = x[0] * x[0] + 4*x[1]*x[1];
    fCon[1] = (x[0] - 2)*(x[0] - 2) + x[1] * x[1];
  }

  if ( *mode == 0 || *mode == 2 ) {
    gCon[0] = 2.0 * x[0];
    gCon[1] = 2.0 * (x[0] - 2);

    gCon[2] = 8.0 * x[1];
    gCon[3] = 2.0 * x[1];
  }
}

/*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/

int main( int argc , char* argv[] )
{
  snProblem toy;

  int    i, info;
  int    Cold   =  0;
  int    n      =  2;
  int    m      =  3;
  int    ne     =  5;
  int    nnCon  =  2;
  int    nnObj  =  0;
  int    nnJac  =  2;
  int    iObj   =  2;
  double ObjAdd =  0;

  int    leniw  = 2000;
  int    lenrw  = 2000;
  int    iw[2000];
  double rw[2000];

  int    hs[5], locJ[3], indJ[5];
  double bl[5], bu[5], x[5], pi[3], rc[5], valJ[5];

  int    nS, nInf;
  double objective, sInf;


  // snInit must be called first.
  //   last entry -> summary on if > 0, else summary off
  snInitW( &toy, "ToyC", "ToyC.out", 1, iw, leniw, rw, lenrw );

  // User workspace allocated and set.
  //   May be accesed in the user-defined function toyusrfun.
  toy.leniu = 2;
  toy.iu = (int*)malloc( sizeof(int) * toy.leniu );
  toy.iu[0] = 0;
  toy.iu[1] = 1;


  // Set bounds
  bl[0] =     0;  bu[0] = 1e20;
  bl[1] = -1e20;  bu[1] = 1e20;
  bl[2] = -1e20;  bu[2] =    4;
  bl[3] = -1e20;  bu[3] =    5;
  bl[4] = -1e20;  bu[4] = 1e20;


  // Initialize states, x and multipliers
  for ( i = 0; i < n+m; i++ ) {
    hs[i] = 0;
    x[i]  = 0;
    rc[i] = 0;
  }

  for ( i = 0; i < m; i++ ) {
    pi[i] = 0;
  }

  x[0] = 1.0;
  x[1] = 1.0;

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

  // Read options, set options.
  info = setSpecsfile( &toy, "sntoy.spc" );
  setIntParameter( &toy, "Verify level", 3 );
  setIntParameter( &toy, "Derivative option", 3 );

  nS   = 0;

  info = solveC( &toy, Cold, m, n, ne, nnCon, nnObj, nnJac,
		 iObj, ObjAdd,
		 toyusrfun,
		 valJ, indJ, locJ,
		 bl, bu, hs, x, pi, rc,
		 &objective, &nS, &nInf, &sInf );

  printf ("SNOPT Exit, INFO = %d\n", info );
  printf ("Final nonlinear objective = %16.8f\n", objective );

  // Deallocate space.
  free( toy.iu );
  deleteSNOPT( &toy );

  return 0;
}
