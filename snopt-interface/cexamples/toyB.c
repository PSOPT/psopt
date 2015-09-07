#include <stdlib.h>
#include <stdio.h>
#include "snopt_cwrap.h"

/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
/*  toycon  toyobj  main                                                       */
/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */

void toycon
( int *mode, int *nnCon, int *nnJac, int  *negCon,  double x[],
  double  fCon[],  double gCon[],  int *nState,
  char   cu[], int   *lencu,
  int    iu[], int   *leniu,
  double ru[], int   *lenru )
{
  if ( *mode == 0 || *mode == 2 ) {
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

void toyobj
( int *mode, int *nnObj,  double x[],
  double *fObj,  double gObj[], int *nState,
  char   cu[], int   *lencu,
  int    iu[], int   *leniu,
  double ru[], int   *lenru )
{
  if ( *mode == 0 || *mode == 2 ) {
    *fObj = 0;
  }

  if ( *mode == 0 || *mode == 2 ) {
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
  double objective;

  // snInit must be called first.
  //   9, 6 are print and summary unit numbers (for Fortran).
  //   6 == standard out
  snInit ( &toy, "ToyB", "ToyB.out", 9, 6 );

  // Set the problem size and other data.
  // This will allocate arrays inside snProblem struct.
  setProblemSize ( &toy, m, n, ne, nnCon, nnJac, nnObj );
  setObjective   ( &toy, iObj, ObjAdd );
  setFuncon      ( &toy, (snConB) toycon );
  setFunobj      ( &toy, (snObjB) toyobj );

  // User workspace allocated and set.
  //   May be accesed in the user-defined function toyusrfun.
  toy.leniu = 2;
  toy.iu = (int*)malloc( sizeof(int) * toy.leniu );
  toy.iu[0] = 0;
  toy.iu[1] = 1;


  // Set bounds
  toy.bl[0] =     0;  toy.bu[0] = 1e20;
  toy.bl[1] = -1e20;  toy.bu[1] = 1e20;
  toy.bl[2] = -1e20;  toy.bu[2] =    4;
  toy.bl[3] = -1e20;  toy.bu[3] =    5;
  toy.bl[4] = -1e20;  toy.bu[4] = 1e20;


  // Initialize states, x and multipliers
  for ( i = 0; i < n+m; i++ ) {
    toy.hs[i] = 0;
    toy.x[i]  = 0;
    toy.rc[i] = 0;
  }

  for ( i = 0; i < m; i++ ) {
    toy.pi[i] = 0;
  }

  toy.x[0] = 1.0;
  toy.x[1] = 1.0;

  // Set up the Jacobian matrix
  // Column 1
  toy.locJ[0] = 0;

  toy.indJ[0] = 0;
  toy.valJ[0] = 0;

  toy.indJ[1] = 1;
  toy.valJ[1] = 0;

  // Column 2
  toy.locJ[1] = 2;

  toy.indJ[2] = 0;
  toy.valJ[2] = 0;

  toy.indJ[3] = 1;
  toy.valJ[3] = 0;

  toy.indJ[4] = 2;
  toy.valJ[4] = 1;

  toy.locJ[2] = 5;

  // Read options, set options.
  info = setSpecsfile   ( &toy, "sntoy.spc" );
  setIntParameter( &toy, "Verify level", 3 );
  setIntParameter( &toy, "Derivative option", 3 );

  info = solveB( &toy, Cold, &objective );

  // Deallocate space.
  free ( toy.iu );
  deleteSNOPT ( &toy );

  return 0;
}
