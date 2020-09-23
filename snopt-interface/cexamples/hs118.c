#include <stdlib.h>
#include <stdio.h>
#include "sqopt_cwrap.h"

void hs118Hx
(int *nnH, double x[], double Hx[], int *nState,
  char   cu[], int   *lencu,
  int    iu[], int   *leniu,
  double ru[], int   *lenru) {

  int i;
  for (i = 0; i < 5; i++) {
    Hx[3*i]   = .0002*x[3*i];
    Hx[3*i+1] = .0002*x[3*i+1];
    Hx[3*i+2] = .0003*x[3*i+2];
  }
}

/*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/

int main(int argc , char* argv[])
{
  sqProblem hs118prob;

  int    i, info;
  int    Cold   = 0;
  int    n      = 15;
  int    m      = 18;
  int    nnH    = 15;
  int    ncObj  =  0;
  int    lenA   = 54, neA;
  int    iObj   = 17;
  double ObjAdd =  0;

  int    indA[54], locA[16];
  double valA[54];

  int    eType[33], hs[33];
  double x[33], bl[33], bu[33], pi[18], rc[33], cObj[1];

  int    nS, nInf;
  double objective, sInf;

  double infBnd = 1.0e20;

  // sqInit must be called first.
  sqInit(&hs118prob, "hs118", "hs118.out", 1);


  // User workspace allocated and set.
  //   May be accesed in the user-defined function hs118hx.
  hs118prob.leniu = 2;
  hs118prob.iu = (int*)malloc(sizeof(int) * hs118prob.leniu);
  hs118prob.iu[0] = 0;
  hs118prob.iu[1] = 1;


  // Set bounds
  for (i = n; i < n+m; i++) {
    bl[i] = 0;
    bu[i] = infBnd;
  }
  bl[n+iObj] = -infBnd;

  bl[n]  =  -7;
  bu[n]  =   6;

  bl[n+1]  =  -7;
  bu[n+1]  =   6;

  bl[n+2]  =  -7;
  bu[n+2]  =   6;

  bl[n+3]  =  -7;
  bu[n+3]  =   6;

  bl[n+4]  =  -7;
  bu[n+4]  =   7;

  bl[n+5]  =  -7;
  bu[n+5]  =   7;

  bl[n+6]  =  -7;
  bu[n+6]  =   7;

  bl[n+7]  =  -7;
  bu[n+7]  =   7;

  bl[n+8]  =  -7;
  bu[n+8]  =   6;

  bl[n+9] =  -7;
  bu[n+9] =   6;

  bl[n+10] =  -7;
  bu[n+10] =   6;

  bl[n+11] =  -7;
  bu[n+11] =   6;

  bl[n+12] =  60;
  bl[n+13] =  50;
  bl[n+14] =  70;
  bl[n+15] =  85;
  bl[n+16] = 100;

  for (i = 0; i < n; i++) {
    bl[i] = 0;
    bu[i] = infBnd;
  }

  bl[ 0] =   8;
  bu[ 0] =  21;

  bl[ 1] =  43;
  bu[ 1] =  57;

  bl[ 2] =   3;
  bu[ 2] =  16;

  bu[ 3] =  90;
  bu[ 4] = 120;
  bu[ 5] =  60;
  bu[ 6] =  90;
  bu[ 7] = 120;
  bu[ 8] =  60;
  bu[ 9] =  90;
  bu[10] = 120;
  bu[11] =  60;
  bu[12] =  90;
  bu[13] = 120;
  bu[14] =  60;

  // Initialize states, x and multipliers
  for (i = 0; i < n+m; i++) {
    eType[i] = 0;
    hs[i]    = 0;
    x[i]     = 0;
    rc[i]    = 0;
  }

  for (i = 0; i < m; i++) {
    pi[i]   = 0;
  }

  x[ 0]  =  20;
  x[ 1]  =  55;
  x[ 2]  =  15;
  x[ 3]  =  20;
  x[ 4]  =  60;
  x[ 5]  =  20;
  x[ 6]  =  20;
  x[ 7]  =  60;
  x[ 8]  =  20;
  x[ 9]  =  20;
  x[10]  =  60;
  x[11]  =  20;
  x[12]  =  20;
  x[13]  =  60;
  x[14]  =  20;


  // Set up the Jacobian matrix
  locA[ 0] =  0;

  indA[ 0] =  0;
  indA[ 1] = 12;
  indA[ 2] = 17;

  valA[ 0] = -1.0;
  valA[ 1] =  1.0;
  valA[ 2] =  2.3;

  // Column 2.
  locA[ 1] =  3;

  indA[ 3] =  4;
  indA[ 4] = 12;
  indA[ 5] = 17;

  valA[ 3] = -1.0;
  valA[ 4] =  1.0;
  valA[ 5] =  1.7;

  // Column 3.
  locA[ 2] =  6;

  indA[ 6] =  8;
  indA[ 7] = 12;
  indA[ 8] = 17;

  valA[ 6] = -1.0;
  valA[ 7] =  1.0;
  valA[ 8] =  2.2;

  // Column 4.
  locA[ 3] = 9;

  indA[ 9] =  0;
  indA[10] =  1;
  indA[11] = 13;
  indA[12] = 17;

  valA[ 9] =  1.0;
  valA[10] = -1.0;
  valA[11] =  1.0;
  valA[12] =  2.3;

  // Column 5.
  locA[ 4] = 13;

  indA[13] =  4;
  indA[14] =  5;
  indA[15] = 13;
  indA[16] = 17;

  valA[13]  =  1.0;
  valA[14]  = -1.0;
  valA[15]  =  1.0;
  valA[16]  =  1.7;

  // Column 6.
  locA[5]  = 17;

  indA[17] =  8;
  indA[18] =  9;
  indA[19] = 13;
  indA[20] = 17;

  valA[17] =  1.0;
  valA[18] = -1.0;
  valA[19] =  1.0;
  valA[20] =  2.2;

  // Column 7.
  locA[6]  = 21;

  indA[21] =  1;
  indA[22] =  2;
  indA[23] = 14;
  indA[24] = 17;

  valA[21] =  1.0;
  valA[22] = -1.0;
  valA[23] =  1.0;
  valA[24] =  2.3;

  // Column 8.
  locA[7]  = 25;

  indA[25] =  5;
  indA[26] =  6;
  indA[27] = 14;
  indA[28] = 17;

  valA[25] =  1.0;
  valA[26] = -1.0;
  valA[27] =  1.0;
  valA[28] =  1.7;

  // Column 9.
  locA[8]  = 29;

  indA[29] =  9;
  indA[30] = 10;
  indA[31] = 14;
  indA[32] = 17;

  valA[29] =  1.0;
  valA[30] = -1.0;
  valA[31] =  1.0;
  valA[32] =  2.2;

  // Column 10.
  locA[9] = 33;

  indA[33] =  2;
  indA[34] =  3;
  indA[35] = 15;
  indA[36] = 17;

  valA[33] =  1.0;
  valA[34] = -1.0;
  valA[35] =  1.0;
  valA[36] =  2.3;

  // Column 11.
  locA[10] = 37;

  indA[37] =  6;
  indA[38] =  7;
  indA[39] = 15;
  indA[40] = 17;

  valA[37] =  1.0;
  valA[38] = -1.0;
  valA[39] =  1.0;
  valA[40] =  1.7;

  // Column 12.
  locA[11] = 41;

  indA[41] = 10;
  indA[42] = 11;
  indA[43] = 15;
  indA[44] = 17;

  valA[41] =  1.0;
  valA[42] = -1.0;
  valA[43] =  1.0;
  valA[44] =  2.2;

  // Column 13.
  locA[12] = 45;

  indA[45] =  3;
  indA[46] = 16;
  indA[47] = 17;

  valA[45] =  1.0;
  valA[46] =  1.0;
  valA[47] =  2.3;

  // Column 14.
  locA[13] = 48;

  indA[48] =  7;
  indA[49] = 16;
  indA[50] = 17;

  valA[48] =  1.0;
  valA[49] =  1.0;
  valA[50] =  1.7;

  // Column 15.
  locA[14] = 51;

  indA[51] = 11;
  indA[52] = 16;
  indA[53] = 17;

  valA[51] =  1.0;
  valA[52] =  1.0;
  valA[53] =  2.2;

  // locA[n]-1 points to the last nonzero of the nth column.
  locA[15] = 54;

  neA = 54;

  // Read options, set options.
  info = sqSpecsfile(&hs118prob, "hs118.spc");

  info = solveQ(&hs118prob, Cold, hs118Hx,
		 m, n, neA, ncObj, nnH, iObj, ObjAdd,
		 valA, indA, locA,
		 bl, bu, cObj, eType, hs, x, pi, rc,
		 &objective, &nS, &nInf, &sInf);

  // Deallocate space.
  free(hs118prob.iu);
  deleteSQOPT(&hs118prob);

  return 0;
}
