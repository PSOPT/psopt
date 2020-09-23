#include <stdio.h>
#include <string.h>
#include <iostream>
#include <math.h>

#include "snoptProblem.hpp"

using namespace std;

void usrFG(int    *Status, int *n,    double x[],
	   int    *needF,  int *lenF,  double F[],
	   int    *needG,  int *lenG,  double G[],
	   char      *cu,  int *lencu,
	   int    iu[],    int *leniu,
	   double ru[],    int *lenru) {
  int neG = 0;
  int jx1, jx2, ju, ode1, ode2, ObjRow;
  double alpha = 0.0, tf = 1.0;
  double  Fobj = 0.0, Gobj0 = 0.0;
  double Gobj1;

  int    nH    = *n / 3 - 1;
  double h     = tf / nH;

  if (*Status == 2) {
    std::cout << "Last call to usrFG iu[5]: " << iu[5] << std::endl;
  }

  jx1    = 0;
  jx2    = jx1 + nH + 1;
  ju     = jx2 + nH + 1;

  ode1   = 0;
  ode2   = ode1 + nH;
  ObjRow = ode2 + nH;

  for (int i = 0; i < nH; i++) {
    if (*needF > 0) {
      Fobj += alpha*h*pow(x[ju+i+1] - x[ju+i], 2);

      // subject to ode1 {i in 0..(nh-1]}
      F[ode1+i] = x[jx1+i+1] - x[jx1+i]
	- 0.5*h*( x[ju+i]  *(10.0*x[jx2+i]   - x[jx1+i])
		    + x[ju+i+1]*(10.0*x[jx2+i+1] - x[jx1+i+1]));
      // subject to ode2 {i in 0..[nh-1]}
      F[ode2+i] = x[jx2+i+1] - x[jx2+i]
	- 0.5*h*(x[ju+i]  *(x[jx1+i]   - 10.0*x[jx2+i])
		   - (1.0-x[ju+i])  *x[jx2+i]
		   +x[ju+i+1]*(x[jx1+i+1] - 10.0*x[jx2+i+1])
		   - (1.0-x[ju+i+1])*x[jx2+i+1]);
    }

    if (*needG > 0) {
      // Objective gradient.
      Gobj1  =  2.0*alpha*h*(x[ju+i+1] - x[ju+i]);
      G[neG] =  Gobj0 - Gobj1;
      neG    = neG + 1;
      Gobj0  =  Gobj1;

      // First ode constraint.
      G[neG] =  - 1.0 + 0.5*h*x[ju+i];
      neG    = neG + 1;

      G[neG] =    1.0 + 0.5*h*x[ju+i+1];
      neG    = neG + 1;

      G[neG] =        - 5.0*h*x[ju+i];
      neG    = neG + 1;

      G[neG] =        - 5.0*h*x[ju+i+1];
      neG    = neG + 1;

      G[neG] =        - 0.5*h*(10.0*x[jx2+i]   - x[jx1+i]);
      neG    = neG + 1;

      G[neG] =        - 0.5*h*(10.0*x[jx2+i+1] - x[jx1+i+1]);
      neG    = neG + 1;

      // second ode constraint.
      G[neG] =        - 0.5*h*x[ju+i];
      neG    = neG + 1;

      G[neG] =        - 0.5*h*x[ju+i+1];
      neG = neG + 1;

      G[neG] = - 1.0  + 0.5*h*(9.0*x[ju+i]   + 1.0);
      neG    = neG + 1;

      G[neG] =   1.0  + 0.5*h*(9.0*x[ju+i+1] + 1.0);
      neG    = neG + 1;

      G[neG] = - 0.5*h*(x[jx1+i]   - 9.0*x[jx2+i]);
      neG    = neG + 1;

      G[neG] = - 0.5*h*(x[jx1+i+1] - 9.0*x[jx2+i+1]);
      neG    = neG + 1;
    }
  }

  if (*needG > 0) {
    G[neG] = -Gobj0;
    neG    = neG + 1;
  }

  if (*needF > 0) {
    F[ObjRow] = Fobj;
  }
}

void mySTOP
( int *iAbort, int KTcond[], int *MjrPrt, int *minimz,
  int *m, int *maxS, int *n, int *nb,
  int *nnCon0, int *nnCon, int *nnObj0, int *nnObj, int *nS,
  int *itn, int *nMajor, int *nMinor, int *nSwap,
  double *condHz, int *iObj, double *sclObj, double *ObjAdd,
  double *fObj, double *fMrt, double *PenNrm, double *step,
  double *prInf, double *duInf, double *vimax, double *virel, int hs[],
  int *ne, int *nlocJ, int locJ[], int indJ[], double Jcol[], int *negCon,
  double Ascale[], double bl[], double bu[], double Fx[], double fCon[], double gCon[], double gObj[],
  double yCon[], double pi[], double rc[], double rg[], double x[],
  char cu[], int *lencu, int iu[], int *leniu, double ru[], int *lenru,
  char cw[], int *lencw, int iw[], int *leniw, double rw[], int *lenrw ) {

  std::cout << "hi from mystop fun" << std::endl;
  iu[5] = 15;
  std::cout << "bye from mystop fun" << std::endl;
  iAbort = 0;
}


int main(int argc, char **argv) {
  snoptProblemA catmixa;

  int Cold  = 0;

  int nH    = 1000;
  int n     = 3*(nH+1);
  int nCon  = 2*nH;
  int nF    = nCon + 1;

  int    ObjRow;
  double ObjAdd = -1.0;

  double *x      = new double[n];
  double *xlow   = new double[n];
  double *xupp   = new double[n];
  double *xmul   = new double[n];
  int    *xstate = new    int[n];

  double *F      = new double[nF];
  double *Flow   = new double[nF];
  double *Fupp   = new double[nF];
  double *Fmul   = new double[nF];
  int    *Fstate = new int[nF];

  int lenA   = 2;
  int *iAfun = new int[lenA];
  int *jAvar = new int[lenA];
  double *A  = new double[lenA];

  int lenG   = 14*nH + 1;
  int *iGfun = new int[lenG];
  int *jGvar = new int[lenG];

  int iopt;
  double ropt;

  int nS = 0, nInf = 0, neA = 0, neG = 0;
  int jx1, jx2, ju, ode1, ode2, Obj;
  double sInf;

  double inf = 1.0e20;

  int    leniu = 100, lenru = 100;
  int    iu[leniu];
  double ru[lenru];

  jx1 = 0;
  jx2 = jx1 + nH + 1;
  ju  = jx2 + nH + 1;

  ode1 = 0;
  ode2 = ode1 + nH;
  Obj  = ode2 + nH;

  ObjRow = Obj;

  // Linear terms first
  iAfun[neA] = ObjRow;
  jAvar[neA] = jx1 + nH;
  A[neA]     = 1.0;
  neA        = neA + 1;

  iAfun[neA] = ObjRow;
  jAvar[neA] = jx2 + nH;
  A[neA]     = 1.0;
  neA        = neA + 1;


  for (int i = 0; i < nH; i++) {
    iGfun[neG] = ObjRow;
    jGvar[neG] = ju + i;
    neG        = neG + 1;

    // First ode constraint
    iGfun[neG] = ode1 + i;
    jGvar[neG] = jx1  + i;
    neG        = neG + 1;

    iGfun[neG] = ode1 + i;
    jGvar[neG] = jx1  + i + 1;
    neG        = neG + 1;

    iGfun[neG] = ode1 + i;
    jGvar[neG] = jx2  + i;
    neG        = neG + 1;

    iGfun[neG] = ode1 + i;
    jGvar[neG] = jx2  + i + 1;
    neG        = neG + 1;

    iGfun[neG] = ode1 + i;
    jGvar[neG] = ju   + i;
    neG        = neG + 1;

    iGfun[neG] = ode1 + i;
    jGvar[neG] = ju   + i + 1;
    neG        = neG + 1;

    // Second ode constraint
    iGfun[neG] = ode2 + i;
    jGvar[neG] = jx1  + i;
    neG        = neG + 1;

    iGfun[neG] = ode2 + i;
    jGvar[neG] = jx1  + i + 1;
    neG        = neG + 1;

    iGfun[neG] = ode2 + i;
    jGvar[neG] = jx2  + i;
    neG        = neG + 1;

    iGfun[neG] = ode2 + i;
    jGvar[neG] = jx2  + i + 1;
    neG        = neG + 1;

    iGfun[neG] = ode2 + i;
    jGvar[neG] = ju   + i;
    neG        = neG + 1;

    iGfun[neG] = ode2 + i;
    jGvar[neG] = ju   + i + 1;
    neG        = neG + 1;
  }

  iGfun[neG] = ObjRow;
  jGvar[neG] = ju + nH;
  neG        = neG + 1;

  // Set bounds, states and initial values.
  for (int i = 0; i <= nH; i++) {
    // x1 variable
    xlow[jx1+i]   = -inf;
    xupp[jx1+i]   =  inf;
    x[jx1+i]      =  1.0;
    xstate[jx1+i] =  0;

    // x2 variable
    xlow[jx2+i]   = -inf;
    xupp[jx2+i]   =  inf;
    x[jx2+i]      =  0.0;
    xstate[jx2+i] =  0;

    // u variable
    xlow[ju+i]   = 0.0;
    xupp[ju+i]   = 1.0;
    x[ju+i]      = 0.0;
    xstate[ju+i] = 0;
  }

  xlow[jx1] = 1.0;
  xupp[jx1] = 1.0;
  x[jx1]    = 1.0;

  xlow[jx2] = 0.0;
  xupp[jx2] = 0.0;
  x[jx2]    = 0.0;

  // Bounds on F (all equalities)
  for (int i = 0; i < nCon; i++) {
    Flow[i] = 0.0;
    Fupp[i] = 0.0;
    Fmul[i] = 0.0;
  }

  // Objective row
  Fmul[ObjRow] = 0.0;
  Flow[ObjRow] = -inf;
  Fupp[ObjRow] =  inf;

  catmixa.initialize     ("", 1);  // no print file, summary on
  catmixa.setProbName    ("catmix");

  catmixa.setPrintFile   ("catmix.out");  // ok now add a print file
  catmixa.setIntParameter("Verify level ", 3);
  catmixa.setParameter("Sticky parameters yes");

  catmixa.getIntParameter("Verify level", iopt);
  std::cout << "Verify level set to: " << iopt << std::endl;


  for (int i = 0; i < leniu; i++) {
    iu[i] = 9;
    ru[i] = 0.9;
  }

  catmixa.setUserI(iu, leniu);
  catmixa.setUserR(ru, lenru);

  catmixa.setSTOP(mySTOP);

  catmixa.solve          (Cold, nF, n, ObjAdd, ObjRow, usrFG,
			  iAfun, jAvar, A, neA,
			  iGfun, jGvar, neG,
			  xlow, xupp, Flow, Fupp,
			  x, xstate, xmul,
			  F, Fstate, Fmul,
			  nS, nInf, sInf);

  catmixa.getRealParameter("Major optimality tolerance", ropt);
  std::cout << "Opt tolerance was set to: " << ropt << std::endl;



  delete []iAfun;  delete []jAvar;  delete []A;
  delete []iGfun;  delete []jGvar;

  delete []x;      delete []xlow;   delete []xupp;
  delete []xmul;   delete []xstate;

  delete []F;      delete []Flow;   delete []Fupp;
  delete []Fmul;   delete []Fstate;

}
