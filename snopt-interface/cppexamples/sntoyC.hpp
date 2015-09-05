#ifndef TOYFUNCTIONC
#define TOYFUNCTIONC

#ifdef __cplusplus
extern "C" {
#endif

  void toyusrfun ( int *mode,  int *nnObj, int *nnCon,
		   int *nnJac, int *nnL,   int *negCon, double x[],
		   double *fObj,  double gObj[],
		   double fCon[], double gCon[], int *Status,
		   char    *cu, int *lencu,
		   int    iu[], int *leniu,
		   double ru[], int *lenru );

#ifdef __cplusplus
}
#endif

#endif
