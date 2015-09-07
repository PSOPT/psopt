#ifndef TOYFUNCTIONB
#define TOYFUNCTIONB

#ifdef __cplusplus
extern "C" {
#endif

void toyobjB ( int *mode,  int *nnObj, double x[],
	       double *fObj,  double gObj[], int *nState,
	       char    *cu, int *lencu,
	       int    iu[], int *leniu,
	       double ru[], int *lenru );

void toyconB ( int *mode,  int *nnCon, int *nnJac, int *negCon,
	       double x[], double fCon[], double gCon[], int *nState,
	       char    *cu, int *lencu,
	       int    iu[], int *leniu,
	       double ru[], int *lenru );

#ifdef __cplusplus
}
#endif

#endif
