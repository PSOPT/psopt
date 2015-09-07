#ifndef TOYFUNCTIONA
#define TOYFUNCTIONA

#ifdef __cplusplus
extern "C" {
#endif

  void toyusrf_ ( int    *Status, int *n,    double x[],
		  int    *needF,  int *neF,  double F[],
		  int    *needG,  int *neG,  double G[],
		  char   *cu,     int *lencu,
		  int    iu[],    int *leniu,
		  double ru[],    int *lenru );

  void toyusrfg_( int    *Status, int *n,    double x[],
		  int    *needF,  int *neF,  double F[],
		  int    *needG,  int *neG,  double G[],
		  char   *cu,     int *lencu,
		  int    iu[],    int *leniu,
		  double ru[],    int *lenru);

#ifdef __cplusplus
}
#endif

#endif
