
/*********************************************************************************************

This file is part of the PSOPT library, a software tool for computational optimal control

Copyright (C) 2009 Victor M. Becerra

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

Author:    Dr. Victor M. Becerra
           University of Reading
           School of Systems Engineering
           P.O. Box 225, Reading RG6 6AY
           United Kingdom
           e-mail: v.m.becerra@reading.ac.uk

**********************************************************************************************/

#define PSOPT_RELEASE_STRING  "4 BETA"

/* Define to the C type corresponding to Fortran INTEGER */
#define FORTRAN_INTEGER_TYPE int


#ifdef WIN32
extern "C" {
_CRTIMP  int * __cdecl errno(void) { static int i=0; return &i; };
}
#endif

#define CINDEX( i )    ((i)-1)

#ifndef MAX
#define MAX(a, b) ( (a)>(b)?  (a):(b) )
#endif
#ifndef MIN
#define MIN(a, b) ( (a)<(b)?  (a):(b) )
#endif

#include "dmatrixv.h"


#undef max
#undef min
#undef abs


#define ADOLC_VERSION_2

#define FREE_ARG char*

#ifndef WIN32
#include <stdlib.h>
#endif

#include <math.h>
#include <string.h>
#include <time.h>
#include <assert.h>


using namespace std;


#ifdef WIN32

#ifdef ADOLC_VERSION_1
#include "drivers.h"
#include "interfaces.h"
#include "adalloc.h"
#include "adouble.h"
#include "sparse.h"
#endif

#ifdef ADOLC_VERSION_2

#include <adolc/drivers/drivers.h>
#include <adolc/interfaces.h>
#include <adolc/adalloc.h>
#include <adolc/adouble.h>
#include <adolc/sparse/sparsedrivers.h>
#include <adolc/taping.h>
#endif



#else



#ifdef ADOLC_VERSION_1
#include <adolc/drivers/drivers.h>
#include <adolc/interfaces.h>
#include <adolc/adalloc.h>
#include <adolc/adouble.h>
#include <adolc/sparse/sparse.h>
#endif

#ifdef ADOLC_VERSION_2
#include <adolc/drivers/drivers.h>
#include <adolc/interfaces.h>
#include <adolc/adalloc.h>
#include <adolc/adouble.h>
#include <adolc/sparse/sparsedrivers.h>
#include <adolc/taping.h>
#endif

#endif

#include <string>
using std::string;


class ADMatrix {
  adouble* a;
  int n;
  int m;
  public:
  // Constructor
  ADMatrix( int n, int m);
  ADMatrix( int n );
  ADMatrix();
  // Destructor
  ~ADMatrix();
  adouble elem(int i, int j) const;
  adouble operator()(int i, int j) const;
  adouble operator()(int i) const;
  adouble& elem(int i, int j);
  adouble& operator()(int i, int j);
  adouble& operator()(int i);
  int GetNoRows() const { return n;}
  int GetNoCols() const { return m;}
  adouble* GetPr() { return a;}
  adouble* GetConstPr() const { return a;}
  void Print(const char* text) const;
  void Resize(int nn, int mm);
  void FillWithZeros();

};


struct dual_str {
  DMatrix* Hamiltonian;
  DMatrix* costates;
  DMatrix* path;
  DMatrix* events;
  DMatrix* linkages;
};


typedef struct dual_str Dual;


typedef struct {
int     nnodes;
int     nvars;
int     ncons;
int     n_obj_evals;
int     n_con_evals;
int     n_jacobian_evals;
int     n_hessian_evals;
int     n_ode_rhs_evals;
double  epsilon_max;
double  CPU_time;
string  method;
} MeshStats;




struct guess_str {
  DMatrix controls;
  DMatrix states;
  DMatrix parameters;
  DMatrix time;
};

typedef struct guess_str Guess;


struct alg_str {
  int       nlp_iter_max;
  double    nlp_tolerance;
  string    nlp_method;
  string    scaling;
  string    derivatives;
  string    constraint_scaling;
  string    ps_method;
  string    collocation_method;
  string    hessian;
  string    defect_scaling;
  string    diff_matrix;
  string    parameter_statistics;
  double    jac_sparsity_ratio;
  double    hess_sparsity_ratio;
  int       print_level; // 1: detailed output on screen and files (default), 0: no output
  int       save_sparsity_pattern;
  int       nsteps_error_integration;
  int       parameter_estimation_norm;


  double    ode_tolerance;
  double    mr_max_increment_factor;
  int	    mr_max_iterations;
  int       mr_min_extrapolation_points;
  int       mr_initial_increment;
  double    mr_kappa;
  int       mr_M1;
  string    mesh_refinement;
  int       switch_order;
  double    ipopt_max_cpu_time;


};

typedef struct alg_str Alg;


struct ulbounds_str {

    DMatrix events;
    DMatrix path;
    DMatrix states;
    DMatrix controls;
    DMatrix parameters;

    double StartTime;
    double EndTime;

};


typedef struct ulbounds_str ULBounds;

struct bounds_str {
   ULBounds lower;
   ULBounds upper;
};

struct name_str {
   string* states_;
   string* controls_;
   string* parameters_;

   string& states(int i) {return states_[i-1];}
   string& controls(int i) { return controls_[i-1];}
   string& parameters(int i) {return parameters_[i-1];}
};

typedef struct name_str Name;

struct units_str {
    string* states_;
    string* controls_;
    string* parameters_;
    string   time;
    string& states(int i) {return states_[i-1];}
    string& controls(int i) { return controls_[i-1];}
    string& parameters(int i) {return parameters_[i-1];}
};

typedef struct units_str Units;

typedef struct bounds_str Bounds;


struct scaling_str {
   DMatrix  controls;
   DMatrix  states;
   DMatrix  parameters;
   DMatrix  defects;
   DMatrix  path;
   DMatrix  events;
   double    time;
};

typedef struct scaling_str Scaling;

struct prob_scaling_str {
   DMatrix   linkages;
   double    objective;
};

typedef struct prob_scaling_str ProbScaling;


struct phases_str {

   int nstates;

   int ncontrols;

   int nparameters;

   int nobserved;

   int nsamples;

   int nevents;

   int npath;

   int current_number_of_intervals;

   Bounds bounds;

   Guess  guess;

   DMatrix nodes;

   Scaling scale;

   bool zero_cost_integrand;

   DMatrix observations;

   DMatrix covariance;

   DMatrix residual_weights;

   DMatrix observation_nodes;

   double  regularization_factor;

   Name name;

   Units units;

};

typedef struct phases_str Phases;


struct prob_ul_bounds {
     DMatrix linkage;
     DMatrix times;
};

typedef prob_ul_bounds ProbULBounds;


struct prob_bounds_str {
    ProbULBounds lower;
    ProbULBounds upper;
};

typedef struct prob_bounds_str ProbBounds;

typedef struct work_str Workspace;


struct prob_str {

   int nphases;

   Phases* phase;

   ProbScaling scale;

   int nlinkages;

   bool multi_segment_flag;

   bool continuous_controls_flag;

   ProbBounds bounds;

   string  name;

   string  outfilename;

   void* user_data;

   Phases&   phases(int iphase);

   adouble (*endpoint_cost)(adouble* initial_states, adouble* final_states, adouble* parameters,adouble& t0, adouble& tf, adouble* xad, int iphase, Workspace* workspace);

   adouble (*integrand_cost)(adouble* states, adouble* controls, adouble* parameters, adouble& time,  adouble* xad, int iphase, Workspace* workspace);

   void (*dae)(adouble* derivatives, adouble* path, adouble* states, adouble* controls, adouble* parameters, adouble& time, adouble* xad, int iphase, Workspace* workspace);

   void (*events)(adouble* e, adouble* initial_states, adouble* final_states, adouble* parameters,adouble& t0, adouble& tf, adouble* xad, int iphase, Workspace* workspace);

   void (*linkages)( adouble* linkages, adouble* xad, Workspace* workspace);

   void (*observation_function)(adouble* observed_variable, adouble* states, adouble* controls, adouble* parameters, adouble& time, int k, adouble* xad, int iphase, Workspace* workspace);

};

typedef struct prob_str Prob;




struct sol_str {
   DMatrix *states;
   DMatrix *controls;
   DMatrix *nodes;
   DMatrix *parameters;
   DMatrix *relative_errors;
   Dual     dual;
   DMatrix *integrand_cost;
   double  *endpoint_cost;
   double  *integrated_cost;
   double   cost;
   adouble* xad;
   int      nlp_return_code;
   double   cpu_time;
   bool     error_flag;
   string   error_msg;
   int       mesh_refinement_iterations;
   MeshStats*  mesh_stats;
   string   start_date_and_time;
   string   end_date_and_time;
   Prob* problem;
   DMatrix& get_states_in_phase(int iphase);
   DMatrix& get_controls_in_phase(int iphase);
   DMatrix& get_time_in_phase(int iphase);
   DMatrix& get_parameters_in_phase(int iphase);
   DMatrix& get_dual_costates_in_phase(int iphase);
   DMatrix& get_dual_hamiltonian_in_phase(int iphase);
   DMatrix& get_dual_path_in_phase(int iphase);
   DMatrix& get_dual_events_in_phase(int iphase);
   DMatrix& get_dual_linkages();
   DMatrix& get_relative_local_error_in_phase(int iphase);
   double   get_cost() { return cost; }
};

typedef struct sol_str Sol;




typedef struct {
    DMatrix* dfdx_j;
    DMatrix* F1;
    DMatrix* F2;
    DMatrix* F3;
    DMatrix* F4;
    double*  x;
    double*  g;
    adouble* xad;
    adouble* gad;
} GRWORK;


typedef struct {
  int nsegments;
  int nstates;
  int nparameters;
  int ncontrols;
  int npath;
  int ninitial_events;
  int nfinal_events;
  int nobserved;
  int nsamples;
  bool continuous_controls;
  DMatrix nodes;
} MSdata;


typedef struct {

   int** colindex;
   int*  size;
   int   number;

} IGroup;

struct work_str {

   Sol*      solution;
   Prob*     problem;
   Alg*      algorithm;
   DMatrix*  P;
   DMatrix*  sindex;
   DMatrix*  w;
   DMatrix*  D;
   DMatrix*  D2;
   DMatrix*  snodes;
   DMatrix*  old_snodes;
   DMatrix*  xlb;
   DMatrix*  xub;
   DMatrix*  x0;
   DMatrix*  lambda;
   DMatrix*  dual_costates;
   DMatrix*  dual_path;
   DMatrix*  dual_events;
   DMatrix*  Xsnopt;
   DMatrix*  gsnopt;
   DMatrix*  DerivResid;
   DMatrix*  h;
   DMatrix*  Xdotgg;
   DMatrix*  e;
   DMatrix*  hgg;
   DMatrix*  Xip;
   DMatrix*  JacRow;
   DMatrix*  Gip;
   DMatrix*  GFip;
   SparseMatrix*  Ax;
   SparseMatrix*  As;
   SparseMatrix*  Gsp;
   DMatrix*  control_scaling;
   DMatrix*  control_shift;
   DMatrix*  state_scaling;
   DMatrix*  parameter_scaling;
   DMatrix*  state_shift;
   DMatrix*  deriv_scaling;
   DMatrix*  path_scaling;
   DMatrix*  event_scaling;
   DMatrix*  constraint_scaling;
   DMatrix*  linkage;
   DMatrix*  prev_states;
   DMatrix*  prev_costates;
   DMatrix*  prev_controls;
   DMatrix*  prev_param;
   DMatrix*  prev_path;
   DMatrix*  prev_nodes;
   DMatrix*  prev_t0;
   DMatrix*  prev_tf;
   DMatrix*  Xdot;
   DMatrix*  JacCol1;
   DMatrix*  JacCol2;
   DMatrix*  JacCol3;
   DMatrix*  xp;
   DMatrix*  emax_history;
   DMatrix*  order_reduction;
   DMatrix*  old_relative_errors;
   DMatrix*  error_scaling_weights;


   double    obj_scaling;
   GRWORK*   grw;
   int       nvars;
   int       ncons;
   int       **current_nodes;
   int       jac_done;
   int*      iArow;
   int*      jAcol;
   int*      iGrow;
   int*      jGcol;
   double*   jac_values;
   double*   jac_Aij;
   double*   jac_Gij;
   int       jac_nnz;
   int       jac_nnzA;
   int       jac_nnzG;
   double*   nrm_row;
   unsigned int*      hess_ir;
   unsigned int*      hess_jc;
   unsigned int*      iGfun;
   unsigned int*      jGvar;
   unsigned int*      iGfun1;
   unsigned int*      jGvar1;
   unsigned int*      iGfun2;
   unsigned int*      jGvar2;
   int       use_constraint_scaling;
   int       F_nnz;
   double*   G2;
   adouble*  xad;
   adouble*  gad;
   adouble**  states;
   adouble**  controls;
   adouble**  parameters;
   adouble**  resid;
   adouble**  derivatives;
   adouble**  initial_states;
   adouble**  final_states;
   adouble**  initial_controls;
   adouble**  final_controls;
   adouble**  events;
   adouble**  path;
   adouble**  states_traj;
   adouble**  derivs_traj;
   adouble**  second_derivs_traj;
   adouble*   linkages;
   adouble*   fgad;
   adouble*   time_array_tmp;
   adouble*   single_trajectory_tmp;
   adouble*   L_ad_tmp;
   adouble*   u_spline;
   adouble*   z_spline;
   adouble*   y2a_spline;
   adouble**   states_next;
   adouble**   controls_next;
   adouble**   derivatives_next;
   adouble**   path_next;
   adouble**   path_bar;
   adouble**   states_bar;
   adouble**   controls_bar;
   adouble**   derivatives_bar;
   adouble**   observed_variable;
   adouble**   observed_residual;
   adouble**   lam_resid;
   adouble**   interp_states_pe;
   adouble**   interp_controls_pe;
   double*    lambda_d;
   double*    fg;
   bool       trace_f_done;
   IGroup*    igroup;
   char       text[2000];
   FILE*      psopt_solution_summary_file;
   FILE*      mesh_statistics;
   FILE*      mesh_statistics_tex;
   int        current_mesh_refinement_iteration;
   bool       auto_linked_flag;
   bool       enable_nlp_counters;
   string     differential_defects;
   clock_t    start_ticks;

// tape tags to be used by ADOL_C

  int tag_f     ;
  int tag_g 	;
  int tag_hess 	;
  int tag_fg 	;
  int tag_gc    ;
  void *user_data;

};



struct xad_str {
    adouble *xad;
    Workspace *workspace;
};

typedef xad_str XAD;





// Function prototypes

void ScalarGradientAD( adouble (*fun)(adouble *, Workspace*), DMatrix& x, DMatrix* grad, bool* flag, int itag, Workspace* workspace );

void EfficientlyComputeJacobianNonZeros( void fun(DMatrix& x, DMatrix* f, Workspace* ), DMatrix& x,
                int nf, double *nzvalue, int nnz, int* iArow, int* jAcol, IGroup* igroup, GRWORK* grw, Workspace* workspace );

void compute_jacobian_of_constraints_with_respect_to_variables(DMatrix& Jc, DMatrix& X, DMatrix& XL, DMatrix& XU, Workspace* workspace);

void compute_jacobian_of_residual_vector_with_respect_to_variables(DMatrix& Jr, DMatrix& X, DMatrix& XL, DMatrix& XU, Workspace* workspace);

void getIndexGroups( IGroup* igroup, int nrows, int ncols, int nnz, int* iArow, int* jAcol, Workspace* workspace);

void getIndexGroups( IGroup* igroup, int nrows, int ncols, int nnz, int* iArow, int* jAcol);

void deleteIndexGroups(IGroup* igroup, int ncols );

void psopt(Sol& solution, Prob& problem, Alg& algorithm);

void psopt_level2_setup(Prob& problem, Alg& algorithm);

void initialize_solution(Sol& solution, Prob& problem, Alg& algorithm, Workspace* workspace);

void lglnodes(int N, DMatrix& x, DMatrix& w, DMatrix& P, DMatrix& D, Workspace* workspace);

void cglnodes(int N, DMatrix& x, DMatrix& w,  DMatrix& D, Workspace* workspace);

void copy_decision_variables(Sol& solution, DMatrix& x, Prob& problem, Alg& algorithm, Workspace* workspace);

double ff( DMatrix& x );

void gg( DMatrix& x, DMatrix* g );

void  define_initial_nlp_guess(DMatrix& x0, DMatrix& lambda, Sol& solution, Prob& problem, Alg& algorithm, Workspace* workspace);

void determine_scaling_factors_for_variables(Sol& solution, Prob& problem, Alg& algorithm);

void determine_objective_scaling(DMatrix& X,Sol& solution, Prob& problem, Alg& algorithm, Workspace* workspace );

void determine_constraint_scaling_factors(DMatrix & X, Sol& solution, Prob& problem, Alg& algorithm, Workspace* workspace);

void  define_nlp_bounds(DMatrix& xlb,DMatrix&  xub,Prob& problem, Alg& algorithm, Workspace* workspace);

double convert_to_original_time(double tbar,double t0,double tf);

void resize_solution(Sol& solution, Prob& problem, Alg& algorithm);

void hot_start_nlp_guess(DMatrix& x0,DMatrix& lambda, Sol& solution,Prob& problem,Alg& algorithm, DMatrix* prev_states, DMatrix* prev_controls, DMatrix* prev_costates, DMatrix* prev_path, DMatrix* prev_nodes, DMatrix* prev_param, DMatrix& prev_t0, DMatrix& prev_tf, Workspace* workspace );

void lagrange_interpolation(DMatrix& y, DMatrix& x, DMatrix& pointx, DMatrix& pointy);

double smooth_fabs(double x, double eps);

adouble smooth_fabs(adouble x, double eps);

void Jacobian( void fun(DMatrix& x, DMatrix* f ), DMatrix& x, DMatrix* grad, GRWORK* grw );

void JacobianRow( void fun(DMatrix& x, DMatrix* f, Workspace* ), DMatrix& x, int iRow, int nf,
                  DMatrix* JacRow, GRWORK* grw, Workspace* workspace );

void JacobianColumn( void fun(DMatrix& x, DMatrix* f, Workspace* ), DMatrix& x, DMatrix& xlb, DMatrix& xub, int jCol,
                DMatrix* JacColumn, GRWORK* grw, Workspace* workspace );

void evaluate_matrix_of_integrated_errors_in_phase(DMatrix& eta, int iphase, adouble* xad, int nsteps, Workspace* workspace);

void evaluate_solution(Prob& problem,Alg& algorithm,Sol& solution, Workspace* workspace);

void compute_next_mesh_size( Prob& problem, Alg& algorithm, Sol& solution, Workspace* workspace );

int get_max_nodes(Prob& problem,int iphase, Alg* algorithm);

void estimate_order_reduction(Prob& problem,Alg& algorithm,Sol& solution, Workspace* workspace);

void zero_order_reduction(Prob& problem,Alg& algorithm,Sol& solution, Workspace* workspace);

void construct_new_mesh(Prob& problem,Alg& algorithm,Sol& solution, Workspace* workspace);

bool check_for_equidistributed_error(Prob& problem,Alg& algorithm,Sol& solution);

bool use_local_collocation(Alg & algorithm);

bool use_global_collocation(Alg & algorithm);

void bilinear_interpolation(adouble* z, adouble& x, adouble& y, DMatrix& X, DMatrix& Y, DMatrix& Z);

void spline_2d_interpolation(adouble* z, adouble& x, adouble& y, DMatrix& X, DMatrix& Y, DMatrix& Z, Workspace* workspace);

void smooth_bilinear_interpolation(adouble* z, adouble& x, adouble& y, DMatrix& X, DMatrix& Y, DMatrix& Z);

void smooth_linear_interpolation(adouble* y, adouble& x, adouble* Xdata, adouble* Ydata, int n);

void get_times(adouble *t0, adouble *tf, adouble* xad, int iphase, Workspace* workspace);

void get_states(adouble* states, adouble* xad, int iphase, int k, Workspace* workspace);

void get_controls_bar(adouble* controls_bar, adouble* xad, int iphase, int k, Workspace* workspace);

void rr_num(DMatrix& X, DMatrix* residual_vector, Workspace* workspace);

void extract_parameter_covariance(DMatrix& Cp, DMatrix& C, Workspace* workspace);



#define INF  (1.0e19)
#define inf  (1.0e19)

static const double pi = 3.141592653589793;

int NLP_interface(
         Alg& algorithm,
         DMatrix* x0,
         double (*f)(DMatrix&, Workspace*),
	 void (*g)(DMatrix&,DMatrix*,Workspace*),
	 int m          = 0,
	 int neqcstr    = 0,
	 DMatrix* xlb   = NULL,
	 DMatrix* xub   = NULL,
         DMatrix* lambda= NULL,
         int hotflag    = 0,
         int iprint     = 1,
         Workspace* workspace=NULL,
         void* user_data = NULL  );



void ScalarGradient( double (*fun)(DMatrix& x, Workspace*), DMatrix& x,DMatrix* grad, GRWORK* grw, Workspace* workspace );

void DetectJacobianSparsity(void fun(DMatrix& x, DMatrix* f, Workspace* ), DMatrix& x, int nf,
                           int* nnzA, int* iArow, int* jAcol, double* Aij,
                           int* nnzG, int* jGrow, int* jGcol,
                           GRWORK* grw, Workspace* workspace);

void ComputeJacobianNonZeros( void fun(DMatrix& x, DMatrix* f ), DMatrix& x, int nf, double *nzvalue, int nnz, int* iArow, int* jAcol, GRWORK* grw, Workspace* workspace );

void Jacobian( void fun(DMatrix& x, DMatrix* f ), DMatrix& x,
                DMatrix* grad, GRWORK* grw );

void initialize_workspace_vars(Prob& problem, Alg& algorithm, Sol& solution, Workspace* workspace);

void resize_workspace_vars(Prob& problem, Alg& algorithm, Sol& solution, Workspace* workspace);

int get_number_nlp_vars(Prob& problem, Workspace* workspace);

int get_number_nlp_constraints(Prob& problem, Workspace* workspace);

int get_max_number_nlp_vars(Prob& problem, Alg& algorithm);

int get_max_number_nlp_constraints(Prob& problem, Alg& algorithm);

void psopt_level1_setup(Prob& problem);

int get_max_nodes_in_all_phases(Prob& problem, Alg& algorithm);

adouble convert_to_original_time_ad(double tbar,adouble& t0,adouble& tf);

int get_nvars_phase_i(Prob& problem, int i, Workspace* workspace);

int get_ncons_phase_i(Prob& problem, int i, Workspace* workspace);

void get_constraint_bounds(double* g_l, double* g_u, Workspace* workspace);

void mtrx_mul_trans(adouble* a,double* b,adouble* ab,int na,int ma,int nb,int mb);

void get_initial_states(adouble* states, adouble* xad, int i, Workspace* workspace);

void get_final_states(adouble* states, adouble* xad, int i, Workspace* workspace);

void get_controls(adouble* controls, adouble* xad, int i, int k, Workspace* workspace);

void get_final_controls(adouble* controls, adouble* xad, int i, Workspace* workspace);

void get_initial_controls(adouble* controls, adouble* xad, int i, Workspace* workspace);

void get_individual_control_trajectory(adouble *control_traj, int control_index, int iphase, adouble* xad, Workspace* workspace);

void get_individual_state_trajectory(adouble *state_traj, int state_index, int iphase, adouble* xad, Workspace* workspace);

void get_delayed_control(adouble* delayed_control, int control_index, int iphase, adouble& time, double delay, adouble* xad, Workspace* workspace);

void get_delayed_state(adouble* delayed_state, int state_index, int iphase, adouble& time, double delay, adouble* xad, Workspace* workspace);

void get_interpolated_state(adouble* interp_state, int state_index, int iphase, adouble& time, adouble* xad, Workspace* workspace);

void get_interpolated_control(adouble* interp_control, int control_index, int iphase, adouble& time, adouble* xad, Workspace* workspace);

void get_state_derivative(adouble* state_derivative, int state_index, int iphase, adouble& time, adouble* xad, Workspace* workspace);

void get_control_derivative(adouble* control_derivative, int control_index, int iphase, adouble& time, adouble* xad, Workspace* workspace);

int get_number_of_controls(Prob& problem, int iphase);

int get_number_of_states(Prob& problem,int iphase);

int get_number_of_events(Prob& problem,int iphase);

int get_number_of_path_constraints(Prob& problem,int iphase);

int get_number_of_nodes(Prob& problem,int iphase);

int get_number_of_parameters(Prob& problem,int iphase);

int get_number_of_linkages(Prob& problem);

int get_number_of_phases(Prob& problem);

void lagrange_interpolation_ad(adouble* y, adouble& x, adouble* pointx, adouble* pointy, int npoints, Workspace* workspace);

void linear_interpolation(adouble* y, adouble& x, adouble* pointx, adouble* pointy, int npoints);

void linear_interpolation(DMatrix& y, double x, DMatrix& pointx, DMatrix& pointy, int npoints);

void linear_interpolation(DMatrix& y, DMatrix& x, DMatrix& pointx, DMatrix& pointy, int npoints);

void linear_interpolation(adouble* y, adouble x, DMatrix& pointx, DMatrix& pointy, int npoints);

void smooth_linear_interpolation(adouble* y, adouble& x, DMatrix& Xdata, DMatrix& Ydata, int n);

void spline_interpolation(adouble* y, adouble& x, adouble* pointx, adouble* pointy, int npoints, Workspace* workspace);

void spline_interpolation(adouble* y, adouble& x, DMatrix& Xdata, DMatrix& Ydata, int n);

void zoh_interpolation(adouble* y, adouble x, DMatrix& pointx, DMatrix& pointy, int npoints);


adouble get_initial_time(adouble* xad, int iphase, Workspace* workspace);

adouble get_final_time(adouble* xad, int iphase, Workspace* workspace);

void get_parameters(adouble* parameters, adouble* xad, int iphase, Workspace* workspace);

bool useAutomaticDifferentiation(Alg& algorithm);

void gg_ad( adouble* xad, adouble* gad, Workspace* workspace );

double ff_num(DMatrix& x, Workspace* workspace);

void gg_num( DMatrix& x, DMatrix* g, Workspace* workspace );

void fg_ad( adouble* x, adouble* fg, Workspace* workspace);

void compute_derivatives_trajectory( DMatrix& Xdot, Prob& problem, Sol& solution,  int i, Workspace* workspace );

adouble integrate( adouble (*integrand)(adouble*,adouble*,adouble*,adouble&,adouble*,int, Workspace* workspace), adouble* xad, int i, Workspace* workspace );

void auto_link(adouble* linkages, int* index, adouble* xad, int iphase_a, int iphase_b, Workspace* workspace);

void auto_link_2(adouble* linkages, int* index, adouble* xad, int iphase_a, int iphase_b, Workspace* workspace);

void plot(DMatrix& x, DMatrix& y,const string& title,
          char* xlabel, char* ylabel, char* legend=NULL, char* terminal=NULL, char* output=NULL);

void plot(DMatrix& x1, DMatrix& y1, DMatrix& x2, DMatrix& y2, const string& title,
          char* xlabel, char* ylabel, char* legend=NULL, char* terminal=NULL, char* output=NULL);

void plot(DMatrix& x1, DMatrix& y1, DMatrix& x2, DMatrix& y2, DMatrix& x3, DMatrix& y3,
          const string& title, char* xlabel, char* ylabel, char* legend=NULL, char* terminal=NULL, char* output=NULL);

void multiplot(DMatrix& x, DMatrix& y, const string& title, char* xlabel, char* ylabel, char* legend, int nrows=0, int ncols=0,  char* terminal=NULL, char* output=NULL ) ;

void spplot(DMatrix& x1a, DMatrix& y1a, DMatrix& x2a, DMatrix& y2a, const string& title, char* xlabel, char* ylabel, char* legend, char* terminal=NULL, char* output=NULL);

void polar(DMatrix& theta, DMatrix& r, const string& title,
           char* legend=NULL, char* terminal=NULL, char* output=NULL);

void polar(DMatrix& theta, DMatrix& r, DMatrix& theta2, DMatrix& r2, const string& title,
            char* legend=NULL, char* terminal=NULL, char* output=NULL);

void polar(DMatrix& theta, DMatrix& r, DMatrix& theta2, DMatrix& r2,  DMatrix& theta3, DMatrix& r3, const string& title,
            char* legend=NULL, char* terminal=NULL, char* output=NULL);

void surf(DMatrix& x, DMatrix& y, DMatrix& z, const string& title, char* xlabel, char* ylabel, char* zlabel, char* terminal=NULL, char* output=NULL, char* view=NULL);

void plot3(DMatrix& x, DMatrix& y, DMatrix& z, const string& title, char* xlabel, char* ylabel, char* zlabel, char* terminal=NULL, char* output=NULL, char* view=NULL);

void psopt_error_message(const char *error_text);

adouble dot(adouble* x, adouble* y, int n);

adouble smooth_heaviside(adouble x, double a);

adouble smooth_sign(adouble x, double a);

void cross(adouble* x, adouble* y, adouble* z);

void validate_user_input(Prob& problem, Alg& algorithm, Workspace* workspace);

void print_psopt_summary(Prob& problem, Alg& algorithm, Sol& solution, Workspace* workspace);

void psopt_main(Sol& solution, Prob& problem, Alg& algorithm);

void clip_vector_given_bounds(DMatrix& xp, DMatrix& xlb, DMatrix& xub);

void psopt_print(Workspace* workspace, char* msg);

int auto_link_count(Prob& problem, int nstates);

void multi_segment_setup(Prob& problem, Alg& algorithm, MSdata& msdata);

int auto_link2_count(Prob& problem, int nstates, int ncontrols);

void auto_phase_setup(Prob& problem, int n_final_events,DMatrix& nodes);

void  auto_phase_bounds(Prob& problem);

void  auto_phase_guess(Prob& problem,DMatrix& controls, DMatrix& states, DMatrix& param, DMatrix& time);

void auto_link_multiple(adouble* linkages, adouble* xad,int nphases, Workspace* workspace);

void auto_link2_multiple(adouble* linkages, adouble* xad,int nphases, Workspace* workspace);

void product_ad(const DMatrix& A, const adouble* x, int nx, adouble* y);

void sum_ad(const adouble* a, const adouble*b, int n, adouble* c);

void subtract_ad(const adouble* a, const adouble*b, int n, adouble* c);

void inverse_ad(adouble* minput, int n, adouble* minv);

void product_ad(adouble* Apr,adouble* Bpr, int na, int ma, int nb, int mb, adouble* ABpr);

bool need_midpoint_controls(Alg& algorithm, Workspace* workspace);

void inverse_ad(const ADMatrix& minput, ADMatrix* minv);

void product_ad(const ADMatrix& A,const ADMatrix& B, ADMatrix* AB);

void sum_ad(const ADMatrix& A,const ADMatrix& B, ADMatrix* AB);

void subtract_ad(const ADMatrix& A,const ADMatrix& B, ADMatrix* AB);

void print_constraint_summary(Prob& problem, Sol& solution, Workspace* workspace);

int get_number_of_mesh_refinement_iterations(Prob& problem, Alg& algorithm);

void update_mesh_statistics(Prob& problem,Alg& solution);

void   chronometer_tic(Workspace* workspace);

double chronometer_toc(Workspace* workspace);

void print_iterations_summary(Prob& problem,Alg& algorithm,Sol& solution, Workspace* workspace);

void print_iterations_summary_tex(Prob& problem,Alg& algorithm,Sol& solution, Workspace* workspace);

void get_local_time( string& date_and_time);

void print_algorithm_summary(Prob& problem, Alg& algorithm, Sol& solution, Workspace* workspace);

void print_solution_summary(Prob& problem, Alg& algorithm, Sol& solution,Workspace* workspace);

void transpose_ad(adouble* Apr, int na, int ma,  adouble* Atpr);

void resample_trajectory(DMatrix& Y, DMatrix& X, DMatrix& Ydata, DMatrix& Xdata);

void load_parameter_estimation_data(Prob& problem, int iphase, char* filename);

bool compute_parameter_statistics(DMatrix& Qp, DMatrix& p, DMatrix& plow, DMatrix& phigh, DMatrix& r, Workspace* workspace);

double inverse_twotailed_t_cdf(double A, int ndf);

double nint(double x);

int get_iphase_offset(Prob& problem, int iphase,Workspace* workspace);

adouble ff_ad(adouble* xad, Workspace* workspace);


void rk4_propagate( void (*dae)(adouble* derivatives, adouble* path, adouble* states,
         adouble* controls, adouble* parameters, adouble& time,
        adouble* xad, int iphase, Workspace* workspace),
        DMatrix& control_trajectory,
        DMatrix& time_vector,
        DMatrix& initial_state,
	DMatrix& parameters,
        Prob & problem,
        int iphase,
        DMatrix& state_trajectory, Workspace* workspace);

void rkf_propagate( void (*dae)(adouble* derivatives, adouble* path, adouble* states,
         adouble* controls, adouble* parameters, adouble& time,
        adouble* xad, int iphase, Workspace* workspace),
        DMatrix& control_trajectory,
        DMatrix& time_vector,
        DMatrix& initial_state,
	DMatrix& parameters,
        double tolerance,
        double hmin,
	double hmax,
        Prob & problem,
        int iphase,
        DMatrix& state_trajectory,
        DMatrix& new_time_vector,
	DMatrix& new_control_trajectory, Workspace* workspace);


void auto_split_observations(Prob& problem, DMatrix& observation_nodes, DMatrix& observations);

adouble endpoint_cost_for_parameter_estimation(adouble* initial_states, adouble* final_states, adouble* parameters,adouble& t0, adouble& tf, adouble* xad, int iphase, Workspace* workspace);


extern "C" {
   int dgeqrf_(integer *m, integer *n, doublereal *a, integer *
	lda, doublereal *tau, doublereal *work, integer *lwork, integer *info);
int dormqr_(char *side, char *trans, integer *m, integer *n,
        integer *k, doublereal *a, integer *lda, doublereal *tau, doublereal *
        c__, integer *ldc, doublereal *work, integer *lwork, integer *info,
        ftnlen side_len, ftnlen trans_len);

}



#ifdef USE_SNOPT
extern "C" {
void snPSOPTusrf_
( int    *Status, int *n,    double x[],
  int    *needF,  int *neF,  double F[],
  int    *needG,  int *neG,  double G[],
  char       *cu, int *lencu,
  int    iu[],    int *leniu,
  double ru[],    int *lenru );




}


extern Workspace* tempsnoptworkspace;



#endif // USE_SNOPT



#ifdef USE_IPOPT

#ifdef PACKAGE
#undef PACKAGE
#undef PACKAGE_VERSION
#undef PACKAGE_BUGREPORT
#undef PACKAGE_NAME
#undef PACKAGE_STRING
#undef PACKAGE_TARNAME
#undef VERSION
#endif


#include "IpIpoptApplication.hpp"
#include "ipopt_psopt.h"



using namespace Ipopt;

#endif


