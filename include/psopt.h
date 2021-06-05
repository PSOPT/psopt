
/*********************************************************************************************
This file is part of the PSOPT library, a software tool for computational optimal control
Copyright (C) 2009-2020 Victor M. Becerra
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
Author:    Professor Victor M. Becerra
Address:   University of Portsmouth
           School of Energy and Electronic Engineering
           Portsmouth PO1 3DJ
           United Kingdom
e-mail:    v.m.becerra@ieee.org
**********************************************************************************************/
#ifndef PSOPT_H
#define PSOPT_H

#include <Eigen/Dense>

#include <limits>
using std::numeric_limits;

namespace PSOPT {
    constexpr double inf        =   std::numeric_limits<double>::infinity();
    constexpr double pi         =   4.0*atan(1.0);
}


// using namespace Eigen;
using Eigen::MatrixXd;
using Eigen::RowVectorXi;



// typedef Eigen::SparseMatrix<double> SpMat; // declares a column-major sparse matrix type of double


#define DMatrix MatrixXd  // In case a change of type is forgotten



#ifdef WIN32
extern "C" {
_CRTIMP  int * __cdecl errno(void) { static int i=0; return &i; };
}
#endif

#define CINDEX( i )    ((i)-1)



#define FREE_ARG char*

#ifndef WIN32
#include <stdlib.h>
#endif

#include <math.h>
#include <string.h>
#include <time.h>
#include <assert.h>
#include <memory>



using namespace std;


#include <adolc/drivers/drivers.h>
#include <adolc/interfaces.h>
#include <adolc/adalloc.h>
#include <adolc/adouble.h>
#include <adolc/sparse/sparsedrivers.h>
#include <adolc/taping.h>


#include <string>
using std::string;


class TripletSparseMatrix;

class AutoDiffMatrix {
  adouble* a;
  int n;
  int m;
  public:
  // Constructor
  AutoDiffMatrix( int n, int m);
  AutoDiffMatrix( int n );
  AutoDiffMatrix();
  // Destructor
  ~AutoDiffMatrix();
  adouble elem(int i, int j) const;
  adouble operator()(int i, int j) const;
  adouble operator()(int i) const;
  adouble& elem(int i, int j);
  adouble& operator()(int i, int j);
  adouble& operator()(int i);
  int rows() const { return n;}
  int cols() const { return m;}
  adouble* GetPr() { return a;}
  adouble* GetConstPr() const { return a;}
  void Print(const char* text) const;
  void resize(int nn, int mm);
  void setZero();
  void setOne();

};


class dual_str {
public:
  MatrixXd* Hamiltonian;
  MatrixXd* costates;
  MatrixXd* path;
  MatrixXd* events;
  MatrixXd* linkages;

  dual_str()
  {
    Hamiltonian = NULL;
    costates = NULL;
    path = NULL;
    events = NULL;
    linkages = NULL;
  }

  ~dual_str()
  {
    if (Hamiltonian) delete [] Hamiltonian;
    if (costates) delete [] costates;
    if (path) delete [] path;
    if (events) delete [] events;
    if (linkages) delete linkages;
  }
};


typedef class dual_str Dual;


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
  MatrixXd controls;
  MatrixXd states;
  MatrixXd parameters;
  MatrixXd time;
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
  string    ipopt_linear_solver;
  double    jac_sparsity_ratio;
  double    hess_sparsity_ratio;
  int       print_level; // 1: detailed output on screen and files (default), 0: no output
  int       save_sparsity_pattern;
  int       nsteps_error_integration;
  int       parameter_estimation_norm;



  double    ode_tolerance;
  double    mr_max_increment_factor;
  int	      mr_max_iterations;
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

    MatrixXd events;
    MatrixXd path;
    MatrixXd states;
    MatrixXd controls;
    MatrixXd parameters;

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
   MatrixXd  controls;
   MatrixXd  states;
   MatrixXd  parameters;
   MatrixXd  defects;
   MatrixXd  path;
   MatrixXd  events;
   double    time;
};

typedef struct scaling_str Scaling;

struct prob_scaling_str {
   MatrixXd   linkages;
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

   RowVectorXi nodes;

   Scaling scale;

   bool zero_cost_integrand;

   MatrixXd observations;

   MatrixXd covariance;

   MatrixXd residual_weights;

   MatrixXd observation_nodes;

   double  regularization_factor;

   Name name;

   Units units;

};

typedef struct phases_str Phases;


struct prob_ul_bounds {
     MatrixXd linkage;
     MatrixXd times;
};

typedef prob_ul_bounds ProbULBounds;


struct prob_bounds_str {
    ProbULBounds lower;
    ProbULBounds upper;
};

typedef struct prob_bounds_str ProbBounds;

typedef class work_str Workspace;


class prob_str {
public:

   prob_str()
   {
       phase = NULL;
   }
   ~prob_str()
   {
       if (phase)
       {
         delete [] phase;
       }
   }

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

typedef class prob_str Prob;




class sol_str {
public:
   sol_str()
   {
      states = NULL;
      controls = NULL;
      nodes = NULL;
      parameters = NULL;
      relative_errors = NULL;
      integrand_cost = NULL;
      endpoint_cost = NULL;
      integrated_cost = NULL;
      xad = NULL;
      mesh_stats = NULL;
   }
   ~sol_str()
   {
      if (this->states) delete [] this->states;
      if (this->controls) delete [] this->controls;
      if (this->nodes) delete [] this->nodes;
      if (this->integrand_cost) delete [] this->integrand_cost;
      if (this->parameters) delete [] this->parameters;
      if (this->relative_errors) delete [] this->relative_errors;
      if (this->endpoint_cost) delete [] this->endpoint_cost;
      if (this->integrated_cost) delete [] this->integrated_cost;
      if (this->mesh_stats) delete [] this->mesh_stats;
   }
   MatrixXd *states;

   MatrixXd *controls;
   MatrixXd *nodes;
   MatrixXd *parameters;
   MatrixXd *relative_errors;
   MatrixXd *integrand_cost;
   Dual     dual;

   double  *endpoint_cost;
   double  *integrated_cost;
   double   cost;
   adouble* xad;
   int      nlp_return_code;
   double   cpu_time;
   int      error_flag;
   string   error_msg;
   int       mesh_refinement_iterations;
   MeshStats*  mesh_stats;
   string   start_date_and_time;
   string   end_date_and_time;
   Prob* problem;
   MatrixXd& get_states_in_phase(int iphase);
   MatrixXd& get_controls_in_phase(int iphase);
   MatrixXd& get_time_in_phase(int iphase);
   MatrixXd& get_parameters_in_phase(int iphase);
   MatrixXd& get_dual_costates_in_phase(int iphase);
   MatrixXd& get_dual_hamiltonian_in_phase(int iphase);
   MatrixXd& get_dual_path_in_phase(int iphase);
   MatrixXd& get_dual_events_in_phase(int iphase);
   MatrixXd& get_dual_linkages();
   MatrixXd& get_relative_local_error_in_phase(int iphase);
   double   get_cost() { return cost; }
};

typedef class sol_str Sol;




typedef struct {
    MatrixXd* dfdx_j;
    MatrixXd* F1;
    MatrixXd* F2;
    MatrixXd* F3;
    MatrixXd* F4;
    double*  x;
    double*  g;
    adouble* xad;
    adouble* gad;
} GRWORK;


class MSdata {
public:
  MSdata() {
    nodes.resize(1);  
  }
  MSdata(long n) {
    nodes.resize(n);  
  }
  int nsegments;
  int nstates;
  int nparameters;
  int ncontrols;
  int npath;
  int ninitial_events;
  int nfinal_events;
  bool continuous_controls;
  RowVectorXi nodes;
};


typedef struct {

   int** colindex;
   int*  size;
   int   number;

} IGroup;

class work_str {
public:
   work_str(Prob& problem, Alg& algorithm, Sol& solution);
   ~work_str();
   long unsigned int nphases;

   Sol*      solution;
   Prob*     problem;
   Alg*      algorithm;
   MatrixXd*  P;
   RowVectorXi*  sindex;
   MatrixXd*  w;
   MatrixXd*  D;
   MatrixXd*  D2;
   MatrixXd*  snodes;
   MatrixXd*  old_snodes;
   MatrixXd*  xlb;
   MatrixXd*  xub;
   MatrixXd*  x0;
   MatrixXd*  lambda;
   MatrixXd*  dual_costates;
   MatrixXd*  dual_path;
   MatrixXd*  dual_events;
   MatrixXd*  Xsnopt;
   MatrixXd*  gsnopt;
   MatrixXd*  DerivResid;
   MatrixXd*  h;
   MatrixXd*  Xdotgg;
   MatrixXd*  e;
   MatrixXd*  hgg;
   MatrixXd*  Xip;
   MatrixXd*  JacRow;
   MatrixXd*  Gip;
   MatrixXd*  GFip;
   TripletSparseMatrix*  Ax;
   TripletSparseMatrix*  As;
   TripletSparseMatrix*  Gsp;
   MatrixXd*  control_scaling;
   MatrixXd*  control_shift;
   MatrixXd*  state_scaling;
   MatrixXd*  parameter_scaling;
   MatrixXd*  state_shift;
   MatrixXd*  deriv_scaling;
   MatrixXd*  path_scaling;
   MatrixXd*  event_scaling;
   MatrixXd*  constraint_scaling;
   MatrixXd*  linkage;
   MatrixXd*  prev_states;
   MatrixXd*  prev_costates;
   MatrixXd*  prev_controls;
   MatrixXd*  prev_param;
   MatrixXd*  prev_path;
   MatrixXd*  prev_nodes;
   MatrixXd*  prev_t0;
   MatrixXd*  prev_tf;
   MatrixXd*  Xdot;
   MatrixXd*  JacCol1;
   MatrixXd*  JacCol2;
   MatrixXd*  JacCol3;
   MatrixXd*  xp;
   MatrixXd*  emax_history;
   MatrixXd*  order_reduction;
   MatrixXd*  old_relative_errors;
   MatrixXd*  error_scaling_weights;


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
   int*      iGfun1;
   int*      jGvar1;
   unsigned  int*      iGfun2;
   unsigned  int*      jGvar2;
   int       use_constraint_scaling;
   int       F_nnz;
   double*   G2;
   double*   G3;
   double*   G4;
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
  
// A persistent variable for warm starts in SNOPT7
  int nS;  

};



struct xad_str {
    adouble *xad;
    Workspace *workspace;
};

typedef xad_str XAD;





// Function prototypes

void ScalarGradientAD( adouble (*fun)(adouble *, Workspace*), MatrixXd& x, MatrixXd* grad, bool* flag, int itag, Workspace* workspace );

void EfficientlyComputeJacobianNonZeros( void fun(MatrixXd& x, MatrixXd* f, Workspace* ), MatrixXd& x,
                int nf, double *nzvalue, int nnz, int* iArow, int* jAcol, IGroup* igroup, GRWORK* grw, Workspace* workspace );

void compute_jacobian_of_constraints_with_respect_to_variables(MatrixXd& Jc, MatrixXd& X, MatrixXd& XL, MatrixXd& XU, Workspace* workspace);

void compute_jacobian_of_residual_vector_with_respect_to_variables(MatrixXd& Jr, MatrixXd& X, MatrixXd& XL, MatrixXd& XU, Workspace* workspace);

void getIndexGroups( IGroup* igroup, int nrows, int ncols, int nnz, int* iArow, int* jAcol, Workspace* workspace);

void getIndexGroups( IGroup* igroup, int nrows, int ncols, int nnz, int* iArow, int* jAcol);

void deleteIndexGroups(IGroup* igroup, int ncols );

int psopt(Sol& solution, Prob& problem, Alg& algorithm);

void psopt_level2_setup(Prob& problem, Alg& algorithm);

void initialize_solution(Sol& solution, Prob& problem, Alg& algorithm, Workspace* workspace);

void lglnodes(int N, MatrixXd& x, MatrixXd& w, MatrixXd& P, MatrixXd& D, Workspace* workspace);

void cglnodes(int N, MatrixXd& x, MatrixXd& w,  MatrixXd& D, Workspace* workspace);

void copy_decision_variables(Sol& solution, MatrixXd& x, Prob& problem, Alg& algorithm, Workspace* workspace);

double ff( MatrixXd& x );

void gg( MatrixXd& x, MatrixXd* g );

void fg_num(MatrixXd& x, MatrixXd* fg, Workspace* workspace);

void  define_initial_nlp_guess(MatrixXd& x0, MatrixXd& lambda, Sol& solution, Prob& problem, Alg& algorithm, Workspace* workspace);

void determine_scaling_factors_for_variables(Sol& solution, Prob& problem, Alg& algorithm);

void determine_objective_scaling(MatrixXd& X,Sol& solution, Prob& problem, Alg& algorithm, Workspace* workspace );

void determine_constraint_scaling_factors(MatrixXd & X, Sol& solution, Prob& problem, Alg& algorithm, Workspace* workspace);

void  define_nlp_bounds(MatrixXd& xlb,MatrixXd&  xub,Prob& problem, Alg& algorithm, Workspace* workspace);

double convert_to_original_time(double tbar,double t0,double tf);

void resize_solution(Sol& solution, Prob& problem, Alg& algorithm);

void hot_start_nlp_guess(MatrixXd& x0,MatrixXd& lambda, Sol& solution,Prob& problem,Alg& algorithm, MatrixXd* prev_states, MatrixXd* prev_controls, MatrixXd* prev_costates, MatrixXd* prev_path, MatrixXd* prev_nodes, MatrixXd* prev_param, MatrixXd& prev_t0, MatrixXd& prev_tf, Workspace* workspace );

void lagrange_interpolation(MatrixXd& y, MatrixXd& x, MatrixXd& pointx, MatrixXd& pointy);

double smooth_fabs(double x, double eps);

adouble smooth_fabs(adouble x, double eps);

void Jacobian( void fun(MatrixXd& x, MatrixXd* f ), MatrixXd& x, MatrixXd* grad, GRWORK* grw );

void JacobianRow( void fun(MatrixXd& x, MatrixXd* f, Workspace* ), MatrixXd& x, int iRow, int nf,
                  MatrixXd* JacRow, GRWORK* grw, Workspace* workspace );

void JacobianColumn( void fun(MatrixXd& x, MatrixXd* f, Workspace* ), MatrixXd& x, MatrixXd& xlb, MatrixXd& xub, int jCol,
                MatrixXd* JacColumn, GRWORK* grw, Workspace* workspace );

void evaluate_matrix_of_integrated_errors_in_phase(MatrixXd& eta, int iphase, adouble* xad, int nsteps, Workspace* workspace);

void evaluate_solution(Prob& problem,Alg& algorithm,Sol& solution, Workspace* workspace);

void compute_next_mesh_size( Prob& problem, Alg& algorithm, Sol& solution, Workspace* workspace );

int get_max_nodes(Prob& problem,int iphase, Alg* algorithm);

void estimate_order_reduction(Prob& problem,Alg& algorithm,Sol& solution, Workspace* workspace);

void zero_order_reduction(Prob& problem,Alg& algorithm,Sol& solution, Workspace* workspace);

void construct_new_mesh(Prob& problem,Alg& algorithm,Sol& solution, Workspace* workspace);

bool check_for_equidistributed_error(Prob& problem,Alg& algorithm,Sol& solution);

bool use_local_collocation(Alg & algorithm);

bool use_global_collocation(Alg & algorithm);

void bilinear_interpolation(adouble* z, adouble& x, adouble& y, MatrixXd& X, MatrixXd& Y, MatrixXd& Z);

void spline_2d_interpolation(adouble* z, adouble& x, adouble& y, MatrixXd& X, MatrixXd& Y, MatrixXd& Z, Workspace* workspace);

void smooth_bilinear_interpolation(adouble* z, adouble& x, adouble& y, MatrixXd& X, MatrixXd& Y, MatrixXd& Z);

void smooth_linear_interpolation(adouble* y, adouble& x, adouble* Xdata, adouble* Ydata, int n);

void get_times(adouble *t0, adouble *tf, adouble* xad, int iphase, Workspace* workspace);

void get_states(adouble* states, adouble* xad, int iphase, int k, Workspace* workspace);

void get_controls_bar(adouble* controls_bar, adouble* xad, int iphase, int k, Workspace* workspace);

void rr_num(MatrixXd& X, MatrixXd* residual_vector, Workspace* workspace);

void extract_parameter_covariance(MatrixXd& Cp, MatrixXd& C, Workspace* workspace);

void get_scaled_decision_variables_and_bounds(MatrixXd& x, MatrixXd& xlb, MatrixXd& xub, Workspace* workspace);


int NLP_interface(
         Alg& algorithm,
         MatrixXd* x0,
         double (*f)(MatrixXd&, Workspace*),
	 void (*g)(MatrixXd&,MatrixXd*,Workspace*),
	 int m          = 0,
	 int neqcstr    = 0,
	 MatrixXd* xlb   = NULL,
	 MatrixXd* xub   = NULL,
         MatrixXd* lambda= NULL,
         int hotflag    = 0,
         int iprint     = 1,
         Workspace* workspace=NULL,
         void* user_data = NULL  );



void ScalarGradient( double (*fun)(MatrixXd& x, Workspace*), MatrixXd& x,MatrixXd* grad, GRWORK* grw, Workspace* workspace );

void DetectJacobianSparsity(void fun(MatrixXd& x, MatrixXd* f, Workspace* ), MatrixXd& x, int nf,
                           int* nnzA, int* iArow, int* jAcol, double* Aij,
                           int* nnzG, int* jGrow, int* jGcol,
                           GRWORK* grw, Workspace* workspace);
                           
void DetectJacobianSparsityAD(void fun(MatrixXd& x, MatrixXd* f, Workspace* ), MatrixXd& x, int nf,
                           int* nnzA, int* iArow, int* jAcol, double* Aij,
                           int* nnzG, int* jGrow, int* jGcol,
                           GRWORK* grw, Workspace* workspace);

void ComputeJacobianNonZeros( void fun(MatrixXd& x, MatrixXd* f ), MatrixXd& x, int nf, double *nzvalue, int nnz, int* iArow, int* jAcol, GRWORK* grw, Workspace* workspace );

void Jacobian( void fun(MatrixXd& x, MatrixXd* f ), MatrixXd& x,
                MatrixXd* grad, GRWORK* grw );

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

int get_total_number_of_parameters(Prob& problem);

int get_number_of_linkages(Prob& problem);

int get_number_of_phases(Prob& problem);

void lagrange_interpolation_ad(adouble* y, adouble& x, adouble* pointx, adouble* pointy, int npoints, Workspace* workspace);

void linear_interpolation(adouble* y, adouble& x, adouble* pointx, adouble* pointy, int npoints);

void linear_interpolation(MatrixXd& y, double x, MatrixXd& pointx, MatrixXd& pointy, int npoints);

void linear_interpolation(MatrixXd& y, MatrixXd& x, MatrixXd& pointx, MatrixXd& pointy, int npoints);

void linear_interpolation(adouble* y, adouble x, MatrixXd& pointx, MatrixXd& pointy, int npoints);

void smooth_linear_interpolation(adouble* y, adouble& x, MatrixXd& Xdata, MatrixXd& Ydata, int n);

void spline_interpolation(adouble* y, adouble& x, adouble* pointx, adouble* pointy, int npoints, Workspace* workspace);

void spline_interpolation(adouble* y, adouble& x, MatrixXd& Xdata, MatrixXd& Ydata, int n);

void zoh_interpolation(adouble* y, adouble x, MatrixXd& pointx, MatrixXd& pointy, int npoints);


adouble get_initial_time(adouble* xad, int iphase, Workspace* workspace);

adouble get_final_time(adouble* xad, int iphase, Workspace* workspace);

void get_parameters(adouble* parameters, adouble* xad, int iphase, Workspace* workspace);

bool useAutomaticDifferentiation(Alg& algorithm);

void gg_ad( adouble* xad, adouble* gad, Workspace* workspace );

double ff_num(MatrixXd& x, Workspace* workspace);

void gg_num( MatrixXd& x, MatrixXd* g, Workspace* workspace );

void fg_ad( adouble* x, adouble* fg, Workspace* workspace);

void compute_derivatives_trajectory( MatrixXd& Xdot, Prob& problem, Sol& solution,  int i, Workspace* workspace );

adouble integrate( adouble (*integrand)(adouble*,adouble*,adouble*,adouble&,adouble*,int, Workspace* workspace), adouble* xad, int i, Workspace* workspace );

void auto_link(adouble* linkages, int* index, adouble* xad, int iphase_a, int iphase_b, Workspace* workspace);

void auto_link_2(adouble* linkages, int* index, adouble* xad, int iphase_a, int iphase_b, Workspace* workspace);

void plot(const MatrixXd& x, const MatrixXd& y,const string& title,
          const char* xlabel, const char* ylabel, const char* legend=NULL, const char* terminal=NULL, const char* output=NULL);

void plot(const MatrixXd& x1, const MatrixXd& y1, const MatrixXd& x2, const MatrixXd& y2, const string& title,
          const char* xlabel, const char* ylabel, const char* legend=NULL, const char* terminal=NULL, const char* output=NULL);

void plot(const MatrixXd& x1, const MatrixXd& y1, const MatrixXd& x2, const MatrixXd& y2, const MatrixXd& x3, const MatrixXd& y3,
          const string& title, const char* xlabel, const char* ylabel, const char* legend=NULL, const char* terminal=NULL, const char* output=NULL);

void multiplot(const MatrixXd& x, const MatrixXd& y, const string& title, const char* xlabel, const char* ylabel, const char* legend, int nrows=0, int ncols=0, const char* terminal=NULL, const char* output=NULL ) ;

void spplot(const MatrixXd& x1a, const MatrixXd& y1a, const MatrixXd& x2a, const MatrixXd& y2a, const string& title, const char* xlabel, const char* ylabel, const char* legend, const char* terminal=NULL, const char* output=NULL);

void polar(const MatrixXd& theta, const MatrixXd& r, const string& title,
           const char* legend=NULL, const char* terminal=NULL, const char* output=NULL);

void polar(const MatrixXd& theta, const MatrixXd& r, const MatrixXd& theta2, const MatrixXd& r2, const string& title,
            const char* legend=NULL, const char* terminal=NULL, const char* output=NULL);

void polar(const MatrixXd& theta, const MatrixXd& r, const MatrixXd& theta2, const MatrixXd& r2,  const MatrixXd& theta3, const MatrixXd& r3, const string& title,
            const char* legend=NULL, const char* terminal=NULL, const char* output=NULL);

void surf(const MatrixXd& x, const MatrixXd& y, const MatrixXd& z, const string& title, const char* xlabel, const char* ylabel, const char* zlabel, const char* terminal=NULL, const char* output=NULL, const char* view=NULL);

void plot3(const MatrixXd& x, const MatrixXd& y, const MatrixXd& z, const string& title, const char* xlabel, const char* ylabel, const char* zlabel, const char* terminal=NULL, const char* output=NULL, const char* view=NULL);

void psopt_error_message(const char *error_text);

adouble dot(adouble* x, adouble* y, int n);

adouble smooth_heaviside(adouble x, double a);

adouble smooth_sign(adouble x, double a);

void cross(adouble* x, adouble* y, adouble* z);

void validate_user_input(Prob& problem, Alg& algorithm, Workspace* workspace);

void print_psopt_summary(Prob& problem, Alg& algorithm, Sol& solution, Workspace* workspace);

void psopt_main(Sol& solution, Prob& problem, Alg& algorithm, unique_ptr<Workspace>& workspace_up);

void clip_vector_given_bounds(MatrixXd& xp, MatrixXd& xlb, MatrixXd& xub);

void psopt_print(Workspace* workspace, const char* msg);

int auto_link_count(Prob& problem, int nstates);

void multi_segment_setup(Prob& problem, Alg& algorithm, MSdata& msdata);

int auto_link2_count(Prob& problem, int nstates, int ncontrols);

void auto_phase_setup(Prob& problem, int n_final_events,RowVectorXi& nodes);

void  auto_phase_bounds(Prob& problem);

void  auto_phase_guess(Prob& problem,MatrixXd& controls, MatrixXd& states, MatrixXd& param, MatrixXd& time);

void auto_link_multiple(adouble* linkages, adouble* xad,int nphases, Workspace* workspace);

void auto_link2_multiple(adouble* linkages, adouble* xad,int nphases, Workspace* workspace);

void product_ad(const MatrixXd& A, const adouble* x, int nx, adouble* y);

void sum_ad(const adouble* a, const adouble*b, int n, adouble* c);

void subtract_ad(const adouble* a, const adouble*b, int n, adouble* c);

void inverse_ad(adouble* minput, int n, adouble* minv);

void product_ad(adouble* Apr,adouble* Bpr, int na, int ma, int nb, int mb, adouble* ABpr);

bool need_midpoint_controls(Alg& algorithm, Workspace* workspace);

void inverse_ad(const AutoDiffMatrix& minput, AutoDiffMatrix* minv);

void product_ad(const AutoDiffMatrix& A,const AutoDiffMatrix& B, AutoDiffMatrix* AB);

void sum_ad(const AutoDiffMatrix& A,const AutoDiffMatrix& B, AutoDiffMatrix* AB);

void subtract_ad(const AutoDiffMatrix& A,const AutoDiffMatrix& B, AutoDiffMatrix* AB);

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

void resample_trajectory(MatrixXd& Y, MatrixXd& X, MatrixXd& Ydata, MatrixXd& Xdata);

void load_parameter_estimation_data(Prob& problem, int iphase, const char* filename);

bool compute_parameter_statistics(MatrixXd& Qp, MatrixXd& p, MatrixXd& plow, MatrixXd& phigh, MatrixXd& r, Workspace* workspace);

double inverse_twotailed_t_cdf(double A, int ndf);

double nint(double x);

int get_iphase_offset(Prob& problem, int iphase,Workspace* workspace);

adouble ff_ad(adouble* xad, Workspace* workspace);


void rk4_propagate( void (*dae)(adouble* derivatives, adouble* path, adouble* states,
         adouble* controls, adouble* parameters, adouble& time,
        adouble* xad, int iphase, Workspace* workspace),
        MatrixXd& control_trajectory,
        MatrixXd& time_vector,
        MatrixXd& initial_state,
	MatrixXd& parameters,
        Prob & problem,
        int iphase,
        MatrixXd& state_trajectory, Workspace* workspace);

void rkf_propagate( void (*dae)(adouble* derivatives, adouble* path, adouble* states,
         adouble* controls, adouble* parameters, adouble& time,
        adouble* xad, int iphase, Workspace* workspace),
        MatrixXd& control_trajectory,
        MatrixXd& time_vector,
        MatrixXd& initial_state,
	MatrixXd& parameters,
        double tolerance,
        double hmin,
	double hmax,
        Prob & problem,
        int iphase,
        MatrixXd& state_trajectory,
        MatrixXd& new_time_vector,
	MatrixXd& new_control_trajectory, Workspace* workspace);


void auto_split_observations(Prob& problem, MatrixXd& observation_nodes, MatrixXd& observations);

adouble endpoint_cost_for_parameter_estimation(adouble* initial_states, adouble* final_states, adouble* parameters,adouble& t0, adouble& tf, adouble* xad, int iphase, Workspace* workspace);


//extern "C" {
//   int dgeqrf_(integer *m, integer *n, doublereal *a, integer *
//	lda, doublereal *tau, doublereal *work, integer *lwork, integer *info);
//int dormqr_(char *side, char *trans, integer *m, integer *n,
//        integer *k, doublereal *a, integer *lda, doublereal *tau, doublereal *
//        c__, integer *ldc, doublereal *work, integer *lwork, integer *info,
//        ftnlen side_len, ftnlen trans_len);
//
//}



//! ErrorHandler class
/**
   This is a C++ class intended to handle error conditions.
*/
class ErrorHandler
{
    public:
        //! A string of characters which contains the error message
        string error_message;
        //! A constructor which takes the error message as an argument and assigns it to error_message
        /**
          \param m is the error message string.
          \sa function error_message().
        */
        ErrorHandler(const string m);
};


#ifdef WIN32
#define DEC_THREAD __declspec( thread )
#else
#define DEC_THREAD __thread
#endif


class PSOPT_extras {

//! Flag to indicate error condition
   static DEC_THREAD int      errorFlag;
//! Print level flag,  1: output sent to sderr, 0: no output sent
   static DEC_THREAD int      print_level;
//! variable to store start time after tic() call.
   static DEC_THREAD time_t  start_time;
 //! clock_t variable
   static DEC_THREAD clock_t start_clock;


public:

  //! Sets the print level
  /**
      \param  plevel desired print level
      \return void
  */
  
   static void     tic(void);

   static double   toc(void);

   static clock_t GetStartTicks(void) { return start_clock; }

   static void SetStartTicks(clock_t st) { start_clock=st; }

   static double GetEPS() { return std::numeric_limits<double>::epsilon();  }

   static void SetPrintLevel( int plevel );

   static int PrintLevel();
   
   static void   RiseErrorFlag()   { errorFlag = true; }


   friend void error_message(const char *input_text);
  //! This function, which is to be used in conjunction with function toc(), starts counting elapsed CPU time.
   friend void tic(void);
  //! This function, which is to be used in conjunction with function tic(), stops counting CPU time, and it prints and returns the elapsed time in seconds since the function tic() was called.
  /**
      \return the elapsed time in seconds.
  */
   friend double toc();



};

void error_message(const char *input_text);
void tic(void);
double toc();

void sort_vector(MatrixXd& A, RowVectorXi& sindex);
void sort(MatrixXd& m);
void rearrange_vector(MatrixXd& A, RowVectorXi& sindex);
void Save(const MatrixXd& m, const char* filename);
void Print(const MatrixXd& m, const char* text);
MatrixXd linspace(double X1, double X2, long N);
MatrixXd reshape(MatrixXd& A, int n, int m);
MatrixXd zeros(long nrows, long ncols);
MatrixXd eye(long nrows);
MatrixXd ones(long nrows, long ncols);
MatrixXd elemProduct(const MatrixXd& m1,const MatrixXd& m2);
MatrixXd elemDivision(const MatrixXd& m1,const MatrixXd& m2);
MatrixXd load_data(const char* filename, long nrows, long ncols);
MatrixXd Abs(const MatrixXd& m);
MatrixXd tra(const MatrixXd& m);
MatrixXd stack_columns(const MatrixXd& m);
MatrixXd sum_columns(const MatrixXd& A);
MatrixXd GaussianRandom(long nrows, long ncols);
MatrixXd sin(const MatrixXd& m);
MatrixXd cos(const MatrixXd& m);
MatrixXd tan(const MatrixXd& m);
MatrixXd asin(const MatrixXd& m);
MatrixXd acos(const MatrixXd& m);
MatrixXd atan(const MatrixXd& m);
MatrixXd sinh(const MatrixXd& m);
MatrixXd cosh(const MatrixXd& m);
MatrixXd tanh(const MatrixXd& m);
MatrixXd exp(const MatrixXd& m);
MatrixXd log(const MatrixXd& m);
MatrixXd log10(const MatrixXd& m);
MatrixXd sqrt(const MatrixXd& m);
long length(const MatrixXd& m);
double Max(const MatrixXd& m);
double Max(const MatrixXd& m, long* i);
double Min(const MatrixXd& m);
double Min(const MatrixXd& m, long* i);
double mean(const MatrixXd& m);
double MaxAbs(const MatrixXd& m);
double sum(MatrixXd& A);
bool any(const MatrixXd& m);


bool isEmpty(const MatrixXd& m);     
bool isSymmetric(const MatrixXd& m); 



#ifdef PACKAGE
#undef PACKAGE
#undef PACKAGE_VERSION
#undef PACKAGE_BUGREPORT
#undef PACKAGE_NAME
#undef PACKAGE_STRING
#undef PACKAGE_TARNAME
#undef VERSION
#endif


#include <IpIpoptApplication.hpp>
#include "ipopt_psopt.h"



using namespace Ipopt;



//! TripletSparseMatrix class
/**
*/
class TripletSparseMatrix {
protected:
//! Array to store matrix elements
   double *a;
//! Number of matrix rows
   int n;
//! Number of matrix columns
   int m;
//! Number of allocated elements in a
   int asize;
//! Number of non-zero elements
   int nz;
//! Array of nonzero row indices
   int *RowIndx;
//! Array of nonzero column indices
   int *ColIndx;

public:

   // Interface functions
  //! Gets the number of rows from the calling object
  /**
      \return integer value with the number of rows.
  */
   int rows() const  { return n; }
  //! Gets the number of columns from the calling object
  /**
      \return integer value with the number of columns.
  */
   int cols() const { return m; }
  //! Gets the number of non-zero elements from the calling object
  /**
      \return integer value with the number of non-zero elements.
  */
   int GetNonZero() const {return nz;}
  //! Returns a int pointer to the start of an array of row indices for the non-zero elements
  /**
      \return int pointer to the start of the array
  */
   int* GetRowIndx_C_Array() const {return RowIndx;}
  //! Returns a int pointer to the start of an array of column indices for the non-zero elements
  /**
      \return int pointer to the start of the array
  */
   int* GetColIndx_C_Array() const {return ColIndx;}
  //! Returns a double pointer to the start of an array where the non-zero elements are stored
  /**
      \return double pointer to the start of the array
  */
   double *GetPr()  const  { return a; }
  //! Inserts a non-zero value at a specified location of the sparse matrix. The function re-allocates storage as required.
  /**
      \param i row index (starting from 1)
      \param j column index (starting from 1)
      \param val non-zero double value to be inserted.
      \return void
  */
   void InsertNonZero(int i, int j, double val);
  //! Saves a sparse matrix in triplet form. The first row of the saved file contains the number of rows, number of columns and the number of non-zeros of the matrix. Each subsequent row contains the row index, the column index and the corresponding nonzero value.
  /**
      \param fname name of the file to be created
      \return void
  */
   void Save(const char* fname) const;
  //! Saves a sparse sparsity pattern for the sparse matrix, such that the saved matrix contains only zeros and asterisks.
  /**
      \param fname name of the file to be created
      \return void
  */
   void SaveSparsityPattern(const char* fname) const;
  //! Loads a sparse matrix in triplet form. The first row of the specified file must contain the number of rows, number of columns and the number of non-zeros of the matrix. Each subsequent row must contain the row index, the column index and the corresponding nonzero value, separated by spaces. An error is thrown is the file cannot be opened.
  /**
      \param fname name of the file to be read
      \return void
  */
   void Load(const char* fname);
  //! This function transposes a sparse matrix and returns the transposed sparse matrix.
  /**
      \param  A  sparse matrix to be transposed
      \return a temporary TripletSparseMatrix object with the result of the operation
  */
   friend TripletSparseMatrix tra(const TripletSparseMatrix& A);

  //! Extracts a specified column from a sparse matrix and returns a MatrixXd object.
  /**
      \param  j: column index
      \return Temporary MatrixXd object with the specified column
  */
   MatrixXd col(int j) const;
  //! Extracts a specified row from a sparse matrix and returns a MatrixXd object.
  /**
      \param  i: row index
      \return Temporary MatrixXd object with the specified row
  */
   MatrixXd row(int i) const;

   // Display functions
  //! Prints a sparse matrix in triplet format.
  /**
      \param  text: label to identify the sparse matrix to be printed
      \return void
  */
   void Print(const char* text);

   // Constructors
  //! Default constructor. Creates an empty sparse matrix with zero rows and columns.
   TripletSparseMatrix(void); // Default constructor
  //! Constructor using triplet arrays.
  /**
      \param  aa: array with double non-zero elements
      \param  nn: number of rows
      \param  mm: number of columns
      \param  nnz: number of non-zero elements
      \param  RowIndxArg: array of int values with the row indices of the nonzero elements (indices start from 1).
      \param  ColIndxArg: array of int values with the column indices of the nonzero elements (indices start from 1).
  */
   TripletSparseMatrix(double* aa, int nn, int mm, int nnz, int* RowIndxArg, int* ColIndxArg); // Constructor using arrays


  //! Copy constructor. Creates a copy of the argument.
  /**
        \param  A: TripletSparseMatrix object
  */
   TripletSparseMatrix( const TripletSparseMatrix& A); // copy constructor
  //! Constructor without element asignment. Creates a TripletSparseMatrix object with storage for a specified number of non-zero values.
  /**
        \param  nn: number of rows.
        \param  mm: number of columns.
        \param  nnz: number of non-zero values.
  */
   TripletSparseMatrix( int nn, int mm, int nnz); // Constructor without element assignment

   // Destructor
  //! Desctructor.
   ~TripletSparseMatrix();

   // Other functions
  //! Resizes a TripletSparseMatrix object
  /**
      \param  nnew: new number of rows.
      \param  mnew: new number of columns.
      \param  nznew: new number of non-zero values.
      \return void
  */
   void resize(int nnew, int mnew, int nznew);


  //! Transposes a TripletSparseMatrix object.
  /**
      \return void
  */
   void Transpose();
  
   //!  This function computes and returns the element-wise product of two sparse matrices of the same dimensions. If the dimensions of the two input matrices are not the same, an error is thrown.
  /**
      \param  A is a TripletSparseMatrix object.
      \param B is a TripletSparseMatrix object.
      \return a temporary TripletSparseMatrix object with the result of the operation
  */
   friend TripletSparseMatrix elemProduct(const TripletSparseMatrix A, const TripletSparseMatrix& B);
  
  //! Creates a sparse identity matrix of specified dimension.
  /**
      \param  n: number of rows and columns of the sparse matrix to be created.
      \return A temporary TripletSparseMatrix object with the result of the operation
  */
   friend TripletSparseMatrix speye(int n);

  
   friend double enorm(const TripletSparseMatrix& A);


   //!  This function computes and returns the element-wise absolute value of a sparse matrix A.
  /**
      \param  A is a TripletSparseMatrix object.
      \return a reference to a temporary TripletSparseMatrix object with the result of the operation.
  */
   friend TripletSparseMatrix Abs(const TripletSparseMatrix& A);
 

   //! This function eliminates zero elements from a sparse matrix and deletes unnecessary storage.
  /**
      \return void
  */
   void Compress();

   // Operators
  //! Sparse matrix addition and substitution operator. The sizes of the matrices being added must be the same, otherwise an error is thrown.The left hand side object is replaced with the result of the operation.
  /**
      \param rval:  TripletSparseMatrix object located right hand side of the operator.
      \return Reference to the calling TripletSparseMatrix object
  */
   TripletSparseMatrix& operator += (const TripletSparseMatrix &rval);
  //! Sparse matrix subtraction and substitution operator. The sizes of the matrices being subtracted must be the same, otherwise an error is thrown. The left hand side object is replaced with the result of the operation.
  /**
      \param rval:  TripletSparseMatrix object located right hand side of the operator.
      \return Reference to the calling TripletSparseMatrix object
  */
   TripletSparseMatrix& operator -= (const TripletSparseMatrix &rval);
  
  //! Computes the product of a sparse matrix (left hand side of the operator) times a real scalar (right hand side value) and replaces the left hand side object with the result of the operation.
  /**
      \param Arg: double value that will multiply each non-zero element of the sparse matrix.
      \return Reference the calling TripletSparseMatrix object
  */
   TripletSparseMatrix& operator*= (double Arg);
  //! Computes the division of a sparse matrix (left hand side of the operator) by a real scalar (right hand side value) and replaces the left hand side object with the result of the operation.
  /**
      \param Arg: double value that will divide each non-zero element of the sparse matrix.
      \return Reference the calling TripletSparseMatrix object
  */
   TripletSparseMatrix& operator/= (double Arg);
  //! Sparse matrix addition operator. The row and column sizes of the matrices being added must be the same, otherwise an error is thrown.
  /**
      \param rval:  sparse matrix located at the right hand side of the operator.
      \return Reference to a temporary TripletSparseMatrix object with the result of the operation
  */
   TripletSparseMatrix operator+ (const TripletSparseMatrix& rval) const;
  //! Sparse matrix subtraction operator. The row and column sizes of the matrices being subtracted must be the same, otherwise an error is thrown.
  /**
      \param rval:  sparse matrix located at the right hand side of the operator.
      \return Reference to a temporary TripletSparseMatrix object with the result of the operation
  */
   TripletSparseMatrix operator - (const TripletSparseMatrix& rval) const;
  
  //! Computes the product of a sparse matrix (left hand side of the operator) times a real scalar (right hand side value).
  /**
      \param Arg: double value that will multiply each non-zero element of the sparse matrix.
      \return Reference to a temporary object with the result of the operation.
  */
   TripletSparseMatrix operator* (double Arg) const;
//! Computes the product of a triplet sparse matrix (left hand side of the operator) by a MatrixXd object (right hand side of the operator).
  /**
      \param A:  MatrixXd object to be multiplied
      \return Reference to a temporary object with the result of the operation.   
  */    
   TripletSparseMatrix operator* (const MatrixXd& A) const;
  //! Computes the product of a real value (left hand side of the operator) by a sparse matrix (right hand side of the operator).
  /**
      \param Arg: double value that will multiply each non-zero element of the sparse matrix.
      \param A:  TripletSparseMatrix object to be multiplied by a real value.
      \return Reference to a temporary object with the result of the operation.
  */
   friend TripletSparseMatrix operator *(double Arg, const TripletSparseMatrix& A);
  //! Computes the division of a sparse matrix (left hand side of the operator) by a real scalar (right hand side value).
  /**
      \param Arg: double value that will divide each non-zero element of the sparse matrix.
      \return Reference to a temporary object with the result of the operation.
  */
   TripletSparseMatrix operator/ (double Arg) const;

  //! Sparse matrix assignment. The size of the left hand side object is modified if necessary, and the values of all non-zero real elements of the right hand side object and copied to the left hand side object.
  /**
      \param rval: TripletSparseMatrix object at the right hand side of the operator
      \return Reference to the calling TripletSparseMatrix object
  */
   TripletSparseMatrix& operator= (const TripletSparseMatrix& rval);

  //! Matrix indexing. Returns a reference to the matrix element located at the position indicated by the row and column indices. Indices start from 1. If the row and column indices refer to a zero (non-allocated) element, a new element with zero value is allocated so that it can be returned as a reference.
  /**
      \param row: Row index starting from 1.
      \param col: Column index starting from 1.
      \return Reference to the indexed matrix element.
  */
   double& operator() (int row,int col);
  //! Matrix indexing. Returns the value of the matrix element located at the position indicated by the row and column indices. Indices start from 1.
  /**
      \param row: Row index starting from 1.
      \param col: Column index starting from 1.
      \return double value of the element at the indicated position.
  */
   double operator() (int row,int col) const;
  //! Computes the left division of a sparse matrix (left hand side of the operator) by another sparse matrix (right hand side value). This is conceptually equivalent to multiplying the inverse of the left object by the right hand side object but it is computed in a more efficient way. The dimensions of the matrices must be consistent, otherwise an error is returned. The left hand side object must be a square matrix.
  /**
      \param B: TripletSparseMatrix object at the right hand side of the operator.
      \return Reference to a temporary TripletSparseMatrix object with the result of the operation
   */

   friend void sp_error_message(const char *error_text);


};

// Declaration of all friend functions of TripletSparseMatrix class


void sp_error_message(const char *error_text);
TripletSparseMatrix elemProduct(const TripletSparseMatrix A, const TripletSparseMatrix& B);
TripletSparseMatrix tra(const TripletSparseMatrix& A);
TripletSparseMatrix speye(int n);
double enorm(const TripletSparseMatrix& A);
TripletSparseMatrix Abs(const TripletSparseMatrix& A);
TripletSparseMatrix operator *(double Arg, const TripletSparseMatrix& A);

#endif
