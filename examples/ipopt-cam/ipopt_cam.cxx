//////////////////////////////////////////////////////////////////////////
////////////////            ipopt_cam.cxx            /////////////////////
//////////////////////////////////////////////////////////////////////////
//////// Title:  Shape optimization of a cam                 /////////////
//////// Reference: Dolan, More and Munson, "Benchmarking    /////////////
////////            Optimization Software with COPS 3.0",   /////////////
////////            ANL/MCS-TM-273, 2004                     /////////////
//////////////////////////////////////////////////////////////////////////
////////     Copyright (c) Victor M. Becerra, 2020-2026       ////////////
//////////////////////////////////////////////////////////////////////////
//////// This is part of the PSOPT software library, which    ////////////
//////// is distributed under the terms of the GNU Lesser     ////////////
//////// General Public License (LGPL)                        ////////////
//////////////////////////////////////////////////////////////////////////
//
// This is the worked example of Chapter 2 of "Computational Optimal Control
// and Estimation". It does NOT use PSOPT: it calls IPOPT directly, through
// IPOPT's TNLP interface, and obtains sparse first and second derivatives by
// automatic differentiation with CppAD. It is shipped with PSOPT so that a
// reader of the book can build and run it, and it is built by the same CMake
// configuration, but nothing in PSOPT depends on it.
//
// The problem. The shape of a cam is circular over an angle of 6*pi/5 of its
// circumference, with radius r_min. The design variables r_i, i = 1..n are the
// radius at n equally spaced angles over an angle of 2*pi/5; the shape is
// symmetric, so only half of the design angle of 4*pi/5 is represented. The
// objective is the area of the valve opening over one rotation, which is
// maximised. The convexity constraints are nonlinear; the objective and the
// curvature constraints are linear. The result is n variables and 2n+3
// inequality constraints.
//
// On automatic differentiation. An earlier version of this program used ADOL-C
// with graph colouring by ColPack. PSOPT now uses CppAD throughout, and this
// example follows, so that a reader needs one AD tool rather than two. The
// structure below is the idiomatic CppAD one and is worth reading for its own
// sake: a SINGLE tape is recorded of the combined function
//
//     fg(x) = [ f(x) ; g(x) ]   from R^n to R^(1+m),
//
// and everything IPOPT asks for is read off that one tape. The gradient of the
// objective is one reverse sweep with weight e_0. The Jacobian of the
// constraints is the rows 1..m of the Jacobian of fg. The Hessian of the
// Lagrangian, sigma*Hess f + sum_i lambda_i * Hess g_i, is the weighted Hessian
// of fg with weight vector [sigma; lambda] -- which is why taping f and g
// together, rather than separately, is the natural thing to do.
//
// One wrinkle, recorded because it will bite anyone porting their own problem.
// A structural sparsity pattern is allowed to be a superset of the true one:
// it is a statement about which second derivatives *might* be nonzero. In
// CppAD 2024 the reverse Hessian sparsity of the unary negation operator is not
// tight -- it treats -u as though it had a nonzero second derivative in u -- so
// writing a bilinear constraint as
//
//     g = -r[i-1]*r[i] - r[i]*r[i+1] + 2*c*r[i-1]*r[i+1];
//
// leads the leading term to be taped as (-r[i-1])*r[i] and puts a spurious
// (i-1,i-1) on the diagonal. Over the whole problem that is 999 structurally
// zero entries out of 2997, half again as many as are really there. Ordering
// each row so that it begins with a positive term, as below, costs nothing,
// changes no value, and recovers the exact 1998.
//
//////////////////////////////////////////////////////////////////////////

#include <cppad/cppad.hpp>
#include <IpTNLP.hpp>
#include <IpIpoptApplication.hpp>

#include <cassert>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <string>
#include <vector>

using namespace Ipopt;

// ---------------------------------------------------------------------------
// Problem size and data.
// ---------------------------------------------------------------------------
// The book reports this example for n = 1000. Another size can be had from
// the command line: "ipopt_cam 500". The default is 1000.
static int          NP    = 1000;      // number of design radii
static const double alpha = 1.5;       // curvature limit
static const double R_min = 1.0;
static const double R_max = 2.0;
static const double R_v   = 1.0;       // valve radius
static const double Pi    = 3.14159265358979323846;

typedef CppAD::AD<double>       adouble;
typedef CppAD::vector<double>   DVector;
typedef CppAD::vector<adouble>  ADVector;
typedef CppAD::vector<size_t>   SizeVector;
typedef CppAD::vector<bool>     BoolVector;
typedef CppAD::sparse_rc<SizeVector>            SparsePattern;
typedef CppAD::sparse_rcv<SizeVector, DVector>  SparseSubset;

// ---------------------------------------------------------------------------
// The problem itself. Templated on the scalar type so that one piece of source
// serves both the double evaluation IPOPT asks for and the AD<double> recording
// the tape is made from. This is the only part of the file a reader wanting to
// solve a different problem needs to change.
// ---------------------------------------------------------------------------
template <class T>
static T cam_objective(const T* r, int n)
{
    // Area of the valve opening, negated because IPOPT minimises.
    T obj = T(0.0);
    for (int i = 0; i < n; i++) obj -= r[i];
    obj *= (Pi * R_v) / double(n);
    return obj;
}

template <class T>
static void cam_constraints(const T* r, int n, T* g)
{
    const double d_theta = 2.0 * Pi / (5.0 * (n + 1));
    const double c_dth   = std::cos(d_theta);

    int j = 0;
    // Convexity, on the interior radii. These are the nonlinear ones. Each row
    // leads with its positive term, for the sparsity reason set out in the
    // header comment; mathematically the ordering is immaterial.
    for (int i = 1; i < n - 1; i++) {
        g[j++] = 2.0*c_dth*r[i-1]*r[i+1] - r[i-1]*r[i] - r[i]*r[i+1];
    }
    // Convexity at the four edges of the design arc.
    g[j++] = 2.0*R_min*c_dth*r[1]   - R_min*r[0]    - r[0]*r[1];
    g[j++] = 2.0*R_min*c_dth*r[0]   - R_min*R_min   - R_min*r[0];
    g[j++] = 2.0*R_max*c_dth*r[n-2] - r[n-2]*r[n-1] - R_max*r[n-1];
    g[j++] = 2.0*c_dth*r[n-1]*r[n-1] - 2.0*R_max*r[n-1];

    // Curvature. Linear.
    for (int i = 0; i < n - 1; i++) g[j++] = r[i+1] - r[i];
    g[j++] = r[0]   - R_min;
    g[j++] = R_max  - r[n-1];
}

static void cam_bounds(int n, int m, double* x_l, double* x_u,
                       double* g_l, double* g_u)
{
    const double d_theta = 2.0 * Pi / (5.0 * (n + 1));
    const double INFTY   = 2e19;

    for (int i = 0; i < n; i++) { x_l[i] = R_min; x_u[i] = R_max; }

    int j = 0;
    // Convexity: (n-2) interior rows and 4 edge rows, each one-sided.
    for (int i = 0; i < (n - 2) + 4; i++) { g_l[j] = -INFTY; g_u[j] = 0.0; j++; }
    // Curvature: (n-1) rows and 2 more, each two-sided.
    for (int i = 0; i < (n - 1) + 2; i++) {
        g_l[j] = -alpha*d_theta; g_u[j] = alpha*d_theta; j++;
    }
    assert(j == m);
    (void) m;
}

static void cam_starting_point(int n, double* r)
{
    for (int i = 0; i < n; i++) r[i] = 0.5 * (R_min + R_max);
}

// ---------------------------------------------------------------------------
// Writes a sparse matrix in the three-column form the MATLAB script plotspy.m
// reads: a header line of "rows cols nnz", then one "row col value" per entry
// with zero-based indices.
// ---------------------------------------------------------------------------
static void save_sparse(const char* filename, int nrow, int ncol,
                        const SizeVector& row, const SizeVector& col,
                        const double* val, size_t nnz)
{
    FILE* f = fopen(filename, "w");
    if (!f) return;
    fprintf(f, "%d\t%d\t%zu\n", nrow, ncol, nnz);
    for (size_t k = 0; k < nnz; k++)
        fprintf(f, "%zu\t%zu\t%e\n", row[k], col[k], val[k]);
    fclose(f);
}

// ---------------------------------------------------------------------------
// The IPOPT problem class. Everything below is generic: it would serve any
// problem whose objective and constraints are written as the two templates
// above.
// ---------------------------------------------------------------------------
class CamNLP : public TNLP
{
public:
    CamNLP() : n_(NP), m_(2*NP + 3) {}
    virtual ~CamNLP() {}

    bool get_nlp_info(Index& n, Index& m, Index& nnz_jac_g,
                      Index& nnz_h_lag, IndexStyleEnum& index_style) override
    {
        record_tape();
        n = n_;  m = m_;
        nnz_jac_g = (Index) jac_subset_.nnz();
        nnz_h_lag = (Index) hes_subset_.nnz();
        index_style = C_STYLE;
        return true;
    }

    bool get_bounds_info(Index n, Number* x_l, Number* x_u,
                         Index m, Number* g_l, Number* g_u) override
    {
        cam_bounds(n, m, x_l, x_u, g_l, g_u);
        return true;
    }

    bool get_starting_point(Index n, bool init_x, Number* x,
                            bool init_z, Number*, Number*,
                            Index, bool init_lambda, Number*) override
    {
        assert(init_x && !init_z && !init_lambda);
        (void) init_x; (void) init_z; (void) init_lambda;
        cam_starting_point(n, x);
        return true;
    }

    bool eval_f(Index n, const Number* x, bool, Number& obj_value) override
    {
        obj_value = cam_objective<double>(x, n);
        return true;
    }

    bool eval_g(Index n, const Number* x, bool, Index, Number* g) override
    {
        cam_constraints<double>(x, n, g);
        return true;
    }

    // One reverse sweep on the combined tape, with weight e_0, is the gradient
    // of the objective and nothing else.
    bool eval_grad_f(Index n, const Number* x, bool, Number* grad_f) override
    {
        DVector xv(n);
        for (Index i = 0; i < n; i++) xv[i] = x[i];
        fg_.Forward(0, xv);

        DVector w(1 + m_);
        for (int i = 0; i < 1 + m_; i++) w[i] = 0.0;
        w[0] = 1.0;

        DVector dw = fg_.Reverse(1, w);
        for (Index i = 0; i < n; i++) grad_f[i] = dw[i];
        return true;
    }

    bool eval_jac_g(Index n, const Number* x, bool, Index m, Index,
                    Index* iRow, Index* jCol, Number* values) override
    {
        if (values == NULL) {
            const SizeVector& row = jac_subset_.row();
            const SizeVector& col = jac_subset_.col();
            for (size_t k = 0; k < jac_subset_.nnz(); k++) {
                iRow[k] = (Index)(row[k] - 1);   // row 0 of fg is the objective
                jCol[k] = (Index) col[k];
            }
            return true;
        }
        DVector xv(n);
        for (Index i = 0; i < n; i++) xv[i] = x[i];
        fg_.sparse_jac_rev(xv, jac_subset_, jac_pattern_, "cppad", jac_work_);

        const DVector& v = jac_subset_.val();
        for (size_t k = 0; k < jac_subset_.nnz(); k++) values[k] = v[k];

        static bool saved = false;
        if (!saved) {
            SizeVector r(jac_subset_.nnz());
            for (size_t k = 0; k < jac_subset_.nnz(); k++)
                r[k] = jac_subset_.row()[k] - 1;
            save_sparse("jacobian.dat", m, n, r, jac_subset_.col(),
                        &v[0], jac_subset_.nnz());
            saved = true;
        }
        return true;
    }

    // The Hessian of the Lagrangian is the Hessian of fg weighted by
    // [obj_factor; lambda]. IPOPT wants the lower-left triangle, which is what
    // the subset was built to hold.
    bool eval_h(Index n, const Number* x, bool, Number obj_factor,
                Index m, const Number* lambda, bool, Index,
                Index* iRow, Index* jCol, Number* values) override
    {
        if (values == NULL) {
            const SizeVector& row = hes_subset_.row();
            const SizeVector& col = hes_subset_.col();
            for (size_t k = 0; k < hes_subset_.nnz(); k++) {
                iRow[k] = (Index) row[k];
                jCol[k] = (Index) col[k];
            }
            return true;
        }
        DVector xv(n), w(1 + m);
        for (Index i = 0; i < n; i++) xv[i] = x[i];
        w[0] = obj_factor;
        for (Index i = 0; i < m; i++) w[1 + i] = lambda[i];

        fg_.sparse_hes(xv, w, hes_subset_, hes_pattern_, "cppad.symmetric", hes_work_);

        const DVector& v = hes_subset_.val();
        for (size_t k = 0; k < hes_subset_.nnz(); k++) values[k] = v[k];

        // Saved for the spy plot with the indices transposed, so that the file
        // holds the UPPER triangle. IPOPT wants the lower one and that is what
        // the subset holds; which triangle is plotted is a presentation choice,
        // and this is the one the figures in the book were drawn from.
        static bool saved = false;
        if (!saved) {
            save_sparse("Hessian.dat", n, n, hes_subset_.col(), hes_subset_.row(),
                        &v[0], hes_subset_.nnz());
            saved = true;
        }
        return true;
    }

    void finalize_solution(SolverReturn, Index n, const Number* x,
                           const Number*, const Number*, Index m,
                           const Number* g, const Number*, Number obj_value,
                           const IpoptData*, IpoptCalculatedQuantities*) override
    {
        printf("\n\nObjective value\n");
        printf("f(x*) = %e\n", obj_value);
        printf("Non-zeros: Jacobian %zu, Hessian of the Lagrangian %zu\n",
               jac_subset_.nnz(), hes_subset_.nnz());

        FILE* f = fopen("ipopt.m", "w");
        if (f) {
            fprintf(f, "\n x = [");
            for (Index i = 0; i < n; i++) fprintf(f, "\n %lf", x[i]);
            fprintf(f, "\n];");
            fprintf(f, "\n\n g = [");
            for (Index i = 0; i < m; i++) fprintf(f, "\n %lf", g[i]);
            fprintf(f, "\n];");
            fclose(f);
        }
    }

private:
    CamNLP(const CamNLP&);
    CamNLP& operator=(const CamNLP&);

    // Record fg(x) = [f(x); g(x)] once, then derive from it the sparsity of the
    // constraint Jacobian and of the Hessian of the Lagrangian.
    void record_tape()
    {
        DVector x0(n_);
        cam_starting_point(n_, &x0[0]);

        ADVector ax(n_);
        for (int i = 0; i < n_; i++) ax[i] = x0[i];

        CppAD::Independent(ax);
        ADVector afg(1 + m_);
        afg[0] = cam_objective<adouble>(&ax[0], n_);
        cam_constraints<adouble>(&ax[0], n_, &afg[1]);
        fg_.Dependent(ax, afg);
        fg_.optimize();

        // Jacobian sparsity of the whole of fg, by a forward sweep from the
        // identity pattern on the independent variables.
        SparsePattern eye(n_, n_, n_);
        for (int i = 0; i < n_; i++) eye.set(i, i, i);
        SparsePattern jac_all;
        fg_.for_jac_sparsity(eye, false, false, false, jac_all);

        // Keep the constraint rows, 1..m. Row 0 is the objective, whose
        // gradient is obtained by a reverse sweep instead.
        {
            size_t keep = 0;
            for (size_t k = 0; k < jac_all.nnz(); k++)
                if (jac_all.row()[k] > 0) keep++;
            SparsePattern p(1 + m_, n_, keep);
            size_t j = 0;
            for (size_t k = 0; k < jac_all.nnz(); k++)
                if (jac_all.row()[k] > 0)
                    p.set(j++, jac_all.row()[k], jac_all.col()[k]);
            jac_pattern_ = p;
            jac_subset_  = SparseSubset(p);
        }

        // Hessian sparsity of the Lagrangian: every component of fg may carry a
        // weight, so all of them are selected.
        BoolVector select_range(1 + m_);
        for (int i = 0; i < 1 + m_; i++) select_range[i] = true;
        SparsePattern hes_all;
        fg_.rev_hes_sparsity(select_range, false, false, hes_all);

        // sparse_hes is given the FULL symmetric pattern -- its symmetric
        // colouring is computed from it and is wrong if handed one triangle --
        // while the subset it fills is the lower-left triangle IPOPT asks for.
        hes_pattern_ = hes_all;
        {
            size_t keep = 0;
            for (size_t k = 0; k < hes_all.nnz(); k++)
                if (hes_all.row()[k] >= hes_all.col()[k]) keep++;
            SparsePattern p(n_, n_, keep);
            size_t j = 0;
            for (size_t k = 0; k < hes_all.nnz(); k++)
                if (hes_all.row()[k] >= hes_all.col()[k])
                    p.set(j++, hes_all.row()[k], hes_all.col()[k]);
            hes_subset_ = SparseSubset(p);
        }
    }

    int n_, m_;
    CppAD::ADFun<double>   fg_;
    SparsePattern          jac_pattern_, hes_pattern_;
    SparseSubset           jac_subset_,  hes_subset_;
    CppAD::sparse_jac_work jac_work_;
    CppAD::sparse_hes_work hes_work_;
};

// ---------------------------------------------------------------------------
int main(int argc, char* argv[])
{
    if (argc > 1) {
        int np = atoi(argv[1]);
        if (np < 4) {
            printf("usage: %s [number of design radii >= 4]\n", argv[0]);
            return 1;
        }
        NP = np;
    }
    printf("Cam design with %d radii.\n", NP);

    SmartPtr<TNLP> nlp = new CamNLP();
    SmartPtr<IpoptApplication> app = IpoptApplicationFactory();

    app->Options()->SetNumericValue("tol", 1.0e-9);
    app->Options()->SetStringValue("output_file", "ipopt.out");
    app->Options()->SetStringValue("nlp_scaling_method", "gradient-based");
    app->Options()->SetStringValue("hessian_approximation", "exact");
    app->Options()->SetIntegerValue("print_level", 5);

    ApplicationReturnStatus status = app->Initialize();
    if (status != Solve_Succeeded) {
        printf("\n\n*** Error during initialization!\n");
        return (int) status;
    }

    status = app->OptimizeTNLP(nlp);

    if (status == Solve_Succeeded)
        printf("\n\n*** The problem has been solved!\n");
    else
        printf("\n\n*** The problem FAILED!\n");

    return (int) status;
}
