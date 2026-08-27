//////////////////////////////////////////////////////////////////////////
////////////////          ipopt_polygon.cxx          /////////////////////
//////////////////////////////////////////////////////////////////////////
//////// Title:  Largest small polygon                       /////////////
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
// The problem. Find the polygon of maximal area, LSP(nv), among polygons with
// nv sides and diameter at most 1. The vertices are given in polar coordinates
// (r_i, theta_i) about one of them, so there are 2*nv variables. The area is a
// sum of triangle areas and is what is maximised. The constraints are the
// ordering of the angles, the diameter bound on every pair of vertices -- which
// is where the (nv-2)(nv-1)/2 nonlinear rows come from, and why the problem
// grows quadratically -- and two equalities fixing the last vertex.
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
// The book reports this example for nv = 16 and nv = 50. Either can be had from
// the command line: "ipopt_polygon 16". The default is 50.
static int          NV = 50;           // number of vertices
static const double Pi = 3.14159265358979323846;

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
static T poly_objective(const T* x, int nv)
{
    // Area of the polygon, negated because IPOPT minimises.
    const T* r     = &x[0];
    const T* theta = &x[nv];

    T obj = T(0.0);
    for (int i = 0; i < nv - 1; i++)
        obj -= 0.5 * r[i+1] * r[i] * sin(theta[i+1] - theta[i]);
    return obj;
}

template <class T>
static void poly_constraints(const T* x, int nv, T* g)
{
    const T* r     = &x[0];
    const T* theta = &x[nv];

    // The angles must be ordered.
    int k = 0;
    for (int i = 0; i < nv - 2; i++) g[k++] = theta[i] - theta[i+1];

    // Diameter: the squared distance between every pair of vertices, bounded
    // by one. These are the nonlinear rows, and there are (nv-2)(nv-1)/2.
    for (int i = 0; i < nv - 2; i++)
        for (int j = i + 1; j < nv - 1; j++)
            g[k++] = r[i]*r[i] + r[j]*r[j] - 2.0*r[i]*r[j]*cos(theta[i] - theta[j]);

    // Two equalities pinning the last vertex.
    g[k++] = r[nv-1];
    g[k++] = theta[nv-1];
}

static void poly_bounds(int nv, int m, double* x_l, double* x_u,
                        double* g_l, double* g_u)
{
    const double INFTY = 2e19;

    for (int i = 0; i < nv; i++)      { x_l[i] = 0.0; x_u[i] = 1.0; }
    for (int i = nv; i < 2*nv; i++)   { x_l[i] = 0.0; x_u[i] = Pi;  }

    int k = 0;
    for (int i = 0; i < nv - 2; i++) { g_l[k] = -INFTY; g_u[k] = 0.0; k++; }
    for (int i = 0; i < nv - 2; i++)
        for (int j = i + 1; j < nv - 1; j++) { g_l[k] = -INFTY; g_u[k] = 1.0; k++; }
    g_l[k] = g_u[k] = 0.0; k++;      // r[nv-1] = 0
    g_l[k] = g_u[k] = Pi;  k++;      // theta[nv-1] = pi
    assert(k == m);
    (void) m;
}

static void poly_starting_point(int nv, double* x)
{
    double* r     = &x[0];
    double* theta = &x[nv];
    for (int i = 0; i < nv; i++) {
        r[i]     = sin(i * Pi / (nv - 1));
        theta[i] = i * Pi / (nv - 1);
    }
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
class PolygonNLP : public TNLP
{
public:
    PolygonNLP() : n_(2*NV), m_((NV-2) + ((NV-2)*(NV-1)/2 + 2)) {}
    virtual ~PolygonNLP() {}

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
        poly_bounds(NV, m, x_l, x_u, g_l, g_u);
        (void) n;
        return true;
    }

    bool get_starting_point(Index n, bool init_x, Number* x,
                            bool init_z, Number*, Number*,
                            Index, bool init_lambda, Number*) override
    {
        assert(init_x && !init_z && !init_lambda);
        (void) init_x; (void) init_z; (void) init_lambda;
        poly_starting_point(NV, x);
        (void) n;
        return true;
    }

    bool eval_f(Index n, const Number* x, bool, Number& obj_value) override
    {
        obj_value = poly_objective<double>(x, NV);
        (void) n;
        return true;
    }

    bool eval_g(Index n, const Number* x, bool, Index, Number* g) override
    {
        poly_constraints<double>(x, NV, g);
        (void) n;
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
        printf("\n\nObjective value (area = -f)\n");
        printf("f(x*) = %e\n", obj_value);
        printf("Non-zeros: Jacobian %zu, Hessian of the Lagrangian %zu\n",
               jac_subset_.nnz(), hes_subset_.nnz());

        // Named for the vertex count, because plotpolygon.m draws the two cases
        // side by side and so needs both files present.
        char name[64];
        snprintf(name, sizeof(name), "ipopt%d.m", NV);
        FILE* f = fopen(name, "w");
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
    PolygonNLP(const PolygonNLP&);
    PolygonNLP& operator=(const PolygonNLP&);

    // Record fg(x) = [f(x); g(x)] once, then derive from it the sparsity of the
    // constraint Jacobian and of the Hessian of the Lagrangian.
    void record_tape()
    {
        DVector x0(n_);
        poly_starting_point(NV, &x0[0]);

        ADVector ax(n_);
        for (int i = 0; i < n_; i++) ax[i] = x0[i];

        CppAD::Independent(ax);
        ADVector afg(1 + m_);
        afg[0] = poly_objective<adouble>(&ax[0], NV);
        poly_constraints<adouble>(&ax[0], NV, &afg[1]);
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
        int nv = atoi(argv[1]);
        if (nv < 4) {
            printf("usage: %s [number of vertices >= 4]\n", argv[0]);
            return 1;
        }
        NV = nv;
    }
    printf("Largest small polygon with %d vertices.\n", NV);

    SmartPtr<TNLP> nlp = new PolygonNLP();
    SmartPtr<IpoptApplication> app = IpoptApplicationFactory();

    // The adaptive barrier-parameter strategy is not a detail here. The problem
    // has many local maxima -- a polygon can be locally optimal without being
    // the largest -- and the two strategies follow different central paths and
    // stop at different ones. With the monotone default this program converges
    // in 54 iterations to a polygon of area 0.7268687, which is a genuine local
    // solution but not the best one; with the adaptive strategy it takes 11
    // iterations to reach 0.7840773, the value reported in the book and the
    // known optimum for nv = 50.
    app->Options()->SetStringValue("mu_strategy", "adaptive");
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
