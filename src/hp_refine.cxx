// hp_refine.cxx -- hp-adaptive (Route B) automatic mesh refinement driver.
//
// Liu-Hager-Rao ph scheme for the multi-interval Radau (LGR) discretisation built by
// lgr_nodes_multi. After each solve PSOPT fills solution.relative_errors[i] with the Betts
// per-collocation-node relative ODE error; this driver groups those into the K mesh
// intervals and, for every interval whose error exceeds algorithm.ode_tolerance, decides
// how to refine. The smoothness of each interval is measured from the decay of the
// Legendre coefficients of the state on that interval (the principled indicator), and the
// predicted polynomial degree follows from the exponential-convergence model
//        error ~ exp(-sigma * N).
//
// INCREMENT 3a: pure p-refinement. Every interval needing refinement has its order raised
// to the predicted degree; intervals (and hence breakpoints) are never split here. The
// h-vs-p branch (split when the predicted degree is too large / the interval is non-smooth)
// is added in increment 3b and will reuse the same sigma computed here.
//
// Generated with AI assistance.

#include "psopt.h"
#include <cmath>
#include <vector>
#include <algorithm>
using namespace Eigen;

// Upper limit on N_eff for the hp-auto workspace ceiling. Staying below it keeps the
// single-block LGR node generators in range and bounds the worst-case storage the
// automatic driver may request.
static const int HP_MAX_NODES = 600;

// ---------------------------------------------------------------------------------------
// Workspace-sizing ceiling. The hp mesh grows across refinement iterations, but the
// Workspace is allocated once before the loop, so get_max_nodes must report an a-priori
// upper bound on N_eff. Each of the mr_max_iterations rounds adds at most
// ceil(mr_max_growth_factor * nodes0) collocation points, bounded by the pseudospectral node
// limit so the single-block node generators stay in range.
int hp_node_ceiling( Prob& problem, Alg& algorithm, int iphase )
{
    int nodes0 = (int) problem.phase[iphase-1].nodes(0);
    int per_iter = (int) std::ceil(algorithm.mr_max_growth_factor * nodes0);
    if ( per_iter < 1 ) per_iter = 1;
    int R = nodes0 + algorithm.mr_max_iterations * per_iter;
    if ( R > HP_MAX_NODES - 1 ) R = HP_MAX_NODES - 1;
    if ( R < nodes0 )                    R = nodes0;
    return R;
}

// ---------------------------------------------------------------------------------------
// Legendre polynomials P_0..P_m at xi, by the standard three-term recurrence.
static void legendre_column(double xi, int m, VectorXd& P)
{
    P.resize(m+1);
    P(0) = 1.0;
    if (m >= 1) P(1) = xi;
    for (int k=1; k<m; k++)
        P(k+1) = ( (2*k+1)*xi*P(k) - k*P(k-1) ) / (double)(k+1);
}

// Exponential decay rate sigma of the Legendre-coefficient envelope of a degree-m
// polynomial sampled at m+1 points (xi in [-1,1], local interval coordinate). Larger
// sigma = faster spectral decay = smoother. Returns sigma clamped to [sigma_min, +inf).
// The coefficients are obtained by exact interpolation (square Legendre Vandermonde solve);
// the slope of log|a_k| over the upper half of the spectrum is the decay estimate.
static double legendre_decay_rate(const VectorXd& xi, const VectorXd& vals, double sigma_min)
{
    int np = (int) xi.size();      // m+1 points
    int m  = np - 1;               // polynomial degree
    if (m < 2) return 1.0;         // too few points to estimate; treat as smooth

    MatrixXd V(np, np);
    VectorXd P;
    for (int i=0;i<np;i++){ legendre_column(xi(i), m, P); V.row(i) = P.transpose(); }
    VectorXd a = V.colPivHouseholderQr().solve(vals);   // Legendre coefficients a_0..a_m

    // log|a_k| with a floor to avoid log(0); least-squares slope over upper-half degrees.
    int k0 = std::max(1, (int)std::ceil(m/2.0));
    int n  = m - k0 + 1;
    if (n < 2) k0 = 1, n = m;      // fall back to all degrees >=1 for small m
    double Sk=0, Sy=0, Skk=0, Sky=0;
    for (int k=k0;k<=m;k++){
        double ak = std::fabs(a(k));
        double y  = std::log( std::max(ak, 1e-16) );
        Sk+=k; Sy+=y; Skk+=(double)k*k; Sky+=(double)k*y;
    }
    double denom = n*Skk - Sk*Sk;
    double slope = (std::fabs(denom) > 1e-30) ? (n*Sky - Sk*Sy)/denom : 0.0;
    double sigma = -slope;
    if (!(sigma > sigma_min)) sigma = sigma_min;   // also catches NaN
    return sigma;
}

// ---------------------------------------------------------------------------------------
// The driver. Rewrites problem.phase[i].hp_orders (and, from 3b, hp_breakpoints) in place.
void hp_refine_driver( Prob& problem, Alg& algorithm, Sol& solution, Workspace* workspace )
{
    const double eps_tol   = algorithm.ode_tolerance;
    const double sigma_min = 0.05;

    // Per-interval polynomial-degree bounds (literature N_min / N_max, Patterson-Hager-Rao
    // 2015; Darby-Hager-Rao 2011), now exposed as Alg fields.
    const int    N_min     = algorithm.mr_min_order;  // floor order for a freshly split sub-interval
    const int    N_max     = algorithm.mr_max_order;  // p->h switch: never p-refine beyond this
    const double sigma_lo  = 0.5;      // below this an interval is non-smooth -> prefer h
    const double min_width = 1.0e-4;   // smallest sub-interval (in (0,1)) the driver will create

    for (int i=0; i<problem.nphases; i++)
    {
        if ( !hp_mesh_active(problem.phase[i]) ) continue;   // driver seeds this at iter 1

        int    nstates = problem.phase[i].nstates;
        int    K       = hp_num_intervals(problem.phase[i]);
        int    N_eff   = problem.phase[i].current_number_of_intervals;   // = sum(hp_orders)
        int    R       = hp_node_ceiling(problem, algorithm, i+1);

        const MatrixXd& err = solution.relative_errors[i];        // 1 x N_eff (per node)
        const MatrixXd  X   = solution.get_states_in_phase(i+1);   // nstates x (N_eff+1)
        const MatrixXd& sn  = workspace->snodes[i];               // (N_eff+1) x 1, tau in [-1,1]

        RowVectorXi old_orders = problem.phase[i].hp_orders;
        RowVectorXd old_breaks = problem.phase[i].hp_breakpoints;   // size K-1, in (0,1)

        bool gauss = ( algorithm.collocation_method == "Gauss" );

        // Map each interval to its error columns / storage nodes. Radau shares breakpoints, so
        // storage = sum(orders) and the per-interval stride is the order. Gauss stores every
        // non-collocated breakpoint separately, so storage = sum(orders) + K - 1 and the stride
        // is order+1; the last interval's trailing column is clamped by the err/N_eff bound.
        std::vector<int> start(K+1, 0);
        for (int j=0;j<K;j++) start[j+1] = start[j] + old_orders(j) + (gauss ? 1 : 0);

        // per-iteration node budget: bound mesh growth to ~mr_max_growth_factor and never
        // exceed the workspace ceiling R. This keeps the interval count from exploding (worst
        // intervals are served first; the rest wait for a later iteration) and guarantees the
        // once-allocated workspace is never overrun.
        int N_target = (int) std::ceil( N_eff * (1.0 + algorithm.mr_max_growth_factor) );
        if (N_target < N_eff + 1) N_target = N_eff + 1;
        if (N_target > R)         N_target = R;

        // interval error and smoothness
        std::vector<double> e_j(K, 0.0), sig(K, 0.0);
        for (int j=0;j<K;j++) {
            int s=start[j], eend=start[j+1];
            for (int c=s; c<eend && c<N_eff; c++) e_j[j] = std::max(e_j[j], err(0,c));
            if (e_j[j] > eps_tol) {
                // Legendre-decay smoothness over the interval's own nodes. Radau uses its
                // order+1 nodes (right endpoint = shared breakpoint). Gauss adds the stored
                // right breakpoint (order+2 nodes) so the fit spans the true interval; the
                // last Gauss interval has no stored right breakpoint (the appended x_f), so
                // it falls back to order+1 nodes.
                int np   = old_orders(j) + 1 + ((gauss && j<K-1) ? 1 : 0);
                int last = s + np - 1;
                double a=sn(s), b=sn(last), half=0.5*(b-a);
                VectorXd xi(np);
                for (int c=0;c<np;c++) xi(c)=(sn(s+c)-0.5*(a+b))/half;
                double sg=1.0e30;
                for (int l=0;l<nstates;l++){ VectorXd v(np);
                    for (int c=0;c<np;c++) v(c)=X(l,s+c);
                    sg=std::min(sg, legendre_decay_rate(xi,v,sigma_min)); }
                sig[j]=sg;
            }
        }

        // decisions: 0 keep, 1 p (to p_new[j]), 2 h (binary split, sub-order h_sub[j]).
        std::vector<int> action(K,0), p_new(K,0), h_sub(K,0);

        // serve intervals worst-error first, spending the budget.
        std::vector<int> idx;
        for (int j=0;j<K;j++) if (e_j[j] > eps_tol) idx.push_back(j);
        std::sort(idx.begin(), idx.end(), [&](int p,int q){ return e_j[p] > e_j[q]; });

        int proj = N_eff;
        for (size_t r=0; r<idx.size() && proj < N_target; r++) {
            int j  = idx[r];
            int nj = old_orders(j);
            double left  = (j==0)   ? 0.0 : old_breaks(j-1);
            double right = (j==K-1) ? 1.0 : old_breaks(j);
            bool smooth  = ( sig[j] >= sigma_lo );

            // p-refinement when the interval is smooth and still below the order cap.
            if ( smooth && nj < N_max ) {
                int dN_raw = std::max(1, (int)std::ceil( std::log(e_j[j]/eps_tol)/sig[j] ));
                int dN_cap = std::max(2, (int)std::ceil(algorithm.mr_max_growth_factor*nj)+1);
                int N_new  = std::min( nj + std::min(dN_raw, dN_cap), N_max );
                int cost   = N_new - nj;
                if (cost >= 1 && proj + cost <= N_target) { action[j]=1; p_new[j]=N_new; proj+=cost; }
                else if (proj + 1 <= N_target)            { action[j]=1; p_new[j]=nj+1;  proj+=1;   }
                continue;
            }

            // otherwise h-refinement: binary split, sub-order chosen so the split is affordable
            // (roughly node-preserving) while still adding resolution. Too-narrow intervals fall
            // back to a unit p-bump if any order headroom remains.
            if ( (right-left) >= 2.0*min_width ) {
                int sub  = std::max(N_min, std::min(N_max, nj/2 + 1));
                int cost = 2*sub - nj + (gauss ? 1 : 0);   // Gauss: the new interface breakpoint is an extra stored node
                if (proj + cost <= N_target) { action[j]=2; h_sub[j]=sub; proj+=cost; }
                else if (nj < N_max && proj + 1 <= N_target) { action[j]=1; p_new[j]=nj+1; proj+=1; }
            } else if (nj < N_max && proj + 1 <= N_target) {
                action[j]=1; p_new[j]=nj+1; proj+=1;
            }
        }

        // assemble the new mesh from the decisions
        std::vector<int> no; std::vector<double> nb;
        for (int j=0;j<K;j++) {
            double left  = (j==0)   ? 0.0 : old_breaks(j-1);
            double right = (j==K-1) ? 1.0 : old_breaks(j);
            if (action[j]==2) {                          // binary h-split
                nb.push_back( 0.5*(left+right) );
                no.push_back( h_sub[j] ); no.push_back( h_sub[j] );
            } else if (action[j]==1) {                   // p
                no.push_back( p_new[j] );
            } else {                                     // keep
                no.push_back( old_orders(j) );
            }
            if (j < K-1) nb.push_back(right);
        }

        // proj was held <= N_target <= R, so storage <= R by construction; assert defensively.
        // Storage is sum(orders) for Radau and sum(orders)+K-1 for Gauss (the interface nodes).
        int tot=0; for (size_t q=0;q<no.size();q++) tot+=no[q];
        int storage = tot + (gauss ? (int)no.size() - 1 : 0);
        if (storage > R) {                               // should not happen; clamp if it does
            while (storage > R) {
                int qm=-1, om=N_min;
                for (size_t q=0;q<no.size();q++) {
                    if (no[q]>om) { om=no[q]; qm=(int)q; }
                }
                if (qm<0) break;
                no[qm]-=1; tot--; storage--;
            }
        }

        RowVectorXi new_orders(no.size());
        for (size_t q=0;q<no.size();q++) new_orders((int)q)=no[q];
        RowVectorXd new_breaks(nb.size());
        for (size_t q=0;q<nb.size();q++) new_breaks((int)q)=nb[q];
        problem.phase[i].hp_orders      = new_orders;
        problem.phase[i].hp_breakpoints = new_breaks;
    }
}
