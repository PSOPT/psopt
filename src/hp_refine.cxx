/*********************************************************************************************

This file is part of the PSOPT library, a software tool for computational optimal control

Copyright (C) 2009-2026 Victor M. Becerra

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
           School of Electrical and Mechanical Engineering
           Portsmouth PO1 3DJ
           United Kingdom
e-mail:    v.m.becerra@ieee.org

**********************************************************************************************/

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
// Estimate where the control switches inside one mesh interval, in phase coordinate (0,1).
//
// A bang-bang control sampled on a polynomial basis does not look like a step: it rings
// about the two plateaus it is trying to join. The largest node-to-node jump is therefore a
// poor estimate of the switch, because the ringing produces jumps of its own near the
// plateaus. The crossings of the level midway between the plateaus are far more stable --
// the ringing sits about the plateaus, not about the middle -- so those are what is used.
//
// An interval is only considered to hold a switch if the control varies over a good
// fraction of its range across the whole phase; that keeps the detector off intervals where
// the control is merely curved.
static void detect_control_switches(const MatrixXd& U, const MatrixXd& sn, int ncontrols,
                                    int s, int np, const VectorXd& urange_phase,
                                    double min_width, std::vector<double>& out)
{
    out.clear();
    if ( np < 3 ) return;
    for (int l=0;l<ncontrols;l++) {
        if ( urange_phase(l) <= 0.0 ) continue;
        double umin=U(l,s), umax=U(l,s);
        for (int c=s;c<s+np;c++){ umin=std::min(umin,U(l,c)); umax=std::max(umax,U(l,c)); }
        if ( (umax-umin) < 0.25*urange_phase(l) ) continue;      // no real transition here
        double umid = 0.5*(umin+umax);
        for (int c=s;c<s+np-1;c++) {
            double a=U(l,c)-umid, b=U(l,c+1)-umid;
            if ( (a<0.0 && b>0.0) || (a>0.0 && b<0.0) ) {
                double th = sn(c) + (sn(c+1)-sn(c))*(-a)/(b-a);  // tau of the crossing
                out.push_back( 0.5*(th+1.0) );                   // -> phase coordinate (0,1)
            }
        }
    }
    std::sort(out.begin(), out.end());
    std::vector<double> merged;                                  // one switch seen in two
    for (size_t q=0;q<out.size();q++)                            // controls is one switch
        if ( merged.empty() || out[q]-merged.back() > 2.0*min_width ) merged.push_back(out[q]);
    out.swap(merged);
}

// The m-1 interior breakpoints for a split of (left,right) into m pieces: the detected
// switches first, then bisection of whatever piece is widest until there are enough. With no
// switches detected this reproduces the uniform split exactly.
static void breakpoints_for_split(double left, double right, int m,
                                  const std::vector<double>& sw, double min_width,
                                  std::vector<double>& bk)
{
    bk.clear();
    for (size_t q=0;q<sw.size() && (int)bk.size()<m-1;q++) {
        double pnt = sw[q];
        if ( pnt-left < min_width || right-pnt < min_width ) continue;
        if ( !bk.empty() && pnt-bk.back() < min_width )      continue;
        bk.push_back(pnt);
    }
    if ( bk.empty() ) {                                          // uniform, as before
        for (int q=1;q<m;q++) bk.push_back( left + (right-left)*(double)q/(double)m );
        return;
    }
    while ( (int)bk.size() < m-1 ) {
        std::vector<double> pts; pts.push_back(left);
        for (size_t q=0;q<bk.size();q++) pts.push_back(bk[q]);
        pts.push_back(right);
        int qm=-1; double wm=-1.0;
        for (size_t q=0;q+1<pts.size();q++){ double w=pts[q+1]-pts[q]; if (w>wm){wm=w;qm=(int)q;} }
        if ( qm<0 || wm < 2.0*min_width ) break;
        bk.push_back( 0.5*(pts[qm]+pts[qm+1]) );
        std::sort(bk.begin(), bk.end());
    }
    std::sort(bk.begin(), bk.end());
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

        int    nstates  = problem.phase[i].nstates;
        int    ncontrols= problem.phase[i].ncontrols;
        int    K       = hp_num_intervals(problem.phase[i]);
        int    N_eff   = problem.phase[i].current_number_of_intervals;   // = sum(hp_orders)
        int    R       = hp_node_ceiling(problem, algorithm, i+1);

        const MatrixXd& err = solution.relative_errors[i];        // 1 x N_eff (per node)
        const MatrixXd  X   = solution.get_states_in_phase(i+1);   // nstates x (N_eff+1)
        const MatrixXd  U   = solution.get_controls_in_phase(i+1); // ncontrols x (N_eff+1)
        const MatrixXd& sn  = workspace->snodes[i];               // (N_eff+1) x 1, tau in [-1,1]

        // range of each control over the whole phase: the scale against which an interval's
        // own control variation is judged large enough to be a switch rather than curvature.
        VectorXd urange_phase(ncontrols>0 ? ncontrols : 1); urange_phase.setZero();
        for (int l=0;l<ncontrols;l++) {
            double lo=U(l,0), hi=U(l,0);
            for (int c=0;c<=N_eff;c++){ lo=std::min(lo,U(l,c)); hi=std::max(hi,U(l,c)); }
            urange_phase(l) = hi-lo;
        }

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
        std::vector< std::vector<double> > sw_pos(K);   // detected switch positions, in (0,1)
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

                // The control has to be looked at as well as the state. Where the Hamiltonian
                // is linear in the control the solution is bang-bang: the states stay
                // continuous and, over a short interval, read as perfectly smooth, while the
                // control jumps. Judging from the states alone then p-refines an interval that
                // straddles a switch, and no polynomial degree approximates a discontinuity, so
                // the error stops falling and the refinement stalls. Measured on examples/wheat:
                // on every interval whose error exceeded the tolerance the state decay rate was
                // between 2.0 and 5.5 -- comfortably "smooth" -- while the control decay rate on
                // the same interval sat on the sigma_min floor of 0.05.
                //
                // Gauss stores each interval's left breakpoint as a non-collocated node, at
                // which the control is not a variable of the NLP; the control fit skips it and
                // uses the interval's own collocation points.
                if ( ncontrols > 0 ) {
                    int su  = s + (gauss ? 1 : 0);
                    int npu = old_orders(j) + (gauss ? 0 : 1);
                    if ( npu >= 3 && su + npu - 1 <= N_eff ) {
                        double au=sn(su), bu=sn(su+npu-1), halfu=0.5*(bu-au);
                        VectorXd xiu(npu);
                        for (int c=0;c<npu;c++) xiu(c)=(sn(su+c)-0.5*(au+bu))/halfu;
                        for (int l=0;l<ncontrols;l++){ VectorXd v(npu);
                            for (int c=0;c<npu;c++) v(c)=U(l,su+c);
                            sg=std::min(sg, legendre_decay_rate(xiu,v,sigma_min)); }
                    }
                }
                sig[j]=sg;

                // and where, inside this interval, the control switches
                if ( algorithm.mr_switch_detection ) {
                    int su  = s + (gauss ? 1 : 0);
                    int npu = old_orders(j) + (gauss ? 0 : 1);
                    if ( ncontrols > 0 && su + npu - 1 <= N_eff )
                        detect_control_switches(U, sn, ncontrols, su, npu, urange_phase,
                                                min_width, sw_pos[j]);
                }
            }
        }

        // decisions: 0 keep, 1 p (to p_new[j]), 2 h (split into h_m[j] pieces of h_sub[j]).
        std::vector<int> action(K,0), p_new(K,0), h_sub(K,0), h_m(K,0);
        std::vector< std::vector<double> > h_bk(K);     // the split's interior breakpoints

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
                // Split into m pieces, m being the smallest number for which m sub-intervals
                // capped at N_max can still carry the order the interval already had. A binary
                // split cannot do that once nj > 2*N_max: two pieces of order N_max total less
                // than nj, so the "refined" mesh is a strictly coarser discretisation than the
                // one it replaces and the error is free to rise. That is not hypothetical. An
                // initial mesh of 32 Gauss nodes is one interval of order 31, and with the
                // default N_max of 12 the first refinement replaced it by two intervals of
                // order 12, taking the total order from 31 down to 24 and the stored nodes from
                // 32 down to 26. The budget test admitted it because the cost was negative.
                // How many pieces. Two demands: enough that the pieces can carry the order
                // the interval already had (m_ord), and one piece per arc when the control's
                // switches have been located. A cut placed on a switch puts the discontinuity
                // on a breakpoint, which is the one place the discretisation is allowed to be
                // discontinuous, and leaves the pieces either side smooth enough for
                // p-refinement to work on them afterwards. A cut at the midpoint leaves the
                // switch inside a piece, where no polynomial degree resolves it.
                int m_ord = std::max(2, (int)std::ceil( (double)(nj+1)/(double)N_max ));
                int m     = std::max(m_ord, (int)sw_pos[j].size() + 1);

                // What order each piece gets. Dividing the parent's order among the pieces is
                // right only when the split was forced by the order cap: a piece created to
                // isolate an arc still has to resolve that arc to the same tolerance, so it
                // wants the parent's order, not a share of it. Take the largest order the
                // budget allows, and only then give up pieces.
                int sub=0, cost=0;
                bool placed=false;
                while ( m >= 2 && !placed ) {
                    int sub_lo = std::max(N_min, std::min(N_max, (int)std::ceil((double)(nj+1)/(double)m)));
                    // Only a piece created to isolate an arc asks for more than its share of
                    // the parent's order. When the split was forced by the order cap alone,
                    // sub_hi == sub_lo and this reduces exactly to dividing the order up.
                    int sub_hi = (m > m_ord) ? std::max(N_min, std::min(N_max, nj)) : sub_lo;
                    for (sub=sub_hi; sub>=sub_lo; sub--) {
                        cost = m*sub - nj + (gauss ? m-1 : 0);   // Gauss: each new interface breakpoint is an extra stored node
                        if (cost >= 1 && proj + cost <= N_target) { placed=true; break; }
                    }
                    if (placed) break;
                    if (m == m_ord) break;                       // cannot afford even the minimum
                    m--; if (m < m_ord) m = m_ord;
                }
                if (placed) {
                    action[j]=2; h_sub[j]=sub; h_m[j]=m; proj+=cost;
                    breakpoints_for_split(left, right, m, sw_pos[j], min_width, h_bk[j]);
                }
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
            if (action[j]==2) {                          // h-split at h_bk[j]
                int m = (int)h_bk[j].size() + 1;
                for (size_t q=0;q<h_bk[j].size();q++) nb.push_back( h_bk[j][q] );
                for (int q=0;q<m;q++) no.push_back( h_sub[j] );
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

        // PSOPT_HP_DEBUG traces the decision the driver took on each interval. The mesh
        // schedule is otherwise invisible in the output, and a refinement that stalls is very
        // hard to diagnose without it.
        if ( getenv("PSOPT_HP_DEBUG") ) {
            fprintf(stderr,"[hp] phase %d  K=%d N_eff=%d N_target=%d R=%d\n", i+1,K,N_eff,N_target,R);
            for (int j=0;j<K;j++)
                fprintf(stderr,"[hp]   old j=%d order=%d e=%.3e sigma=%.3f action=%d p_new=%d h_sub=%d switches=%d\n",
                        j, old_orders(j), e_j[j], sig[j], action[j], p_new[j], h_sub[j],
                        (int)sw_pos[j].size());
            int dbg_tot=0; for (size_t q=0;q<no.size();q++) dbg_tot+=no[q];
            fprintf(stderr,"[hp]   new K=%d orders=[", (int)no.size());
            for (size_t q=0;q<no.size();q++) fprintf(stderr,"%d%s",no[q],q+1<no.size()?",":"");
            fprintf(stderr,"] sum_orders %d -> %d  storage %d -> %d\n",
                    (int)old_orders.sum(), dbg_tot, N_eff+1, dbg_tot+(gauss?(int)no.size()-1:0)+1);
        }

        RowVectorXi new_orders(no.size());
        for (size_t q=0;q<no.size();q++) new_orders((int)q)=no[q];
        RowVectorXd new_breaks(nb.size());
        for (size_t q=0;q<nb.size();q++) new_breaks((int)q)=nb[q];
        problem.phase[i].hp_orders      = new_orders;
        problem.phase[i].hp_breakpoints = new_breaks;
    }
}
