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


#include "psopt.h"
#include <vector>
#include <algorithm>
#include <cmath>
#include <map>
#include <set>
#include <unordered_map>
#include <chrono>

// Bring std names into this translation unit (formerly leaked via psopt.h).
using namespace std;
using namespace Ipopt;   // Ipopt names no longer leak from psopt.h



// constructor
IPOPT_PSOPT::IPOPT_PSOPT(Workspace *pr, void *user_data)
{
    workspace       = pr;
    _user_data      = user_data;
}

//destructor
IPOPT_PSOPT::~IPOPT_PSOPT()
{}

adouble Lagrangian_ad(adouble* xad, double* lambda, double& obj_factor, Index& m, Workspace* workspace)
{
	adouble L;
	adouble f;
	adouble *g = workspace->gad.get();
	Index i;

	L = obj_factor*ff_ad(xad, workspace);

	gg_ad(xad, g, workspace);

	for(i=0; i<m ; i++) {
		L += lambda[ i ]*g[ i ];
	}

	return L;
}


// ----------------------------------------------------------------------------
// Numerical sparse Lagrangian Hessian (algorithm.hessian == "numerical").
// Lets IPOPT use an exact-Hessian step under numerical derivatives instead of
// falling back to limited-memory L-BFGS. The Hessian of the Lagrangian
// H = obj_factor*grad^2 f + sum_k lambda_k grad^2 g_k is symmetric; we supply
// its lower triangle. Sparsity is taken structurally from the already-detected
// constraint-Jacobian pattern (every variable pair the objective couples is, in
// PSOPT's collocation NLP, also coupled by some constraint); values are formed
// by finite differences of the scalar Lagrangian. See ComputeHessianNonZeros...
// ----------------------------------------------------------------------------

// Scalar Lagrangian evaluated numerically, reusing the same ff_num/gg_num that
// feed IPOPT's eval_f/eval_g (so the Hessian is consistent with the problem
// IPOPT actually sees). G uses the shared constraint buffer workspace->Gip.
static double Lagrangian_num(MatrixXd& X, const double* lambda, double obj_factor,
                             int m, Workspace* workspace)
{
    MatrixXd& G = *workspace->Gip;
    double L = obj_factor * ff_num(X, workspace);
    gg_num(X, &G, workspace);
    for (int i = 0; i < m; i++) L += lambda[i] * G(i);
    return L;
}

// One probe of the Hessian sparsity. At point X with multipliers lambda, estimate each
// candidate (superset) entry by a cheap FORWARD mixed difference of the Lagrangian and
// accumulate, per entry, the largest magnitude seen across probes (maxAbsE) together with the
// global maximum (globalMaxE). Thresholding is deferred to the caller so that a single GLOBAL
// threshold is applied: a probe where the Hessian happens to be small cannot lower the bar and
// admit round-off as a false positive. Forward accuracy is irrelevant here - a structurally-zero
// entry has every mixed a-b derivative vanishing, so it reads as round-off. The 2n single-variable
// perturbations are reused, so the probe costs ~ (2n + #off-diagonals) Lagrangian evaluations.
static void ProbeHessianPattern(MatrixXd& X, const double* lambda, double obj_factor, int m,
                                const std::vector<long>& keys, int n,
                                std::vector<double>& maxAbsE, double& globalMaxE,
                                double& maxLmag, Workspace* workspace)
{
    const double qrt_eps = 1.2243397568e-04;          // eps^(1/4)
    const double L0 = Lagrangian_num(X, lambda, obj_factor, m, workspace);
    if (std::fabs(L0) > maxLmag) maxLmag = std::fabs(L0);

    std::vector<double> Lp((size_t) n), h((size_t) n);
    for (int i = 0; i < n; i++) {
        const double xi = X(i);
        const double hi = qrt_eps * std::max(1.0, std::fabs(xi));
        h[(size_t) i] = hi;
        X(i) = xi + hi; Lp[(size_t) i] = Lagrangian_num(X, lambda, obj_factor, m, workspace);
        X(i) = xi;
    }

    const size_t N = keys.size();
    for (size_t e = 0; e < N; e++) {
        const int a = (int) (keys[e] / n), b = (int) (keys[e] % n);
        if (a == b) continue;                          // diagonals kept unconditionally
        const double xa = X(a), xb = X(b), ha = h[(size_t) a], hb = h[(size_t) b];
        X(a) = xa + ha; X(b) = xb + hb;
        const double Lpp = Lagrangian_num(X, lambda, obj_factor, m, workspace);
        X(a) = xa; X(b) = xb;
        const double Eab = std::fabs((Lpp - Lp[(size_t) a] - Lp[(size_t) b] + L0) / (ha * hb));
        if (Eab > maxAbsE[e]) maxAbsE[e] = Eab;
        if (Eab > globalMaxE) globalMaxE = Eab;
    }
}

// Build the lower-triangular Lagrangian-Hessian sparsity pattern (row >= col). A structural
// superset is first formed from the (non-constant) constraint-Jacobian pattern: union over
// each constraint row of (variables in that row) x (variables in that row), plus every
// diagonal. For pseudospectral / free-final-time discretisations this superset can be far
// denser than the true Hessian (a t_f-scaled linear term registers as non-constant yet has
// zero second derivative), so it is then pruned by probing: the FD Hessian is evaluated at
// several non-degenerate points (results unioned) and entries that are zero at every probe are
// dropped. Up to two (rows,cols,count)
// source blocks may be supplied (the second is optional; pass count2=0). Fills
// workspace->hess_ir/hess_jc with the pruned pattern and returns the number of nonzeros.
static int DetectHessianSparsityNumerical(int n, int m,
                                          const int* iR1, const int* jC1, int cnt1,
                                          const int* iR2, const int* jC2, int cnt2,
                                          Workspace* workspace)
{
    // ---- 1. structural superset ------------------------------------------------------------
    std::vector< std::vector<int> > rowcols((size_t) m);
    for (int e = 0; e < cnt1; e++) { int r = iR1[e]; if (r >= 0 && r < m) rowcols[r].push_back(jC1[e]); }
    for (int e = 0; e < cnt2; e++) { int r = iR2[e]; if (r >= 0 && r < m) rowcols[r].push_back(jC2[e]); }

    std::vector<long> keys;
    for (int r = 0; r < m; r++) {
        std::vector<int>& cols = rowcols[(size_t) r];
        std::sort(cols.begin(), cols.end());
        cols.erase(std::unique(cols.begin(), cols.end()), cols.end());
        const size_t d = cols.size();
        for (size_t p = 0; p < d; p++)
            for (size_t q = 0; q <= p; q++)            // cols sorted ascending => cols[p] >= cols[q]
                keys.push_back((long) cols[p] * (long) n + (long) cols[q]);
    }
    for (int i = 0; i < n; i++) keys.push_back((long) i * (long) n + (long) i);   // all diagonals

    std::sort(keys.begin(), keys.end());
    keys.erase(std::unique(keys.begin(), keys.end()), keys.end());

    // ---- 2. probe-and-prune (union over several non-degenerate probes) ---------------------
    const size_t N = keys.size();
    std::vector<char> keep(N, 0);
    for (size_t e = 0; e < N; e++) if (keys[e] / n == keys[e] % n) keep[e] = 1;   // keep all diagonals

    MatrixXd& X   = *workspace->Xip;
    MatrixXd& x0  = *workspace->x0;
    MatrixXd& xlb = *workspace->xlb;
    MatrixXd& xub = *workspace->xub;

    // Several probes at distinct points and multiplier vectors are unioned: an entry is kept if
    // it is nonzero at ANY probe. This guards against an entry that is incidentally zero at one
    // point or one lambda (state- or multiplier-dependent sparsity). Each probe nudges every
    // component away from zero by a distinct amount (clipped to the variable bounds) and uses
    // distinct, all-nonzero multipliers.
    const int    NPROBE   = 3;
    const double pAmp[3]  = { 0.05, 0.17, 0.31 };   // point nudge amplitudes (fractions of scale)
    const double pFrq[3]  = { 0.70, 1.30, 0.40 };   // point phase frequencies
    const double pPhs[3]  = { 0.30, 1.70, 3.10 };
    const double lFrq[3]  = { 2.30, 1.10, 3.70 };   // multiplier frequencies
    const double lPhs[3]  = { 1.10, 2.40, 0.60 };

    std::vector<double> lam((size_t) m);
    std::vector<double> maxAbsE(N, 0.0);
    double globalMaxE = 0.0;
    double maxLmag    = 0.0;
    for (int p = 0; p < NPROBE; p++) {
        for (int i = 0; i < n; i++) {
            const double xi = x0(i);
            double xp = xi + pAmp[p] * std::max(1.0, std::fabs(xi))
                             * (1.0 + 0.5 * std::sin(pFrq[p] * (double) i + pPhs[p]));
            const double lo = xlb(i), hi = xub(i);
            if (lo > -1.0e19 && xp < lo) xp = lo;
            if (hi <  1.0e19 && xp > hi) xp = hi;
            X(i) = xp;
        }
        for (int k = 0; k < m; k++)
            lam[(size_t) k] = 1.0 + 0.5 * std::sin(lFrq[p] * (double) k + lPhs[p]);
        ProbeHessianPattern(X, lam.data(), 1.0, m, keys, n, maxAbsE, globalMaxE, maxLmag, workspace);
    }
    // Keep an off-diagonal entry if its largest magnitude over all probes clears a single global
    // threshold (diagonals already kept). The threshold is the larger of (i) a relative bar,
    // 1e-6 of the global maximum second derivative, and (ii) a round-off floor: a forward mixed
    // difference of the Lagrangian cannot resolve entries below ~ eps*|L|/delta^2, and with
    // delta ~ eps^(1/4) that floor is ~ 1e-6*|L|. The floor stops probes with a large Lagrangian
    // magnitude from admitting round-off as false positives.
    const double thresh = std::max(1.0e-6 * globalMaxE, 1.0e-6 * maxLmag);
    for (size_t e = 0; e < N; e++) if (maxAbsE[e] > thresh) keep[e] = 1;

    // ---- 3. write the pruned pattern (guarded against the fixed buffer) ---------------------
    int nnz = 0;
    for (size_t e = 0; e < N; e++) if (keep[e]) nnz++;

    if (nnz > workspace->hess_nnz_capacity) {
        snprintf(workspace->text, sizeof(workspace->text),
                 "numerical Hessian needs algorithm.hess_sparsity_ratio just above %f",
                 (double) nnz / ((double) n * (double) n));
        error_message(workspace->text);
    }

    int w = 0;
    for (size_t e = 0; e < N; e++) if (keep[e]) {
        workspace->hess_ir[w] = (unsigned int) (keys[e] / n);
        workspace->hess_jc[w] = (unsigned int) (keys[e] % n);
        w++;
    }
    return nnz;
}

// Fill the Hessian values over the detected lower-triangular pattern by central
// second differences of the scalar Lagrangian: diagonal entries reuse the 2n
// single-variable perturbations; off-diagonal entries use the central mixed
// stencil (four evaluations each). Step h ~ eps^(1/4) (near-optimal for second
// differences). H2 will reduce the off-diagonal cost via symmetric colouring.
static void ComputeHessianNonZerosNumerical(MatrixXd& X, const double* lambda, double obj_factor,
                                            int m, double* values, int nnz_hess, Workspace* workspace)
{
    const int n = workspace->nvars;
    const double qrt_eps = 1.2243397568e-04;   // eps^(1/4); near-optimal step for 2nd differences

    const double L0 = Lagrangian_num(X, lambda, obj_factor, m, workspace);

    // Single-variable +/- perturbations, reused for the diagonal entries.
    std::vector<double> Lp((size_t) n), Lm((size_t) n), h((size_t) n);
    for (int i = 0; i < n; i++) {
        const double xi = X(i);
        const double hi = qrt_eps * std::max(1.0, std::fabs(xi));
        h[(size_t) i] = hi;
        X(i) = xi + hi; Lp[(size_t) i] = Lagrangian_num(X, lambda, obj_factor, m, workspace);
        X(i) = xi - hi; Lm[(size_t) i] = Lagrangian_num(X, lambda, obj_factor, m, workspace);
        X(i) = xi;
    }

    for (int e = 0; e < nnz_hess; e++) {
        const int a = (int) workspace->hess_ir[e];
        const int b = (int) workspace->hess_jc[e];
        if (a == b) {
            // central second difference, O(h^2)
            values[e] = (Lp[(size_t) a] - 2.0 * L0 + Lm[(size_t) a]) / (h[(size_t) a] * h[(size_t) a]);
        } else {
            // central mixed second difference, O(h^2): (L++ - L+- - L-+ + L--)/(4 ha hb)
            const double xa = X(a), xb = X(b);
            const double ha = h[(size_t) a], hb = h[(size_t) b];
            X(a) = xa + ha; X(b) = xb + hb; const double Lpp = Lagrangian_num(X, lambda, obj_factor, m, workspace);
            X(a) = xa + ha; X(b) = xb - hb; const double Lpm = Lagrangian_num(X, lambda, obj_factor, m, workspace);
            X(a) = xa - ha; X(b) = xb + hb; const double Lmp = Lagrangian_num(X, lambda, obj_factor, m, workspace);
            X(a) = xa - ha; X(b) = xb - hb; const double Lmm = Lagrangian_num(X, lambda, obj_factor, m, workspace);
            X(a) = xa; X(b) = xb;
            values[e] = (Lpp - Lpm - Lmp + Lmm) / (4.0 * ha * hb);
        }
    }
}


// ---- H2a: index-set Hessian -- colouring map and verification ------------------------------
//
// The index-set ("colouring") Hessian computes the sparse Lagrangian Hessian in O(gamma^2)
// group perturbations rather than H1's O(2n + 4*nnz_offdiag) per-entry differences (gamma is the
// CPR colour count, ~ the per-node state+control coupling, which stays small and mesh-independent
// as n grows). It reuses the constraint-Jacobian colouring already built for the sparse Jacobian
// (workspace->igroup): two columns share a group only if they never co-occur in a constraint row,
// so each constraint row has at most one variable per group -- the prerequisite for recovering a
// constraint's Hessian entry from a per-row second difference (H2b).
//
// H2a installs only the colour map and a verification harness; the Hessian values are still
// produced by H1's direct method, so the solve is unchanged.

// Build hess_col_group[j] = colour group of variable j, from workspace->igroup.
static void BuildHessianColourMap(int n, Workspace* workspace)
{
    IGroup* ig = workspace->igroup.get();
    int* grp = workspace->hess_col_group.get();
    for (int j = 0; j < n; j++) grp[j] = -1;
    for (int k = 0; k < ig->number; k++)
        for (int t = 0; t < ig->size[k]; t++)
            grp[ ig->colindex[k][t] ] = k;
}

// Verify the colouring supports per-row Hessian recovery, and quantify objective-only coupling.
// Reports: gamma (group count); per-row group collisions (must be 0 -- two variables of one
// constraint row in the same group would be inseparable); and how many pruned off-diagonal
// Hessian entries are NOT covered by any constraint row (objective-only couplings, which H2b
// must obtain from a direct objective-Hessian computation rather than the constraint colouring).
static void VerifyHessianColouring(int n, int m, int nele_hess, Workspace* workspace)
{
    IGroup* ig = workspace->igroup.get();
    const int* grp  = workspace->hess_col_group.get();
    const int  nnzG = workspace->jac_nnzG;
    const int* iG   = workspace->iGrow.get();
    const int* jG   = workspace->jGcol.get();

    std::vector< std::vector<int> > rowvars((size_t) m);
    for (int e = 0; e < nnzG; e++) { int r = iG[e]; if (r >= 0 && r < m) rowvars[(size_t) r].push_back(jG[e]); }

    int per_row_collisions = 0;
    std::set<long> ccoup;                                   // constraint-coupled pairs (lower triangle)
    for (int r = 0; r < m; r++) {
        std::vector<int>& vs = rowvars[(size_t) r];
        std::vector<int> gs; gs.reserve(vs.size());
        for (size_t t = 0; t < vs.size(); t++) gs.push_back(grp[vs[t]]);
        std::sort(gs.begin(), gs.end());
        for (size_t t = 1; t < gs.size(); t++) if (gs[t] == gs[t-1]) { per_row_collisions++; break; }
        for (size_t p = 0; p < vs.size(); p++)
            for (size_t q = 0; q < p; q++) {
                long a = vs[p], b = vs[q]; if (a < b) { long t = a; a = b; b = t; }
                ccoup.insert(a * (long) n + b);
            }
    }

    int total_off = 0, obj_only = 0;
    for (int e = 0; e < nele_hess; e++) {
        const int a = (int) workspace->hess_ir[e], b = (int) workspace->hess_jc[e];
        if (a == b) continue;
        total_off++;
        const long key = (long) a * (long) n + (long) b;
        if (ccoup.find(key) == ccoup.end()) obj_only++;
    }

    // Structural diagnostics for the gamma blow-up: the Jacobian colouring gamma is bounded
    // below by the densest non-constant Jacobian row; the *Hessian* colouring is bounded below
    // by the densest Hessian column. A single variable coupled to almost everything (e.g. the
    // free final time, whose time-scaling multiplies the differentiation matrix) forces gamma
    // up regardless of the colouring, and is a candidate for direct (un-coloured) handling.
    int max_jac_row = 0;
    for (int r = 0; r < m; r++) if ((int) rowvars[(size_t) r].size() > max_jac_row) max_jac_row = (int) rowvars[(size_t) r].size();

    std::vector<int> hdeg((size_t) n, 0);
    for (int e = 0; e < nele_hess; e++) {
        const int a = (int) workspace->hess_ir[e], b = (int) workspace->hess_jc[e];
        if (a == b) continue;
        hdeg[(size_t) a]++; hdeg[(size_t) b]++;
    }
    int max_hess_deg = 0, dense_cols = 0;
    for (int j = 0; j < n; j++) { if (hdeg[(size_t) j] > max_hess_deg) max_hess_deg = hdeg[(size_t) j]; if (hdeg[(size_t) j] > n / 4) dense_cols++; }

    snprintf(workspace->text, sizeof(workspace->text),
        "\nH2a index-set colouring: groups gamma = %d (n = %d); per-row group collisions = %d "
        "(must be 0); off-diagonal Hessian entries = %d, constraint-covered = %d, objective-only = %d\n"
        "  diagnostics: densest non-constant Jacobian row = %d; densest Hessian column = %d; "
        "Hessian columns with degree > n/4 = %d\n",
        ig->number, n, per_row_collisions, total_off, total_off - obj_only, obj_only,
        max_jac_row, max_hess_deg, dense_cols);
    psopt_print(workspace, workspace->text);
}


// ---- H2b: index-set Hessian value recovery (Option A) --------------------------------------
//
// Computes the sparse Lagrangian Hessian E = obj_factor*grad^2 f + sum_r lambda_r grad^2 g_r into
// out[] (parallel to the pattern hess_ir/hess_jc) by the index-set method:
//
//   * diagonal entries:           direct central 2nd differences of the Lagrangian (2n evals),
//                                 always correct, and side-steps the same-group diagonal-collision
//                                 case (two variables in one group cannot share a constraint row,
//                                 so the per-row scheme never separates same-group diagonals).
//   * off-diagonal constraint part: per-constraint recovery from group-perturbed constraint
//                                 vectors. Each group pair (k,l) that shares a constraint row costs
//                                 four vector g evaluations; the CPR property (<=1 variable per row
//                                 per group) makes each constraint's four-corner difference isolate
//                                 a single Hessian entry (grad^2 g_r)_{a,b}, which is accumulated
//                                 with weight lambda_r. Group pairs that share no row are skipped.
//   * off-diagonal objective part:  direct central mixed 2nd differences of the objective, added
//                                 with weight obj_factor. The objective is one global function whose
//                                 group second differences would mix entries across nodes, so it is
//                                 handled directly (its Hessian is block-diagonal + the dense time
//                                 columns, so this is comparatively cheap).
//
// The split is exact: central-mixed(L) = obj_factor*central-mixed(f) + sum_r lambda_r*central-mixed(g_r),
// and a whole-group perturbation reproduces the single-variable perturbation each constraint sees
// (the other group members do not appear in that row), so out[] matches the H1 direct values to
// round-off. H2b uses this for validation only (the solve still uses the direct values); H2c wires
// it in as the value path with an auto-select on gamma.
static void ComputeHessianIndexSetNumerical(MatrixXd& X, const double* lambda, double obj_factor,
                                            int n, int m, double* out, int nele_hess, Workspace* workspace)
{
    const double qrt_eps = 1.2243397568e-04;   // eps^(1/4)
    IGroup* ig  = workspace->igroup.get();
    const int* grp = workspace->hess_col_group.get();

    std::vector<double> h((size_t) n);
    for (int j = 0; j < n; j++) h[(size_t) j] = qrt_eps * std::max(1.0, std::fabs(X(j)));

    std::unordered_map<long,int> emap;
    emap.reserve((size_t) nele_hess * 2);
    for (int e = 0; e < nele_hess; e++)
        emap[(long) workspace->hess_ir[e] * (long) n + (long) workspace->hess_jc[e]] = e;
    for (int e = 0; e < nele_hess; e++) out[e] = 0.0;

    // ---- diagonals: direct central Lagrangian second differences ----
    const double L0 = Lagrangian_num(X, lambda, obj_factor, m, workspace);
    for (int a = 0; a < n; a++) {
        std::unordered_map<long,int>::iterator it = emap.find((long) a * (long) n + (long) a);
        if (it == emap.end()) continue;
        const double xa = X(a), ha = h[(size_t) a];
        X(a) = xa + ha; const double Lp = Lagrangian_num(X, lambda, obj_factor, m, workspace);
        X(a) = xa - ha; const double Lm = Lagrangian_num(X, lambda, obj_factor, m, workspace);
        X(a) = xa;
        out[it->second] = (Lp - 2.0 * L0 + Lm) / (ha * ha);
    }

    // ---- off-diagonal constraint part: index-set recovery ----
    const int nnzG = workspace->jac_nnzG;
    const int* iG  = workspace->iGrow.get();
    const int* jG  = workspace->jGcol.get();
    std::vector< std::unordered_map<int,int> > rgv((size_t) m);          // row -> (group -> variable)
    std::vector< std::vector<int> >            group_rows((size_t) ig->number);
    for (int e = 0; e < nnzG; e++) {
        const int r = iG[e], j = jG[e];
        if (r < 0 || r >= m) continue;
        const int g = grp[j];
        if (g < 0) continue;
        rgv[(size_t) r][g] = j;
        group_rows[(size_t) g].push_back(r);
    }

    MatrixXd xbase = X;
    MatrixXd Gpp(m,1), Gpm(m,1), Gmp(m,1), Gmm(m,1);
    std::vector<int> sr, sa, sb;                                         // shared (row, var_k, var_l)
    for (int k = 0; k < ig->number; k++) {
        for (int l = k + 1; l < ig->number; l++) {
            sr.clear(); sa.clear(); sb.clear();
            for (size_t idx = 0; idx < group_rows[(size_t) k].size(); idx++) {
                const int r = group_rows[(size_t) k][idx];
                std::unordered_map<int,int>::iterator itl = rgv[(size_t) r].find(l);
                if (itl == rgv[(size_t) r].end()) continue;
                sr.push_back(r); sa.push_back(rgv[(size_t) r][k]); sb.push_back(itl->second);
            }
            if (sr.empty()) continue;

            for (int corner = 0; corner < 4; corner++) {
                const double sk = (corner & 2) ? -1.0 : 1.0;
                const double sl = (corner & 1) ? -1.0 : 1.0;
                for (int t = 0; t < ig->size[k]; t++) { const int j = ig->colindex[k][t]; X(j) = xbase(j) + sk * h[(size_t) j]; }
                for (int t = 0; t < ig->size[l]; t++) { const int j = ig->colindex[l][t]; X(j) = xbase(j) + sl * h[(size_t) j]; }
                MatrixXd& G = (corner == 0) ? Gpp : (corner == 1) ? Gpm : (corner == 2) ? Gmp : Gmm;
                gg_num(X, &G, workspace);
                for (int t = 0; t < ig->size[k]; t++) X(ig->colindex[k][t]) = xbase(ig->colindex[k][t]);
                for (int t = 0; t < ig->size[l]; t++) X(ig->colindex[l][t]) = xbase(ig->colindex[l][t]);
            }

            for (size_t s = 0; s < sr.size(); s++) {
                const int r = sr[s], a = sa[s], b = sb[s];
                int A = a, B = b; if (A < B) { int tt = A; A = B; B = tt; }
                std::unordered_map<long,int>::iterator it = emap.find((long) A * (long) n + (long) B);
                if (it == emap.end()) continue;
                const double dgr = Gpp(r) - Gpm(r) - Gmp(r) + Gmm(r);
                out[it->second] += lambda[r] * dgr / (4.0 * h[(size_t) a] * h[(size_t) b]);
            }
        }
    }

    // ---- off-diagonal objective part: direct central mixed differences of the objective ----
    // When the objective off-diagonal pattern has been detected (H2d), compute it only where the
    // objective contributes; otherwise (not yet detected) compute over the full pattern (safe).
    const bool obj_restricted = workspace->hess_obj_detected;
    const char* obj_mark = workspace->hess_obj_offdiag.get();
    for (int e = 0; e < nele_hess; e++) {
        const int a = (int) workspace->hess_ir[e], b = (int) workspace->hess_jc[e];
        if (a == b) continue;
        if (obj_restricted && !obj_mark[e]) continue;
        const double xa = X(a), xb = X(b), ha = h[(size_t) a], hb = h[(size_t) b];
        X(a) = xa + ha; X(b) = xb + hb; const double fpp = ff_num(X, workspace);
        X(a) = xa + ha; X(b) = xb - hb; const double fpm = ff_num(X, workspace);
        X(a) = xa - ha; X(b) = xb + hb; const double fmp = ff_num(X, workspace);
        X(a) = xa - ha; X(b) = xb - hb; const double fmm = ff_num(X, workspace);
        X(a) = xa; X(b) = xb;
        out[e] += obj_factor * (fpp - fpm - fmp + fmm) / (4.0 * ha * hb);
    }
}


// H2c: auto-select index-set vs H1-direct for the off-diagonal Hessian values.
// Both methods share the 2n diagonal cost and evaluate every constraint per perturbation, so the
// gg work is what differs: H1-direct evaluates 4 per off-diagonal entry (4*nnz_off), the index-set
// method evaluates 4 per *distinct interacting colour-group pair* (4*P). Each off-diagonal entry
// (a,b) maps to exactly one group pair (group(a),group(b)), so P <= nnz_off always; the index-set
// method strictly compresses the gg work precisely when group pairs are reused across entries
// (P < nnz_off), which is the mesh-independent local/hp structure. Returns the decision and reports
// P and nnz_off through the out-parameters.
static bool DecideHessianIndexSet(int n, int m, int nele_hess, Workspace* workspace,
                                  int& P_out, int& noff_out)
{
    IGroup* ig = workspace->igroup.get();
    const int* grp  = workspace->hess_col_group.get();
    const int  nnzG = workspace->jac_nnzG;
    const int* iG   = workspace->iGrow.get();
    const int* jG   = workspace->jGcol.get();

    std::vector< std::vector<int> > rowvars((size_t) m);
    for (int e = 0; e < nnzG; e++) { const int r = iG[e]; if (r >= 0 && r < m) rowvars[(size_t) r].push_back(jG[e]); }

    std::set<long> pairs;
    for (int r = 0; r < m; r++) {
        std::vector<int>& vs = rowvars[(size_t) r];
        for (size_t p = 0; p < vs.size(); p++)
            for (size_t q = 0; q < p; q++) {
                int k = grp[vs[p]], l = grp[vs[q]];
                if (k < 0 || l < 0 || k == l) continue;
                if (k > l) { int t = k; k = l; l = t; }
                pairs.insert((long) k * (long) ig->number + (long) l);
            }
    }

    int noff = 0;
    for (int e = 0; e < nele_hess; e++) if (workspace->hess_ir[e] != workspace->hess_jc[e]) noff++;

    P_out = (int) pairs.size();
    noff_out = noff;
    return (P_out < noff);
}


// H2d: detect the objective's own off-diagonal Hessian pattern, so the index-set value path
// computes the objective contribution only where grad^2 f is nonzero instead of over the full
// off-diagonal pattern. Robustness follows H1: the objective second derivative is probed at
// several non-degenerate nudged points (reusing ProbeHessianPattern with lambda = 0, which makes
// it probe the objective alone) and an entry is marked if it clears the relative / round-off-floor
// threshold at ANY probe. Marks workspace->hess_obj_offdiag and reports the count.
static int DetectObjectiveHessianPattern(int n, int m, int nele_hess, Workspace* workspace)
{
    char* mark = workspace->hess_obj_offdiag.get();
    for (int e = 0; e < nele_hess; e++) mark[e] = 0;

    std::vector<long> keys; std::vector<int> eidx;
    keys.reserve((size_t) nele_hess); eidx.reserve((size_t) nele_hess);
    for (int e = 0; e < nele_hess; e++) {
        const int a = (int) workspace->hess_ir[e], b = (int) workspace->hess_jc[e];
        if (a == b) continue;
        keys.push_back((long) a * (long) n + (long) b);
        eidx.push_back(e);
    }
    if (keys.empty()) return 0;

    MatrixXd& X   = *workspace->Xip;
    MatrixXd& x0  = *workspace->x0;
    MatrixXd& xlb = *workspace->xlb;
    MatrixXd& xub = *workspace->xub;
    MatrixXd  Xsave = X;                                   // restore the real iterate afterwards

    std::vector<double> zero((size_t) m, 0.0);
    std::vector<double> maxAbsE(keys.size(), 0.0);
    double gmax = 0.0, fmax = 0.0;

    const int    NPROBE  = 3;
    const double pAmp[3] = { 0.05, 0.17, 0.31 };
    const double pFrq[3] = { 0.70, 1.30, 0.40 };
    const double pPhs[3] = { 0.30, 1.70, 3.10 };
    for (int p = 0; p < NPROBE; p++) {
        for (int i = 0; i < n; i++) {
            const double xi = x0(i);
            double xp = xi + pAmp[p] * std::max(1.0, std::fabs(xi))
                             * (1.0 + 0.5 * std::sin(pFrq[p] * (double) i + pPhs[p]));
            const double lo = xlb(i), hi = xub(i);
            if (lo > -1.0e19 && xp < lo) xp = lo;
            if (hi <  1.0e19 && xp > hi) xp = hi;
            X(i) = xp;
        }
        ProbeHessianPattern(X, zero.data(), 1.0, m, keys, n, maxAbsE, gmax, fmax, workspace);
    }
    X = Xsave;

    const double thresh = std::max(1.0e-6 * gmax, 1.0e-6 * fmax);
    int count = 0;
    for (size_t j = 0; j < keys.size(); j++)
        if (maxAbsE[j] > thresh) { mark[eidx[j]] = 1; count++; }
    return count;
}


bool check_no_cancel(void *user_data)
{
#ifdef WIN32
    if (user_data)
    {
        HANDLE *handle = (HANDLE *)user_data;
        int no_cancel = (WAIT_OBJECT_0!=WaitForSingleObject(*handle, 0));
        if (!no_cancel)
            fprintf(stderr, "\n --- User cancel event received! ---");
        return no_cancel;
    }
#endif
    return true;
}


bool IPOPT_PSOPT::intermediate_callback(AlgorithmMode mode,
                                       Index iter, Number obj_value,
                                       Number inf_pr, Number inf_du,
                                       Number mu, Number d_norm,
                                       Number regularization_size,
                                       Number alpha_du, Number alpha_pr,
                                       Index ls_trials,
                                       const IpoptData* ip_data,
                                       IpoptCalculatedQuantities* ip_cq)
{
    return check_no_cancel(_user_data);
}

// returns the size of the problem
bool IPOPT_PSOPT::get_nlp_info(Index& n, Index& m, Index& nnz_jac_g,
                             Index& nnz_h_lag, IndexStyleEnum& index_style)
{
  int nnz;
  int nnzA = 0;
  int nnzG = 0;
  int i;
  double jsratio;


  // Number of variables
  n = workspace->nvars;

  // Number of constraints in g(x)
  m = workspace->ncons;

  MatrixXd *X0 = workspace->x0.get();

  double  *x  = &(*X0)(0);

  if( !useAutomaticDifferentiation(*workspace->algorithm) ) {


     DetectJacobianSparsity(gg_num, *X0, m,  &nnzA,  workspace->iArow.get(), workspace->jAcol.get(), workspace->jac_Aij.get(),
                                             &nnzG,  workspace->iGrow.get(), workspace->jGcol.get(),
                                             workspace->grw.get(), workspace );

     nnz = nnzA+nnzG;

     jsratio = (double) ((double)  nnz/((double) (n*m)));

     if (jsratio > workspace->algorithm->jac_sparsity_ratio)
     {
           snprintf(workspace->text,sizeof(workspace->text), "increase algorithm.jac_sparsity_ratio to just above %f", jsratio);
           error_message(workspace->text);
     }

     snprintf(workspace->text,sizeof(workspace->text),"\nJacobian sparsity detected numerically:");
     psopt_print(workspace,workspace->text);
     snprintf(workspace->text,sizeof(workspace->text),"\n*** %i nonzero elements out of %i [ratio=%f]", nnz, n*m, jsratio );
     psopt_print(workspace,workspace->text);
     snprintf(workspace->text,sizeof(workspace->text),"\n*** %i nonzero elements are constant", nnzA );
     psopt_print(workspace,workspace->text);
     snprintf(workspace->text,sizeof(workspace->text),"\n*** %i nonzero elements are not constant", nnzG );
     psopt_print(workspace,workspace->text);


  }


  if( useAutomaticDifferentiation(*workspace->algorithm) ) {

	psopt_ad::ad_record(workspace->ad_g, n, m, x,
		[&](const adouble* xin, adouble* yout){ gg_ad(const_cast<adouble*>(xin), yout, workspace); });
	psopt_ad::SparseTriplet J = psopt_ad::ad_sparse_jacobian(workspace->ad_g, x, /*reuse=*/false);
	nnz = J.nnz();
	if (nnz > workspace->jac_nnz_capacity) {
		snprintf(workspace->text,sizeof(workspace->text),
			"detected Jacobian non-zeros (%d) exceed the allocated buffer (%d); increase algorithm.jac_sparsity_ratio to just above %f",
			nnz, workspace->jac_nnz_capacity, (double) nnz/((double)(n*m)) );
		error_message(workspace->text);
	}
	for(i=0;i<nnz;i++)
	{
		workspace->jGcol[i] = J.col[i];
		workspace->iGrow[i] = J.row[i];
	}

        snprintf(workspace->text,sizeof(workspace->text),"\nJacobian sparsity detected using ADOLC:");
        psopt_print(workspace,workspace->text);

        jsratio = (double) ((double)  nnz/((double) (n*m)));

        if (jsratio > workspace->algorithm->jac_sparsity_ratio) {
           snprintf(workspace->text,sizeof(workspace->text), "increase algorithm.jac_sparsity_ratio to just above %f", jsratio);
           error_message(workspace->text);
        }

        snprintf(workspace->text,sizeof(workspace->text),"\n%i nonzero elements out of %i [ratio=%f]\n", nnz, n*m, jsratio);
        psopt_print(workspace,workspace->text);

  } // end if (autoderiv)

  int activate_hess;

  if (workspace->algorithm->hessian=="exact")
      activate_hess = 1;
  else
      activate_hess = 0;

  if( activate_hess && useAutomaticDifferentiation(*workspace->algorithm)  ) {

	double  obj_factor = 1.0;

   std::vector<double> lambda_full((size_t)m, 1.0);   // all-ones -> full, lambda-independent Hessian pattern
   double *lambda = lambda_full.data();
   int nnz_hess;
	/* Lagrangian Hessian structure via AD backend (step 2b) */
	psopt_ad::ad_record(workspace->ad_hess, n, 1, x,
		[&](const adouble* xin, adouble* yout){
			yout[0] = Lagrangian_ad(const_cast<adouble*>(xin), lambda, obj_factor, m, workspace); });
	psopt_ad::SparseTriplet Hs = psopt_ad::ad_sparse_hessian(workspace->ad_hess, x, /*reuse=*/false);
	nnz_hess = Hs.nnz();
	if (nnz_hess > workspace->hess_nnz_capacity) {
		snprintf(workspace->text,sizeof(workspace->text),
			"detected Hessian non-zeros (%d) exceed the allocated buffer (%d); increase algorithm.hess_sparsity_ratio to just above %f",
			nnz_hess, workspace->hess_nnz_capacity, (double) nnz_hess/((double)(n*n)) );
		error_message(workspace->text);
	}
	for (i=0; i< nnz_hess; i++) {
		workspace->hess_ir[i] = Hs.row[i];
		workspace->hess_jc[i] = Hs.col[i];
	}

       snprintf(workspace->text,sizeof(workspace->text),"\nHessian sparsity detected using ADOLC:");
       psopt_print(workspace,workspace->text);
       double hsratio = (double) ((double)  nnz_hess/((double) (n*n)));
       if (hsratio > workspace->algorithm->hess_sparsity_ratio) {
            snprintf(workspace->text,sizeof(workspace->text), "increase algorithm.hess_sparsity_ratio to just above %f", hsratio);
            error_message(workspace->text);
       }

       snprintf(workspace->text,sizeof(workspace->text),"\n%i nonzero elements out of %i [ratio = %f] \n", nnz_hess, n*n, hsratio );
       psopt_print(workspace,workspace->text);

       nnz_h_lag = nnz_hess;

  } // end if (autoderiv)

  else if (workspace->algorithm->hessian=="numerical") {

     // Lagrangian-Hessian sparsity built structurally from the constraint-Jacobian
     // pattern. Only the NON-CONSTANT Jacobian entries (iGrow/jGcol) can contribute
     // to the Hessian: a constant dg_k/dx_j means g_k is linear in x_j, so its second
     // derivative vanishes. Excluding the constant block (iArow) keeps the pattern
     // tight -- crucially, it drops the dense but linear differentiation-matrix coupling
     // of pseudospectral methods. Any genuinely nonlinear pair has both variables in the
     // non-constant set, so nothing is lost. Automatic derivatives expose only the full
     // pattern (nnz), with no constant/non-constant split.
     int nnz_hess;
     if ( useAutomaticDifferentiation(*workspace->algorithm) )
        nnz_hess = DetectHessianSparsityNumerical(n, m,
                       workspace->iGrow.get(), workspace->jGcol.get(), nnz,
                       (const int*)NULL, (const int*)NULL, 0, workspace);
     else
        nnz_hess = DetectHessianSparsityNumerical(n, m,
                       workspace->iGrow.get(), workspace->jGcol.get(), nnzG,
                       (const int*)NULL, (const int*)NULL, 0, workspace);

     snprintf(workspace->text,sizeof(workspace->text),"\nHessian sparsity built numerically from the Jacobian pattern:");
     psopt_print(workspace,workspace->text);
     double hsratio = (double) ((double) nnz_hess/((double) (n*n)));
     snprintf(workspace->text,sizeof(workspace->text),"\n%i nonzero elements out of %i [ratio = %f] \n", nnz_hess, n*n, hsratio );
     psopt_print(workspace,workspace->text);

     nnz_h_lag = nnz_hess;

  } // end numerical Hessian

    nnz_jac_g = nnz;


  // the hessian is in this case assumed to be a square dense matrix but we
  // only need the lower left corner (since it is symmetric)
  if( ( !useAutomaticDifferentiation(*workspace->algorithm) || workspace->algorithm->hessian!="exact" )
        && workspace->algorithm->hessian!="numerical" )
        nnz_h_lag = (int) ((n*n)+n)/2;

/*   *
     * *
     * * *
     * * * *
     * * * * *
*/

  // use the C style indexing (0-based)
  index_style = TNLP::C_STYLE;

  return true;
}

// returns the variable bounds
bool IPOPT_PSOPT::get_bounds_info(Index n, Number* x_l, Number* x_u,
                                Index m, Number* g_l, Number* g_u)
{

  // here, the n and m we gave IPOPT in get_nlp_info are passed back to us.
  // If desired, we could assert to make sure they are what we think they are.
  assert(n == workspace->nvars);
  assert(m == workspace->ncons);


    double *xlb = &(*workspace->xlb)(0);

    double *xub = &(*workspace->xub)(0);


  int  j;


  // lower bounds on x
  for (j=0; j<workspace->nvars; j++) {
    x_l[j] = xlb[j];
  }

  // upper bounds on x
  for (j=0; j<workspace->nvars; j++) {
    x_u[j] = xub[j];
  }

  get_constraint_bounds(g_l, g_u, workspace);

  return true;
}

// returns the initial point for the problem
bool IPOPT_PSOPT::get_starting_point(Index n, bool init_x, Number* x,
                                   bool init_z, Number* z_L, Number* z_U,
                                   Index m, bool init_lambda,
                                   Number* lambda)
{
  // Here, we assume we only have starting values for x, if you code
  // your own NLP, you can provide starting values for the dual variables
  // if you wish
  assert(init_x == true);
//  assert(init_z == false);
//  assert(init_lambda == false);

  Index i;

//double *x0 = (workspace->x0)->GetPr();
  double *x0 = &(*workspace->x0)(0);

  // initialize to the given starting point

  for (i=0; i<workspace->nvars;i++)
  {
	  x[i] = x0[i];
  }


  return true;
}

// returns the value of the objective function
bool IPOPT_PSOPT::eval_f(Index n, const Number* x, bool new_x, Number& obj_value)
{
  assert(n == workspace->nvars);

  MatrixXd& X = *workspace->Xip;

  memcpy( &X(0), x, workspace->nvars*sizeof(double) );

  obj_value = ff_num(X, workspace);

  return true;
}

// return the gradient of the objective function grad_{x} f(x)
bool IPOPT_PSOPT::eval_grad_f(Index n, const Number* x, bool new_x, Number* grad_f)
{
  assert(n == workspace->nvars);

  MatrixXd& X = *workspace->Xip;

  MatrixXd& GF = *workspace->GFip;


  memcpy( &X(0), x, workspace->nvars*sizeof(double) );

  if(!useAutomaticDifferentiation(*workspace->algorithm))
     ScalarGradient( ff_num, X, &GF , workspace->grw.get(), workspace );
  else
     ScalarGradientAD( ff_ad, X, &GF, &workspace->trace_f_done, workspace->ad_f, workspace );

  memcpy( grad_f, &GF(0), workspace->nvars*sizeof(double));

  return true;
}

// return the value of the constraints: g(x)
bool IPOPT_PSOPT::eval_g(Index n, const Number* x, bool new_x, Index m, Number* g)
{
  assert(n == workspace->nvars);
  assert(m == workspace->ncons);

  MatrixXd& X = *workspace->Xip;

  MatrixXd& G  = *workspace->Gip;

  memcpy( &X(0), x, workspace->nvars*sizeof(double) );

  gg_num(X, &G, workspace);

  memcpy( g, &G(0), workspace->ncons*sizeof(double) );

  return true;
}




void save_jacobian_sparsity_pattern(Index* rindex, Index* cindex, long nvars, long ncols, long nnz, Workspace* workspace)
{
    // Saves the jacobian structure in triplet form with 1 at the non-zeros.

    FILE* jac_file;

    jac_file = fopen("jacobian_pattern.dat", "w");

    if (jac_file == NULL) {
         error_message("save_jacobian_sparsity_pattern(): error opening jacobian pattern file");
    }

    double one  = 1.0;
    double zero = 0.0;

    fprintf(jac_file,"% li %li  %f", ncols, nvars, zero );

    for (int i = 1; i< nnz; i++ ) {
//        fprintf(jac_file,"\n%li %li  %f", rindex[i]+1, cindex[i]+1, one ); // SPARSITY PATTERN SAVED USING 1-BASED INDICES (MATLAB-STYLE)
         fprintf(jac_file,"\n%i %i  %f", rindex[i]+1, cindex[i]+1, one ); // SPARSITY PATTERN SAVED USING 1-BASED INDICES (MATLAB-STYLE)
    }

    fclose(jac_file);

    workspace->algorithm->save_sparsity_pattern = 0;

}


// THE UPDATE OF THIS FUNCTION NEEDS TO BE FINALISED - EIGEN_UPDATE.

// return the structure or values of the jacobian
bool IPOPT_PSOPT::eval_jac_g(Index n, const Number* x, bool new_x,
                           Index m, Index nele_jac, Index* iRow, Index *jCol,
                           Number* values)
{
  MatrixXd& X = *workspace->Xip;

  double *xpr = &(*workspace->xp)(0); 

  int nnzA, nnzG, i;


  if (values == NULL) {
  // return the structure of the jacobian
    X = *workspace->x0;
    if (!useAutomaticDifferentiation(*workspace->algorithm))
    {


        	nnzA = workspace->jac_nnzA;
        	nnzG = workspace->jac_nnzG;


        	getIndexGroups( workspace->igroup.get(), m, n, nnzG, workspace->iGrow.get(), workspace->jGcol.get(), workspace);

        	// H2: the colour map and auto-select decision are derived from this colouring, so they
        	// must be rebuilt whenever the colouring is (i.e. at every mesh refinement). The
        	// workspace persists across meshes, so reset the build-once flags here.
        	workspace->hess_maps_built = false;
        	workspace->hess_obj_detected = false;

	     	for (i=0;i<nnzG;i++)
      	{
				iRow[i] = workspace->iGrow[i]; // EIGEN_UPDATE
				jCol[i] = workspace->jGcol[i]; // EIGEN_UPDATE
			}

			for (i=0;i<nnzA;i++)
			{
				iRow[i+nnzG] =workspace->iArow[i]; // EIGEN_UPDATE
				jCol[i+nnzG] =workspace->jAcol[i]; // EIGEN_UPDATE
			}

    }


    if (useAutomaticDifferentiation(*workspace->algorithm)) {


			for(i=0;i<nele_jac;i++)
			{
				iRow[i] = workspace->iGrow[i];
				jCol[i] = workspace->jGcol[i];
			}

    } // End if (autoderiv)
    
  }
  else {
    // return the values of the jacobian of the constraints
    if (!useAutomaticDifferentiation(*workspace->algorithm)) {

          memcpy( &X(0), x, workspace->nvars*sizeof(double) );
          
          // Compute by sparse finite differences only the non-constant Jacobian elements...
          EfficientlyComputeJacobianNonZeros(gg_num, X, m, values, workspace->jac_nnzG, workspace->iGrow.get(),workspace->jGcol.get(), workspace->igroup.get(), workspace->grw.get(), workspace );

          // Now include in array values[] the constant Jacobian elements calculated previously
          for(i=0;i<workspace->jac_nnzA;i++) {
                  values[workspace->jac_nnzG+i] = workspace->jac_Aij[i];
          }

    }

    /* find the jacobian values */

    if (useAutomaticDifferentiation(*workspace->algorithm)) {

	   for (i=0;i<n;i++) {
		   xpr[i] = x[i];
	   }
	   psopt_ad::SparseTriplet J = psopt_ad::ad_sparse_jacobian(workspace->ad_g, xpr, /*reuse=*/false);
	   for(i=0;i<J.nnz();i++) {
                values[i] = J.val[i];
	   }



    }

    if (workspace->enable_nlp_counters) {
        workspace->solution->mesh_stats[ workspace->current_mesh_refinement_iteration-1 ].n_jacobian_evals++;
    }

  }

  if (workspace->algorithm->save_sparsity_pattern == 1) {
     save_jacobian_sparsity_pattern(iRow , jCol, n, m, nele_jac, workspace );
  }

  return true;
}


void IPOPT_PSOPT::finalize_solution(SolverReturn status,
                                  Index n, const Number* x, const Number* z_L, const Number* z_U,
                                  Index m, const Number* g, const Number* lambda,
                                  Number obj_value,
				  const IpoptData* ip_data,
				  IpoptCalculatedQuantities* ip_cq)
{
  // here is where we would store the solution to variables, or write to a file, etc
  // so we could use the solution.

  Sol* solution = workspace->solution;


    memcpy( &(*workspace->x0)(0), x, n*sizeof(double) );


    memcpy( &(*workspace->lambda)(0), lambda, m*sizeof(double) );

  for(int ii=0;ii<n;ii++) solution->xad[ii]=x[ii];

}






//return the structure or values of the hessian
bool IPOPT_PSOPT::eval_h(Index n, const Number* x, bool new_x,
                       Number obj_factor, Index m, const Number* lambda,
                       bool new_lambda, Index nele_hess, Index* iRow,
                       Index* jCol, Number* values)
{

  int i;

  const bool numH    = (workspace->algorithm->hessian=="numerical");
  const bool exactAD = (workspace->algorithm->hessian=="exact")
                       && useAutomaticDifferentiation(*workspace->algorithm);

  if (!numH && !exactAD) return false;

  if (values == NULL) {
    // structure phase: the lower-triangular pattern (both paths) is in hess_ir/hess_jc
    for(i=0;i<nele_hess;i++)
    {
        iRow[i] = workspace->hess_ir[i];
        jCol[i] = workspace->hess_jc[i];
    }
  }
  else {
    // return the values of the Hessian

    if (numH) {

       MatrixXd& X = *workspace->Xip;
       memcpy( &X(0), x, workspace->nvars*sizeof(double) );
       double  obj_factor_d = obj_factor;
       double* lambda_d     = workspace->lambda_d.get();
       for(i=0;i<m;i++) lambda_d[i] = lambda[i];

       // H2c: one-time-per-mesh setup. Build the colour map, count the distinct interacting
       // group pairs, and auto-select the value path (index-set when it compresses the gg work,
       // else H1-direct). hess_maps_built is reset in eval_jac_g whenever the colouring rebuilds,
       // so this re-runs at every mesh refinement.
       if (!workspace->hess_maps_built && workspace->igroup->number > 0) {
           BuildHessianColourMap(n, workspace);
           int P = 0, noff = 0;
           workspace->hess_use_indexset = DecideHessianIndexSet(n, m, nele_hess, workspace, P, noff);
           int objn = -1;
           if (workspace->hess_use_indexset) {
               // H2d: restrict the objective off-diagonal work to its own pattern. Allocate the
               // mark buffer lazily here (only when the index-set method is actually selected) so
               // that problems on the H1-direct path keep an unperturbed heap and stay
               // bit-identical.
               if (!workspace->hess_obj_offdiag)
                   workspace->hess_obj_offdiag = make_unique<char[]>(workspace->hess_nnz_capacity);
               objn = DetectObjectiveHessianPattern(n, m, nele_hess, workspace);
               workspace->hess_obj_detected = true;
           }
           workspace->hess_maps_built = true;
           if (workspace->algorithm->hessian_verify) {
               VerifyHessianColouring(n, m, nele_hess, workspace);
               snprintf(workspace->text, sizeof(workspace->text),
                   "H2c auto-select: interacting group pairs P = %d, off-diagonal entries = %d -> "
                   "%s (gg-work factor %.2f)%s", P, noff,
                   workspace->hess_use_indexset ? "index-set" : "H1-direct",
                   (P > 0) ? (double) noff / (double) P : 1.0,
                   "\n");
               psopt_print(workspace, workspace->text);
               if (objn >= 0) {
                   snprintf(workspace->text, sizeof(workspace->text),
                       "H2d objective off-diagonal pattern: %d of %d entries carry the objective "
                       "(objective ff work reduced %.1fx)\n",
                       objn, noff, (objn > 0) ? (double) noff / (double) objn : (double) noff);
                   psopt_print(workspace, workspace->text);
               }
           }
       }

       // value path: the auto-selected method writes the Hessian values.
       if (workspace->hess_maps_built && workspace->hess_use_indexset)
           ComputeHessianIndexSetNumerical(X, lambda_d, obj_factor_d, n, m, values, nele_hess, workspace);
       else
           ComputeHessianNonZerosNumerical(X, lambda_d, obj_factor_d, m, values, nele_hess, workspace);

       // Optional one-shot check: compare the FD Hessian to the CppAD Hessian at this same
       // (x, lambda, obj_factor). Deferred until the Lagrangian Hessian is non-trivial (early
       // iterates can have zero multipliers and a zero objective Hessian, which would compare
       // zeros). Guarded - if the problem's functions are not AD-traceable (the usual reason for
       // numerical derivatives) the AD evaluation throws and is skipped rather than failing the
       // solve. The AD work buffers (gad, ad_hess) are always allocated.
       if (workspace->algorithm->hessian_verify && !workspace->hess_verify_done) {
           double vmax = 0.0;
           for (i = 0; i < nele_hess; i++) if (std::fabs(values[i]) > vmax) vmax = std::fabs(values[i]);
           if (vmax > 1.0e-10) {
             workspace->hess_verify_done = true;        // attempt only once, at a meaningful iterate

             // H2b/H2c: cross-check the two finite-difference methods against each other and time
             // both at this iterate (cost study). Pure finite differences - no AD needed.
             {
                 std::vector<double> vis((size_t) nele_hess, 0.0), vdir((size_t) nele_hess, 0.0);
                 const int REP = 5;
                 std::chrono::steady_clock::time_point ta = std::chrono::steady_clock::now();
                 for (int rep = 0; rep < REP; rep++)
                     ComputeHessianIndexSetNumerical(X, lambda_d, obj_factor_d, n, m, vis.data(), nele_hess, workspace);
                 std::chrono::steady_clock::time_point tb = std::chrono::steady_clock::now();
                 for (int rep = 0; rep < REP; rep++)
                     ComputeHessianNonZerosNumerical(X, lambda_d, obj_factor_d, m, vdir.data(), nele_hess, workspace);
                 std::chrono::steady_clock::time_point tc = std::chrono::steady_clock::now();
                 const double t_is  = std::chrono::duration<double, std::micro>(tb - ta).count() / REP;
                 const double t_dir = std::chrono::duration<double, std::micro>(tc - tb).count() / REP;

                 double maxd = 0.0, scale = 0.0, maxd_off = 0.0;
                 for (int e = 0; e < nele_hess; e++) {
                     if (std::fabs(vdir[e]) > scale) scale = std::fabs(vdir[e]);
                     const double d = std::fabs(vis[e] - vdir[e]);
                     if (d > maxd) maxd = d;
                     if (workspace->hess_ir[e] != workspace->hess_jc[e] && d > maxd_off) maxd_off = d;
                 }
                 snprintf(workspace->text, sizeof(workspace->text),
                     "\nH2 index-set vs H1-direct: max|IS-direct| = %.3e (off-diagonal %.3e; scale %.3e)\n"
                     "H2c cost study: index-set %.1f us/eval_h, H1-direct %.1f us/eval_h, speed-up %.2fx "
                     "(selected: %s)\n",
                     maxd, maxd_off, scale, t_is, t_dir, (t_is > 0.0) ? t_dir / t_is : 1.0,
                     workspace->hess_use_indexset ? "index-set" : "H1-direct");
                 psopt_print(workspace, workspace->text);
             }

             try {
               std::vector<double> xloc(x, x + n);
               psopt_ad::ad_record(workspace->ad_hess, n, 1, x,
                   [&](const adouble* xin, adouble* yout){
                       yout[0] = Lagrangian_ad(const_cast<adouble*>(xin), lambda_d, obj_factor_d, m, workspace); });
               psopt_ad::SparseTriplet Had =
                   psopt_ad::ad_sparse_hessian(workspace->ad_hess, xloc.data(), /*reuse=*/false);

               std::map<long,double> admap;
               for (int e = 0; e < Had.nnz(); e++) {
                   long r = Had.row[e], c = Had.col[e];
                   if (r < c) { long t = r; r = c; c = t; }     // normalise to lower triangle
                   admap[r * (long) n + c] = Had.val[e];
               }
               double maxabs = 0.0, hscale = 0.0;
               int num_extra = 0;
               for (int e = 0; e < nele_hess; e++) {
                   const long key = (long) workspace->hess_ir[e] * (long) n + (long) workspace->hess_jc[e];
                   const double vnum = values[e];
                   std::map<long,double>::iterator it = admap.find(key);
                   const double vad = (it == admap.end()) ? 0.0 : it->second;
                   if (it != admap.end()) { if (std::fabs(vad) > hscale) hscale = std::fabs(vad); admap.erase(it); }
                   else if (std::fabs(vnum) > 0.0) num_extra++;
                   const double d = std::fabs(vnum - vad);
                   if (d > maxabs) maxabs = d;
               }
               int n_missed = 0; double max_missed = 0.0;       // AD entries absent from the FD pattern
               for (std::map<long,double>::iterator it = admap.begin(); it != admap.end(); ++it) {
                   n_missed++;
                   if (std::fabs(it->second) > max_missed) max_missed = std::fabs(it->second);
                   if (std::fabs(it->second) > hscale)     hscale     = std::fabs(it->second);
               }
               snprintf(workspace->text, sizeof(workspace->text),
                   "\nNumerical-Hessian check vs CppAD: max|FD-AD| = %.3e (Hessian scale %.3e); "
                   "AD entries missed by pruning = %d (max |missed| = %.3e); extra FD entries = %d\n",
                   maxabs, hscale, n_missed, max_missed, num_extra);
               psopt_print(workspace, workspace->text);
           } catch (...) {
               snprintf(workspace->text, sizeof(workspace->text),
                   "\nNumerical-Hessian check: CppAD evaluation unavailable for this problem; check skipped\n");
               psopt_print(workspace, workspace->text);
           }
           }
       }

       if (workspace->enable_nlp_counters) {
           workspace->solution->mesh_stats[ workspace->current_mesh_refinement_iteration-1 ].n_hessian_evals++;
       }

    }
    else if (nele_hess>0) {     // exactAD
//    	double *xpr = workspace->Xsnopt->GetPr();
       double *xpr = &(*workspace->Xsnopt)(0);


// *******************************************************************
	double  obj_factor_d = obj_factor;
	double*  lambda_d     = workspace->lambda_d.get();
	for(i=0;i<m;i++)
		lambda_d[i] = lambda[i];
	psopt_ad::ad_record(workspace->ad_hess, n, 1, x,
		[&](const adouble* xin, adouble* yout){ yout[0] = Lagrangian_ad(const_cast<adouble*>(xin), lambda_d, obj_factor_d, m, workspace); });
	for (i=0;i<n;i++) {
		xpr[i] = x[i];
	}
	psopt_ad::SparseTriplet H = psopt_ad::ad_sparse_hessian(workspace->ad_hess, xpr, /*reuse=*/false);
	nele_hess = H.nnz();
	for(i=0;i<nele_hess;i++) {
		values[i] = H.val[i];
	}

	if (workspace->enable_nlp_counters) {
	    workspace->solution->mesh_stats[ workspace->current_mesh_refinement_iteration-1 ].n_hessian_evals++;
	}

    }

  }

  return true;

}



