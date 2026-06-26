// pseudospectral_rg.cxx
//
// Legendre-Gauss-Radau (LGR) and Legendre-Gauss (LG) primitives for PSOPT,
// written in the same Eigen idiom as src/pseudospectral.cxx and designed to
// drop into the existing square (norder+1) storage scaffold:
//
//   * snodes are returned ASCENDING, length norder+1, matching the post-sort
//     ordering that psopt.cxx relies on (so no descending/negate/sort dance).
//   * the rectangular differentiation matrix occupies the top rows of a
//     (norder+1) x (norder+1) matrix; the trailing row(s) are zero, to be
//     handled as zero-padded / closure rows by the defect assembly.
//
// Radau (the method wired first): N = norder collocation points including the
// initial point tau = -1, plus the non-collocated terminal tau = +1 as the
// last storage node.  This keeps PSOPT's final_states / events semantics valid.
//
// Node generation uses the von Winckel Newton iterations (capped), cross-checked
// to machine precision against an independent scipy reference.
//
// Build (standalone self-test):
//   g++ -O2 -DRG_STANDALONE_TEST -I/usr/include/eigen3 pseudospectral_rg.cxx -o rg && ./rg
//
// Inside PSOPT this file is compiled WITHOUT RG_STANDALONE_TEST and the two
// generators are declared in psopt.h (see wiring notes).

#include <Eigen/Dense>
#include <cmath>
#include <algorithm>
#include <vector>

using Eigen::MatrixXd;
using Eigen::RowVectorXd;
using Eigen::RowVectorXi;

#ifndef PSOPT_PI
#define PSOPT_PI 3.14159265358979323846
#endif

// ---------------------------------------------------------------------------
// Rectangular differentiation matrix on a set of state points `s`, evaluated
// at collocation points `c`:  Drect(k,m) = L_m'(c_k), with L_m the Lagrange
// basis polynomial built on `s`.  Shapes: s is M, c is K, Drect is K x M.
// ---------------------------------------------------------------------------
static MatrixXd rg_diffmat(const MatrixXd& s, const MatrixXd& c)
{
    const int M = (int)s.size();
    const int K = (int)c.size();
    MatrixXd D(K, M);
    for (int i = 0; i < K; ++i) {
        const double xc = c(i);
        for (int j = 0; j < M; ++j) {
            double total = 0.0;
            for (int m = 0; m < M; ++m) {
                if (m == j) continue;
                double term = 1.0 / (s(j) - s(m));
                for (int l = 0; l < M; ++l) {
                    if (l == j || l == m) continue;
                    term *= (xc - s(l)) / (s(j) - s(l));
                }
                total += term;
            }
            D(i, j) = total;
        }
    }
    return D;
}

// ---------------------------------------------------------------------------
// Legendre-Gauss-Radau generator.
//   N  = number of collocation points (= current_number_of_intervals).
//   x  : (N+1) x 1 ascending storage nodes  [ N Radau pts incl. -1 ; +1 ].
//   w  : (N+1) x 1 quadrature weights, w(N) = 0 (terminal slot, no quadrature).
//   D  : (N+1) x (N+1); top N rows = rectangular collocation matrix, last row 0.
// ---------------------------------------------------------------------------
void lgr_nodes(int N, MatrixXd& x, MatrixXd& w, MatrixXd& D)
{
    const int Kc = N;             // number of Radau collocation points
    const int M  = N + 1;         // number of storage / state points
    MatrixXd col(Kc, 1), wc(Kc, 1), xold(Kc, 1), P(Kc, Kc + 1);

    // Chebyshev-Gauss-Radau initial guess
    for (int i = 0; i < Kc; ++i)
        col(i) = -std::cos(2.0 * PSOPT_PI * i / (2.0 * (Kc - 1) + 1.0));
    xold.setConstant(2.0);

    int iter = 0;
    while ((col - xold).cwiseAbs().maxCoeff() > 1e-15 && iter++ < 100) {
        xold = col;
        for (int j = 0; j <= Kc; ++j) P(0, j) = std::pow(-1.0, (double)j);
        for (int i = 1; i < Kc; ++i) { P(i, 0) = 1.0; P(i, 1) = col(i); }
        for (int k = 2; k <= Kc; ++k)
            for (int i = 1; i < Kc; ++i)
                P(i, k) = ((2.0 * k - 1.0) * col(i) * P(i, k - 1)
                           - (k - 1.0) * P(i, k - 2)) / (double)k;
        for (int i = 1; i < Kc; ++i)
            col(i) = xold(i) - ((1.0 - xold(i)) / Kc)
                     * (P(i, Kc - 1) + P(i, Kc)) / (P(i, Kc - 1) - P(i, Kc));
        col(0) = -1.0;
    }

    // collocation weights
    wc(0) = 2.0 / ((double)Kc * Kc);
    for (int i = 1; i < Kc; ++i)
        wc(i) = (1.0 - col(i)) / ((double)Kc * Kc * P(i, Kc - 1) * P(i, Kc - 1));

    // assemble storage nodes (ascending): Radau collocation pts, then terminal +1
    x.resize(M, 1);
    for (int i = 0; i < Kc; ++i) x(i) = col(i);
    x(M - 1) = 1.0;

    w.resize(M, 1);
    for (int i = 0; i < Kc; ++i) w(i) = wc(i);
    w(M - 1) = 0.0;               // terminal node carries no quadrature weight

    // rectangular differentiation matrix in the top Kc rows of a square store
    MatrixXd Drect = rg_diffmat(x, col);   // Kc x M
    D.resize(M, M);
    D.setZero();
    D.topRows(Kc) = Drect;        // last row stays zero (terminal closure / pad)
}

// ---------------------------------------------------------------------------
// Multi-interval (hp) Legendre-Gauss-Radau assembler with SHARED breakpoints.
//
// A K-interval LGR mesh with shared interior breakpoints is structurally
// identical to a single LGR block of effective order N_eff = sum(orders): it has
// N_eff collocation points plus one non-collocated global terminal node, so the
// storage layout produced here has exactly the same shape as lgr_nodes(N_eff),
// and every downstream consumer (defects, cost quadrature, events, variable and
// constraint counts, offsets) reuses unchanged. The only structural differences,
// all confined to the three returned arrays, are:
//   * x : per-interval LGR collocation points affine-mapped onto their
//         sub-interval of tau in [-1,1], concatenated, each shared breakpoint
//         counted once, with the global terminal tau=+1 last (ascending overall).
//   * w : each interval's local weights scaled by Dtau_j/2 so the composite
//         quadrature on [-1,1] equals the sum of the per-interval quadratures;
//         the terminal slot w(N_eff) = 0.
//   * D : block structure. Each interval's local rectangular block (the top n_j
//         rows of lgr_nodes(n_j)) is scaled by 2/Dtau_j (the affine d/dtau
//         factor) and placed so consecutive blocks overlap by exactly one
//         column - the shared breakpoint - which imposes interior C0 continuity
//         implicitly. The terminal (last) row is zero, the single-block closure.
//
//   breakpoints : (K-1) interior breaks, strictly increasing, in (0,1) of the
//                 normalised phase time (0 -> tau=-1, 1 -> tau=+1).
//   orders      : (K) per-interval Radau orders n_j (collocation points each).
//
// AI provenance: this composite assembler was generated with AI assistance
// (increment 1 of the Route-B hp-adaptive work) and is validated by a standalone
// regression check against lgr_nodes plus polynomial differentiation/quadrature
// exactness tests; see the increment-1 notes.
// ---------------------------------------------------------------------------
void lgr_nodes_multi(const RowVectorXd& breakpoints, const RowVectorXi& orders,
                     MatrixXd& x, MatrixXd& w, MatrixXd& D)
{
    const int K = (int) orders.size();
    int N_eff = 0;
    for (int j = 0; j < K; ++j) N_eff += orders(j);
    const int M = N_eff + 1;

    // breakpoints in (0,1) -> sub-interval endpoints s_j on [-1,1].
    std::vector<double> s(K + 1);
    s[0] = -1.0;  s[K] = 1.0;
    for (int j = 1; j < K; ++j) s[j] = -1.0 + 2.0 * breakpoints(j - 1);

    x.resize(M, 1);  x.setZero();
    w.resize(M, 1);  w.setZero();
    D.resize(M, M);  D.setZero();

    int c = 0;                                    // running composite collocation index
    for (int j = 0; j < K; ++j) {
        const int    n  = orders(j);
        const double a  = s[j], b = s[j + 1];
        const double dt = b - a;                  // sub-interval width in tau (> 0)

        MatrixXd xj, wj, Dj;
        lgr_nodes(n, xj, wj, Dj);                 // local block on [-1,1], order n

        // collocation nodes of this interval are the first n local nodes; the
        // local terminal node (xi=+1) is the next interval's first node (shared
        // breakpoint) or, for the last interval, the global terminal handled below.
        for (int r = 0; r < n; ++r) {
            x(c + r) = a + (xj(r) + 1.0) * 0.5 * dt;   // local xi -> tau in [a,b]
            w(c + r) = wj(r) * dt * 0.5;               // composite quadrature weight
        }
        // local rectangular block (top n rows, n+1 columns) scaled by 2/dt; local
        // column n (terminal) lands on composite column c+n = the shared breakpoint
        // (first node of the next interval) or the global terminal (last interval).
        const double dscale = 2.0 / dt;
        for (int row = 0; row < n; ++row)
            for (int col = 0; col <= n; ++col)
                D(c + row, c + col) = dscale * Dj(row, col);
        c += n;
    }
    x(M - 1) = 1.0;     // global terminal tau=+1 (= s[K]); w(M-1)=0, D last row 0
}

// ---------------------------------------------------------------------------
// Legendre-Gauss generator (node primitives validated; NLP wiring deferred
// pending the final_states design decision - see notes).
//   N  = number of interior collocation points.
//   x  : (N+1) x 1 ascending storage nodes  [ -1 ; N Gauss pts ].
//   w  : (N+1) x 1 weights, w(0) = 0 (initial node carries no Gauss weight).
//   D  : (N+1) x (N+1); rows 1..N = rectangular collocation matrix, row 0 zero.
// ---------------------------------------------------------------------------
// Pure m-point Gauss-Legendre nodes and weights mapped to the unit interval [0,1].
// Used by the integrated-residual transcription to integrate ||r||^2 over each mesh
// interval at points distinct from the collocation knots (where the representing
// polynomial makes the residual vanish by construction). The Newton iteration is the
// same one validated in lg_nodes; here we return only the m interior nodes/weights,
// rescaled so that sum(w01) = 1 (i.e. integral_0^1 g dtau ~ sum_j w01_j g(tau_j)).
void gauss_legendre_unit(int m, MatrixXd& nodes01, MatrixXd& w01)
{
    const int Kc = m;
    const int n  = Kc - 1;
    const int N2 = Kc + 1;
    MatrixXd y(Kc,1), y0(Kc,1), Lp(Kc,1), L(Kc, N2+1), wc(Kc,1);

    for (int i = 0; i < Kc; ++i) {
        double xu = (Kc>1) ? (-1.0 + 2.0*i/(double)(Kc-1)) : 0.0;
        y(i) = std::cos((2.0*i + 1.0)*PSOPT_PI/(2.0*n + 2.0))
             + (0.27/Kc)*std::sin(PSOPT_PI*xu*n/(double)N2);
    }
    y0.setConstant(2.0);
    int iter = 0;
    while ((y - y0).cwiseAbs().maxCoeff() > 1e-15 && iter++ < 100) {
        L.col(0).setOnes();
        L.col(1) = y;
        for (int k = 2; k <= Kc; ++k)
            L.col(k) = ((2.0*k - 1.0)*y.cwiseProduct(L.col(k-1)) - (k - 1.0)*L.col(k-2)) / (double)k;
        for (int i = 0; i < Kc; ++i)
            Lp(i) = N2*(L(i, Kc-1) - y(i)*L(i, N2-1)) / (1.0 - y(i)*y(i));
        y0 = y;
        for (int i = 0; i < Kc; ++i) y(i) = y0(i) - L(i, N2-1)/Lp(i);
    }
    for (int i = 0; i < Kc; ++i)
        wc(i) = 2.0/((1.0 - y(i)*y(i))*Lp(i)*Lp(i)) * std::pow((double)N2/Kc, 2.0);

    // ascending sort
    for (int i = 0; i < Kc; ++i)
        for (int j = i+1; j < Kc; ++j)
            if (y(j) < y(i)) { std::swap(y(i), y(j)); std::swap(wc(i), wc(j)); }

    // map (-1,1) -> (0,1):  tau = (xi+1)/2,  w01 = wc/2
    nodes01.resize(Kc,1);
    w01.resize(Kc,1);
    for (int i = 0; i < Kc; ++i) {
        nodes01(i) = 0.5*(y(i) + 1.0);
        w01(i)     = 0.5*wc(i);
    }
}

void lg_nodes(int N, MatrixXd& x, MatrixXd& w, MatrixXd& D)
{
    const int Kc = N;             // interior Gauss collocation points
    const int M  = N + 1;         // storage points: initial + collocation
    const int n  = Kc - 1;
    const int N2 = Kc + 1;
    MatrixXd y(Kc, 1), y0(Kc, 1), Lp(Kc, 1), L(Kc, N2 + 1), wc(Kc, 1);

    for (int i = 0; i < Kc; ++i) {
        double xu = -1.0 + 2.0 * i / (double)(Kc - 1);
        y(i) = std::cos((2.0 * i + 1.0) * PSOPT_PI / (2.0 * n + 2.0))
             + (0.27 / Kc) * std::sin(PSOPT_PI * xu * n / (double)N2);
    }
    y0.setConstant(2.0);

    int iter = 0;
    while ((y - y0).cwiseAbs().maxCoeff() > 1e-15 && iter++ < 100) {
        L.col(0).setOnes();
        L.col(1) = y;
        for (int k = 2; k <= Kc; ++k)
            L.col(k) = ((2.0 * k - 1.0) * y.cwiseProduct(L.col(k - 1))
                        - (k - 1.0) * L.col(k - 2)) / (double)k;
        for (int i = 0; i < Kc; ++i)
            Lp(i) = N2 * (L(i, Kc - 1) - y(i) * L(i, N2 - 1)) / (1.0 - y(i) * y(i));
        y0 = y;
        for (int i = 0; i < Kc; ++i) y(i) = y0(i) - L(i, N2 - 1) / Lp(i);
    }
    for (int i = 0; i < Kc; ++i)
        wc(i) = 2.0 / ((1.0 - y(i) * y(i)) * Lp(i) * Lp(i)) * std::pow((double)N2 / Kc, 2.0);

    // ascending sort of collocation pts/weights
    for (int i = 0; i < Kc; ++i)
        for (int j = i + 1; j < Kc; ++j)
            if (y(j) < y(i)) { std::swap(y(i), y(j)); std::swap(wc(i), wc(j)); }

    x.resize(M, 1);
    x(0) = -1.0;
    for (int i = 0; i < Kc; ++i) x(i + 1) = y(i);

    w.resize(M, 1);
    w(0) = 0.0;                   // initial node carries no Gauss weight
    for (int i = 0; i < Kc; ++i) w(i + 1) = wc(i);

    MatrixXd Drect = rg_diffmat(x, y);    // Kc x M, collocation at the Gauss pts
    D.resize(M, M);
    D.setZero();
    D.bottomRows(Kc) = Drect;     // row 0 (initial, non-collocated) stays zero
}

// ===========================================================================
#ifdef RG_STANDALONE_TEST
#include <cstdio>

static int pass(const char* nm, bool ok, double d) {
    std::printf("  [%s] %s  (%.2e)\n", ok ? "PASS" : "FAIL", nm, d);
    return ok ? 1 : 0;
}

int main()
{
    int all = 1;
    for (int N : {3, 5, 8, 12}) {
        std::printf("\n========  N = %d  ========\n", N);
        MatrixXd x, w, D;

        std::printf("Radau (collocation rows 0..N-1, terminal row N):\n");
        lgr_nodes(N, x, w, D);
        all &= pass("x ascending, x(0)=-1", std::abs(x(0) + 1.0) < 1e-14 &&
                    (x.bottomRows(N) - x.topRows(N)).minCoeff() > 0, 0.0);
        all &= pass("x(N) = +1 (terminal)", std::abs(x(N) - 1.0) < 1e-14, std::abs(x(N) - 1.0));
        all &= pass("weights sum to 2", std::abs(w.sum() - 2.0) < 1e-12, std::abs(w.sum() - 2.0));
        {
            // top N rows annihilate constants and differentiate exactly to deg N
            MatrixXd Dc = D.topRows(N);
            all &= pass("D rows sum 0", Dc.rowwise().sum().cwiseAbs().maxCoeff() < 1e-9,
                        Dc.rowwise().sum().cwiseAbs().maxCoeff());
            double me = 0;
            MatrixXd col = x.topRows(N);
            for (int deg = 1; deg <= N; ++deg) {
                MatrixXd xs = x.array().pow((double)deg);
                MatrixXd dn = Dc * xs;
                MatrixXd de = (double)deg * col.array().pow((double)(deg - 1));
                me = std::max(me, (dn - de).cwiseAbs().maxCoeff());
            }
            all &= pass("D exact poly deg<=N", me < 1e-9, me);
            all &= pass("terminal row is zero", D.row(N).cwiseAbs().maxCoeff() == 0.0, 0.0);
        }

        std::printf("Gauss (collocation rows 1..N, initial row 0):\n");
        lg_nodes(N, x, w, D);
        all &= pass("x(0) = -1 (initial)", std::abs(x(0) + 1.0) < 1e-14, 0.0);
        all &= pass("weights sum to 2", std::abs(w.sum() - 2.0) < 1e-12, std::abs(w.sum() - 2.0));
        {
            MatrixXd Dc = D.bottomRows(N);
            all &= pass("D rows sum 0", Dc.rowwise().sum().cwiseAbs().maxCoeff() < 1e-9,
                        Dc.rowwise().sum().cwiseAbs().maxCoeff());
            double me = 0;
            MatrixXd col = x.bottomRows(N);
            for (int deg = 1; deg <= N; ++deg) {
                MatrixXd xs = x.array().pow((double)deg);
                MatrixXd dn = Dc * xs;
                MatrixXd de = (double)deg * col.array().pow((double)(deg - 1));
                me = std::max(me, (dn - de).cwiseAbs().maxCoeff());
            }
            all &= pass("D exact poly deg<=N", me < 1e-9, me);
            all &= pass("initial row is zero", D.row(0).cwiseAbs().maxCoeff() == 0.0, 0.0);
        }
    }
    std::printf("\n%s\n", all ? "ALL CHECKS PASSED" : "SOME CHECKS FAILED");
    return all ? 0 : 1;
}
#endif
