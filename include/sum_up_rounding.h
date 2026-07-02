#ifndef PSOPT_SUM_UP_ROUNDING_H
#define PSOPT_SUM_UP_ROUNDING_H

// ---------------------------------------------------------------------------
// Sum-up rounding (header-only utility)
//
// References:
//   S. Sager, G. Reinelt, H.G. Bock, "Direct methods with maximal lower bound
//     for mixed-integer optimal control problems", Math. Program. 118 (2009)
//     109-149.
//   S. Sager, H.G. Bock, M. Diehl, "The integer approximation error in
//     mixed-integer optimal control", Math. Program. 133 (2012) 1-23.
//
// Reconstruct a single-active-mode integer control from the relaxed weights of
// an outer-convexified (SOS1) integer control, by the accumulator sweep. The
// routine is value-agnostic: it operates on the relaxed weights and the mesh
// spacing only, so it serves the binary and the general M-value cases alike.
//
// Given relaxed weights omega_i on a control mesh of N intervals for M modes
// (each column ideally a convex combination summing to one), the sweep keeps,
// for every mode, the running integral of the difference between the relaxed
// weight and the emitted binary indicator. On each interval it activates the
// mode whose accumulated relaxed integral is largest, so the emitted integer
// control tracks the relaxed one in an integral sense. The maximum accumulated
// difference (integral_gap) is the quantity bounded, to first order in the mesh
// spacing h, by the Sager-Bock-Diehl theorem; refining the mesh drives it to
// zero linearly in h.
//
// This is an approximation, not a certified integer optimum, and the emitted
// control is not guaranteed feasible for the original path/terminal constraints
// (the bound is on the trajectory, not on constraint satisfaction); callers
// should re-simulate the rounded control and check the violation.
// ---------------------------------------------------------------------------

#include <Eigen/Dense>
#include <algorithm>

// weights      : M x N, relaxed omega_i per interval (M modes, N intervals).
// h            : 1 x N, interval widths (may be non-uniform).
// mode_index   : (out) 1 x N, activated mode 0..M-1 on each interval.
// integral_gap : (out) max over modes and prefixes of |sum_{k<=n} (omega_i - b_i) h_k|.
// n_switches   : (out) number of mode changes along the mesh.
inline void sum_up_rounding(const Eigen::MatrixXd&    weights,
                            const Eigen::RowVectorXd& h,
                            Eigen::RowVectorXi&       mode_index,
                            double&                   integral_gap,
                            int&                      n_switches)
{
    const int M = static_cast<int>(weights.rows());
    const int N = static_cast<int>(weights.cols());

    mode_index.resize(N);
    integral_gap = 0.0;
    n_switches   = 0;
    if (M <= 0 || N <= 0) return;

    // Per-mode accumulated (relaxed - emitted) integral.
    Eigen::VectorXd phi = Eigen::VectorXd::Zero(M);

    int prev = -1;
    for (int n = 0; n < N; ++n) {
        const double hn = h(n);

        // Accumulate the relaxed integral over this interval.
        phi += weights.col(n) * hn;

        // Activate the most-owed mode; lowest index on ties keeps it deterministic.
        int i = 0;
        double best = phi(0);
        for (int k = 1; k < M; ++k) {
            if (phi(k) > best) { best = phi(k); i = k; }
        }
        mode_index(n) = i;

        // Subtract the emitted (binary) integral of the activated mode.
        phi(i) -= hn;

        // Diagnostics: accumulated deviation and switch count.
        integral_gap = std::max(integral_gap, phi.cwiseAbs().maxCoeff());
        if (prev >= 0 && i != prev) ++n_switches;
        prev = i;
    }
}

// Convenience overload: also map the activated modes to their control values.
// values          : 1 x M, the discrete value associated with each mode.
// rounded_control : (out) 1 x N reconstructed integer control values.
inline void sum_up_rounding(const Eigen::MatrixXd&    weights,
                            const Eigen::RowVectorXd& h,
                            const Eigen::RowVectorXd& values,
                            Eigen::RowVectorXi&       mode_index,
                            Eigen::RowVectorXd&       rounded_control,
                            double&                   integral_gap,
                            int&                      n_switches)
{
    sum_up_rounding(weights, h, mode_index, integral_gap, n_switches);
    const int N = static_cast<int>(mode_index.size());
    rounded_control.resize(N);
    for (int n = 0; n < N; ++n) rounded_control(n) = values(mode_index(n));
}

#endif // PSOPT_SUM_UP_ROUNDING_H
