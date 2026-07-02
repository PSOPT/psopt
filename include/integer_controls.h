#ifndef PSOPT_INTEGER_CONTROLS_H
#define PSOPT_INTEGER_CONTROLS_H

// ---------------------------------------------------------------------------
// Integer-control helper building blocks (header-only)
//
// These support the two-stage integer-control workflow (outer convexification
// plus sum-up rounding). They are the reusable numeric primitives underlying
// the ergonomic declaration API: the convex-combination assembly and the SOS1
// residual are used by the outer-convexification set-up, and the guess helper
// expands a mode-index sequence into one-hot weight guesses.
//
// In v1 a phase may declare at most one integer control. Declaring it records
// the control index and its admissible value set on the phase; the expansion
// itself (extra weight-controls, the SOS1 equality path constraint, and the
// convex-combination dynamics) is performed by psopt_level2_setup.
//
// See also include/sum_up_rounding.h for reconstruction of the integer control
// from the relaxed weights.
// ---------------------------------------------------------------------------

#include "psopt.h"
#include "sum_up_rounding.h"

// Declare the single integer control of a phase (v1). control_index is the
// index of the control in the user's control layout; values holds the M
// admissible discrete values. Recorded on the phase and consumed by
// psopt_level2_setup.
inline void declare_integer_control(Prob& problem, int iphase,
                                    int control_index, const RowVectorXd& values)
{
    problem.phases(iphase).integer_control.control_index = control_index;
    problem.phases(iphase).integer_control.values        = values;
}

// Convex combination of per-mode derivative arrays:
//     derivatives[k] = sum_{i=0}^{m-1} weights[i] * mode_derivs[i][k],  k < nstates.
// Templated on the scalar type so it serves both the adouble dynamics (and the
// outer-convexification wrapper) and plain double (unit tests). mode_derivs is
// an array of m pointers, each to an nstates-length per-mode derivative array.
template<class T>
inline void convex_combine(T* derivatives, int nstates, int m,
                           const T* weights, const T* const* mode_derivs)
{
    for (int k = 0; k < nstates; ++k) {
        T s = T(0);
        for (int i = 0; i < m; ++i)
            s += weights[i] * mode_derivs[i][k];
        derivatives[k] = s;
    }
}

// SOS1 residual sum_{i} weights[i], constrained to equal 1 as an equality path
// constraint in the convexified problem.
template<class T>
inline T sos1_residual(const T* weights, int m)
{
    T s = T(0);
    for (int i = 0; i < m; ++i) s += weights[i];
    return s;
}

// Expand a mode-index sequence (each entry in 0..m-1, length N) into an m x N
// one-hot weight-guess matrix: weights_out(i, n) = 1 if mode_sequence(n) == i,
// else 0. Out-of-range indices contribute an all-zero column.
inline void integer_guess_to_weights(int m, const RowVectorXi& mode_sequence,
                                     MatrixXd& weights_out)
{
    const int N = (int)mode_sequence.size();
    weights_out = MatrixXd::Zero(m, N);
    for (int n = 0; n < N; ++n) {
        const int mode = mode_sequence(n);
        if (mode >= 0 && mode < m)
            weights_out(mode, n) = 1.0;
    }
}

// Result of reconstructing a phase's integer control from the relaxed weights.
struct IntegerControlReconstruction {
    RowVectorXd control;          // 1 x P: rounded integer control value per interval
    RowVectorXi mode_index;       // 1 x P: index of the active admissible value per interval
    RowVectorXd interval_widths;  // 1 x P: interval widths (handy for forward simulation)
    double      integral_gap;     // accumulated sum-up-rounding gap
    int         n_switches;       // number of mode switches
};

// One-call reconstruction of a phase's integer control from a solved problem.
// After psopt(), the solution controls are in the weights layout; this extracts
// the M trailing weight rows, forms the interval widths from the phase time, and
// applies sum-up rounding using the admissible values recorded by
// declare_integer_control. Returns the rounded control together with the mode
// sequence, the interval widths, the integral gap, and the switch count. If the
// phase declares no integer control, the returned structure is empty.
inline IntegerControlReconstruction
reconstruct_integer_control(Sol& solution, Prob& problem, int iphase)
{
    const RowVectorXd values = problem.phases(iphase).integer_control.values;
    const int M = static_cast<int>(values.size());

    IntegerControlReconstruction r;
    r.integral_gap = 0.0;
    r.n_switches   = 0;
    if (M <= 0) return r;   // no integer control declared on this phase

    MatrixXd u = solution.get_controls_in_phase(iphase);   // weights layout
    MatrixXd t = solution.get_time_in_phase(iphase);
    const int N = static_cast<int>(t.cols());
    const int P = N - 1;

    MatrixXd    W = u.bottomRows(M).leftCols(P);   // M x P: weights at interval left nodes
    RowVectorXd h(P);
    for (int i = 0; i < P; ++i) h(i) = t(0, i + 1) - t(0, i);

    r.interval_widths = h;
    sum_up_rounding(W, h, values, r.mode_index, r.control, r.integral_gap, r.n_switches);
    return r;
}

#endif // PSOPT_INTEGER_CONTROLS_H
