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
// A phase may declare any number of integer controls. Each declaration records a
// control index and its admissible value set on the phase; the expansion itself
// (the weight-controls, the SOS1 equality path constraint, and the convex-
// combination dynamics) is performed by the IntegerControlExpansionGuard at the
// top of psopt().
//
// With K integer controls of sizes M_1..M_K the convexification is over the
// CARTESIAN PRODUCT of the admissible sets, so the phase carries P = M_1*...*M_K
// weight-controls and one SOS1 constraint over them. The product is what makes the
// relaxation tight: a separate weight vector per control would enter the dynamics
// multiplicatively, the relaxation would no longer be convex in the weights, and
// the equality between the relaxed optimum and the integer infimum would be lost.
// Two binary controls therefore cost four weights per node, three cost eight, and
// the cost grows as the product -- which is the price of an exact relaxation and
// the reason to keep the admissible sets small.
//
// See also include/sum_up_rounding.h for reconstruction of the integer control
// from the relaxed weights.
// ---------------------------------------------------------------------------

#include "psopt.h"
#include "sum_up_rounding.h"

// Declare an integer control of a phase. control_index is the index of the control
// in the user's control layout; values holds its M admissible discrete values. Call
// once per integer control; declarations accumulate, and re-declaring the same
// control index replaces the earlier entry rather than adding a second one. The
// order of declaration fixes the order in which the controls are reported by
// reconstruct_integer_controls.
inline void declare_integer_control(Prob& problem, int iphase,
                                    int control_index, const RowVectorXd& values)
{
    std::vector<IntegerControl>& v = problem.phases(iphase).integer_controls;
    for (size_t k = 0; k < v.size(); ++k) {
        if (v[k].control_index == control_index) { v[k].values = values; return; }
    }
    IntegerControl ic;
    ic.control_index = control_index;
    ic.values        = values;
    v.push_back(ic);
}

// Number of product modes P = M_1 * ... * M_K for a phase, and the decoding of a
// product index into the per-control mode indices. The first declared control is
// the most slowly varying digit, so mode 0 is "every control at its first
// admissible value".
inline int integer_control_product_size(const std::vector<IntegerControl>& v)
{
    int P = 1;
    for (size_t k = 0; k < v.size(); ++k) P *= (int) v[k].values.size();
    return (v.empty() ? 0 : P);
}

inline void integer_control_decode(const std::vector<IntegerControl>& v, int p,
                                   std::vector<int>& mode_of_control)
{
    const int K = (int) v.size();
    mode_of_control.assign(K, 0);
    for (int k = K - 1; k >= 0; --k) {
        const int Mk = (int) v[k].values.size();
        mode_of_control[k] = p % Mk;
        p /= Mk;
    }
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

// One-call reconstruction of a phase's integer controls from a solved problem.
// After psopt(), the solution controls are in the weights layout; this extracts the
// P trailing weight rows (P the product of the admissible set sizes), forms the
// interval widths from the phase time, applies sum-up rounding over the product
// modes, and then decodes each chosen product mode into a value for each declared
// control. Returns one reconstruction per declared integer control, in declaration
// order. The integral gap and the switch count are properties of the product mode
// sequence and are therefore the same in every entry; n_switches counts the
// instants at which any control changes.
inline std::vector<IntegerControlReconstruction>
reconstruct_integer_controls(Sol& solution, Prob& problem, int iphase)
{
    const std::vector<IntegerControl>& v = problem.phases(iphase).integer_controls;
    const int K = (int) v.size();
    std::vector<IntegerControlReconstruction> out;
    if (K <= 0) return out;                 // no integer control declared on this phase

    const int P = integer_control_product_size(v);

    MatrixXd u = solution.get_controls_in_phase(iphase);   // weights layout
    MatrixXd t = solution.get_time_in_phase(iphase);
    const int N  = static_cast<int>(t.cols());
    const int Ni = N - 1;

    MatrixXd    W = u.bottomRows(P).leftCols(Ni);   // P x Ni: weights at interval left nodes
    RowVectorXd h(Ni);
    for (int i = 0; i < Ni; ++i) h(i) = t(0, i + 1) - t(0, i);

    // Round over the product modes: the "values" are the product indices themselves,
    // which are then decoded into one value per declared control.
    RowVectorXd product_ids(P);
    for (int p = 0; p < P; ++p) product_ids(p) = (double) p;

    RowVectorXi product_mode;
    RowVectorXd product_value;
    double gap = 0.0;
    int    nsw = 0;
    sum_up_rounding(W, h, product_ids, product_mode, product_value, gap, nsw);

    out.resize(K);
    std::vector<int> modes;
    for (int k = 0; k < K; ++k) {
        out[k].control.resize(Ni);
        out[k].mode_index.resize(Ni);
        out[k].interval_widths = h;
        out[k].integral_gap    = gap;
        out[k].n_switches      = nsw;
    }
    for (int i = 0; i < Ni; ++i) {
        integer_control_decode(v, product_mode(i), modes);
        for (int k = 0; k < K; ++k) {
            out[k].mode_index(i) = modes[k];
            out[k].control(i)    = v[k].values(modes[k]);
        }
    }
    return out;
}

// Convenience overload for the common case of a single integer control, and for
// picking one control out of several. `which` is the declaration ordinal. Returns
// an empty structure if the phase declares no integer control, or if `which` is out
// of range.
inline IntegerControlReconstruction
reconstruct_integer_control(Sol& solution, Prob& problem, int iphase, int which = 0)
{
    std::vector<IntegerControlReconstruction> all =
        reconstruct_integer_controls(solution, problem, iphase);
    if (which < 0 || which >= (int) all.size()) {
        IntegerControlReconstruction empty;
        empty.integral_gap = 0.0;
        empty.n_switches   = 0;
        return empty;
    }
    return all[which];
}

#endif // PSOPT_INTEGER_CONTROLS_H
