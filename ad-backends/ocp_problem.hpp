#pragma once
// Representative PSOPT-style optimal-control kernel, templated on the scalar type
// so any operator-overloading AD tool (ADOL-C, CppAD, Sacado) — and, at the
// double level, Enzyme/CppADCodeGen — can differentiate the *same* source.
//
//  ocp_dae : R^n -> R^n   nonlinear coupled dynamics residual (the DAE)
//  ocp_cost: R^n -> R     scalar Lagrange cost
//
// The expressions use the transcendentals/products typical of PSOPT DAEs
// (sin/cos/exp, state coupling) so the derivative work is representative.
#include <cmath>

template <class T>
inline void ocp_dae(const T* x, T* r, int n) {
    using std::sin; using std::cos; using std::exp;
    for (int i = 0; i < n; ++i) {
        const T xi = x[i];
        const T xm = x[(i - 1 + n) % n];
        const T xp = x[(i + 1) % n];
        r[i] = sin(xi) + T(0.5) * xi * xp - T(0.3) * xm * xm
             + exp(-xi * xi) - T(0.1) * cos(xp);
    }
}

template <class T>
inline T ocp_cost(const T* x, int n) {
    using std::sin;
    T s = T(0);
    for (int i = 0; i < n; ++i)
        s += (x[i] - T(1)) * (x[i] - T(1)) + T(0.1) * sin(x[i] * x[(i + 1) % n]);
    return s;
}
