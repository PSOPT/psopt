#pragma once
// Shared reference derivatives (central finite differences), validation, and timing.
#include "ocp_problem.hpp"
#include <vector>
#include <cstdio>
#include <cmath>
#include <chrono>
#include <string>

using std::vector;

// Deterministic test point (no RNG — reproducible).
inline vector<double> test_point(int n) {
    vector<double> x(n);
    for (int i = 0; i < n; ++i) x[i] = 0.3 + 0.1 * std::sin(0.7 * i) + 0.01 * i;
    return x;
}

// Central-difference reference Jacobian (row-major n*n) of ocp_dae.
inline vector<double> fd_jacobian(const vector<double>& x) {
    int n = x.size();
    vector<double> J(n * n, 0.0), xp = x, rp(n), rm(n);
    const double h = 1e-6;
    for (int j = 0; j < n; ++j) {
        xp = x; xp[j] += h; ocp_dae(xp.data(), rp.data(), n);
        xp = x; xp[j] -= h; ocp_dae(xp.data(), rm.data(), n);
        for (int i = 0; i < n; ++i) J[i * n + j] = (rp[i] - rm[i]) / (2 * h);
    }
    return J;
}

// Central-difference reference gradient of ocp_cost.
inline vector<double> fd_gradient(const vector<double>& x) {
    int n = x.size();
    vector<double> g(n), xp;
    const double h = 1e-6;
    for (int j = 0; j < n; ++j) {
        xp = x; xp[j] += h; double fp = ocp_cost(xp.data(), n);
        xp = x; xp[j] -= h; double fm = ocp_cost(xp.data(), n);
        g[j] = (fp - fm) / (2 * h);
    }
    return g;
}

inline double max_abs_diff(const vector<double>& a, const vector<double>& b) {
    double m = 0; for (size_t i = 0; i < a.size(); ++i) m = std::max(m, std::fabs(a[i] - b[i]));
    return m;
}

// Report: backend name, max |J-Jref|, max |g-gref|, and per-derivative-eval time (us).
inline int report(const std::string& backend, int n,
                  const vector<double>& J, const vector<double>& g,
                  double us_per_eval) {
    auto Jref = fd_jacobian(test_point(n));
    auto gref = fd_gradient(test_point(n));
    double dJ = max_abs_diff(J, Jref), dg = max_abs_diff(g, gref);
    bool ok = (dJ < 1e-5) && (dg < 1e-5);
    std::printf("%-14s n=%-4d  maxErr(J)=%.2e  maxErr(grad)=%.2e  %.1f us/eval  [%s]\n",
                backend.c_str(), n, dJ, dg, us_per_eval, ok ? "OK" : "MISMATCH");
    return ok ? 0 : 1;
}

// time a callable over `reps`, return microseconds per rep
template <class F>
inline double time_us(F&& f, int reps) {
    auto t0 = std::chrono::high_resolution_clock::now();
    for (int r = 0; r < reps; ++r) f();
    auto t1 = std::chrono::high_resolution_clock::now();
    return std::chrono::duration<double, std::micro>(t1 - t0).count() / reps;
}
