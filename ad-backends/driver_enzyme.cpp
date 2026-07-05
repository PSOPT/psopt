// Enzyme backend: LLVM-level AD. Enzyme differentiates the *compiled* double
// functions (no retemplating of the source). Forward mode gives Jacobian columns
// (seed one input), reverse mode gives the cost gradient. Built with the
// llvm-20 clang + ClangEnzyme-20 plugin.
#include "common.hpp"
#include <algorithm>

extern "C" {
    int enzyme_dup, enzyme_const, enzyme_dupnoneed;
    void   __enzyme_fwddiff(void*, ...);
    double __enzyme_autodiff(void*, ...);
}

// Clear call targets for Enzyme (avoid full inlining before the pass runs).
__attribute__((noinline)) static void   dae_wrap (const double* x, double* r, int n) { ocp_dae(x, r, n); }
__attribute__((noinline)) static double cost_wrap(const double* x, int n)            { return ocp_cost(x, n); }

int main(int argc, char** argv) {
    int n = (argc > 1) ? std::atoi(argv[1]) : 20;
    int reps = (argc > 2) ? std::atoi(argv[2]) : 2000;
    auto x = test_point(n);
    vector<double> J(n * n), g(n);

    double us = time_us([&]{
        vector<double> dx(n, 0.0), r(n, 0.0), dr(n, 0.0);
        for (int j = 0; j < n; ++j) {          // one forward sweep per input column
            std::fill(dx.begin(), dx.end(), 0.0); dx[j] = 1.0;
            std::fill(dr.begin(), dr.end(), 0.0);
            __enzyme_fwddiff((void*)dae_wrap,
                enzyme_dup,       x.data(), dx.data(),
                enzyme_dupnoneed, r.data(), dr.data(),
                enzyme_const,     n);
            for (int i = 0; i < n; ++i) J[i * n + j] = dr[i];   // column j
        }
        std::fill(g.begin(), g.end(), 0.0);    // reverse mode for the gradient
        __enzyme_autodiff((void*)cost_wrap,
            enzyme_dup,   x.data(), g.data(),
            enzyme_const, n);
    }, reps);

    return report("Enzyme", n, J, g, us);
}
