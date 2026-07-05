// CppAD backend: operator-overloading AD (coin-or/CppAD). Tape once, Jacobian().
#include "common.hpp"
#include <cppad/cppad.hpp>

int main(int argc, char** argv) {
    int n = (argc > 1) ? std::atoi(argv[1]) : 20;
    int reps = (argc > 2) ? std::atoi(argv[2]) : 2000;
    using CppAD::AD;
    auto x = test_point(n);
    vector<double> J, g;

    // Tape the DAE.
    vector<AD<double>> ax(n), ar(n);
    for (int i = 0; i < n; ++i) ax[i] = x[i];
    CppAD::Independent(ax);
    ocp_dae(ax.data(), ar.data(), n);
    CppAD::ADFun<double> fdae(ax, ar);

    // Tape the cost (scalar output).
    vector<AD<double>> ax2(n), ay(1);
    for (int i = 0; i < n; ++i) ax2[i] = x[i];
    CppAD::Independent(ax2);
    ay[0] = ocp_cost(ax2.data(), n);
    CppAD::ADFun<double> fcost(ax2, ay);

    double us = time_us([&]{
        J = fdae.Jacobian(x);   // length n*n, row-major
        g = fcost.Jacobian(x);  // length n
    }, reps);

    return report("CppAD", n, J, g, us);
}
