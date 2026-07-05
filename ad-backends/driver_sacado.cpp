// Sacado backend: forward-mode AD (Trilinos/Sacado, Sacado::Fad::DFad). Vector
// forward mode — one sweep seeds all n directions, giving the full Jacobian.
#include "common.hpp"
#include <Sacado_No_Kokkos.hpp>

int main(int argc, char** argv) {
    int n = (argc > 1) ? std::atoi(argv[1]) : 20;
    int reps = (argc > 2) ? std::atoi(argv[2]) : 2000;
    typedef Sacado::Fad::DFad<double> Fad;
    auto x = test_point(n);
    vector<double> J(n * n), g(n);

    double us = time_us([&]{
        // Seed: xf[i] carries an n-vector derivative with unit in slot i.
        vector<Fad> xf(n), rf(n);
        for (int i = 0; i < n; ++i) xf[i] = Fad(n, i, x[i]);
        ocp_dae(xf.data(), rf.data(), n);
        for (int i = 0; i < n; ++i)
            for (int j = 0; j < n; ++j) J[i * n + j] = rf[i].dx(j);
        Fad c = ocp_cost(xf.data(), n);
        for (int j = 0; j < n; ++j) g[j] = c.dx(j);
    }, reps);

    return report("Sacado", n, J, g, us);
}
