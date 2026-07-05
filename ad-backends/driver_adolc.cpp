// ADOL-C backend (PSOPT's native AD): tape once, then jacobian()/gradient().
#include "common.hpp"
#include <adolc/adolc.h>

int main(int argc, char** argv) {
    int n = (argc > 1) ? std::atoi(argv[1]) : 20;
    int reps = (argc > 2) ? std::atoi(argv[2]) : 2000;
    auto x = test_point(n);
    vector<double> J(n * n), g(n);

    // Tape the DAE (tag 0) and the cost (tag 1).
    { trace_on(0);
      vector<adouble> ax(n), ar(n);
      for (int i = 0; i < n; ++i) ax[i] <<= x[i];
      ocp_dae(ax.data(), ar.data(), n);
      double dummy; for (int i = 0; i < n; ++i) ar[i] >>= dummy;
      trace_off(); }
    { trace_on(1);
      vector<adouble> ax(n); adouble ac;
      for (int i = 0; i < n; ++i) ax[i] <<= x[i];
      ac = ocp_cost(ax.data(), n);
      double dummy; ac >>= dummy;
      trace_off(); }

    // Jacobian driver wants double**; use a row-pointer view over J.
    vector<double*> Jrows(n); for (int i = 0; i < n; ++i) Jrows[i] = &J[i * n];

    double us = time_us([&]{
        jacobian(0, n, n, x.data(), Jrows.data());
        gradient(1, n, x.data(), g.data());
    }, reps);

    return report("ADOL-C", n, J, g, us);
}
