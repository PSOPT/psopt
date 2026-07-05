// CppADCodeGen backend: record with CppAD, generate C for the Jacobian, compile
// it to a dynamic library at runtime, then call the compiled derivative. The
// timed section is the *compiled* Jacobian eval (codegen/compile is one-time).
#include "common.hpp"
#include <cppad/cg.hpp>

using namespace CppAD;
using namespace CppAD::cg;

int main(int argc, char** argv) {
    int n = (argc > 1) ? std::atoi(argv[1]) : 20;
    int reps = (argc > 2) ? std::atoi(argv[2]) : 2000;
    typedef CG<double> CGD;
    typedef AD<CGD> ADCG;
    auto x = test_point(n);
    std::vector<double> J, g;

    // Build a combined model with outputs [ dae(0..n-1), cost(n) ].
    std::vector<ADCG> ax(n);
    for (int i = 0; i < n; ++i) ax[i] = x[i];
    Independent(ax);
    std::vector<ADCG> ay(n + 1);
    ocp_dae(ax.data(), ay.data(), n);
    ay[n] = ocp_cost(ax.data(), n);
    ADFun<CGD> fun(ax, ay);

    // Generate + compile the Jacobian.
    ModelCSourceGen<double> cgen(fun, "ocp");
    cgen.setCreateJacobian(true);
    ModelLibraryCSourceGen<double> libcgen(cgen);
    DynamicModelLibraryProcessor<double> proc(libcgen, "ocp_lib");
    GccCompiler<double> compiler("/opt/local/bin/gcc-mp-15");
    std::unique_ptr<DynamicLib<double>> dynlib = proc.createDynamicLibrary(compiler);
    std::unique_ptr<GenericModel<double>> model = dynlib->model("ocp");

    std::vector<double> Jfull;  // ((n+1) x n) row-major
    double us = time_us([&]{ Jfull = model->Jacobian(x); }, reps);

    // Split: first n rows -> DAE Jacobian; last row -> cost gradient.
    J.assign(Jfull.begin(), Jfull.begin() + n * n);
    g.assign(Jfull.begin() + n * n, Jfull.begin() + (n + 1) * n);

    return report("CppADCodeGen", n, J, g, us);
}
