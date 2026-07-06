// Independent native-CasADi solve of the Bryson-Denham problem (NOT via PSOPT),
// for cross-validation against PSOPT and Wolfram IPOPTLink.
//
//   min (1/2) integral_0^1 u^2 dt
//   s.t. x1' = x2,  x2' = u,  x1 <= 1/9
//        x1(0)=0, x2(0)=1, x1(1)=0, x2(1)=-1
// Known optimal cost = 4. Trapezoidal direct collocation via casadi::Opti.
#include <casadi/casadi.hpp>
#include <cstdio>
using namespace casadi;

int main(int argc, char** argv) {
    int N = (argc > 1) ? std::atoi(argv[1]) : 200;
    double dt = 1.0 / N;
    Opti opti;
    MX X = opti.variable(2, N + 1);   // states
    MX U = opti.variable(1, N + 1);   // control
    MX J = 0;
    for (int k = 0; k < N; ++k) {
        MX xk = X(Slice(), k), xk1 = X(Slice(), k + 1);
        MX fk  = vertcat(xk(1),  U(0, k));
        MX fk1 = vertcat(xk1(1), U(0, k + 1));
        opti.subject_to(xk1 == xk + dt / 2 * (fk + fk1));            // trapezoidal dynamics
        J += dt / 2 * (0.5 * U(0, k) * U(0, k) + 0.5 * U(0, k + 1) * U(0, k + 1));
    }
    opti.subject_to(X(0, Slice()) <= 1.0 / 9.0);                    // state path constraint
    opti.subject_to(X(Slice(), 0)  == DM(std::vector<double>{0, 1}));
    opti.subject_to(X(Slice(), N)  == DM(std::vector<double>{0, -1}));
    opti.minimize(J);
    Dict opts; opts["ipopt.print_level"] = 0; opts["print_time"] = false;
    opti.solver("ipopt", opts);
    try {
        OptiSol sol = opti.solve();
        std::printf("native-CasADi  Bryson-Denham  N=%d  objective = %.6f\n",
                    N, (double)sol.value(J));
    } catch (std::exception& e) {
        std::printf("native-CasADi FAILED: %s\n", e.what());
        return 1;
    }
    return 0;
}
