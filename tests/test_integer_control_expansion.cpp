#include "gtest/gtest.h"
#include <Eigen/Dense>
#include "integer_controls.h"

using Eigen::MatrixXd;
using Eigen::RowVectorXd;

namespace {

// Build a two-control phase (one continuous, one integer) ready for expansion.
// integer_index selects which control (0 or 1) is the integer one; the other is
// continuous with bounds [-5,5]. Guess: continuous row = 1..N, integer row =
// alternating admissible values (10,20,...).
void build_two_control_problem(Prob& problem, int integer_index)
{
    const int N = 5;
    problem.nphases   = 1;
    problem.nlinkages = 0;
    psopt_level1_setup(problem);

    problem.phases(1).nstates   = 2;
    problem.phases(1).ncontrols = 2;
    problem.phases(1).nevents   = 0;
    problem.phases(1).npath     = 0;
    problem.phases(1).nodes     << N;

    Alg algorithm;
    psopt_level2_setup(problem, algorithm);

    RowVectorXd values(2); values << 10.0, 20.0;
    declare_integer_control(problem, 1, integer_index, values);

    const int cont = (integer_index == 0) ? 1 : 0;   // the continuous control index
    problem.phases(1).bounds.lower.controls(cont, 0)          = -5.0;
    problem.phases(1).bounds.upper.controls(cont, 0)          =  5.0;
    problem.phases(1).bounds.lower.controls(integer_index, 0) = 10.0;  // ignored
    problem.phases(1).bounds.upper.controls(integer_index, 0) = 20.0;  // ignored

    MatrixXd g(2, N);
    for (int n = 0; n < N; ++n) {
        g(cont, n)          = (double)(n + 1);              // continuous guess 1..N
        g(integer_index, n) = (n % 2 == 0) ? 10.0 : 20.0;  // -> modes 0,1,0,1,0
    }
    problem.phases(1).guess.controls = g;

    problem.dae            = nullptr;   // not solved here
    problem.integrand_cost = nullptr;
}

} // namespace

// Integer control at index 1: the continuous control (index 0) stays at internal
// index 0; two weights are appended; one SOS1 path row is added; the guess is
// remapped to one-hot weights.
TEST(IntegerControlExpansion, IntegerAtIndexOne)
{
    Prob problem;
    build_two_control_problem(problem, /*integer_index=*/1);
    const int N = 5;

    {
        IntegerControlExpansionGuard guard(problem);

        EXPECT_EQ(problem.phases(1).ncontrols, 3);   // (2-1)+2
        EXPECT_EQ(problem.phases(1).npath, 1);       // SOS1 appended

        // continuous control (was index 0) remains at internal index 0
        EXPECT_DOUBLE_EQ(problem.phases(1).bounds.lower.controls(0, 0), -5.0);
        EXPECT_DOUBLE_EQ(problem.phases(1).bounds.upper.controls(0, 0),  5.0);
        // weights appended at internal indices 1,2 with [0,1]
        for (int k = 1; k <= 2; ++k) {
            EXPECT_DOUBLE_EQ(problem.phases(1).bounds.lower.controls(k, 0), 0.0);
            EXPECT_DOUBLE_EQ(problem.phases(1).bounds.upper.controls(k, 0), 1.0);
        }
        // SOS1 equality path bound
        EXPECT_DOUBLE_EQ(problem.phases(1).bounds.lower.path(0, 0), 1.0);
        EXPECT_DOUBLE_EQ(problem.phases(1).bounds.upper.path(0, 0), 1.0);

        // guess: row 0 continuous 1..N; rows 1,2 one-hot from modes 0,1,0,1,0
        ASSERT_EQ(problem.phases(1).guess.controls.rows(), 3);
        ASSERT_EQ(problem.phases(1).guess.controls.cols(), N);
        for (int n = 0; n < N; ++n) {
            EXPECT_DOUBLE_EQ(problem.phases(1).guess.controls(0, n), (double)(n + 1));
            double w0 = (n % 2 == 0) ? 1.0 : 0.0;
            EXPECT_DOUBLE_EQ(problem.phases(1).guess.controls(1, n), w0);
            EXPECT_DOUBLE_EQ(problem.phases(1).guess.controls(2, n), 1.0 - w0);
        }
    } // guard destroyed -> restore

    EXPECT_EQ(problem.phases(1).ncontrols, 2);
    EXPECT_EQ(problem.phases(1).npath, 0);
    EXPECT_EQ(problem.phases(1).bounds.lower.controls.rows(), 2);
    EXPECT_EQ(problem.phases(1).guess.controls.rows(), 2);
}

// Integer control at index 0: the continuous control (user index 1) is remapped
// down to internal index 0 (the k>c => k-1 case), with weights appended after.
TEST(IntegerControlExpansion, IntegerAtIndexZeroRemapsContinuous)
{
    Prob problem;
    build_two_control_problem(problem, /*integer_index=*/0);

    {
        IntegerControlExpansionGuard guard(problem);

        EXPECT_EQ(problem.phases(1).ncontrols, 3);
        // continuous control (was user index 1) now at internal index 0
        EXPECT_DOUBLE_EQ(problem.phases(1).bounds.lower.controls(0, 0), -5.0);
        EXPECT_DOUBLE_EQ(problem.phases(1).bounds.upper.controls(0, 0),  5.0);
        // weights at internal indices 1,2
        for (int k = 1; k <= 2; ++k) {
            EXPECT_DOUBLE_EQ(problem.phases(1).bounds.lower.controls(k, 0), 0.0);
            EXPECT_DOUBLE_EQ(problem.phases(1).bounds.upper.controls(k, 0), 1.0);
        }
        // guess row 0 is the continuous guess (was user row 1) = 1..N
        for (int n = 0; n < 5; ++n)
            EXPECT_DOUBLE_EQ(problem.phases(1).guess.controls(0, n), (double)(n + 1));
    }

    EXPECT_EQ(problem.phases(1).ncontrols, 2);   // restored
}
