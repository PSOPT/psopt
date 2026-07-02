#include "gtest/gtest.h"
#include <Eigen/Dense>
#include "integer_controls.h"

using Eigen::MatrixXd;
using Eigen::RowVectorXd;
using Eigen::RowVectorXi;

// convex_combine forms the weighted sum of the per-mode derivative arrays.
TEST(IntegerControls, ConvexCombineWeightedSum)
{
    double f0[3] = {1.0, 2.0, 3.0};
    double f1[3] = {10.0, 20.0, 30.0};
    const double* md[2] = {f0, f1};
    double w[2] = {0.3, 0.7};

    double d[3];
    convex_combine<double>(d, 3, 2, w, md);

    EXPECT_NEAR(d[0], 0.3*1.0 + 0.7*10.0, 1e-12);
    EXPECT_NEAR(d[1], 0.3*2.0 + 0.7*20.0, 1e-12);
    EXPECT_NEAR(d[2], 0.3*3.0 + 0.7*30.0, 1e-12);
}

// With weights summing to one and identical modes, the result reproduces the
// common value (a partition-of-unity sanity check), for m = 3.
TEST(IntegerControls, ConvexCombinePartitionOfUnity)
{
    double g0[2] = {5.0, -1.0};
    double g1[2] = {5.0, -1.0};
    double g2[2] = {5.0, -1.0};
    const double* md[3] = {g0, g1, g2};
    double w[3] = {0.2, 0.3, 0.5};

    double d[2];
    convex_combine<double>(d, 2, 3, w, md);

    EXPECT_NEAR(d[0], 5.0, 1e-12);
    EXPECT_NEAR(d[1], -1.0, 1e-12);
}

// sos1_residual returns the sum of the weights.
TEST(IntegerControls, Sos1Residual)
{
    double w2[2] = {0.3, 0.7};
    double w3[3] = {0.2, 0.3, 0.5};
    EXPECT_NEAR(sos1_residual<double>(w2, 2), 1.0, 1e-12);
    EXPECT_NEAR(sos1_residual<double>(w3, 3), 1.0, 1e-12);
}

// integer_guess_to_weights builds a one-hot m x N weight matrix; each in-range
// column sums to one and activates exactly the indicated mode.
TEST(IntegerControls, GuessToWeightsOneHot)
{
    RowVectorXi seq(4); seq << 0, 1, 1, 0;
    MatrixXd W;
    integer_guess_to_weights(2, seq, W);

    ASSERT_EQ(W.rows(), 2);
    ASSERT_EQ(W.cols(), 4);
    EXPECT_DOUBLE_EQ(W(0,0), 1.0); EXPECT_DOUBLE_EQ(W(1,0), 0.0);
    EXPECT_DOUBLE_EQ(W(0,1), 0.0); EXPECT_DOUBLE_EQ(W(1,1), 1.0);
    EXPECT_DOUBLE_EQ(W(0,2), 0.0); EXPECT_DOUBLE_EQ(W(1,2), 1.0);
    EXPECT_DOUBLE_EQ(W(0,3), 1.0); EXPECT_DOUBLE_EQ(W(1,3), 0.0);
    for (int n = 0; n < 4; ++n) EXPECT_DOUBLE_EQ(W.col(n).sum(), 1.0);
}

// An out-of-range mode index yields an all-zero column.
TEST(IntegerControls, GuessToWeightsOutOfRange)
{
    RowVectorXi seq(2); seq << 0, 5;   // 5 is out of range for m = 2
    MatrixXd W;
    integer_guess_to_weights(2, seq, W);

    EXPECT_DOUBLE_EQ(W.col(0).sum(), 1.0);
    EXPECT_DOUBLE_EQ(W.col(1).sum(), 0.0);
}

// declare_integer_control records the control index and value set on the phase;
// an undeclared phase is dormant (empty value set, index -1).
TEST(IntegerControls, DeclareIntegerControl)
{
    Prob problem;
    problem.nphases   = 1;
    problem.nlinkages = 0;
    psopt_level1_setup(problem);

    // Dormant by default.
    EXPECT_EQ(problem.phases(1).integer_control.control_index, -1);
    EXPECT_EQ(problem.phases(1).integer_control.values.size(), 0);

    RowVectorXd values(2); values << -0.05236, 0.05236;
    declare_integer_control(problem, 1, 0, values);

    EXPECT_EQ(problem.phases(1).integer_control.control_index, 0);
    ASSERT_EQ(problem.phases(1).integer_control.values.size(), 2);
    EXPECT_DOUBLE_EQ(problem.phases(1).integer_control.values(0), -0.05236);
    EXPECT_DOUBLE_EQ(problem.phases(1).integer_control.values(1),  0.05236);
    EXPECT_GT(problem.phases(1).integer_control.values.size(), 0);   // now active
}
