#include "gtest/gtest.h"
#include <Eigen/Dense>
#include "integer_controls.h"

using Eigen::MatrixXd;
using Eigen::RowVectorXd;

// reconstruct_integer_control extracts the trailing M weight rows from the
// solution (ignoring any leading continuous controls), builds the interval
// widths from the phase time, and applies sum-up rounding using the admissible
// values recorded by declare_integer_control. This test hand-builds a solution
// with one continuous control (row 0, to be ignored) and two mode-weight rows,
// so it also confirms the trailing-rows extraction.
TEST(ReconstructIntegerControl, WeightsToRoundedControl)
{
    Prob problem;
    problem.nphases   = 1;
    problem.nlinkages = 0;
    psopt_level1_setup(problem);

    RowVectorXd values(2); values << 0.0, 1.0;
    declare_integer_control(problem, 1, /*control_index=*/1, values);

    // Hand-built solution in the weights layout: 3 rows x 5 nodes.
    //   row 0: a continuous control (must be ignored by the helper)
    //   row 1: weight on mode 0
    //   row 2: weight on mode 1
    Sol solution;
    solution.problem  = &problem;
    solution.controls = new MatrixXd[1];
    solution.nodes    = new MatrixXd[1];

    MatrixXd C(3, 5);
    C.row(0) << 5.0, 5.0, 5.0, 5.0, 5.0;   // continuous, ignored
    C.row(1) << 0.0, 0.0, 1.0, 1.0, 1.0;   // mode 0 weight
    C.row(2) << 1.0, 1.0, 0.0, 0.0, 0.0;   // mode 1 weight
    solution.controls[0] = C;

    MatrixXd T(1, 5);
    T << 0.0, 1.0, 2.0, 3.0, 4.0;          // uniform width h = 1
    solution.nodes[0] = T;

    IntegerControlReconstruction rec = reconstruct_integer_control(solution, problem, 1);

    // Four intervals; the trailing two rows are the weights (mode1 active on the
    // first two intervals, mode0 on the last two).
    ASSERT_EQ(rec.control.size(), 4);
    ASSERT_EQ(rec.mode_index.size(), 4);
    ASSERT_EQ(rec.interval_widths.size(), 4);

    EXPECT_EQ(rec.mode_index(0), 1);
    EXPECT_EQ(rec.mode_index(1), 1);
    EXPECT_EQ(rec.mode_index(2), 0);
    EXPECT_EQ(rec.mode_index(3), 0);

    EXPECT_DOUBLE_EQ(rec.control(0), 1.0);
    EXPECT_DOUBLE_EQ(rec.control(1), 1.0);
    EXPECT_DOUBLE_EQ(rec.control(2), 0.0);
    EXPECT_DOUBLE_EQ(rec.control(3), 0.0);

    for (int i = 0; i < 4; ++i) EXPECT_DOUBLE_EQ(rec.interval_widths(i), 1.0);

    EXPECT_EQ(rec.n_switches, 1);
    EXPECT_GE(rec.integral_gap, 0.0);
    EXPECT_LE(rec.integral_gap, rec.interval_widths.maxCoeff());   // gap bound <= h_max
}

// An undeclared phase yields an empty reconstruction.
TEST(ReconstructIntegerControl, NoIntegerControlIsEmpty)
{
    Prob problem;
    problem.nphases   = 1;
    problem.nlinkages = 0;
    psopt_level1_setup(problem);   // no declare_integer_control

    Sol solution;
    solution.problem  = &problem;
    solution.controls = new MatrixXd[1];
    solution.nodes    = new MatrixXd[1];
    solution.controls[0] = MatrixXd::Zero(1, 5);
    solution.nodes[0]    = MatrixXd::Zero(1, 5);

    IntegerControlReconstruction rec = reconstruct_integer_control(solution, problem, 1);

    EXPECT_EQ(rec.control.size(), 0);
    EXPECT_EQ(rec.mode_index.size(), 0);
    EXPECT_EQ(rec.n_switches, 0);
}
