#include "gtest/gtest.h"
#include <Eigen/Dense>
#include <vector>
#include "integer_parameters.h"

using Eigen::MatrixXd;
using Eigen::RowVectorXd;

namespace {

// integer_cartesian_advance enumerates the full mixed-radix product, least-
// significant digit first, and reports exactly prod(sizes) combinations.
TEST(IntegerParameters, CartesianAdvanceEnumeratesAllInOrder)
{
    std::vector<int> sizes = {2, 3};
    std::vector<int> idx   = {0, 0};

    std::vector<std::pair<int,int>> seen;
    do {
        seen.push_back(std::make_pair(idx[0], idx[1]));
    } while (integer_cartesian_advance(idx, sizes));

    ASSERT_EQ(seen.size(), 6u);
    // least-significant digit (idx[0]) varies fastest
    EXPECT_EQ(seen[0], std::make_pair(0,0));
    EXPECT_EQ(seen[1], std::make_pair(1,0));
    EXPECT_EQ(seen[2], std::make_pair(0,1));
    EXPECT_EQ(seen[3], std::make_pair(1,1));
    EXPECT_EQ(seen[4], std::make_pair(0,2));
    EXPECT_EQ(seen[5], std::make_pair(1,2));
    // wrapped back to all zeros
    EXPECT_EQ(idx[0], 0);
    EXPECT_EQ(idx[1], 0);
}

// A zero-length product is a single (empty) combination: advance returns false
// immediately.
TEST(IntegerParameters, CartesianAdvanceSingleEmptyCombination)
{
    std::vector<int> sizes;
    std::vector<int> idx;
    EXPECT_FALSE(integer_cartesian_advance(idx, sizes));
}

// nearest_admissible_index snaps to the closest value; ties resolve to the lower
// index; an empty set returns -1.
TEST(IntegerParameters, NearestAdmissibleIndex)
{
    RowVectorXd v(4); v << 0.0, 1.0, 2.0, 3.0;
    EXPECT_EQ(nearest_admissible_index(2.3, v), 2);
    EXPECT_EQ(nearest_admissible_index(-5.0, v), 0);
    EXPECT_EQ(nearest_admissible_index(100.0, v), 3);
    EXPECT_EQ(nearest_admissible_index(0.5, v), 0);   // tie -> lower index
    EXPECT_EQ(nearest_admissible_index(1.4999, v), 1);

    RowVectorXd empty;
    EXPECT_EQ(nearest_admissible_index(1.0, empty), -1);
}

// declare_integer_parameter appends records, and the combination count is the
// product of all declared set sizes across the phase.
TEST(IntegerParameters, DeclareAndCombinationCount)
{
    Prob problem;
    problem.nphases   = 1;
    problem.nlinkages = 0;
    psopt_level1_setup(problem);

    // no declarations yet -> single empty combination
    EXPECT_EQ(integer_parameter_combination_count(problem), 1L);

    RowVectorXd a(4); a << 0.0, 1.0, 2.0, 3.0;   // |a| = 4
    RowVectorXd b(3); b << 5.0, 6.0, 7.0;        // |b| = 3
    declare_integer_parameter(problem, 1, 0, a);
    declare_integer_parameter(problem, 1, 1, b);

    ASSERT_EQ(problem.phases(1).integer_parameters.size(), 2u);
    EXPECT_EQ(problem.phases(1).integer_parameters[0].parameter_index, 0);
    EXPECT_EQ(problem.phases(1).integer_parameters[1].parameter_index, 1);
    EXPECT_EQ(integer_parameter_combination_count(problem), 12L);   // 4 * 3
}

// The combination count saturates (rather than overflowing) for very large
// products, so callers can safely compare it against a threshold.
TEST(IntegerParameters, CombinationCountSaturates)
{
    Prob problem;
    problem.nphases   = 1;
    problem.nlinkages = 0;
    psopt_level1_setup(problem);

    RowVectorXd big(64);
    for (int i = 0; i < 64; ++i) big(i) = i;
    for (int k = 0; k < 7; ++k)            // 64^7 = 2^42 > 2^40 saturation bound
        declare_integer_parameter(problem, 1, k, big);

    EXPECT_EQ(integer_parameter_combination_count(problem),
              PSOPT_INTEGER_COMBINATIONS_SATURATE);
}

} // namespace
