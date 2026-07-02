#include "gtest/gtest.h"
#include <Eigen/Dense>
#include "sum_up_rounding.h"

using Eigen::MatrixXd;
using Eigen::RowVectorXd;
using Eigen::RowVectorXi;

namespace {

// Deterministic, varied SOS1 weight matrix (columns sum to one), no RNG so the
// tests are fully reproducible.
MatrixXd make_sos1(int M, int N)
{
    MatrixXd w(M, N);
    for (int n = 0; n < N; ++n) {
        double s = 0.0;
        for (int k = 0; k < M; ++k) {
            w(k, n) = 1.0 + static_cast<double>(((k + 1) * (n + 3)) % 7);
            s += w(k, n);
        }
        w.col(n) /= s;
    }
    return w;
}

} // namespace

// Constant 1/2 binary on a uniform mesh: the rounding must alternate, switch on
// every interior boundary, and its accumulated gap must equal exactly half an
// interval.
TEST(SumUpRounding, ConstantHalfBinaryAlternates)
{
    const int N = 8; const double hval = 0.25;
    MatrixXd w = MatrixXd::Constant(2, N, 0.5);
    RowVectorXd h = RowVectorXd::Constant(N, hval);

    RowVectorXi idx; double gap; int sw;
    sum_up_rounding(w, h, idx, gap, sw);

    for (int n = 0; n < N; ++n) EXPECT_EQ(idx(n), n % 2);
    EXPECT_EQ(sw, N - 1);
    EXPECT_NEAR(gap, 0.5 * hval, 1e-12);
}

// An input that is already integer must be reproduced exactly, with zero gap.
TEST(SumUpRounding, AlreadyIntegerReproduced)
{
    const int N = 4;
    MatrixXd w(2, N);
    w << 1, 0, 1, 0,
         0, 1, 0, 1;
    RowVectorXd h = RowVectorXd::Constant(N, 0.5);

    RowVectorXi idx; double gap; int sw;
    sum_up_rounding(w, h, idx, gap, sw);

    for (int n = 0; n < N; ++n) EXPECT_EQ(idx(n), n % 2);
    EXPECT_NEAR(gap, 0.0, 1e-15);
}

// The accumulated gap is linear in the mesh spacing: doubling h doubles the gap.
// This is the numerical signature of the O(h) approximation property.
TEST(SumUpRounding, GapIsLinearInH)
{
    const int M = 3, N = 20;
    MatrixXd w = make_sos1(M, N);

    RowVectorXi idx; int sw; double g1, g2;
    sum_up_rounding(w, RowVectorXd::Constant(N, 0.1), idx, g1, sw);
    sum_up_rounding(w, RowVectorXd::Constant(N, 0.2), idx, g2, sw);

    EXPECT_GT(g1, 0.0);
    EXPECT_NEAR(g2, 2.0 * g1, 1e-12);
}

// The accumulated gap never exceeds one (maximum) interval width, for any number
// of modes. This is the tight discrete bound underlying the O(h) guarantee.
TEST(SumUpRounding, IntegralGapBoundedByHmax)
{
    const int N = 30; const double hval = 0.37;
    RowVectorXd h = RowVectorXd::Constant(N, hval);
    for (int M = 2; M <= 5; ++M) {
        MatrixXd w = make_sos1(M, N);
        RowVectorXi idx; double gap; int sw;
        sum_up_rounding(w, h, idx, gap, sw);
        EXPECT_LE(gap, hval + 1e-12) << "failed for M=" << M;
    }
}

// A non-uniform mesh genuinely enters the decision: the same weights on a
// non-uniform mesh round differently from the uniform case.
TEST(SumUpRounding, NonUniformMeshEntersDecision)
{
    const int N = 3;
    MatrixXd w = MatrixXd::Constant(2, N, 0.5);

    RowVectorXi idx_nu, idx_u; double gap; int sw;

    RowVectorXd h_nu(N); h_nu << 1.0, 0.1, 1.0;
    sum_up_rounding(w, h_nu, idx_nu, gap, sw);
    EXPECT_EQ(idx_nu(0), 0); EXPECT_EQ(idx_nu(1), 1); EXPECT_EQ(idx_nu(2), 1);

    RowVectorXd h_u = RowVectorXd::Constant(N, 1.0);
    sum_up_rounding(w, h_u, idx_u, gap, sw);
    EXPECT_EQ(idx_u(0), 0); EXPECT_EQ(idx_u(1), 1); EXPECT_EQ(idx_u(2), 0);

    EXPECT_NE(idx_nu(2), idx_u(2));
}

// Ties resolve to the lowest index, deterministically and reproducibly.
TEST(SumUpRounding, DeterministicTieBreak)
{
    const int M = 3, N = 6;
    MatrixXd w = MatrixXd::Constant(M, N, 1.0 / 3.0);
    RowVectorXd h = RowVectorXd::Constant(N, 1.0);

    RowVectorXi a, b; double ga, gb; int sa, sb;
    sum_up_rounding(w, h, a, ga, sa);
    sum_up_rounding(w, h, b, gb, sb);

    EXPECT_EQ(a(0), 0);          // first all-tie interval picks mode 0
    EXPECT_EQ(a, b);             // fully reproducible
    for (int n = 0; n < N; ++n) { EXPECT_GE(a(n), 0); EXPECT_LT(a(n), M); }
}

// The value-mapping overload assigns each interval the discrete value of its
// activated mode.
TEST(SumUpRounding, ValueMappingOverload)
{
    const int M = 3, N = 6;
    MatrixXd w = MatrixXd::Constant(M, N, 1.0 / 3.0);
    RowVectorXd h = RowVectorXd::Constant(N, 1.0);
    RowVectorXd values(M); values << 10.0, 20.0, 30.0;

    RowVectorXi idx; RowVectorXd rounded; double gap; int sw;
    sum_up_rounding(w, h, values, idx, rounded, gap, sw);

    for (int n = 0; n < N; ++n) EXPECT_DOUBLE_EQ(rounded(n), values(idx(n)));
}
