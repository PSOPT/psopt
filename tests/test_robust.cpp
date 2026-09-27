#include "gtest/gtest.h"
#include <Eigen/Dense>
#include <cmath>
#include <vector>
#include "robust.h"

using Eigen::MatrixXd;
using Eigen::RowVectorXd;

namespace {

RowVectorXd row2(double a, double b)
{
    RowVectorXd v(2); v << a, b; return v;
}

MatrixXd cov2(double sa, double sb, double rho)
{
    MatrixXd C(2, 2);
    C << sa*sa, rho*sa*sb, rho*sa*sb, sb*sb;
    return C;
}

//////////////////////////////////////////////////////////////////////////
//  Special functions
//////////////////////////////////////////////////////////////////////////

// The inverse normal CDF against quantiles everyone knows.
TEST(Robust, NormPpfKnownQuantiles)
{
    EXPECT_NEAR(robust_norm_ppf(0.5),      0.0,               1e-12);
    EXPECT_NEAR(robust_norm_ppf(0.975),    1.959963984540054, 1e-10);
    EXPECT_NEAR(robust_norm_ppf(0.995),    2.575829303548901, 1e-10);
    EXPECT_NEAR(robust_norm_ppf(0.0013499), -3.0,             1e-4);
    // Symmetry, which the two branches of the approximation must not break.
    for (double p = 0.001; p < 0.5; p += 0.037)
        EXPECT_NEAR(robust_norm_ppf(p), -robust_norm_ppf(1.0 - p), 1e-11);
}

// Coverage of a k-sigma set, which falls away with dimension in a way that is easy
// to get wrong -- the Python documentation of this same driver once reported the
// THREE-dimensional figure as the planar one. These are the numbers that settle it.
TEST(Robust, SetCoverageFallsAwayWithDimension)
{
    EXPECT_NEAR(robust_set_coverage(3.0, 1), 0.9973002039367398, 1e-12);
    EXPECT_NEAR(robust_set_coverage(3.0, 2), 0.9888910034617577, 1e-12);
    EXPECT_NEAR(robust_set_coverage(3.0, 3), 0.9707091134651118, 1e-10);
    EXPECT_NEAR(robust_set_coverage(2.5, 2), 0.9560630663765926, 1e-12);
    EXPECT_NEAR(robust_set_coverage(3.5, 2), 0.9978125088818172, 1e-12);
    // In two dimensions the closed form is 1 - exp(-k^2/2).
    for (double k = 0.5; k < 5.0; k += 0.25)
        EXPECT_NEAR(robust_set_coverage(k, 2), 1.0 - std::exp(-0.5*k*k), 1e-12);
    // Monotone in k, decreasing in n.
    EXPECT_LT(robust_set_coverage(3.0, 3), robust_set_coverage(3.0, 2));
    EXPECT_LT(robust_set_coverage(2.0, 2), robust_set_coverage(3.0, 2));
}

//////////////////////////////////////////////////////////////////////////
//  Scenario sets
//////////////////////////////////////////////////////////////////////////

// For one Gaussian parameter the unscented set IS the three-point Gauss-Hermite
// rule: {mu, mu +- sqrt(3) sigma} with weights {2/3, 1/6, 1/6}. The sigma points of
// the estimation literature and the quadrature nodes are the same three numbers,
// which is worth pinning down because the two literatures name them differently.
TEST(Robust, SigmaPointsAreGaussHermiteInOneDimension)
{
    RowVectorXd mu(1); mu << 0.5;
    MatrixXd    C(1, 1); C << 0.15*0.15;

    std::vector<RowVectorXd> pts;
    RowVectorXd              w;
    ASSERT_TRUE(robust_sigma_points(robust_gaussian(mu, C, 3.0), pts, w));

    ASSERT_EQ(pts.size(), 3u);
    EXPECT_NEAR(pts[0](0), 0.5, 1e-12);
    EXPECT_NEAR(pts[1](0), 0.5 + std::sqrt(3.0)*0.15, 1e-12);
    EXPECT_NEAR(pts[2](0), 0.5 - std::sqrt(3.0)*0.15, 1e-12);
    EXPECT_NEAR(w(0), 2.0/3.0, 1e-12);
    EXPECT_NEAR(w(1), 1.0/6.0, 1e-12);
    EXPECT_NEAR(w(2), 1.0/6.0, 1e-12);
}

// In two dimensions with a correlated covariance the set must reproduce both
// moments exactly. That is the whole claim the unscented rule makes, and a
// Cholesky used the wrong way round would still look plausible without it.
TEST(Robust, SigmaPointsReproduceMeanAndCovariance)
{
    const RowVectorXd mu = row2(1.0, -2.0);
    const MatrixXd    C  = cov2(0.2, 0.35, -0.35);

    std::vector<RowVectorXd> pts;
    RowVectorXd              w;
    ASSERT_TRUE(robust_sigma_points(robust_gaussian(mu, C, 2.5), pts, w));
    ASSERT_EQ(pts.size(), 5u);
    EXPECT_NEAR(w.sum(), 1.0, 1e-12);

    RowVectorXd m = RowVectorXd::Zero(2);
    for (size_t k = 0; k < pts.size(); ++k) m += w(0, (int) k)*pts[k];
    EXPECT_NEAR(m(0), mu(0), 1e-12);
    EXPECT_NEAR(m(1), mu(1), 1e-12);

    MatrixXd S = MatrixXd::Zero(2, 2);
    for (size_t k = 0; k < pts.size(); ++k) {
        const RowVectorXd d = pts[k] - mu;
        S += w(0, (int) k)*(d.transpose()*d);
    }
    for (int i = 0; i < 2; ++i)
        for (int j = 0; j < 2; ++j)
            EXPECT_NEAR(S(i, j), C(i, j), 1e-12);
}

// A uniform set matches the second moment too, which is why the spread uses
// (hi-lo)/sqrt(12) rather than the half-width.
TEST(Robust, SigmaPointsMatchUniformVariance)
{
    const RobustUncertainty U = robust_uniform(row2(0.0, -1.0), row2(1.0, 3.0));
    std::vector<RowVectorXd> pts;
    RowVectorXd              w;
    ASSERT_TRUE(robust_sigma_points(U, pts, w));

    RowVectorXd m = RowVectorXd::Zero(2);
    for (size_t k = 0; k < pts.size(); ++k) m += w(0, (int) k)*pts[k];
    EXPECT_NEAR(m(0), 0.5, 1e-12);
    EXPECT_NEAR(m(1), 1.0, 1e-12);

    double v0 = 0.0;
    for (size_t k = 0; k < pts.size(); ++k)
        v0 += w(0, (int) k)*(pts[k](0) - 0.5)*(pts[k](0) - 0.5);
    EXPECT_NEAR(v0, 1.0/12.0, 1e-12);
}

// Beyond three dimensions kappa = 3 - n turns the central weight negative. That is
// admissible in a quadrature and not in a scenario set, every member of which is a
// constraint, so the rule must refuse rather than hand back a negative weight.
TEST(Robust, SigmaPointsRefuseNegativeCentralWeight)
{
    RowVectorXd mu = RowVectorXd::Zero(4);
    MatrixXd    C  = MatrixXd::Identity(4, 4);
    std::vector<RowVectorXd> pts;
    RowVectorXd              w;
    EXPECT_FALSE(robust_sigma_points(robust_gaussian(mu, C, 3.0), pts, w));
    EXPECT_TRUE(pts.empty());
}

//////////////////////////////////////////////////////////////////////////
//  Set geometry
//////////////////////////////////////////////////////////////////////////

TEST(Robust, ContainsAndClipOnACorrelatedEllipse)
{
    const RobustUncertainty U =
        robust_gaussian(row2(1.0, -2.0), cov2(0.2, 0.35, -0.35), 2.5);

    EXPECT_TRUE(robust_contains(U, row2(1.0, -2.0)));
    EXPECT_NEAR(robust_mahalanobis(U, row2(1.0, -2.0)), 0.0, 1e-12);

    // A point far out along the ellipse's short direction is outside even though
    // each component is well within its own marginal interval -- which is the
    // whole reason the correlation may not be discarded.
    const RowVectorXd off = row2(1.0 + 0.45, -2.0 + 0.75);
    EXPECT_GT(robust_mahalanobis(U, off), 2.5);
    EXPECT_FALSE(robust_contains(U, off));

    const RowVectorXd back = robust_clip(U, off);
    EXPECT_TRUE(robust_contains(U, back));
    EXPECT_NEAR(robust_mahalanobis(U, back), 2.5, 1e-9);
    // Clipping leaves an interior point alone.
    const RowVectorXd inside = row2(1.05, -2.05);
    ASSERT_TRUE(robust_contains(U, inside));
    EXPECT_NEAR((robust_clip(U, inside) - inside).norm(), 0.0, 1e-14);
}

TEST(Robust, BoundaryPointsLieOnTheBoundary)
{
    const RobustUncertainty U =
        robust_gaussian(row2(1.0, -2.0), cov2(0.2, 0.35, -0.35), 2.5);
    std::vector<RowVectorXd> b;
    robust_boundary_points(U, b);

    ASSERT_EQ(b.size(), 4u);
    for (size_t k = 0; k < b.size(); ++k) {
        EXPECT_NEAR(robust_mahalanobis(U, b[k]), 2.5, 1e-9);
        EXPECT_TRUE(robust_contains(U, b[k]));
    }
}

TEST(Robust, BoundaryPointsOfABoxAreItsCorners)
{
    std::vector<RowVectorXd> b;
    robust_boundary_points(robust_uniform(row2(-1.0, 0.0), row2(2.0, 5.0)), b);
    ASSERT_EQ(b.size(), 4u);
    double xsum = 0.0, ysum = 0.0;
    for (size_t k = 0; k < b.size(); ++k) { xsum += b[k](0); ysum += b[k](1); }
    EXPECT_NEAR(xsum, 2.0*(-1.0 + 2.0), 1e-12);
    EXPECT_NEAR(ysum, 2.0*(0.0 + 5.0), 1e-12);
}

// The low-discrepancy seeding must stay inside the set and must REACH its
// boundary: the worst parameter very often sits there, and a sequence that
// arrives only by accident would miss exactly the point a certificate is about.
TEST(Robust, LowDiscrepancyFillsTheSetAndReachesItsEdge)
{
    const RobustUncertainty U =
        robust_gaussian(row2(1.0, -2.0), cov2(0.2, 0.35, -0.35), 2.5);
    std::vector<RowVectorXd> pts;
    robust_low_discrepancy(U, 200, 0, pts);

    ASSERT_EQ(pts.size(), 200u);
    double furthest = 0.0;
    for (size_t k = 0; k < pts.size(); ++k) {
        EXPECT_TRUE(robust_contains(U, pts[k]));
        furthest = std::max(furthest, robust_mahalanobis(U, pts[k]));
    }
    EXPECT_NEAR(furthest, 2.5, 1e-9);
}

//////////////////////////////////////////////////////////////////////////
//  The worst-case oracle
//////////////////////////////////////////////////////////////////////////

// A smooth interior maximum the search must find to several figures.
double quadratic_bump(const RowVectorXd& t, void* /*user_data*/)
{
    const double a = t(0) - 1.10, b = t(1) + 1.80;
    return 1.0 - (a*a + 0.5*b*b);
}

TEST(Robust, WorstCaseFindsAnInteriorMaximum)
{
    const RobustUncertainty U =
        robust_gaussian(row2(1.0, -2.0), cov2(0.2, 0.35, -0.35), 2.5);
    RowVectorXd  at;
    const double v = robust_worst_case(U, &quadratic_bump, 0, 128, 3, 1u, at);

    EXPECT_NEAR(at(0), 1.10, 2e-3);
    EXPECT_NEAR(at(1), -1.80, 2e-3);
    EXPECT_NEAR(v, 1.0, 1e-5);
}

// A maximum ON the boundary, which is the case that matters in practice and the
// one a purely interior sampler would miss.
double linear_ramp(const RowVectorXd& t, void* /*user_data*/)
{
    return t(0);
}

TEST(Robust, WorstCaseFindsABoundaryMaximum)
{
    const RobustUncertainty U =
        robust_gaussian(row2(1.0, -2.0), cov2(0.2, 0.35, -0.35), 2.5);
    RowVectorXd  at;
    const double v = robust_worst_case(U, &linear_ramp, 0, 128, 3, 1u, at);

    // The largest first component on the ellipse is mu0 + k*sqrt(C00).
    EXPECT_NEAR(v, 1.0 + 2.5*0.2, 1e-6);
    EXPECT_NEAR(robust_mahalanobis(U, at), 2.5, 1e-6);
}

// Over an explicit set the search is exhaustive rather than a search, and it must
// say so by evaluating each member exactly once.
TEST(Robust, WorstCaseOverAnExplicitSetIsExhaustive)
{
    std::vector<RowVectorXd> pts;
    pts.push_back(row2(0.0, 0.0));
    pts.push_back(row2(1.0, 0.0));
    pts.push_back(row2(2.0, 0.0));
    const RobustUncertainty U = robust_explicit(pts);

    RowVectorXd  at;
    const double v = robust_worst_case(U, &linear_ramp, 0, 128, 3, 1u, at);
    EXPECT_NEAR(v, 2.0, 1e-12);
    EXPECT_NEAR(at(0), 2.0, 1e-12);
    EXPECT_EQ(robust_worst_case_evaluations(U, 128, 3), 3L);
}

// The evaluation count is a promise a caller can budget against, so it must be
// exact rather than indicative.
TEST(Robust, WorstCaseEvaluationCountIsExact)
{
    const RobustUncertainty U =
        robust_gaussian(row2(1.0, -2.0), cov2(0.2, 0.35, -0.35), 2.5);

    struct Counter { static double f(const RowVectorXd&, void* n)
                     { ++*(long*) n; return 0.0; } };
    long n = 0;
    RowVectorXd at;
    (void) robust_worst_case(U, &Counter::f, &n, 64, 3, 1u, at);
    EXPECT_EQ(n, robust_worst_case_evaluations(U, 64, 3));
}

// Two searches with the same seed must agree exactly: a certificate that moved
// from run to run would not be a certificate.
TEST(Robust, WorstCaseIsReproducible)
{
    const RobustUncertainty U =
        robust_gaussian(row2(1.0, -2.0), cov2(0.2, 0.35, -0.35), 2.5);
    RowVectorXd a, b;
    const double va = robust_worst_case(U, &quadratic_bump, 0, 64, 3, 7u, a);
    const double vb = robust_worst_case(U, &quadratic_bump, 0, 64, 3, 7u, b);
    EXPECT_DOUBLE_EQ(va, vb);
    EXPECT_NEAR((a - b).norm(), 0.0, 0.0);
}

}  // namespace
