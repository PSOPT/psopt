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

//////////////////////////////////////////////////////////////////////////
//  The polish rule
//////////////////////////////////////////////////////////////////////////

// The polish step re-solves the final scenario set from the caller's own guess, because
// the warm chain conditions the answer: on the two-link arm the chain finishes at
// t_f = 8.9713 and a cold solve of its own twelve scenarios reaches 7.6387 with no
// violation anywhere in the set. Which of the two is kept is this one function, and
// these are its four cases.
TEST(Robust, PolishPrefersTheCheaperCertifiedDesign)
{
    const double slack = 1.0e-3;
    // Both certify: take the cheaper.
    EXPECT_TRUE (robust_prefer_cold(7.6387, 0.0,      8.9713, 6.7e-05, slack));
    EXPECT_FALSE(robust_prefer_cold(9.5000, 0.0,      8.9713, 6.7e-05, slack));
    // A tie on the objective keeps the warm design, which is the one already reported.
    EXPECT_FALSE(robust_prefer_cold(8.9713, 0.0,      8.9713, 6.7e-05, slack));
}

// The case this rule exists for. A cold design that undercuts the warm one and violates
// the constraints somewhere in the set is not a design, and no size of saving buys it in.
TEST(Robust, PolishNeverTradesACertificateForACheaperObjective)
{
    const double slack = 1.0e-3;
    EXPECT_FALSE(robust_prefer_cold(8.6063, 7.717e-03, 8.9670, 1.418e-04, slack));
    EXPECT_FALSE(robust_prefer_cold(0.0010, 1.0e+02,   8.9670, 1.418e-04, slack));
    // Only the warm design fails: the cold one is taken even though it costs more.
    EXPECT_TRUE (robust_prefer_cold(9.9000, 0.0,       8.9670, 5.0e-02,   slack));
}

// With neither certified there is no certificate to protect, so the rule falls back to
// whichever design comes closer to having one. The objective is ignored here on purpose:
// a cheaper infeasible design is not progress.
TEST(Robust, PolishFallsBackToTheSmallerViolation)
{
    const double slack = 1.0e-3;
    EXPECT_TRUE (robust_prefer_cold(9.9, 2.0e-03, 8.0, 5.0e-02, slack));
    EXPECT_FALSE(robust_prefer_cold(8.0, 5.0e-02, 9.9, 2.0e-03, slack));
    // A violation exactly at the slack counts as certified, matching the loop's own
    // stopping test, which is `worst <= slack`.
    EXPECT_TRUE (robust_prefer_cold(7.0, slack,   8.0, slack,   slack));
}

//////////////////////////////////////////////////////////////////////////
//  The model interface
//////////////////////////////////////////////////////////////////////////

// The inward tightening, which is what lets the generation loop terminate. A two-sided
// bound gives up a fraction of its half width; a pinned bound is an equality and must be
// left exactly alone, since tightening it would make it infeasible rather than tight.
TEST(Robust, MarginsTightenOnlyTwoSidedBounds)
{
    RowVectorXd lo(4), hi(4), mg;
    int one_sided = -1;
    lo << 0.0,  -0.02,  1.0,  0.5;
    hi << 0.0,   0.02,  3.0,  0.5;
    robust_margins(lo, hi, 0.9, mg, one_sided);
    EXPECT_EQ(one_sided, 0);
    EXPECT_DOUBLE_EQ(mg(0), 0.0);                    // pinned
    EXPECT_NEAR(mg(1), 0.1*0.02, 1.0e-15);           // ball of 0.02 keeps nine tenths
    EXPECT_NEAR(mg(2), 0.1*1.0,  1.0e-15);           // half width 1
    EXPECT_DOUBLE_EQ(mg(3), 0.0);                    // pinned away from zero
}

// A one-sided bound has no half width to take a fraction of. It gets no margin and is
// counted, because a design left sitting on such a bound is exactly where the
// between-scenario overshoot appears and the loop will not converge without being told.
TEST(Robust, MarginsCountOneSidedBounds)
{
    RowVectorXd lo(3), hi(3), mg;
    int one_sided = 0;
    lo << -1.0e30, 0.0, -2.0;
    hi <<  1.0,    1.0e30, 2.0;
    robust_margins(lo, hi, 0.5, mg, one_sided);
    EXPECT_EQ(one_sided, 2);
    EXPECT_DOUBLE_EQ(mg(0), 0.0);
    EXPECT_DOUBLE_EQ(mg(1), 0.0);
    EXPECT_NEAR(mg(2), 0.5*2.0, 1.0e-15);
}

// A tighten of 1 asks for no margin at all, which is the setting to use when the events
// already carry an inward margin the user put there.
TEST(Robust, MarginsVanishAtTightenOne)
{
    RowVectorXd lo(2), hi(2), mg;
    int one_sided = 0;
    lo << -1.0, 0.0;
    hi <<  1.0, 4.0;
    robust_margins(lo, hi, 1.0, mg, one_sided);
    EXPECT_DOUBLE_EQ(mg(0), 0.0);
    EXPECT_DOUBLE_EQ(mg(1), 0.0);
}

//////////////////////////////////////////////////////////////////////////
//  Reading a design's control
//////////////////////////////////////////////////////////////////////////

// The reading the transcription means, and every route that reproduces a designed control
// has to use this one: the verifier, the warm start, and a verification integrator of the
// caller's own. The parabola needs a complete table, so a design that reports none is read
// as a chord whatever the algorithm said.
TEST(Robust, TheControlShapeAsksTheAlgorithmAndTheDesignTogether)
{
    Alg alg;
    MatrixXd uf = zeros(1, 5), tf = zeros(1, 5), none;

    alg.transcription_method = "automatic";
    EXPECT_EQ(robust_control_shape(alg, uf, tf), ROBUST_PARABOLA);
    EXPECT_EQ(robust_control_shape(alg, none, none), ROBUST_LINEAR);
    EXPECT_EQ(robust_control_shape(alg, uf, none), ROBUST_LINEAR);   // half a table is none

    alg.transcription_method = "multiple-shooting";
    alg.ms_control_parameterisation = "constant";
    EXPECT_EQ(robust_control_shape(alg, none, none), ROBUST_HELD);
    EXPECT_EQ(robust_control_shape(alg, uf, tf), ROBUST_PARABOLA);
    alg.ms_control_parameterisation = "linear";
    EXPECT_EQ(robust_control_shape(alg, none, none), ROBUST_LINEAR);
}

// At the nodes all three readings agree with the nodal table, which is the invariant the
// warm start rests on: it writes its guess at the nodes and must not move the control
// there. Between them they differ, and the parabola is the only one that sees the midpoint.
TEST(Robust, TheThreeReadingsAgreeAtTheNodesAndDifferBetweenThem)
{
    // One interval, nodes 0 and 1, with a midpoint well off the chord.
    MatrixXd un(1, 2), uf(1, 3);
    un << 0.0, 1.0;
    uf << 0.0, 0.9, 1.0;                 // the chord's midpoint would be 0.5
    double u = -1.0;

    for (int s = 0; s < 3; ++s) {
        const RobustControlShape shape = (s == 0) ? ROBUST_HELD
                                       : (s == 1) ? ROBUST_LINEAR : ROBUST_PARABOLA;
        robust_control_at(un, uf, shape, 0, 0.0, 1, &u);
        EXPECT_DOUBLE_EQ(u, 0.0) << "shape " << s << " moved node 0";
        if (shape == ROBUST_HELD) continue;             // held does not reach node 1
        robust_control_at(un, uf, shape, 0, 1.0, 1, &u);
        EXPECT_DOUBLE_EQ(u, 1.0) << "shape " << s << " moved node 1";
    }

    robust_control_at(un, uf, ROBUST_HELD, 0, 0.5, 1, &u);
    EXPECT_DOUBLE_EQ(u, 0.0);
    robust_control_at(un, uf, ROBUST_LINEAR, 0, 0.5, 1, &u);
    EXPECT_DOUBLE_EQ(u, 0.5);
    robust_control_at(un, uf, ROBUST_PARABOLA, 0, 0.5, 1, &u);
    EXPECT_DOUBLE_EQ(u, 0.9);            // the midpoint value itself, not the chord's
}

// The parabola is the Lagrange interpolant through the three stored values, so a control
// that IS a quadratic is reproduced exactly and the chord is not. Checked against the
// polynomial written out by hand, at points that are not the three nodes.
TEST(Robust, TheParabolaReproducesAQuadraticControlExactly)
{
    MatrixXd un(1, 2), uf(1, 3);
    // u(w) = 2 + 3w - 5w^2, whose values at 0, 1/2 and 1 are 2, 2.25 and 0.
    un << 2.0, 0.0;
    uf << 2.0, 2.25, 0.0;
    const double w[4] = { 0.1, 0.25, 0.7, 0.95 };
    for (int k = 0; k < 4; ++k) {
        double u = 0.0, chord = 0.0;
        robust_control_at(un, uf, ROBUST_PARABOLA, 0, w[k], 1, &u);
        robust_control_at(un, uf, ROBUST_LINEAR, 0, w[k], 1, &chord);
        EXPECT_NEAR(u, 2.0 + 3.0*w[k] - 5.0*w[k]*w[k], 1.0e-14);
        EXPECT_GT(fabs(u - chord), 1.0e-3);          // and the chord is not the same curve
    }
}

// Every row, not only the user's controls. A co-designed gain schedule rides in the
// trailing rows and a caller reading only the first nc would reproduce a design whose gain
// never varies, which is a different controller.
TEST(Robust, TheReadingCoversEveryRowOfTheTable)
{
    MatrixXd un(3, 2), uf(3, 3);
    un << 0.0, 1.0,
          1.0, 3.0,
         -2.0, 2.0;
    uf << 0.0, 0.5, 1.0,
          1.0, 2.0, 3.0,
         -2.0, 0.5, 2.0;
    double u[3] = { 0.0, 0.0, 0.0 };
    robust_control_at(un, uf, ROBUST_PARABOLA, 0, 0.5, 3, u);
    EXPECT_DOUBLE_EQ(u[0], 0.5);
    EXPECT_DOUBLE_EQ(u[1], 2.0);
    EXPECT_DOUBLE_EQ(u[2], 0.5);
}

// The derivation is symmetric and a bound's two ends often do not mean the same thing.
// An end given a margin of its own takes that one; an end left empty still derives its
// own from the half width, so one end can be overridden without disturbing the other.
TEST(Robust, MarginsTakePerEndOverridesAndDeriveTheRest)
{
    RowVectorXd lo(2), hi(2), ovl(2), ovu(0), ml, mu;
    int one_sided = -1;
    lo <<  0.0,  0.0;
    hi << 10.0,  4.0;
    ovl << 0.0, 0.25;
    robust_margins(lo, hi, 0.9, ovl, ovu, ml, mu, one_sided);
    EXPECT_EQ(one_sided, 0);
    EXPECT_DOUBLE_EQ(ml(0), 0.0);                    // asked for nothing on this end
    EXPECT_DOUBLE_EQ(ml(1), 0.25);                   // asked for this much
    EXPECT_NEAR(mu(0), 0.1*5.0, 1.0e-15);            // and both upper ends derived
    EXPECT_NEAR(mu(1), 0.1*2.0, 1.0e-15);
}

// A single entry stands for all of them, which is how a user says "a tenth of a unit on
// every upper end" without writing the length of the vector out.
TEST(Robust, MarginsBroadcastASingleOverrideEntry)
{
    RowVectorXd lo(3), hi(3), ovl(1), ovu(0), ml, mu;
    int one_sided = 0;
    lo << 0.0, -1.0, 2.0;
    hi << 1.0,  1.0, 6.0;
    ovl << 0.05;
    robust_margins(lo, hi, 0.9, ovl, ovu, ml, mu, one_sided);
    EXPECT_DOUBLE_EQ(ml(0), 0.05);
    EXPECT_DOUBLE_EQ(ml(1), 0.05);
    EXPECT_DOUBLE_EQ(ml(2), 0.05);
}

// The case the overrides exist for. A one-sided bound has no half width, so the
// derivation gives it nothing and the driver says so; an explicit margin is what it can
// be given instead, and once it has one there is nothing left to report.
TEST(Robust, AnExplicitMarginReachesAOneSidedBoundAndEndsTheReport)
{
    RowVectorXd lo(1), hi(1), none(0), ovu(1), ml, mu;
    int one_sided = -1;
    lo << -1.0e30;
    hi <<  2.0;

    robust_margins(lo, hi, 0.9, none, none, ml, mu, one_sided);
    EXPECT_EQ(one_sided, 1);                         // nothing to take a fraction of
    EXPECT_DOUBLE_EQ(mu(0), 0.0);

    ovu << 0.01;
    robust_margins(lo, hi, 0.9, none, ovu, ml, mu, one_sided);
    EXPECT_EQ(one_sided, 0);                         // the end that binds has a margin
    EXPECT_DOUBLE_EQ(mu(0), 0.01);
    EXPECT_DOUBLE_EQ(ml(0), 0.0);                    // and the vacuous end still none
}

// A pinned bound is an equality, and an override may not touch it either. A margin there
// does not make the constraint tight, it makes it empty, and the user who wrote a margin
// for a whole vector of constraints did not mean that for the one that happens to be an
// equality.
TEST(Robust, MarginsLeaveAPinnedBoundAloneWhateverTheOverrideSays)
{
    RowVectorXd lo(2), hi(2), ov(1), ml, mu;
    int one_sided = 0;
    lo << 0.5, 0.0;
    hi << 0.5, 1.0;
    ov << 0.2;
    robust_margins(lo, hi, 0.9, ov, ov, ml, mu, one_sided);
    EXPECT_DOUBLE_EQ(ml(0), 0.0);
    EXPECT_DOUBLE_EQ(mu(0), 0.0);
    EXPECT_DOUBLE_EQ(ml(1), 0.2);
    EXPECT_DOUBLE_EQ(mu(1), 0.2);
}

// Zero, one or one per constraint, and nothing else: a vector of some other length
// cannot be matched to the constraints without guessing which one each entry meant.
TEST(Robust, OverrideLengthsAreZeroOneOrOnePerConstraint)
{
    EXPECT_TRUE(robust_margin_override_is_sized(zeros(1, 0), 3));
    EXPECT_TRUE(robust_margin_override_is_sized(zeros(1, 1), 3));
    EXPECT_TRUE(robust_margin_override_is_sized(zeros(1, 3), 3));
    EXPECT_FALSE(robust_margin_override_is_sized(zeros(1, 2), 3));
    EXPECT_FALSE(robust_margin_override_is_sized(zeros(1, 4), 3));
}

// The mechanism the whole model interface rests on: a nominal dae written once for the
// derivative tape, called numerically on plain doubles with the scenario as data. If this
// were not exact the library could not build a warm start out of the user's own equations
// and would have to ask for a second, hand-written copy of the physics.
namespace {
void spring_dae(adouble* d, adouble* path, adouble* x, adouble* u, adouble* p,
                adouble& t, const double* theta, int /*ntheta*/,
                adouble* /*xad*/, int /*iphase*/, Workspace* /*ws*/)
{
    const double k = theta[0], c = theta[1];
    d[0] = x[1];
    d[1] = -k*x[0] - c*x[1] + u[0] + p[0]*sin(t);
    path[0] = x[0]*x[0] + x[1]*x[1];
}
}  // namespace

TEST(Robust, DaeValueAgreesWithTheHandWrittenArithmetic)
{
    RobustModel model;
    model.nstates = 2; model.ncontrols = 1; model.npath = 1; model.nparameters = 1;
    model.dae = &spring_dae;

    const double th[2] = { 3.5, 0.25 };
    const double x[2]  = { 0.7, -1.3 }, u[1] = { 0.4 }, p[1] = { 2.0 }, t = 1.1;
    double d[2], g[1];
    robust_dae_value(model, th, 2, x, u, p, t, d, g);

    EXPECT_DOUBLE_EQ(d[0], x[1]);
    EXPECT_DOUBLE_EQ(d[1], -th[0]*x[0] - th[1]*x[1] + u[0] + p[0]*std::sin(t));
    EXPECT_DOUBLE_EQ(g[0], x[0]*x[0] + x[1]*x[1]);
}

// The scenario reaches the equations as data, so two scenarios give two different plants
// from one function. This is the augmentation, in the small.
TEST(Robust, DaeValueSeparatesTheScenarios)
{
    RobustModel model;
    model.nstates = 2; model.ncontrols = 1; model.npath = 1; model.nparameters = 1;
    model.dae = &spring_dae;

    const double x[2] = { 1.0, 0.0 }, u[1] = { 0.0 }, p[1] = { 0.0 };
    const double soft[2] = { 1.0, 0.0 }, stiff[2] = { 9.0, 0.0 };
    double ds[2], dh[2], g[1];
    robust_dae_value(model, soft,  2, x, u, p, 0.0, ds, g);
    robust_dae_value(model, stiff, 2, x, u, p, 0.0, dh, g);
    EXPECT_DOUBLE_EQ(ds[1], -1.0);
    EXPECT_DOUBLE_EQ(dh[1], -9.0);
}

// A path pointer is optional, because the warm start wants the derivatives and has no use
// for the path residual. Passing null must not write through it.
TEST(Robust, DaeValueToleratesANullPathPointer)
{
    RobustModel model;
    model.nstates = 2; model.ncontrols = 1; model.npath = 1; model.nparameters = 1;
    model.dae = &spring_dae;
    const double th[2] = { 1.0, 0.0 }, x[2] = { 1.0, 2.0 }, u[1] = { 0.0 }, p[1] = { 0.0 };
    double d[2] = { 0.0, 0.0 };
    robust_dae_value(model, th, 2, x, u, p, 0.0, d, 0);
    EXPECT_DOUBLE_EQ(d[0], 2.0);
}


//////////////////////////////////////////////////////////////////////////
//  The default verification integrator
//////////////////////////////////////////////////////////////////////////

namespace {
// x' = theta*u, with the state also reported as a path quantity so that the between-node
// sampling can be tested. theta scales the control, so one design serves several plants.
void ramp_dae(adouble* d, adouble* path, adouble* x, adouble* u, adouble* /*p*/,
              adouble& /*t*/, const double* theta, int /*ntheta*/,
              adouble* /*xad*/, int /*iphase*/, Workspace* /*ws*/)
{
    d[0]    = theta[0]*u[0];
    path[0] = x[0];
}

void terminal_event(adouble* e, adouble* /*xi*/, adouble* xf, adouble* /*p*/,
                    adouble& /*t0*/, adouble& /*tf*/,
                    const double* /*theta*/, int /*ntheta*/,
                    adouble* /*xad*/, int /*iphase*/, Workspace* /*ws*/)
{
    e[0] = xf[0];
}

// An event on the INITIAL state, for the test that the verifier hands the scenario's own
// starting point to the events and not the model's fixed vector.
void initial_state_event(adouble* e, adouble* xi, adouble* /*xf*/, adouble* /*p*/,
                         adouble& /*t0*/, adouble& /*tf*/,
                         const double* /*theta*/, int /*ntheta*/,
                         adouble* /*xad*/, int /*iphase*/, Workspace* /*ws*/)
{
    e[0] = xi[0];
}

RobustModel ramp_model(double e_lo, double e_hi, double p_lo, double p_hi)
{
    RobustModel m;
    m.nstates = 1; m.ncontrols = 1; m.nevents = 1; m.npath = 1;
    m.dae = &ramp_dae; m.events = &terminal_event;
    m.initial_state = zeros(1, 1);
    m.events_lower  = e_lo*ones(1, 1);  m.events_upper = e_hi*ones(1, 1);
    m.path_lower    = p_lo*ones(1, 1);  m.path_upper   = p_hi*ones(1, 1);
    m.verify_substeps = 16;
    return m;
}

RobustDesign ramp_design(double u0, double u1)
{
    RobustDesign d;
    d.time = zeros(1, 2);      d.time(0, 0) = 0.0; d.time(0, 1) = 1.0;
    d.controls = zeros(1, 2);  d.controls(0, 0) = u0; d.controls(0, 1) = u1;
    d.valid = true;
    return d;
}

Alg ms_alg(const char* parameterisation)
{
    Alg a;
    a.transcription_method        = "multiple-shooting";
    a.ms_control_parameterisation = parameterisation;
    return a;
}

RowVectorXd one(double v) { RowVectorXd r(1); r << v; return r; }
}  // namespace

// The reading of the control between the nodes, which is the one thing a verification
// integrator must not get wrong. One design, three readings, three answers computed by
// hand: with u going 0 to 2 over a unit interval and x' = u,
//
//   held      x(1) = 0                            (the control never leaves its first value)
//   linear    x(1) = (0 + 2)/2               = 1
//   parabola  x(1) = integral of -2w^2 + 4w  = 4/3
//
// the parabola being the one through (0, 0), (1/2, 3/2) and (1, 2). The terminal bound is
// [0, 1/2], so the excesses are 0, 1/2 and 5/6. A verifier reading that parabola as a chord
// would report 1/2 for a design that misses by 5/6, and the error is not in the safe
// direction either way, which is what this pins down.
TEST(Robust, VerifierReadsTheControlTheWayTheTranscriptionMeansIt)
{
    RobustModel model = ramp_model(0.0, 0.5, -1.0e30, 1.0e30);
    RobustDesign d = ramp_design(0.0, 2.0);

    Alg held = ms_alg("constant");
    EXPECT_NEAR(robust_model_violation(model, held, one(1.0), d), 0.0, 1.0e-10);

    Alg linear = ms_alg("linear");
    EXPECT_NEAR(robust_model_violation(model, linear, one(1.0), d), 0.5, 1.0e-10);

    // The same nodal table, now with the midpoint the design actually carries.
    d.controls_full = zeros(1, 3);
    d.controls_full(0, 0) = 0.0; d.controls_full(0, 1) = 1.5; d.controls_full(0, 2) = 2.0;
    d.time_full = zeros(1, 3);
    d.time_full(0, 0) = 0.0; d.time_full(0, 1) = 0.5; d.time_full(0, 2) = 1.0;
    Alg quad = ms_alg("quadratic");
    EXPECT_NEAR(robust_model_violation(model, quad, one(1.0), d), 4.0/3.0 - 0.5, 1.0e-9);
}

// A path constraint has to be sampled BETWEEN the nodes. Here the trajectory leaves and
// returns within a single interval: u going 2 to -2 gives x = 2w - 2w^2, whose peak is 1/2
// at the midpoint and which ends at 0. The terminal bound is met and the path bound of 1/5
// is exceeded by 3/10, which only an interior sample can see.
TEST(Robust, VerifierSamplesPathConstraintsBetweenTheNodes)
{
    RobustModel model = ramp_model(-0.01, 0.01, -0.2, 0.2);
    RobustDesign d = ramp_design(2.0, -2.0);
    Alg alg = ms_alg("linear");
    EXPECT_NEAR(robust_model_violation(model, alg, one(1.0), d), 0.3, 1.0e-9);
}

// The scales exist so that a metre and a radian are not added together, and they divide.
TEST(Robust, VerifierAppliesTheConstraintScales)
{
    RobustModel model = ramp_model(-0.01, 0.01, -0.2, 0.2);
    model.path_scale = 2.0*ones(1, 1);
    RobustDesign d = ramp_design(2.0, -2.0);
    Alg alg = ms_alg("linear");
    EXPECT_NEAR(robust_model_violation(model, alg, one(1.0), d), 0.15, 1.0e-9);
}

// The scenario reaches the verifier as data too, so one design is scored against several
// plants. Scaling theta scales the terminal state and so the excess beyond 1/2.
TEST(Robust, VerifierScoresOneDesignAgainstSeveralPlants)
{
    RobustModel model = ramp_model(0.0, 0.5, -1.0e30, 1.0e30);
    RobustDesign d = ramp_design(0.0, 2.0);
    Alg alg = ms_alg("linear");
    EXPECT_NEAR(robust_model_violation(model, alg, one(0.25), d), 0.0, 1.0e-10);
    EXPECT_NEAR(robust_model_violation(model, alg, one(1.0),  d), 0.5, 1.0e-10);
    EXPECT_NEAR(robust_model_violation(model, alg, one(2.0),  d), 1.5, 1.0e-10);
}

// A design whose integration leaves the finite numbers is infinitely bad, and must not come
// back as a quiet zero or as a NaN. NaN does not compare, so a diverged trajectory reported
// as one would be invisible to the worst-case search's running maximum, which is a failure
// this project has already had once.
TEST(Robust, VerifierReportsADivergedTrajectoryAsInfinite)
{
    RobustModel model = ramp_model(0.0, 0.5, -1.0e30, 1.0e30);
    RobustDesign d = ramp_design(0.0, 2.0);
    Alg alg = ms_alg("linear");
    EXPECT_TRUE(std::isinf(robust_model_violation(model, alg, one(1.0e308), d)));
}

// A quadratic reading needs a control history of the right width, and a mismatch has to be
// refused rather than read as whatever happens to be in the array.
TEST(Robust, VerifierRefusesAMalformedControlHistory)
{
    RobustModel model = ramp_model(0.0, 0.5, -1.0e30, 1.0e30);
    RobustDesign d = ramp_design(0.0, 2.0);
    d.controls_full = zeros(1, 2);      // should be 2N-1 = 3
    d.time_full     = zeros(1, 2);
    Alg quad = ms_alg("quadratic");
    EXPECT_TRUE(std::isinf(robust_model_violation(model, quad, one(1.0), d)));
}

// The events evaluated numerically, the companion of robust_dae_value.
TEST(Robust, EventsValueAgreesWithTheHandWrittenArithmetic)
{
    RobustModel model = ramp_model(0.0, 0.5, -1.0, 1.0);
    const double th = 1.0, xi = 0.0, xf = 0.37, par = 0.0;
    double e = -1.0;
    robust_events_value(model, &th, 1, &xi, &xf, &par, 0.0, 1.0, &e);
    EXPECT_DOUBLE_EQ(e, 0.37);
}


//////////////////////////////////////////////////////////////////////////
//  Risk measures
//////////////////////////////////////////////////////////////////////////

namespace {
RobustModel sized_model(int ns, int nc, int ne, int np, int npar)
{
    RobustModel m;
    m.nstates = ns; m.ncontrols = nc; m.nevents = ne; m.npath = np;
    m.nparameters = npar;
    return m;
}
}  // namespace

// The two measures that are weighted sums of per-scenario costs go straight into the
// integrand and the endpoint cost, so the augmented problem is the plain replication and
// nothing more.
TEST(Robust, SumRiskMeasuresAddNothingToTheProblem)
{
    const RobustModel m = sized_model(4, 2, 8, 1, 3);
    int nx = 0, nu = 0, ne = 0, np = 0, npar = 0;
    for (int r = 0; r < 2; ++r) {
        const RobustRisk risk = r ? ROBUST_EXPECTATION : ROBUST_NOMINAL;
        robust_augmented_sizes(m, 6, risk, nx, nu, ne, np, npar);
        EXPECT_EQ(nx,   4*6);
        EXPECT_EQ(nu,   2);
        EXPECT_EQ(ne,   8*6);
        EXPECT_EQ(np,   1*6);
        EXPECT_EQ(npar, 3);
    }
}

// Mean-variance needs each scenario's cost as a quantity of its own, because the variance
// of a Lagrange cost across scenarios is not the integral of anything. That is one extra
// state per scenario, and one extra pinned event per scenario to start it at zero.
TEST(Robust, MeanVarianceCarriesACostStatePerScenario)
{
    const RobustModel m = sized_model(4, 2, 8, 1, 3);
    int nx = 0, nu = 0, ne = 0, np = 0, npar = 0;
    robust_augmented_sizes(m, 6, ROBUST_MEAN_VARIANCE, nx, nu, ne, np, npar);
    EXPECT_EQ(nx,   4*6 + 6);
    EXPECT_EQ(ne,   8*6 + 6);
    EXPECT_EQ(npar, 3);          // no extra parameters: the objective is the endpoint
    EXPECT_EQ(np,   1*6);
}

// CVaR adds the Rockafellar-Uryasev device on top: one eta and one slack parameter per
// scenario, and one inequality row per scenario saying s_k >= J_k - eta. The slacks are
// static parameters rather than a smoothed hinge, so the constraints are exact.
TEST(Robust, CvarAddsEtaAndASlackPerScenario)
{
    const RobustModel m = sized_model(4, 2, 8, 1, 3);
    int nx = 0, nu = 0, ne = 0, np = 0, npar = 0;
    robust_augmented_sizes(m, 6, ROBUST_CVAR, nx, nu, ne, np, npar);
    EXPECT_EQ(nx,   4*6 + 6);
    EXPECT_EQ(ne,   8*6 + 6 + 6);
    EXPECT_EQ(npar, 3 + 6 + 1);
    EXPECT_EQ(np,   1*6);
}

// A problem with no events of its own still gets the rows its risk measure needs, which is
// the case where an off-by-one in the layout would otherwise go unnoticed.
TEST(Robust, RiskRowsSurviveAProblemWithNoEventsOfItsOwn)
{
    const RobustModel m = sized_model(2, 1, 0, 0, 0);
    int nx = 0, nu = 0, ne = 0, np = 0, npar = 0;
    robust_augmented_sizes(m, 5, ROBUST_CVAR, nx, nu, ne, np, npar);
    EXPECT_EQ(nx,   2*5 + 5);
    EXPECT_EQ(ne,   5 + 5);
    EXPECT_EQ(npar, 5 + 1);
    robust_augmented_sizes(m, 5, ROBUST_NOMINAL, nx, nu, ne, np, npar);
    EXPECT_EQ(ne,   0);
    EXPECT_EQ(npar, 0);
}

namespace {
adouble quadratic_integrand(adouble* x, adouble* u, adouble* p, adouble& t,
                            adouble* /*xad*/, int /*iphase*/, Workspace* /*ws*/)
{
    return 0.5*(x[0]*x[0] + u[0]*u[0]) + p[0]*t;
}
}  // namespace

// The integrand evaluated numerically, which is what seeds the cost states of a
// measure that carries them. Starting those at zero would start the objective at a value
// the guessed trajectory contradicts.
TEST(Robust, IntegrandValueAgreesWithTheHandWrittenArithmetic)
{
    RobustModel m = sized_model(1, 1, 0, 0, 1);
    m.integrand_cost = &quadratic_integrand;
    const double x = 0.7, u = -0.3, p = 2.0, t = 1.5;
    EXPECT_DOUBLE_EQ(robust_integrand_value(m, &x, &u, &p, t),
                     0.5*(x*x + u*u) + p*t);
}

// A model with no integrand contributes nothing, which the cost-state seeding relies on
// rather than guarding against separately.
TEST(Robust, IntegrandValueIsZeroWithoutAnIntegrand)
{
    const RobustModel m = sized_model(1, 1, 0, 0, 1);
    const double x = 0.7, u = -0.3, p = 2.0;
    EXPECT_DOUBLE_EQ(robust_integrand_value(m, &x, &u, &p, 1.5), 0.0);
}


//////////////////////////////////////////////////////////////////////////
//  The ancillary feedback
//////////////////////////////////////////////////////////////////////////

// Under feedback the realised control differs from scenario to scenario, so the control
// bounds are no longer the decision variable's bounds and become path rows: one per
// control per CORRECTED scenario. Scenario 0 is the reference and runs open loop, so it
// needs none, and a single-scenario problem has no deviation to correct at all.
TEST(Robust, FeedbackAddsPathRowsForEveryCorrectedScenario)
{
    RobustModel m = sized_model(4, 2, 8, 1, 3);
    m.feedback_kind = ROBUST_FEEDBACK_CONSTANT;
    int nx = 0, nu = 0, ne = 0, np = 0, npar = 0;
    robust_augmented_sizes(m, 6, ROBUST_NOMINAL, nx, nu, ne, np, npar);
    EXPECT_EQ(np,   1*6 + 2*5);
    EXPECT_EQ(nu,   2);               // u_bar is still the only control decision
    EXPECT_EQ(npar, 3);
    robust_augmented_sizes(m, 1, ROBUST_NOMINAL, nx, nu, ne, np, npar);
    EXPECT_EQ(np,   1);               // one scenario, which is its own reference
}

// A given gain costs no decision variables. A co-designed one does, and which kind it is
// decides whether they are static parameters or controls: a constant gain is one number
// per entry for the whole horizon, a schedule is one per entry per node, which is what a
// control already is.
TEST(Robust, CodesignedGainsTakeParametersOrControls)
{
    RobustModel m = sized_model(4, 2, 8, 1, 3);
    int nx = 0, nu = 0, ne = 0, np = 0, npar = 0;

    m.feedback_kind = ROBUST_FEEDBACK_SCHEDULED;
    robust_augmented_sizes(m, 6, ROBUST_NOMINAL, nx, nu, ne, np, npar);
    EXPECT_EQ(nu,   2);
    EXPECT_EQ(npar, 3);

    m.feedback_kind = ROBUST_FEEDBACK_CODESIGN;
    robust_augmented_sizes(m, 6, ROBUST_NOMINAL, nx, nu, ne, np, npar);
    EXPECT_EQ(nu,   2);
    EXPECT_EQ(npar, 3 + 2*4);

    m.feedback_kind = ROBUST_FEEDBACK_CODESIGN_SCHEDULE;
    robust_augmented_sizes(m, 6, ROBUST_NOMINAL, nx, nu, ne, np, npar);
    EXPECT_EQ(nu,   2 + 2*4);
    EXPECT_EQ(npar, 3);
}

namespace {
// x' = u + theta: a constant disturbance the control has to work against, with the state
// reported as the terminal event. One scalar state and one scalar control keep the closed
// loop solvable in closed form, which is the point: the answers below are analytic.
void disturbed_dae(adouble* d, adouble* /*path*/, adouble* /*x*/, adouble* u,
                   adouble* /*p*/, adouble& /*t*/, const double* theta, int /*ntheta*/,
                   adouble* /*xad*/, int /*iphase*/, Workspace* /*ws*/)
{
    d[0] = u[0] + theta[0];
}

void disturbed_event(adouble* e, adouble* /*xi*/, adouble* xf, adouble* /*p*/,
                     adouble& /*t0*/, adouble& /*tf*/, const double* /*theta*/,
                     int /*ntheta*/, adouble* /*xad*/, int /*iphase*/, Workspace* /*ws*/)
{
    e[0] = xf[0];
}

// The terminal bound is [0, 0.2] and the control bound is [-0.2, 0.2], so both the
// violation and the realised control's excess are numbers this file can predict.
RobustModel disturbed_model(double gain)
{
    RobustModel m;
    m.nstates = 1; m.ncontrols = 1; m.nevents = 1; m.npath = 0;
    m.dae = &disturbed_dae; m.events = &disturbed_event;
    m.initial_state  = zeros(1, 1);
    m.events_lower   = zeros(1, 1);     m.events_upper   = 0.2*ones(1, 1);
    m.controls_lower = -0.2*ones(1, 1); m.controls_upper = 0.2*ones(1, 1);
    m.verify_substeps = 64;
    if (gain != 0.0) {
        m.feedback_kind = ROBUST_FEEDBACK_CONSTANT;
        m.feedback_gain = gain*ones(1, 1);
    }
    return m;
}

RobustDesign flat_design(void)
{
    RobustDesign d;
    d.time     = zeros(1, 2);  d.time(0, 1) = 1.0;
    d.controls = zeros(1, 2);                      // u_bar identically zero
    d.valid    = true;
    return d;
}
}  // namespace

// The closed loop, against its own solution. With u_bar = 0, a reference at theta = 0 and
// a plant at theta = d, the deviation e = x - x_ref obeys e' = K e + d from e(0) = 0, so
// at K = -1 and t = 1
//
//     e(1) = d (1 - exp(-1)) = 0.6321205588 d,
//
// against d for the open loop. At d = 0.5 the terminal state is 0.3160602794 closed loop
// and 0.5 open loop, and against the bound [0, 0.2] those are excesses of 0.1160602794
// and 0.3. A verifier that integrated the open loop would report the second for a
// controller that achieves the first, which is not an error in the safe direction.
TEST(Robust, VerifierIntegratesTheClosedLoop)
{
    const double d = 0.5;
    const double e_closed = d*(1.0 - exp(-1.0));
    RobustDesign design = flat_design();
    Alg alg = ms_alg("linear");

    RobustModel open = disturbed_model(0.0);
    EXPECT_NEAR(robust_model_violation(open, alg, one(d), design), d - 0.2, 1.0e-9);

    RobustModel closed = disturbed_model(-1.0);
    const RowVectorXd ref = one(0.0);
    EXPECT_NEAR(robust_model_violation(closed, alg, one(d), design, &ref),
                e_closed - 0.2, 1.0e-9);
}

// The reference is an argument and not an assumption, and leaving it out with a gain in
// place makes the plant its own reference: the deviation is then identically zero, the
// correction with it, and what is verified is the open-loop controller. Pinned here
// because it is the one way to use this function that silently measures something other
// than the design that was solved.
TEST(Robust, VerifierWithoutAReferenceVerifiesTheOpenLoop)
{
    RobustModel closed = disturbed_model(-1.0);
    RobustDesign design = flat_design();
    Alg alg = ms_alg("linear");
    EXPECT_NEAR(robust_model_violation(closed, alg, one(0.5), design), 0.3, 1.0e-9);
}

// A gain of zero is the open loop, which is worth pinning because it is the path every
// non-feedback problem takes through the same code after this patch.
TEST(Robust, AZeroGainIsTheOpenLoop)
{
    RobustModel closed = disturbed_model(0.0);
    closed.feedback_kind = ROBUST_FEEDBACK_CONSTANT;
    closed.feedback_gain = zeros(1, 1);
    RobustDesign design = flat_design();
    Alg alg = ms_alg("linear");
    const RowVectorXd ref = one(0.0);
    double excess = -1.0;
    EXPECT_NEAR(robust_model_violation(closed, alg, one(0.5), design, &ref, &excess),
                0.3, 1.0e-9);
    EXPECT_DOUBLE_EQ(excess, 0.0);
}

// What the gain asks of the actuator, reported separately from the violation. Here
// u_bar = 0 and the realised control is K e, which grows monotonically to
// -0.3160602794 at t = 1 and so leaves the bound [-0.2, 0.2] by 0.1160602794. That
// number is not added to the violation: a design that saturates its actuator and a design
// that misses its target are different faults and want different remedies.
//
// The excess is read at the RUNGE-KUTTA STAGES and not only at the nodes, which is the
// whole point of measuring it in the verifier, so the last of them sits a little beyond
// the final state and the agreement here is to the stage's own accuracy and not to the
// integrator's.
TEST(Robust, VerifierReportsTheRealisedControlExcessSeparately)
{
    const double d = 0.5, e_closed = d*(1.0 - exp(-1.0));
    RobustModel closed = disturbed_model(-1.0);
    RobustDesign design = flat_design();
    Alg alg = ms_alg("linear");
    const RowVectorXd ref = one(0.0);
    double excess = -1.0;
    const double v = robust_model_violation(closed, alg, one(d), design, &ref, &excess);
    EXPECT_NEAR(v,      e_closed - 0.2, 1.0e-9);
    EXPECT_NEAR(excess, e_closed - 0.2, 1.0e-6);   // |K| = 1, so the two coincide here
}

// A correction small enough to stay inside the actuator's range leaves no excess at all,
// which separates the two measures: the same design misses its terminal bound and asks
// for nothing it has not got.
TEST(Robust, ASmallCorrectionLeavesNoControlExcess)
{
    RobustModel closed = disturbed_model(-0.2);      // |K| small, so |K e| < 0.2
    RobustDesign design = flat_design();
    Alg alg = ms_alg("linear");
    const RowVectorXd ref = one(0.0);
    double excess = -1.0;
    const double v = robust_model_violation(closed, alg, one(0.5), design, &ref, &excess);
    EXPECT_GT(v, 0.0);
    EXPECT_DOUBLE_EQ(excess, 0.0);
}

namespace {
// x0(theta) = theta, for a model whose uncertainty is in where the plant starts and not in
// how it behaves. Written against the ramp model, whose dae is x' = theta*u.
void theta_start(const double* theta, int /*ntheta*/, double* x0, void* /*user_data*/)
{
    x0[0] = theta[0];
}
}  // namespace

// An initial state that depends on the parameter, which the verifier has to start from.
// The ramp model has x' = theta*u, so with u going 0 to 2 over a unit interval read as a
// ramp, x(1) = x0 + theta. At theta = 0.25 that is 0.25 + 0.25 = 0.5 when the start moves
// with the parameter, against 0.25 when it does not, and the terminal bound [0, 0.2] turns
// those into excesses of 0.3 and 0.05. A verifier that kept the fixed vector would report
// the second for a design that achieves the first, which is the silent error this exists to
// prevent.
TEST(Robust, VerifierStartsWhereTheScenarioSays)
{
    RobustModel model = ramp_model(0.0, 0.2, -1.0e30, 1.0e30);
    RobustDesign d = ramp_design(0.0, 2.0);
    Alg linear = ms_alg("linear");

    EXPECT_NEAR(robust_model_violation(model, linear, one(0.25), d), 0.05, 1.0e-10);

    model.initial_state_fn = &theta_start;
    EXPECT_NEAR(robust_model_violation(model, linear, one(0.25), d), 0.30, 1.0e-10);
}

// The events are evaluated at that same starting state, not at the vector. Here the event
// is the FINAL state, so this pins the terminal reading; the initial-state reading is
// pinned by the next test, where the event is the initial state itself.
TEST(Robust, EventsSeeTheScenariosOwnInitialState)
{
    RobustModel model = ramp_model(-1.0e30, 1.0e30, -1.0e30, 1.0e30);
    model.events        = &initial_state_event;
    model.events_lower  = zeros(1, 1);
    model.events_upper  = zeros(1, 1);
    model.initial_state_fn = &theta_start;
    RobustDesign d = ramp_design(0.0, 0.0);      // no control, so nothing but x0 matters
    Alg linear = ms_alg("linear");
    EXPECT_NEAR(robust_model_violation(model, linear, one(0.25), d), 0.25, 1.0e-12);
}

// With no function the fixed vector is used, which is what every existing model does and
// what the driver's other tests rest on.
TEST(Robust, TheFixedInitialStateIsTheDefault)
{
    RobustModel model = ramp_model(0.0, 0.2, -1.0e30, 1.0e30);
    double x0 = -1.0;
    const double th = 0.25;
    robust_initial_state(model, &th, 1, &x0);
    EXPECT_DOUBLE_EQ(x0, 0.0);

    model.initial_state_fn = &theta_start;
    robust_initial_state(model, &th, 1, &x0);
    EXPECT_DOUBLE_EQ(x0, 0.25);
}

namespace {
// The smallest model that gets past the driver's other checks, so that what a test of the
// feedback validation refuses is the feedback and not something else.
RobustModel checkable_model(void)
{
    RobustModel m;
    m.nstates = 1; m.ncontrols = 1;
    m.dae = &disturbed_dae;
    m.initial_state = zeros(1, 1);
    m.states_lower  = -1.0*ones(1, 1);  m.states_upper = ones(1, 1);
    m.nodes.resize(1); m.nodes << 5;
    m.guess_states   = zeros(1, 5);
    m.guess_controls = zeros(1, 5);
    m.guess_time     = linspace(0.0, 1.0, 5);
    return m;
}

int try_solve_with(RobustModel& model, RobustSpec& spec)
{
    Prob problem; Alg algorithm; Sol solution;
    RowVectorXd mean(1);  mean << 0.0;
    MatrixXd    cov(1,1); cov  << 1.0;
    spec.uncertainty = robust_gaussian(mean, cov, 2.0);
    spec.verbose     = false;
    return psopt_solve_robust(solution, spec, model, problem, algorithm);
}

int try_solve(RobustModel& model)
{
    RobustSpec spec;
    return try_solve_with(model, spec);
}
}  // namespace

// A co-designed gain is a decision variable and nothing else bounds it. An unbounded one
// runs away into saturation, where the realised control is not the control that was
// designed, so the driver refuses the combination before a solve starts rather than
// returning a design whose certificate is about a different controller.
TEST(Robust, ACodesignedGainWithoutBoundsIsRefused)
{
    RobustModel m = checkable_model();
    m.feedback_kind = ROBUST_FEEDBACK_CODESIGN;
    EXPECT_THROW({ const int rc = try_solve(m); (void) rc; }, ErrorHandler);

    RobustModel s = checkable_model();
    s.feedback_kind = ROBUST_FEEDBACK_CODESIGN_SCHEDULE;
    EXPECT_THROW({ const int rc = try_solve(s); (void) rc; }, ErrorHandler);
}

// A given gain of the wrong shape is refused for the same reason, and a scheduled kind
// with nothing to call is refused rather than quietly designing the open loop.
TEST(Robust, AMalformedGivenGainIsRefused)
{
    RobustModel m = checkable_model();
    m.feedback_kind = ROBUST_FEEDBACK_CONSTANT;
    m.feedback_gain = ones(2, 2);                  // should be ncontrols x nstates
    EXPECT_THROW({ const int rc = try_solve(m); (void) rc; }, ErrorHandler);

    RobustModel s = checkable_model();
    s.feedback_kind = ROBUST_FEEDBACK_SCHEDULED;   // and no feedback_schedule
    EXPECT_THROW({ const int rc = try_solve(s); (void) rc; }, ErrorHandler);
}

namespace {
// A path-constrained robust problem whose answer can be written down.
//
//   maximise x(1)   subject to   xdot = u,   x(0) = 0,   u in [0, 10],
//                                g = x + u - theta <= 2,   |theta| <= 1.
//
// The constraint binds throughout, and it binds hardest at theta = -1, which is the
// scenario the generation loop has to find. There the design satisfies x + u <= 1 - m,
// m being the inward margin on the upper end of g, so u = (1 - m) - x and
//
//   x(t) = (1 - m)(1 - exp(-t)),       x(1) = (1 - m)(1 - 1/e).
//
// The bound on g is ONE-SIDED, so the symmetric derivation has no half width to take a
// fraction of and the margin can come only from path_margin_upper. That is what makes
// this a test of the overrides and not only of the path rows: without one the design sits
// exactly on the bound, which is where the between-scenario overshoot appears.
//
// Every scenario's state copy follows the same xdot = u, so the copies are identical and
// it is the worst scenario's bound that shapes the control. That is why the answer is
// available in closed form, and it is the only thing about this problem that is special.
void margin_dae(adouble* d, adouble* path, adouble* x, adouble* u, adouble* /*p*/,
                adouble& /*t*/, const double* theta, int ntheta,
                adouble* /*xad*/, int /*iphase*/, Workspace* /*ws*/)
{
    d[0]    = u[0];
    path[0] = x[0] + u[0] - ((ntheta > 0) ? theta[0] : 0.0);
}

void margin_event(adouble* e, adouble* xi, adouble* /*xf*/, adouble* /*p*/,
                  adouble& /*t0*/, adouble& /*tf*/, const double* /*theta*/,
                  int /*ntheta*/, adouble* /*xad*/, int /*iphase*/, Workspace* /*ws*/)
{
    e[0] = xi[0];
}

adouble margin_cost(adouble* /*xi*/, adouble* xf, adouble* /*p*/, adouble& /*t0*/,
                    adouble& /*tf*/, adouble* /*xad*/, int /*iphase*/,
                    Workspace* /*ws*/)
{
    return -xf[0];
}

RobustModel margin_model(void)
{
    RobustModel m;
    m.nstates = 1; m.ncontrols = 1; m.nevents = 1; m.npath = 1;
    m.dae = &margin_dae; m.events = &margin_event; m.endpoint_cost = &margin_cost;
    m.initial_state  = zeros(1, 1);
    m.events_lower   = zeros(1, 1);      m.events_upper   = zeros(1, 1);
    m.states_lower   = zeros(1, 1);      m.states_upper   = 2.0*ones(1, 1);
    m.controls_lower = zeros(1, 1);      m.controls_upper = 10.0*ones(1, 1);
    m.path_lower     = -1.0e30*ones(1, 1);
    m.path_upper     = 2.0*ones(1, 1);
    m.t0_lower = 0.0; m.t0_upper = 0.0;
    m.tf_lower = 1.0; m.tf_upper = 1.0;
    m.nodes.resize(1); m.nodes << 40;
    m.guess_states   = zeros(1, 40);
    m.guess_controls = ones(1, 40);
    m.guess_time     = linspace(0.0, 1.0, 40);
    m.verify_substeps = 16;
    return m;
}

// |theta| <= 1 as a one-dimensional Gaussian set of two standard deviations.
int solve_margin_model(RobustModel& model, RobustSpec& spec)
{
    Prob problem; Alg algorithm; Sol solution;
    RowVectorXd mean(1);  mean << 0.0;
    MatrixXd    cov(1,1); cov  << 0.25;
    spec.uncertainty   = robust_gaussian(mean, cov, 2.0);
    spec.verbose       = false;
    spec.slack         = 1.0e-4;
    spec.max_iterations = 6;
    return psopt_solve_robust(solution, spec, model, problem, algorithm);
}
}  // namespace

// A path-constrained RobustModel solved end to end to an answer written down in advance.
// Nothing else in this file solves one, and the gap mattered: a defect in the augmented
// path rows would have passed every row-count test here.
TEST(Robust, APathConstrainedModelReachesItsClosedFormAnswer)
{
    const double exact = 1.0 - exp(-1.0);

    RobustModel m = margin_model();
    m.path_margin_upper = 0.01*ones(1, 1);
    RobustSpec spec;
    const int rc = solve_margin_model(m, spec);
    ASSERT_EQ(rc, 0);
    EXPECT_TRUE(spec.converged);
    // The objective is -x(1), and x(1) = (1 - m)(1 - 1/e) at the worst scenario.
    EXPECT_NEAR(-spec.design.objective, 0.99*exact, 2.0e-3);
    // Which requires the edge of the set, and the starting rule does not contain it: for
    // one Gaussian the unscented points are {0, +-sqrt(3) sigma} and sigma is a half, so
    // the thinnest scenario it starts from is -0.866. The loop has to find the rest.
    double thinnest = 0.0;
    for (size_t i = 0; i < spec.scenarios.size(); ++i)
        thinnest = std::min(thinnest, spec.scenarios[i](0));
    EXPECT_LT(thinnest, -0.99);
}

// And the margin is what it says it is: the answer moves by exactly the amount the closed
// form says, which a margin written into the wrong rows, or applied to the wrong end, or
// applied once per scenario instead of once, would not do.
//
// Both margins here are small enough that the edge of the set still binds. Above
// m = 0.134 it stops binding, the starting rule's -0.866 being enough on its own to keep
// the whole set feasible, and the answer becomes (1.134 - m)(1 - 1/e) instead: measured
// 0.590350 at m = 0.2 against 0.590295 predicted. That regime is left out of the
// assertions because it rests on which points the starting rule happens to contain.
TEST(Robust, APathMarginMovesTheAnswerByWhatItPromises)
{
    const double exact = 1.0 - exp(-1.0);

    RobustModel a = margin_model();
    a.path_margin_upper = 0.02*ones(1, 1);
    RobustSpec sa;
    ASSERT_EQ(solve_margin_model(a, sa), 0);
    EXPECT_NEAR(-sa.design.objective, 0.98*exact, 2.0e-3);

    RobustModel b = margin_model();
    b.path_margin_upper = 0.10*ones(1, 1);
    RobustSpec sb;
    ASSERT_EQ(solve_margin_model(b, sb), 0);
    EXPECT_NEAR(-sb.design.objective, 0.90*exact, 2.0e-3);
}

// The refusals, each before a solve starts. A margin the driver cannot match to the
// constraints, one that would move a bound outwards, and a pair that leaves no box: all
// three reach IPOPT as a locally infeasible problem if they are not caught here, reported
// in the solver's language rather than in the model's.
TEST(Robust, AMarginOverrideOfTheWrongLengthIsRefused)
{
    RobustModel m = margin_model();
    m.path_margin_upper = zeros(1, 2);           // the phase declares one path constraint
    RobustSpec spec;
    EXPECT_THROW({ const int rc = solve_margin_model(m, spec); (void) rc; }, ErrorHandler);
}

TEST(Robust, ANegativeMarginIsRefused)
{
    RobustModel m = margin_model();
    m.path_margin_upper = -0.01*ones(1, 1);
    RobustSpec spec;
    EXPECT_THROW({ const int rc = solve_margin_model(m, spec); (void) rc; }, ErrorHandler);
}

TEST(Robust, MarginsThatCloseABoundsBoxAreRefused)
{
    RobustModel m = margin_model();
    m.path_lower        = zeros(1, 1);           // make it two-sided, [0, 2]
    m.path_margin_lower = 1.5*ones(1, 1);
    m.path_margin_upper = 1.5*ones(1, 1);        // which leaves [1.5, 0.5]
    RobustSpec spec;
    EXPECT_THROW({ const int rc = solve_margin_model(m, spec); (void) rc; }, ErrorHandler);
}

//////////////////////////////////////////////////////////////////////////
//  The verifier's recorded history
//////////////////////////////////////////////////////////////////////////

// Asking for the history must not change what the verifier reports. That is the whole basis
// on which a figure may be drawn from it, and it is easy to break: the recording needs the
// realised control at the final node, and forming one raises the control excess as a side
// effect.
TEST(Robust, AskingForTheVerifyHistoryChangesNothingItReports)
{
    const double d = 0.5;
    RobustDesign design = flat_design();
    Alg alg = ms_alg("linear");
    RobustModel closed = disturbed_model(-1.0);
    const RowVectorXd ref = one(0.0);

    double ex_plain = -1.0, ex_rec = -1.0;
    const double v_plain = robust_model_violation(closed, alg, one(d), design, &ref, &ex_plain);
    std::vector<RobustVerifyHistory> h;
    const double v_rec = robust_model_violation(closed, alg, one(d), design, &ref, &ex_rec, &h);

    EXPECT_EQ(v_plain, v_rec);
    EXPECT_EQ(ex_plain, ex_rec);
}

// The recorded trajectory is the one that was integrated, checked against a closed form.
// With u_bar = 0 and no gain the plant is x' = theta from x(0) = 0, so x(t) = theta t, and
// the recorded control is identically zero.
TEST(Robust, TheVerifyHistoryIsTheTrajectoryThatWasIntegrated)
{
    const double d = 0.5;
    RobustDesign design = flat_design();          // two nodes, t in [0, 1], u_bar = 0
    Alg alg = ms_alg("linear");
    RobustModel open = disturbed_model(0.0);      // verify_substeps = 64

    std::vector<RobustVerifyHistory> h;
    robust_model_violation(open, alg, one(d), design, 0, 0, &h);

    ASSERT_EQ((int) h.size(), 1);
    const int M = 1*64 + 1;                       // one interval, 64 steps, plus the end
    ASSERT_EQ((int) h[0].time.cols(), M);
    EXPECT_EQ((int) h[0].states.rows(), 1);
    EXPECT_EQ((int) h[0].controls.rows(), 1);
    EXPECT_EQ((int) h[0].path.size(), 0);         // this model declares no path constraints

    EXPECT_NEAR(h[0].time(0, 0), 0.0, 1.0e-12);
    EXPECT_NEAR(h[0].time(0, M-1), 1.0, 1.0e-12);
    EXPECT_NEAR(h[0].states(0, 0), 0.0, 1.0e-12);
    EXPECT_NEAR(h[0].states(0, M-1), d, 1.0e-12);
    for (int c = 0; c < M; ++c) {
        EXPECT_NEAR(h[0].states(0, c), d*h[0].time(0, c), 1.0e-12);
        EXPECT_NEAR(h[0].controls(0, c), 0.0, 1.0e-12);
    }
}

// Under a gain the recorded control is the REALISED one and not the designed one, which is
// the column a figure of a closed-loop design has to show. With u_bar = 0, K = -1, a
// reference at theta = 0 and a plant at theta = d, the deviation reaches d(1 - 1/e) at t = 1
// and the realised control there is minus that.
TEST(Robust, TheVerifyHistoryRecordsTheRealisedControl)
{
    const double d = 0.5;
    const double e_end = d*(1.0 - exp(-1.0));
    RobustDesign design = flat_design();
    Alg alg = ms_alg("linear");
    RobustModel closed = disturbed_model(-1.0);
    const RowVectorXd ref = one(0.0);

    std::vector<RobustVerifyHistory> h;
    robust_model_violation(closed, alg, one(d), design, &ref, 0, &h);

    ASSERT_EQ((int) h.size(), 1);
    const int M = (int) h[0].time.cols();
    EXPECT_NEAR(h[0].controls(0, 0), 0.0, 1.0e-12);          // no deviation yet
    EXPECT_NEAR(h[0].controls(0, M-1), -e_end, 1.0e-6);      // and the designed control is 0
}

// The claim a figure rests on: the worst recorded path excess IS the certificate. And the
// boundary the header states, pinned so that it cannot drift silently: the final column is
// recorded and not checked, so on a problem whose path value is still rising at the end the
// last column is worse than the number the verifier returned.
//
// The design is u = 10 held over t in [0, 1], so x(t) = 10t and g = x + u - theta rises from
// 10 - theta to 20 - theta against an upper bound of 2. With 40 nodes and 16 substeps the
// last CHECKED point is one substep short of t = 1.
TEST(Robust, TheVerifyHistoryPathAgreesWithTheCertificate)
{
    RobustModel m = margin_model();               // npath = 1, no path_scale, verify 16
    Alg alg = ms_alg("linear");
    const int N = 40;

    RobustDesign design;
    design.time     = linspace(0.0, 1.0, N);
    design.controls = 10.0*ones(1, N);
    design.valid    = true;

    std::vector<RobustVerifyHistory> h;
    const double v = robust_model_violation(m, alg, one(0.0), design, 0, 0, &h);

    ASSERT_EQ((int) h.size(), 1);
    const int M = (N - 1)*16 + 1;
    ASSERT_EQ((int) h[0].time.cols(), M);
    ASSERT_EQ((int) h[0].path.rows(), 1);

    // The worst excess over the CHECKED columns, which is every column but the last.
    double checked = 0.0;
    for (int c = 0; c < M - 1; ++c)
        checked = std::max(checked, h[0].path(0, c) - m.path_upper(0));
    EXPECT_NEAR(checked, v, 1.0e-9);

    // And the last column, recorded only, is beyond it on this problem.
    const double last = h[0].path(0, M-1) - m.path_upper(0);
    EXPECT_GT(last, v);
}

//////////////////////////////////////////////////////////////////////////
//  Several phases
//////////////////////////////////////////////////////////////////////////

namespace {
void two_phase_dae(adouble* d, adouble* /*path*/, adouble* x, adouble* u, adouble*,
                   adouble&, adouble*, int, Workspace*)
{ d[0] = x[0] + u[0]; }

void two_phase_events(adouble* e, adouble* xi, adouble* /*xf*/, adouble*, adouble&,
                      adouble&, adouble*, int, Workspace*)
{ e[0] = xi[0]; }

void two_phase_linkages(adouble*, adouble*, Workspace*) {}

// Two phases of one state and one control, the first with a pinned event and a fixed
// start time, the second with neither. Nothing is solved: the arithmetic under test is a
// count of variables and equalities, and it reads the bounds and the node counts.
void build_two_phase(Prob& problem, Alg& algorithm, int nlinkages)
{
    problem.name      = "dof over two phases";
    problem.nphases   = 2;
    problem.nlinkages = nlinkages;
    psopt_level1_setup(problem);
    for (int p = 1; p <= 2; ++p) {
        problem.phases(p).nstates   = 1;
        problem.phases(p).ncontrols = 1;
        problem.phases(p).nevents   = 1;
        problem.phases(p).npath     = 0;
        problem.phases(p).nodes     << 6;
    }
    psopt_level2_setup(problem, algorithm);
    for (int p = 1; p <= 2; ++p) {
        problem.phases(p).current_number_of_intervals = 5;
        problem.phases(p).bounds.lower.states(0)   = -1.0;
        problem.phases(p).bounds.upper.states(0)   =  1.0;
        problem.phases(p).bounds.lower.controls(0) = -1.0;
        problem.phases(p).bounds.upper.controls(0) =  1.0;
        problem.phases(p).bounds.lower.StartTime   = 0.0;
        problem.phases(p).bounds.upper.StartTime   = 0.0;
        problem.phases(p).bounds.lower.EndTime     = 1.0;
        problem.phases(p).bounds.upper.EndTime     = 1.0;
    }
    // Phase 1's event is an equality and phase 2's is not, so the two phases differ and
    // the tightest of them is the first.
    problem.phases(1).bounds.lower.events(0) = 0.0;
    problem.phases(1).bounds.upper.events(0) = 0.0;
    problem.phases(2).bounds.lower.events(0) = 0.0;
    problem.phases(2).bounds.upper.events(0) = 1.0;
    problem.dae = &two_phase_dae; problem.events = &two_phase_events;
    problem.linkages = &two_phase_linkages;
    algorithm.transcription_method = "multiple-shooting";
}
}  // namespace

// The whole problem's freedom is the phases' summed, less the linkage equalities, which
// belong to the problem and to no phase. Leaving them out would report freedom the problem
// does not have, which is the direction that lets a starved solve through.
TEST(Robust, DegreesOfFreedomSumOverPhasesAndChargeTheLinkages)
{
    Prob problem; Alg algorithm;
    build_two_phase(problem, algorithm, 3);
    const int d1 = robust_degrees_of_freedom(problem, 1, algorithm);
    const int d2 = robust_degrees_of_freedom(problem, 2, algorithm);
    EXPECT_EQ(robust_linkage_equalities(problem), 3);
    EXPECT_EQ(robust_problem_degrees_of_freedom(problem, algorithm), d1 + d2 - 3);
    EXPECT_LT(d1, d2);                       // the pinned event is phase 1's
    EXPECT_EQ(robust_tightest_phase(problem, algorithm), 1);
}

// A linkage the caller bounded as an inequality is not an equality and removes no freedom.
// Unset bounds mean every linkage is an equality, which is what auto_link builds.
TEST(Robust, LinkageEqualitiesReadTheBoundsWhenThereAreAny)
{
    Prob problem; Alg algorithm;
    build_two_phase(problem, algorithm, 3);
    EXPECT_EQ(robust_linkage_equalities(problem), 3);

    problem.bounds.lower.linkage = zeros(3, 1);
    problem.bounds.upper.linkage = zeros(3, 1);
    problem.bounds.upper.linkage(2) = 1.0;           // one of them is a range
    EXPECT_EQ(robust_linkage_equalities(problem), 2);
}

// A design of one phase reports one phase, whether or not anything filled the vector, so
// that code written against a single-phase design reads the same either way.
TEST(Robust, ADesignCountsItsPhases)
{
    RobustDesign d;
    EXPECT_EQ(d.nphases(), 0);
    d.valid = true;
    EXPECT_EQ(d.nphases(), 1);
    d.phase.resize(3);
    EXPECT_EQ(d.nphases(), 3);
}

// A model of several phases: the sizes are the phase's, and the parameters, which PSOPT
// gives to a phase and not to a problem, are replicated in every one.
TEST(Robust, PhaseSizesAreThePhasesOwn)
{
    RobustModel m = sized_model(4, 2, 8, 1, 3);
    m.later_phases.resize(1);
    m.later_phases[0].nstates   = 2;
    m.later_phases[0].ncontrols = 1;
    m.later_phases[0].nevents   = 3;
    m.later_phases[0].npath     = 0;

    int nx = 0, nu = 0, ne = 0, np = 0, npar = 0;
    robust_augmented_phase_sizes(m, m, 5, ROBUST_NOMINAL, nx, nu, ne, np, npar);
    EXPECT_EQ(nx, 4*5); EXPECT_EQ(nu, 2); EXPECT_EQ(ne, 8*5); EXPECT_EQ(np, 1*5);
    EXPECT_EQ(npar, 3);

    robust_augmented_phase_sizes(m, m.later_phases[0], 5, ROBUST_NOMINAL,
                                 nx, nu, ne, np, npar);
    EXPECT_EQ(nx, 2*5); EXPECT_EQ(nu, 1); EXPECT_EQ(ne, 3*5); EXPECT_EQ(np, 0);
    EXPECT_EQ(npar, 3);          // the same parameters, one copy per phase

    // The model's own call is the first phase, which is the whole of a single-phase
    // problem and is what every existing model asks for.
    int bx = 0, bu = 0, be = 0, bp = 0, bpar = 0;
    robust_augmented_sizes(m, 5, ROBUST_NOMINAL, bx, bu, be, bp, bpar);
    EXPECT_EQ(bx, 4*5); EXPECT_EQ(bu, 2); EXPECT_EQ(be, 8*5); EXPECT_EQ(bp, 1*5);
}

namespace {
// A two-phase model that gets past the other checks, so that what the refusals below refuse
// is the thing under test.
RobustModel two_phase_model(void)
{
    RobustModel m = checkable_model();
    m.later_phases.resize(1);
    RobustPhase& q = m.later_phases[0];
    q.nstates = 1; q.ncontrols = 1; q.nevents = 0; q.npath = 0;
    q.dae = &disturbed_dae;
    q.states_lower = -1.0*ones(1, 1); q.states_upper = ones(1, 1);
    q.nodes.resize(1); q.nodes << 5;
    q.guess_states   = zeros(1, 5);
    q.guess_controls = zeros(1, 5);
    q.guess_time     = linspace(1.0, 2.0, 5);
    return m;
}

void a_link(adouble* link, adouble* xf, adouble& /*tf*/, adouble* xi, adouble& /*t0*/,
            adouble*, const double*, int, void*)
{ link[0] = xi[0] - 2.0*xf[0]; }
}  // namespace

// The three things a second phase does not yet reach, each refused in as many words rather
// than assembled into something that looks right. They are increments, not nonsense.
TEST(Robust, WhatASecondPhaseDoesNotYetReachIsRefused)
{
    RobustModel mv = two_phase_model();
    mv.cost_lower = 0.0; mv.cost_upper = 10.0;
    mv.integrand_cost = &quadratic_integrand;
    RobustSpec s1; s1.risk = ROBUST_MEAN_VARIANCE;
    EXPECT_THROW({ const int rc = try_solve_with(mv, s1); (void) rc; }, ErrorHandler);

    RobustModel fb = two_phase_model();
    fb.feedback_kind = ROBUST_FEEDBACK_CONSTANT;
    fb.feedback_gain = zeros(1, 1);
    RobustSpec s2;
    EXPECT_THROW({ const int rc = try_solve_with(fb, s2); (void) rc; }, ErrorHandler);

    // A linkage of the caller's own with no numerical companion: the verifier would have
    // to cross the boundary by guessing.
    RobustModel lk = two_phase_model();
    lk.link = &a_link; lk.nlink = 1;
    RobustSpec s3;
    EXPECT_THROW({ const int rc = try_solve_with(lk, s3); (void) rc; }, ErrorHandler);
}

// And the same verifier says so on its own, without a solve: a violation of infinity for a
// design it cannot follow across the boundary.
TEST(Robust, TheVerifierRefusesALinkageItCannotCross)
{
    RobustModel m = two_phase_model();
    m.link = &a_link; m.nlink = 1;
    RobustDesign d = ramp_design(0.0, 2.0);
    d.phase.resize(2);
    d.phase[0] = d.phase[1] = RobustPhaseTrajectory();
    d.phase[0].time = d.time; d.phase[0].controls = d.controls;
    d.phase[1].time = d.time; d.phase[1].controls = d.controls;
    Alg alg = ms_alg("linear");
    EXPECT_TRUE(std::isinf(robust_model_violation(m, alg, one(0.25), d)));
}

}  // namespace
