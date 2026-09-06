//////////////////////////////////////////////////////////////////////////////
// test_t_quantile.cpp
//
// inverse_twotailed_t_cdf() multiplies every confidence interval PSOPT reports,
// and for the whole life of the code it returned about 1.96 whatever it was
// asked for: 1.9600 for one degree of freedom, where the answer is 12.706, and
// 1.9604 for thirty-nine, where the answer is 2.023.
//
// The cause was a single token. The last row of the table of quantiles read
//
//     Infinity,  1.282, 1.645, 1.960, 2.326, 2.576;
//
// and `Infinity` is not a number here: with `using namespace Eigen` in scope it
// resolves to Eigen::Infinity, the enum constant -1 that names the L-infinity
// norm. The degrees-of-freedom column therefore ended in -1 rather than in a
// large number, the interpolation took the segment running from 500 down to -1
// as the bracketing one for every query, and interpolating between (500, 1.965)
// and (-1, 1.960) gives almost exactly 1.96 for any argument in between.
//
// Nothing announced it. A quantile of 1.96 is the right-looking answer for a
// large sample, the intervals it produces are plausible, and the error is
// largest exactly where samples are smallest and intervals are quoted with the
// most confidence. The covariance test alongside this one was written against
// correct reference half-widths and still passed, because its tolerance had
// been set wide enough to absorb the transcription error of an 80 node mesh and
// that was also wide enough to absorb a systematic error of three per cent.
//
// The values below are the standard two-tailed t table and are exact entries of
// the table PSOPT carries, so no interpolation tolerance is needed for them.
//////////////////////////////////////////////////////////////////////////////

#include "gtest/gtest.h"
#include <psopt.h>

TEST(TQuantile, ExactTableEntries)
{
    // 95 per cent two-tailed, i.e. t^{0.975}_{ndf}
    struct { int ndf; double t; } ref95[] = {
        {  1, 12.706 }, {  2,  4.303 }, {  3,  3.182 }, {  4,  2.776 },
        {  5,  2.571 }, { 10,  2.228 }, { 18,  2.101 }, { 20,  2.086 },
        { 25,  2.060 }, { 30,  2.042 }, { 39,  2.023 }, { 40,  2.021 },
        { 50,  2.009 }, { 60,  2.000 }, {100,  1.984 }, {500,  1.965 },
    };
    for (auto& r : ref95)
        EXPECT_NEAR(inverse_twotailed_t_cdf(0.95, r.ndf), r.t, 1.0e-9)
            << "95 per cent quantile at " << r.ndf << " degrees of freedom";

    // The other confidence levels the function accepts.
    EXPECT_NEAR(inverse_twotailed_t_cdf(0.80,  1),  3.078, 1.0e-9);
    EXPECT_NEAR(inverse_twotailed_t_cdf(0.90,  1),  6.314, 1.0e-9);
    EXPECT_NEAR(inverse_twotailed_t_cdf(0.98,  1), 31.820, 1.0e-9);
    EXPECT_NEAR(inverse_twotailed_t_cdf(0.99,  1), 63.657, 1.0e-9);
    EXPECT_NEAR(inverse_twotailed_t_cdf(0.99, 39),  2.708, 1.0e-9);
}

TEST(TQuantile, DecreasesTowardsTheNormalQuantile)
{
    // The property that the defect destroyed: the quantile has to depend on the
    // sample size, decrease with it, and approach 1.96 only in the limit.
    double prev = inverse_twotailed_t_cdf(0.95, 1);
    EXPECT_GT(prev, 12.0) << "one degree of freedom must give a large quantile";
    for (int ndf = 2; ndf <= 500; ndf++) {
        double t = inverse_twotailed_t_cdf(0.95, ndf);
        EXPECT_LE(t, prev + 1.0e-12) << "not decreasing at " << ndf;
        EXPECT_GT(t, 1.959)          << "below the normal quantile at " << ndf;
        prev = t;
    }
    // Small samples must be far from the normal quantile: this is the assertion
    // that the old behaviour failed by a factor of six.
    EXPECT_GT(inverse_twotailed_t_cdf(0.95,  1)/1.96, 6.0);
    EXPECT_GT(inverse_twotailed_t_cdf(0.95,  5)/1.96, 1.3);
    EXPECT_GT(inverse_twotailed_t_cdf(0.95, 39)/1.96, 1.03);

    // The interpolation is done in 1/ndf, so the infinite-sample limit is
    // approached properly rather than being pinned at the value for 500.
    EXPECT_NEAR(inverse_twotailed_t_cdf(0.95, 2000000), 1.960, 1.0e-5);
    EXPECT_NEAR(inverse_twotailed_t_cdf(0.95, 1000),    1.9624, 2.0e-4);
    EXPECT_NEAR(inverse_twotailed_t_cdf(0.95, 385),     1.9663, 2.0e-4);
}

TEST(TQuantile, InterpolatesBetweenTableRows)
{
    // 45 lies between the rows for 44 (2.015) and 46 (2.013); the true value
    // is 2.01410. Interpolating in 1/ndf does not land on the midpoint of the
    // two rows, and is the more accurate for it.
    double t45 = inverse_twotailed_t_cdf(0.95, 45);
    EXPECT_LT(t45, 2.015);
    EXPECT_GT(t45, 2.013);
    EXPECT_NEAR(t45, 2.0141, 5.0e-4);

    // 39 degrees of freedom is the case the book's catalytic cracking example
    // reports, and the one the defect got wrong by three per cent.
    EXPECT_NEAR(inverse_twotailed_t_cdf(0.95, 39), 2.023, 1.0e-9);
}
