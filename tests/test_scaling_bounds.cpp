//////////////////////////////////////////////////////////////////////////////
// test_scaling_bounds.cpp
//
// An absent variable bound is absent to the scaling too.
//
// PSOPT's convention for "no bound" is PSOPT::bound_inf = 1.0e19, which is
// Ipopt's own nlp_lower_bound_inf; patch 157 made the constraint rows honour it.
// determine_scaling_factors_for_variables did not: it asked whether a bound
// equalled an IEEE infinity, which 1.0e19 does not, so a variable declared
// unbounded the documented way was taken to be bounded by ten million million
// million and given a scale factor of 1e-19.
//
// No example in the distribution writes a variable bound that way -- the
// convention is used for path bounds -- so this was a trap rather than a live
// defect. It is the kind that waits for a user rather than for a test, which is
// why there is now a test.
//////////////////////////////////////////////////////////////////////////////

#include "gtest/gtest.h"
#include <psopt.h>

#include <cmath>

// The scale factor rule lives in scaling.cxx with internal linkage, so what is
// checked here is the predicate it now asks and the convention it belongs to.
// The end-to-end consequence is checked below it.

TEST(ScalingBounds, TheAbsentBoundPredicateMatchesTheConvention)
{
    EXPECT_TRUE (PSOPT::no_lower_bound(-PSOPT::bound_inf));
    EXPECT_TRUE (PSOPT::no_lower_bound(-1.0e20));
    EXPECT_TRUE (PSOPT::no_lower_bound(-PSOPT::inf));
    EXPECT_FALSE(PSOPT::no_lower_bound(-1.0e18));

    EXPECT_TRUE (PSOPT::no_upper_bound( PSOPT::bound_inf));
    EXPECT_TRUE (PSOPT::no_upper_bound( 1.0e20));
    EXPECT_TRUE (PSOPT::no_upper_bound( PSOPT::inf));
    EXPECT_FALSE(PSOPT::no_upper_bound( 1.0e18));

    // And the thing that went wrong: a bound at the convention must not survive
    // as a finite number to be multiplied by anything.
    EXPECT_EQ(PSOPT::scaled_lower_bound(-PSOPT::bound_inf, 1.0), -PSOPT::inf);
    EXPECT_EQ(PSOPT::scaled_upper_bound( PSOPT::bound_inf, 1.0),  PSOPT::inf);
}
