#include "gtest/gtest.h"
#include <Eigen/Dense>
#include <psopt.h>

typedef Eigen::Matrix<adouble, 3, 1> Vector3ad;

TEST(Utils, crossProduct) {
    Vector3ad a, b, c;
    a << 1, 2, 3;
    b << 1, 5, 7;
    
    cross(a.data(), b.data(), c.data());
    
    ASSERT_EQ(c[0], -1);
    ASSERT_EQ(c[1], -4);
    ASSERT_EQ(c[2], 3);
}
