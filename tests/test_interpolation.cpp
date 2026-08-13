//////////////////////////////////////////////////////////////////////////////
// test_interpolation.cpp
//
// The table-lookup interpolation routines, and in particular the way they choose
// which piece of the interpolant applies.
//
// Every one of these routines takes an active variable x and has to locate the
// interval of the table that contains it. Locating it with an ordinary C++
// comparison works when the routine is called for its value, and fails silently
// when it is called while a derivative tape is being recorded: the comparison is
// resolved once, at the argument the tape was recorded with, and the tape then
// holds a single piece of the interpolant extended over the whole real line. The
// value and the derivative are both wrong as soon as the optimizer moves x into a
// different interval, the solver reports nothing unusual, and the answer is simply
// not the answer to the problem that was posed. The routines now select the piece
// with psopt_cond_lt(), which records the comparison itself.
//
// Two things therefore need testing, and the second is the one that matters:
//
//   * that the selection still returns the right value, checked against the
//     all-double implementations, which never had the problem; and
//   * that the selection survives taping, checked by solving a problem whose
//     optimum lies in a different interval from the one the tape is recorded in.
//
// The second test fails on the previous implementation of spline_interpolation()
// under any taped AD backend.
//////////////////////////////////////////////////////////////////////////////

#include "gtest/gtest.h"
#include <psopt.h>
#include <cmath>

namespace interp_test {

// A table symmetric about s = 2. The natural cubic spline through symmetric data
// is itself symmetric, so it is stationary at s = 2 and takes the tabulated value
// 0 there: the minimum over the range of the table is known exactly, without
// having to trust the spline routine to locate it.
static MatrixXd sym_table_x()
{
    MatrixXd X(1,5);
    X << 0.0, 1.0, 2.0, 3.0, 4.0;
    return X;
}

static MatrixXd sym_table_y()
{
    MatrixXd Y(1,5);
    Y << 4.0, 1.0, 0.0, 1.0, 4.0;
    return Y;
}

// A natural cubic spline written out here rather than called from the library, so
// that the routine under test is checked against something independent of it. The
// interior second derivatives solve the usual tridiagonal system with M(0) and
// M(n-1) set to zero; outside the table the end cubics are extended, which is the
// convention the library's binary search had.
static double natural_cubic_spline(const MatrixXd& X, const MatrixXd& Y, double x)
{
    const int n = (int) X.size();
    Eigen::VectorXd M = Eigen::VectorXd::Zero(n);

    if (n > 2) {
        MatrixXd A = MatrixXd::Zero(n-2,n-2);
        Eigen::VectorXd r = Eigen::VectorXd::Zero(n-2);
        for (int i = 1; i <= n-2; i++) {
            double hm = X(i)-X(i-1), hp = X(i+1)-X(i);
            if (i-2 >= 0)   A(i-1,i-2) = hm;
            A(i-1,i-1) = 2.0*(hm+hp);
            if (i <= n-3)   A(i-1,i)   = hp;
            r(i-1) = 6.0*((Y(i+1)-Y(i))/hp - (Y(i)-Y(i-1))/hm);
        }
        M.segment(1,n-2) = A.fullPivLu().solve(r);
    }

    int i = 0;
    while (i < n-2 && x >= X(i+1)) i++;   // the piece the library selects

    double h = X(i+1)-X(i);
    double a = (X(i+1)-x)/h, b = (x-X(i))/h;
    return a*Y(i) + b*Y(i+1) + (a*a*a-a)*h*h/6.0*M(i) + (b*b*b-b)*h*h/6.0*M(i+1);
}

} // namespace interp_test

// --------------------------------------------------------------------------
// Values.
// --------------------------------------------------------------------------

// The adouble overload of spline_interpolation() must agree with a natural cubic
// spline written independently of the library. The comparison covers the whole
// table and runs past both ends, so it also pins the extrapolation convention:
// below the first abscissa the first cubic is used, above the last abscissa the
// last one.
TEST(Interpolation, SplineAgreesWithAnIndependentNaturalSpline)
{
    using namespace interp_test;

    MatrixXd Xd = sym_table_x();
    MatrixXd Yd = sym_table_y();
    const int n = 5;

    const int m = 201;
    MatrixXd Xq = linspace(-1.0, 5.0, m);

    for (int i = 0; i < m; i++) {
        adouble x = Xq(i), y;
        spline_interpolation(&y, x, Xd, Yd, n);
        EXPECT_NEAR(y.value(), natural_cubic_spline(Xd, Yd, Xq(i)), 1.0e-12)
            << "at x = " << Xq(i);
    }

    // The spline interpolates, and is symmetric about the centre of a symmetric
    // table. Both are properties the frozen-interval version lost.
    for (int i = 0; i < n; i++) {
        adouble x = Xd(i), y;
        spline_interpolation(&y, x, Xd, Yd, n);
        EXPECT_NEAR(y.value(), Yd(i), 1.0e-12);
    }
    for (double d = 0.1; d < 2.0; d += 0.1) {
        adouble xl = 2.0-d, xr = 2.0+d, yl, yr;
        spline_interpolation(&yl, xl, Xd, Yd, n);
        spline_interpolation(&yr, xr, Xd, Yd, n);
        EXPECT_NEAR(yl.value(), yr.value(), 1.0e-12);
    }
}

// Linear interpolation reproduces an affine function exactly, everywhere,
// including outside the table where it extrapolates with the end segments.
TEST(Interpolation, LinearReproducesAnAffineFunction)
{
    MatrixXd X(1,4), Y(1,4);
    X << 0.0, 1.0, 2.5, 4.0;
    Y = 3.0*X + MatrixXd::Constant(1,4,-1.0);

    for (double xv = -2.0; xv <= 6.0; xv += 0.25) {
        adouble x = xv, y;
        linear_interpolation(&y, x, X, Y, 4);
        EXPECT_NEAR(y.value(), 3.0*xv-1.0, 1.0e-12) << "at x = " << xv;
    }
}

// Zero order hold: piecewise constant, taking the value at the left end of the
// interval containing x.
TEST(Interpolation, ZohSelectsTheLeftHandValue)
{
    MatrixXd X(1,4), Y(1,4);
    X << 0.0, 1.0, 2.0, 3.0;
    Y << 7.0, -2.0, 5.0, 11.0;

    const double xq[6]   = {0.5, 1.0, 1.5, 2.0, 2.5, 2.9};
    const double yref[6] = {7.0, -2.0, -2.0, 5.0, 5.0, 5.0};

    for (int i = 0; i < 6; i++) {
        adouble x = xq[i], y;
        zoh_interpolation(&y, x, X, Y, 4);
        EXPECT_NEAR(y.value(), yref[i], 1.0e-12) << "at x = " << xq[i];
    }
}

// Bilinear interpolation reproduces a bilinear function exactly. This also covers
// the selection of the second index, which was previously decided by comparing x,
// not y, against the last entry of Y.
TEST(Interpolation, BilinearReproducesABilinearFunction)
{
    MatrixXd X(1,4), Y(1,3), Z(4,3);
    X << 0.0, 1.0, 2.0, 3.0;
    Y << 0.0, 1.5, 3.0;

    for (int i = 0; i < 4; i++)
        for (int j = 0; j < 3; j++)
            Z(i,j) = 1.0 + 2.0*X(i) + 3.0*Y(j) + 4.0*X(i)*Y(j);

    for (double xv = 0.0; xv <= 3.0; xv += 0.3) {
        for (double yv = 0.0; yv <= 3.0; yv += 0.3) {
            adouble x = xv, y = yv, z;
            bilinear_interpolation(&z, x, y, X, Y, Z);
            EXPECT_NEAR(z.value(), 1.0 + 2.0*xv + 3.0*yv + 4.0*xv*yv, 1.0e-10)
                << "at (" << xv << "," << yv << ")";
        }
    }
}

// --------------------------------------------------------------------------
// Taping.
// --------------------------------------------------------------------------

namespace interp_test {

// x' = u on [0,1] with x(0) = 0.5, and a terminal cost read off the spline through
// the symmetric table. The optimizer must drive x(1) to the minimizer of the
// spline, which is s = 2 exactly. The initial guess holds x at 0.5, so the tape is
// recorded in the interval [0,1] while the answer lies two intervals away: an
// implementation that fixes the interval when the tape is recorded minimizes the
// first cubic extrapolated over the whole range instead, and lands on a bound.
//
// The quadratic control penalty is there only to make the control unique; its
// weight is small enough that it moves x(1) by far less than the tolerance below.

adouble endpoint_cost(adouble* i0, adouble* xf, adouble* p, adouble& t0, adouble& tf,
                      adouble* xad, int iphase, Workspace* w)
{
    MatrixXd X = sym_table_x();
    MatrixXd Y = sym_table_y();
    adouble c;
    spline_interpolation(&c, xf[0], X, Y, 5);
    return c;
}

adouble integrand_cost(adouble* s, adouble* u, adouble* p, adouble& t,
                       adouble* xad, int iphase, Workspace* w)
{ return 1.0e-6*u[0]*u[0]; }

void dae(adouble* d, adouble* path, adouble* states, adouble* controls, adouble* p,
         adouble& time, adouble* xad, int iphase, Workspace* w)
{ d[0] = controls[0]; }

void events(adouble* e, adouble* i0, adouble* xf, adouble* p, adouble& t0,
            adouble& tf, adouble* xad, int iphase, Workspace* w)
{ e[0] = i0[0]; }

void linkages(adouble* linkages, adouble* xad, Workspace* w) { }

} // namespace interp_test

TEST(Interpolation, SplineSelectionSurvivesTaping)
{
    using namespace interp_test;

    Alg algorithm; Sol solution; Prob problem;
    const int N = 20;

    problem.name        = "Spline interval selection under AD";
    problem.outfilename = "test_interpolation.txt";
    problem.nphases     = 1;
    problem.nlinkages   = 0;
    psopt_level1_setup(problem);

    problem.phases(1).nstates   = 1;
    problem.phases(1).ncontrols = 1;
    problem.phases(1).nevents   = 1;
    problem.phases(1).npath     = 0;
    problem.phases(1).nodes     << N;
    psopt_level2_setup(problem, algorithm);

    problem.phases(1).bounds.lower.states   << 0.0;    // the range of the table
    problem.phases(1).bounds.upper.states   << 4.0;
    problem.phases(1).bounds.lower.controls << -20.0;
    problem.phases(1).bounds.upper.controls <<  20.0;
    problem.phases(1).bounds.lower.events   << 0.5;
    problem.phases(1).bounds.upper.events   << 0.5;
    problem.phases(1).bounds.lower.StartTime = 0.0;
    problem.phases(1).bounds.upper.StartTime = 0.0;
    problem.phases(1).bounds.lower.EndTime   = 1.0;
    problem.phases(1).bounds.upper.EndTime   = 1.0;

    problem.integrand_cost = &integrand_cost;
    problem.endpoint_cost  = &endpoint_cost;
    problem.dae            = &dae;
    problem.events         = &events;
    problem.linkages       = &linkages;

    problem.phases(1).guess.states   = 0.5*ones(1,N);
    problem.phases(1).guess.controls = zeros(1,N);
    problem.phases(1).guess.time     = linspace(0.0,1.0,N);

    algorithm.nlp_method         = "IPOPT";
    algorithm.scaling            = "automatic";
    algorithm.derivatives        = "automatic";
    algorithm.nlp_iter_max       = 1000;
    algorithm.nlp_tolerance      = 1.e-8;
    algorithm.collocation_method = "trapezoidal";
    algorithm.print_level        = 0;

    psopt(solution, problem, algorithm);

    ASSERT_EQ(solution.error_flag, 0);

    MatrixXd x = solution.get_states_in_phase(1);
    double xf  = x(0, x.cols()-1);

    EXPECT_NEAR(xf, 2.0, 1.0e-3);
    EXPECT_NEAR(solution.cost, 0.0, 1.0e-4);

    // Not at a bound: this is what the frozen-interval implementation produced.
    EXPECT_GT(xf, 0.1);
    EXPECT_LT(xf, 3.9);
}
