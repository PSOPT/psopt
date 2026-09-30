// src/robust.cxx
//
// Robust optimal control by scenario augmentation, v1. See include/robust.h for
// what this is and for the division of labour between the driver and the user.
//
// Copyright (c) Victor M. Becerra, 2026. Part of the PSOPT library (LGPL).

#include "psopt.h"
#include "robust.h"

#include <algorithm>
#include <cstdio>
#include <cstring>
#include <limits>
#include <vector>

//////////////////////////////////////////////////////////////////////////
///////////////////  Constructors  ///////////////////////////////////////
//////////////////////////////////////////////////////////////////////////

RobustUncertainty robust_gaussian(const RowVectorXd& mean, const MatrixXd& cov,
                                  double truncate)
{
    RobustUncertainty U;
    U.kind     = RobustUncertainty::GAUSSIAN;
    U.mean     = mean;
    U.cov      = cov;
    U.truncate = truncate;
    return U;
}

RobustUncertainty robust_uniform(const RowVectorXd& lo, const RowVectorXd& hi)
{
    RobustUncertainty U;
    U.kind = RobustUncertainty::UNIFORM;
    U.lo   = lo;
    U.hi   = hi;
    return U;
}

RobustUncertainty robust_explicit(const std::vector<RowVectorXd>& points)
{
    RobustUncertainty U;
    U.kind   = RobustUncertainty::EXPLICIT;
    U.points = points;
    U.weights.resize(1, (int) points.size());
    for (size_t k = 0; k < points.size(); ++k)
        U.weights(0, (int) k) = points.empty() ? 0.0 : 1.0/(double) points.size();
    return U;
}

int robust_dimension(const RobustUncertainty& U)
{
    switch (U.kind) {
        case RobustUncertainty::GAUSSIAN: return (int) U.mean.size();
        case RobustUncertainty::UNIFORM:  return (int) U.lo.size();
        default: return U.points.empty() ? 0 : (int) U.points[0].size();
    }
}

//////////////////////////////////////////////////////////////////////////
///////////////////  Special functions  //////////////////////////////////
//////////////////////////////////////////////////////////////////////////

// Acklam's rational approximation to the inverse standard normal CDF, followed by
// one Halley refinement against erfc. The approximation alone is good to about
// 1.15e-9 relative; the refinement takes it to machine precision, which costs one
// erfc and is worth having because the same function is exposed to users.
double robust_norm_ppf(double p)
{
    static const double a[6] = { -3.969683028665376e+01,  2.209460984245205e+02,
                                 -2.759285104469687e+02,  1.383577518672690e+02,
                                 -3.066479806614716e+01,  2.506628277459239e+00 };
    static const double b[5] = { -5.447609879822406e+01,  1.615858368580409e+02,
                                 -1.556989798598866e+02,  6.680131188771972e+01,
                                 -1.328068155288572e+01 };
    static const double c[6] = { -7.784894002430293e-03, -3.223964580411365e-01,
                                 -2.400758277161838e+00, -2.549732539343734e+00,
                                  4.374664141464968e+00,  2.938163982698783e+00 };
    static const double d[4] = {  7.784695709041462e-03,  3.224671290700398e-01,
                                  2.445134137142996e+00,  3.754408661907416e+00 };
    const double plow = 0.02425, phigh = 1.0 - plow;

    if (p <= 0.0) return -std::numeric_limits<double>::infinity();
    if (p >= 1.0) return  std::numeric_limits<double>::infinity();

    double x;
    if (p < plow) {
        double q = std::sqrt(-2.0*std::log(p));
        x = (((((c[0]*q + c[1])*q + c[2])*q + c[3])*q + c[4])*q + c[5]) /
            ((((d[0]*q + d[1])*q + d[2])*q + d[3])*q + 1.0);
    } else if (p > phigh) {
        double q = std::sqrt(-2.0*std::log(1.0 - p));
        x = -(((((c[0]*q + c[1])*q + c[2])*q + c[3])*q + c[4])*q + c[5]) /
             ((((d[0]*q + d[1])*q + d[2])*q + d[3])*q + 1.0);
    } else {
        double q = p - 0.5, r = q*q;
        x = (((((a[0]*r + a[1])*r + a[2])*r + a[3])*r + a[4])*r + a[5])*q /
            (((((b[0]*r + b[1])*r + b[2])*r + b[3])*r + b[4])*r + 1.0);
    }

    // One Halley step on Phi(x) - p = 0.
    const double e = 0.5*std::erfc(-x/std::sqrt(2.0)) - p;
    const double u = e*std::sqrt(2.0*M_PI)*std::exp(x*x/2.0);
    return x - u/(1.0 + x*u/2.0);
}

// The regularized lower incomplete gamma P(a, x): series below x = a+1, continued
// fraction above, which is where each converges quickly.
double robust_lower_gamma(double a, double x)
{
    if (x <= 0.0 || a <= 0.0) return 0.0;
    const double lg = std::lgamma(a);

    if (x < a + 1.0) {
        double ap = a, sum = 1.0/a, del = sum;
        for (int n = 0; n < 500; ++n) {
            ap  += 1.0;
            del *= x/ap;
            sum += del;
            if (std::fabs(del) < std::fabs(sum)*1.0e-16) break;
        }
        return sum*std::exp(-x + a*std::log(x) - lg);
    }

    // Lentz's method on the continued fraction for Q(a, x).
    const double tiny = 1.0e-300;
    double b = x + 1.0 - a, c = 1.0/tiny, d = 1.0/b, h = d;
    for (int i = 1; i < 500; ++i) {
        const double an = -i*(i - a);
        b += 2.0;
        d  = an*d + b; if (std::fabs(d) < tiny) d = tiny;
        c  = b + an/c; if (std::fabs(c) < tiny) c = tiny;
        d  = 1.0/d;
        const double del = d*c;
        h *= del;
        if (std::fabs(del - 1.0) < 1.0e-16) break;
    }
    return 1.0 - std::exp(-x + a*std::log(x) - lg)*h;
}

double robust_set_coverage(double k, int n)
{
    if (n <= 0 || k <= 0.0) return 0.0;
    return robust_lower_gamma(0.5*n, 0.5*k*k);
}

//////////////////////////////////////////////////////////////////////////
///////////////////  Set geometry  ///////////////////////////////////////
//////////////////////////////////////////////////////////////////////////

// Lower Cholesky factor of the covariance. The map z -> mean + L z carries the
// unit ball to the set and the standard normal to the distribution, which is what
// every Gaussian operation below is written in terms of.
static MatrixXd cholesky_of(const MatrixXd& cov)
{
    return Eigen::LLT<MatrixXd>(cov).matrixL();
}

static RowVectorXd from_unit(const RobustUncertainty& U, const RowVectorXd& z)
{
    const MatrixXd L = cholesky_of(U.cov);
    return U.mean + (L*z.transpose()).transpose();
}

double robust_mahalanobis(const RobustUncertainty& U, const RowVectorXd& theta)
{
    if (U.kind == RobustUncertainty::GAUSSIAN) {
        const Eigen::VectorXd dx = (theta - U.mean).transpose();
        return std::sqrt(std::max(0.0, dx.dot(U.cov.llt().solve(dx))));
    }
    if (U.kind == RobustUncertainty::UNIFORM) {
        double worst = 0.0;
        for (int i = 0; i < U.lo.size(); ++i) {
            const double c = 0.5*(U.lo(i) + U.hi(i));
            const double h = 0.5*(U.hi(i) - U.lo(i));
            if (h > 0.0) worst = std::max(worst, std::fabs(theta(i) - c)/h);
        }
        return worst;
    }
    double best = std::numeric_limits<double>::infinity();
    for (size_t k = 0; k < U.points.size(); ++k)
        best = std::min(best, (double) (U.points[k] - theta).norm());
    return best;
}

bool robust_contains(const RobustUncertainty& U, const RowVectorXd& theta)
{
    switch (U.kind) {
        case RobustUncertainty::GAUSSIAN:
            return robust_mahalanobis(U, theta) <= U.truncate + 1.0e-12;
        case RobustUncertainty::UNIFORM:
            for (int i = 0; i < U.lo.size(); ++i)
                if (theta(i) < U.lo(i) - 1.0e-12 || theta(i) > U.hi(i) + 1.0e-12)
                    return false;
            return true;
        default:
            for (size_t k = 0; k < U.points.size(); ++k)
                if ((U.points[k] - theta).norm() <= 1.0e-12) return true;
            return false;
    }
}

RowVectorXd robust_clip(const RobustUncertainty& U, const RowVectorXd& theta)
{
    if (U.kind == RobustUncertainty::GAUSSIAN) {
        const double r = robust_mahalanobis(U, theta);
        if (r <= U.truncate || r == 0.0) return theta;
        return U.mean + (theta - U.mean)*(U.truncate/r);
    }
    if (U.kind == RobustUncertainty::UNIFORM) {
        RowVectorXd out = theta;
        for (int i = 0; i < U.lo.size(); ++i)
            out(i) = std::min(std::max(theta(i), U.lo(i)), U.hi(i));
        return out;
    }
    double      best = std::numeric_limits<double>::infinity();
    RowVectorXd out  = theta;
    for (size_t k = 0; k < U.points.size(); ++k) {
        const double d = (U.points[k] - theta).norm();
        if (d < best) { best = d; out = U.points[k]; }
    }
    return out;
}

bool robust_sigma_points(const RobustUncertainty& U,
                         std::vector<RowVectorXd>& points, RowVectorXd& weights)
{
    if (U.kind == RobustUncertainty::EXPLICIT) {
        points  = U.points;
        weights = U.weights;
        return true;
    }

    const int    n     = robust_dimension(U);
    const double kappa = 3.0 - (double) n;
    if (n + kappa <= 0.0 || kappa < 0.0) return false;
    const double c = std::sqrt((double) n + kappa);

    RowVectorXd centre(n);
    MatrixXd    spread(n, n);
    if (U.kind == RobustUncertainty::GAUSSIAN) {
        centre = U.mean;
        spread = cholesky_of(U.cov);
    } else {
        // For a uniform distribution each component's standard deviation is
        // (hi-lo)/sqrt(12), and spreading by sqrt(n+kappa) times THAT is what makes
        // the set match the second moment rather than merely span the box.
        spread.setZero();
        for (int i = 0; i < n; ++i) {
            centre(i)    = 0.5*(U.lo(i) + U.hi(i));
            spread(i, i) = (U.hi(i) - U.lo(i))/std::sqrt(12.0);
        }
    }

    points.clear();
    points.push_back(centre);
    std::vector<double> w;
    w.push_back(kappa/((double) n + kappa));
    for (int k = 0; k < n; ++k) {
        const RowVectorXd col = c*spread.col(k).transpose();
        RowVectorXd plus  = centre + col;
        RowVectorXd minus = centre - col;
        if (U.kind == RobustUncertainty::UNIFORM) {
            plus  = robust_clip(U, plus);
            minus = robust_clip(U, minus);
        }
        points.push_back(plus);
        points.push_back(minus);
        w.push_back(0.5/((double) n + kappa));
        w.push_back(0.5/((double) n + kappa));
    }
    weights.resize(1, (int) w.size());
    for (size_t k = 0; k < w.size(); ++k) weights(0, (int) k) = w[k];
    return true;
}

// The i-th element of the van der Corput sequence in the given base.
static double van_der_corput(long i, int base)
{
    double f = 1.0/(double) base, r = 0.0;
    while (i > 0) {
        r += f*(double) (i % base);
        i /= base;
        f /= (double) base;
    }
    return r;
}

void robust_low_discrepancy(const RobustUncertainty& U, int m, int skip,
                            std::vector<RowVectorXd>& points)
{
    static const int primes[8] = { 2, 3, 5, 7, 11, 13, 17, 19 };
    const int n = robust_dimension(U);
    points.clear();
    if (m <= 0 || n <= 0) return;

    if (U.kind == RobustUncertainty::EXPLICIT) {
        points = U.points;                  // the set IS the list; nothing to fill
        return;
    }

    for (int j = 0; j < m; ++j) {
        RowVectorXd h(n);
        for (int i = 0; i < n; ++i)
            h(i) = van_der_corput((long) (j + skip + 1), primes[i % 8]);

        if (U.kind == RobustUncertainty::UNIFORM) {
            RowVectorXd p(n);
            for (int i = 0; i < n; ++i) p(i) = U.lo(i) + h(i)*(U.hi(i) - U.lo(i));
            points.push_back(p);
        } else {
            // Uniform by volume in the unit ball: a direction from a normal
            // deviate, and a radius stratified over (0, 1] with 1 attained.
            RowVectorXd g(n);
            double      nrm = 0.0;
            for (int i = 0; i < n; ++i) {
                g(i) = robust_norm_ppf(std::min(std::max(h(i), 1.0e-12),
                                                1.0 - 1.0e-12));
                nrm += g(i)*g(i);
            }
            nrm = std::sqrt(nrm);
            if (nrm == 0.0) nrm = 1.0;
            const double r = std::pow((double) (j + 1)/(double) m, 1.0/(double) n);
            points.push_back(from_unit(U, g*(U.truncate*r/nrm)));
        }
    }
}

void robust_boundary_points(const RobustUncertainty& U,
                            std::vector<RowVectorXd>& points)
{
    const int n = robust_dimension(U);
    points.clear();
    if (n <= 0) return;

    if (U.kind == RobustUncertainty::EXPLICIT) { points = U.points; return; }

    if (U.kind == RobustUncertainty::GAUSSIAN) {
        const MatrixXd L = cholesky_of(U.cov);
        for (int k = 0; k < n; ++k) {
            const RowVectorXd col = U.truncate*L.col(k).transpose();
            points.push_back(U.mean + col);
            points.push_back(U.mean - col);
        }
        return;
    }

    if (n <= 8) {                                   // every corner of the box
        const long N = 1L << n;
        for (long mask = 0; mask < N; ++mask) {
            RowVectorXd p(n);
            for (int i = 0; i < n; ++i)
                p(i) = (mask >> i) & 1L ? U.hi(i) : U.lo(i);
            points.push_back(p);
        }
    } else {                                        // face centres beyond
        RowVectorXd mid(n);
        for (int i = 0; i < n; ++i) mid(i) = 0.5*(U.lo(i) + U.hi(i));
        for (int i = 0; i < n; ++i) {
            RowVectorXd a = mid, b = mid;
            a(i) = U.lo(i);
            b(i) = U.hi(i);
            points.push_back(a);
            points.push_back(b);
        }
    }
}

//////////////////////////////////////////////////////////////////////////
///////////////////  The worst-case oracle  //////////////////////////////
//////////////////////////////////////////////////////////////////////////

static const int ROBUST_CLOUD_SWEEPS = 6;

long robust_worst_case_evaluations(const RobustUncertainty& U, int n_seed,
                                   int n_refine)
{
    std::vector<RowVectorXd> b;
    robust_boundary_points(U, b);
    if (U.kind == RobustUncertainty::EXPLICIT)
        return (long) U.points.size();              // exhaustive, once
    const int keep = std::max(1, n_refine);
    const int per  = std::max(8, 64/keep);
    return (long) b.size() + (long) std::max(0, n_seed)
           + (long) ROBUST_CLOUD_SWEEPS*keep*per;
}

namespace {
struct Candidate {
    double      value;
    RowVectorXd at;
};
bool worse_first(const Candidate& a, const Candidate& b) { return a.value > b.value; }
}

double robust_worst_case(const RobustUncertainty& U,
                         RobustViolationFn violation, void* user_data,
                         int n_seed, int n_refine, unsigned seed,
                         RowVectorXd& worst_at)
{
    std::vector<RowVectorXd> seeds, more;
    robust_boundary_points(U, seeds);
    if (U.kind != RobustUncertainty::EXPLICIT) {
        robust_low_discrepancy(U, std::max(0, n_seed), (int) (seed % 1000u), more);
        seeds.insert(seeds.end(), more.begin(), more.end());
    }

    std::vector<Candidate> best;
    for (size_t k = 0; k < seeds.size(); ++k) {
        Candidate c;
        c.at    = seeds[k];
        c.value = violation(seeds[k], user_data);
        best.push_back(c);
    }
    if (best.empty()) { worst_at = RowVectorXd(); return 0.0; }
    std::sort(best.begin(), best.end(), worse_first);

    if (U.kind == RobustUncertainty::EXPLICIT) {     // the set is the list
        worst_at = best[0].at;
        return best[0].value;
    }

    const int keep = std::max(1, n_refine);
    if ((int) best.size() > keep) best.resize(keep);

    // A shrinking cloud around the best few. Deterministic given `seed`: the
    // displacements come from the same van der Corput machinery as the seeding, so
    // a certificate is reproducible, which a rand() would not give.
    const int n   = robust_dimension(U);
    const int per = std::max(8, 64/keep);
    double    rad = 0.25;
    long      tick = (long) (seed % 7919u) + 1;
    for (int sweep = 0; sweep < ROBUST_CLOUD_SWEEPS; ++sweep) {
        std::vector<Candidate> pool = best;
        for (size_t b = 0; b < best.size(); ++b) {
            for (int j = 0; j < per; ++j) {
                RowVectorXd z(n);
                for (int i = 0; i < n; ++i) {
                    static const int bases[4] = { 2, 3, 5, 7 };
                    const double h = van_der_corput(tick, bases[i % 4]);
                    z(i) = robust_norm_ppf(std::min(std::max(h, 1.0e-12),
                                                    1.0 - 1.0e-12));
                }
                ++tick;
                Candidate c;
                c.at = robust_clip(U, best[b].at
                                   + rad*(U.kind == RobustUncertainty::GAUSSIAN
                                          ? (from_unit(U, z) - U.mean)
                                          : RowVectorXd(z.cwiseProduct(
                                                0.5*(U.hi - U.lo)))));
                c.value = violation(c.at, user_data);
                pool.push_back(c);
            }
        }
        std::sort(pool.begin(), pool.end(), worse_first);
        if ((int) pool.size() > keep) pool.resize(keep);
        best = pool;
        rad *= 0.45;
    }

    worst_at = best[0].at;
    return best[0].value;
}

//////////////////////////////////////////////////////////////////////////
///////////////////  The scenario budget  ////////////////////////////////
//////////////////////////////////////////////////////////////////////////

namespace {
struct DofParts { int dof, ncontrols, nodes, free_times, free_params, neq, free_x; };

// Does this transcription impose a defect at every stored node? The global Lobatto
// schemes do and nothing else does, and it is the one question the budget turns on: it
// decides whether a scenario's own initial state costs a degree of freedom.
bool collocates_every_node(Alg& algorithm)
{
    if ( is_multiple_shooting(algorithm) ) return false;
    return ( algorithm.collocation_method == "Legendre"
             || algorithm.collocation_method == "Chebyshev" );
}

DofParts dof_parts(Prob& problem, int iphase, Alg& algorithm)
{
    DofParts d;
    d.ncontrols = problem.phases(iphase).ncontrols;
    d.nodes     = problem.phases(iphase).current_number_of_intervals + 1;

    const Bounds& b = problem.phases(iphase).bounds;
    d.free_times = (b.lower.StartTime != b.upper.StartTime ? 1 : 0)
                 + (b.lower.EndTime   != b.upper.EndTime   ? 1 : 0);

    d.free_params = 0;
    for (int i = 0; i < problem.phases(iphase).nparameters; ++i)
        if (b.upper.parameters(i) > b.lower.parameters(i)) ++d.free_params;

    d.neq = 0;
    for (int i = 0; i < problem.phases(iphase).nevents; ++i)
        if (b.upper.events(i) <= b.lower.events(i)) ++d.neq;

    d.free_x = collocates_every_node(algorithm)
             ? 0 : problem.phases(iphase).nstates;

    d.dof = d.free_x + d.ncontrols*d.nodes + d.free_times + d.free_params - d.neq;
    return d;
}
}

int robust_degrees_of_freedom(Prob& problem, int iphase, Alg& algorithm)
{
    return dof_parts(problem, iphase, algorithm).dof;
}

bool robust_guess_is_finite(Prob& problem, int iphase, const char** offender)
{
    struct Check {
        static bool finite(const MatrixXd& m)
        {
            for (int i = 0; i < m.rows(); ++i)
                for (int j = 0; j < m.cols(); ++j)
                    if (!std::isfinite(m(i, j))) return false;
            return true;
        }
    };
    const Guess& g = problem.phases(iphase).guess;
    if (!Check::finite(g.states))
        { if (offender) *offender = "states";     return false; }
    if (!Check::finite(g.controls))
        { if (offender) *offender = "controls";   return false; }
    if (!Check::finite(g.time))
        { if (offender) *offender = "time";       return false; }
    if (!Check::finite(g.parameters))
        { if (offender) *offender = "parameters"; return false; }
    return true;
}

void robust_dof_message(Prob& problem, int iphase, int nscenarios, Alg& algorithm,
                        char* buffer, size_t buffer_size)
{
    if (!buffer || buffer_size == 0) return;
    const DofParts d   = dof_parts(problem, iphase, algorithm);
    const int      M   = std::max(1, nscenarios);
    const int      per = std::max(1, d.neq/M);

    if ( d.free_x > 0 ) {
        snprintf(buffer, buffer_size,
            "%d scenarios leave the transcription about %d degrees of freedom: "
            "%d free initial state(s), %d control(s) at %d nodes, %d free time "
            "endpoint(s) and %d free static parameter(s), against %d equality "
            "events -- some %d per scenario.\n"
            "  A scenario's own initial condition is FREE: it pins the state values "
            "at t0 that the defect equations leave undetermined, and the two cancel. "
            "What costs a scenario is an equality event beyond that, and a pinned "
            "TERMINAL state is the one to look for -- one open-loop control cannot "
            "steer several different plants to the same point, so relax it to a "
            "tolerance, which costs nothing at all because an inequality event takes "
            "no degree of freedom.",
            nscenarios, d.dof, d.free_x, d.ncontrols, d.nodes, d.free_times,
            d.free_params, d.neq, per);
        return;
    }

    snprintf(buffer, buffer_size,
        "%d scenarios leave the transcription about %d degrees of freedom: "
        "%d control(s) at %d nodes, %d free time endpoint(s) and %d free static "
        "parameter(s), against %d equality events -- some %d per scenario.\n"
        "  This transcription collocates EVERY stored node, so a scenario's initial "
        "condition is not free: the defect block already holds one condition per "
        "stored value and the initial condition is one more. A replicated state "
        "therefore costs %d degrees of freedom per scenario whatever the node count, "
        "and raising the mesh does not help.\n"
        "  Use transcription_method = \"multiple-shooting\", or collocation_method "
        "= \"trapezoidal\", neither of which has this deficit: on the two-link arm "
        "both carry eighty scenarios where this one is refused above twelve.",
        nscenarios, d.dof, d.ncontrols, d.nodes, d.free_times, d.free_params,
        d.neq, per, problem.phases(iphase).nstates/M);
}

//////////////////////////////////////////////////////////////////////////
///////////////////  The driver  /////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////

// A trial solve is accepted only if psopt() reported no PSOPT-level error and the
// NLP itself converged. The same test psopt_solve_integer uses, and for the same
// reason: error_flag == 0 alone is not sufficient, because IPOPT may catch an
// internal exception and return a failure status while it stays 0.
static inline bool robust_trial_succeeded(const Sol& s)
{
    return s.error_flag == 0 && (s.nlp_return_code == 0 || s.nlp_return_code == 1);
}

namespace {
// Binds the user's (theta, solution) violation callback to the oracle's
// (theta, void*) signature for the duration of one search.
struct OracleContext {
    RobustSpec*         spec;
    const RobustDesign* design;
    long                calls;
};

double oracle_bridge(const RowVectorXd& theta, void* user_data)
{
    OracleContext* c = (OracleContext*) user_data;
    ++c->calls;
    return c->spec->violation(theta, *c->design, c->spec->user_data);
}

// Take a copy of what a solve produced. Sol cannot be held; this can.
RobustDesign design_of(Sol& solution, Prob& problem)
{
    RobustDesign d;
    d.time      = solution.get_time_in_phase(1);
    d.controls  = solution.get_controls_in_phase(1);
    d.states    = solution.get_states_in_phase(1);
    if (problem.phases(1).nparameters > 0)
        d.parameters = solution.get_parameters_in_phase(1);
    d.objective = solution.get_cost();
    d.valid     = true;
    return d;
}
}

//////////////////////////////////////////////////////////////////////////
///////////////////  The model interface  ////////////////////////////////
//////////////////////////////////////////////////////////////////////////

void robust_margins(const RowVectorXd& lower, const RowVectorXd& upper,
                    double tighten, RowVectorXd& margin, int& n_one_sided)
{
    const int n = (int) lower.size();
    margin = zeros(1, n);
    n_one_sided = 0;
    if ((int) upper.size() != n) return;
    for (int j = 0; j < n; ++j) {
        const double half = 0.5*(upper(j) - lower(j));
        if (!std::isfinite(half) || fabs(half) > 1.0e29) { ++n_one_sided; continue; }
        if (half <= 0.0) continue;               // pinned: an equality, left alone
        margin(j) = (1.0 - tighten)*half;
    }
}

void robust_dae_value(const RobustModel& model, const double* theta, int ntheta,
                      const double* states, const double* controls,
                      const double* parameters, double time,
                      double* derivatives, double* path)
{
    const int ns = model.nstates, nc = model.ncontrols;
    const int np = model.npath, npar = model.nparameters;

    std::vector<adouble> xa(ns > 0 ? ns : 1), ua(nc > 0 ? nc : 1);
    std::vector<adouble> pa(npar > 0 ? npar : 1), da(ns > 0 ? ns : 1);
    std::vector<adouble> ga(np > 0 ? np : 1);
    adouble ta = time;

    for (int j = 0; j < ns;   ++j) xa[j] = states[j];
    for (int j = 0; j < nc;   ++j) ua[j] = controls[j];
    for (int j = 0; j < npar; ++j) pa[j] = parameters[j];

    model.dae(&da[0], &ga[0], &xa[0], &ua[0], &pa[0], ta, theta, ntheta, 0, 1, 0);

    for (int j = 0; j < ns; ++j) derivatives[j] = da[j].value();
    if (path) for (int j = 0; j < np; ++j) path[j] = ga[j].value();
}

// What the library needs in order to loop the user's nominal equations. Held in
// problem.robust_data, which exists so that user_data stays the user's.
struct RobustAugmentation {
    const RobustModel*       model;
    std::vector<RowVectorXd> scenarios;
    RobustAugmentation() : model(0) {}
};

static void robust_augmented_dae(adouble* derivatives, adouble* path, adouble* states,
                                 adouble* controls, adouble* parameters, adouble& time,
                                 adouble* xad, int iphase, Workspace* workspace)
{
    const RobustAugmentation& A =
        *((const RobustAugmentation*) workspace->problem->robust_data);
    const int ns = A.model->nstates, np = A.model->npath;
    const int M  = (int) A.scenarios.size();
    for (int i = 0; i < M; ++i)
        A.model->dae(derivatives + ns*i, np > 0 ? path + np*i : path,
                     states + ns*i, controls, parameters, time,
                     A.scenarios[i].data(), (int) A.scenarios[i].size(),
                     xad, iphase, workspace);
}

static void robust_augmented_events(adouble* e, adouble* initial_states,
                                    adouble* final_states, adouble* parameters,
                                    adouble& t0, adouble& tf,
                                    adouble* xad, int iphase, Workspace* workspace)
{
    const RobustAugmentation& A =
        *((const RobustAugmentation*) workspace->problem->robust_data);
    const int ns = A.model->nstates, ne = A.model->nevents;
    const int M  = (int) A.scenarios.size();
    if (ne == 0) return;
    for (int i = 0; i < M; ++i)
        A.model->events(e + ne*i, initial_states + ns*i, final_states + ns*i,
                        parameters, t0, tf,
                        A.scenarios[i].data(), (int) A.scenarios[i].size(),
                        xad, iphase, workspace);
}

static void robust_no_linkages(adouble*, adouble*, Workspace*) {}

// The warm start: each scenario's own plant integrated through the previous control.
//
// A guess that is merely the right shape leaves the augmented problem with M copies of an
// infeasible arc and the solver free to wander to a distant local minimum. The control is
// read as held under a constant parameterisation and as a straight line otherwise, which
// is approximate for the two parameterisations that carry a midpoint control; that is
// admissible HERE, because this is a guess, and is not admissible in a verification,
// where it would measure a different controller.
//
// Returns false when the integration did not stay finite even at the largest substep
// count, which is the caller's signal to fall back to the nominal guess.
static bool robust_warm_states(const RobustModel& model,
                               const std::vector<RowVectorXd>& scenarios,
                               const RobustDesign& previous, Alg& algorithm,
                               MatrixXd& x_guess)
{
    const int ns = model.nstates, nc = model.ncontrols, M = (int) scenarios.size();
    const int N  = (int) previous.time.cols();
    if (N < 2 || (int) previous.controls.cols() != N) return false;
    if ((int) model.initial_state.size() != ns)       return false;

    const bool held = is_multiple_shooting(algorithm)
                      && algorithm.ms_control_parameterisation == "constant";
    std::vector<double> par(model.nparameters > 0 ? model.nparameters : 1, 0.0);
    for (int j = 0; j < model.nparameters && j < (int) previous.parameters.size(); ++j)
        par[j] = previous.parameters(j);

    x_guess = zeros(ns*M, N);

    for (int sub = (model.warm_substeps > 0 ? model.warm_substeps : 8);
         sub <= robust_warm_substeps_max; sub *= 2) {
        bool ok = true;
        for (int i = 0; i < M && ok; ++i) {
            std::vector<double> x(ns), y(ns), k1(ns), k2(ns), k3(ns), k4(ns);
            std::vector<double> uA(nc > 0 ? nc : 1), uH(nc > 0 ? nc : 1),
                                uB(nc > 0 ? nc : 1);
            for (int j = 0; j < ns; ++j) {
                x[j] = model.initial_state(j);
                x_guess(ns*i + j, 0) = x[j];
            }
            for (int c = 0; c < N - 1 && ok; ++c) {
                const double ta = previous.time(0, c);
                const double h  = (previous.time(0, c+1) - ta)/sub;
                for (int s = 0; s < sub; ++s) {
                    const double w0 =  s        /(double) sub;
                    const double wh = (s + 0.5) /(double) sub;
                    const double w1 = (s + 1.0) /(double) sub;
                    for (int j = 0; j < nc; ++j) {
                        const double a = previous.controls(j, c);
                        const double b = previous.controls(j, c+1);
                        uA[j] = held ? a : a + w0*(b - a);
                        uH[j] = held ? a : a + wh*(b - a);
                        uB[j] = held ? a : a + w1*(b - a);
                    }
                    const double* th = scenarios[i].data();
                    const int     nt = (int) scenarios[i].size();
                    robust_dae_value(model, th, nt, &x[0], &uA[0], &par[0],
                                     ta + s*h, &k1[0], 0);
                    for (int j = 0; j < ns; ++j) y[j] = x[j] + 0.5*h*k1[j];
                    robust_dae_value(model, th, nt, &y[0], &uH[0], &par[0],
                                     ta + (s + 0.5)*h, &k2[0], 0);
                    for (int j = 0; j < ns; ++j) y[j] = x[j] + 0.5*h*k2[j];
                    robust_dae_value(model, th, nt, &y[0], &uH[0], &par[0],
                                     ta + (s + 0.5)*h, &k3[0], 0);
                    for (int j = 0; j < ns; ++j) y[j] = x[j] + h*k3[j];
                    robust_dae_value(model, th, nt, &y[0], &uB[0], &par[0],
                                     ta + (s + 1.0)*h, &k4[0], 0);
                    for (int j = 0; j < ns; ++j)
                        x[j] += (h/6.0)*(k1[j] + 2.0*k2[j] + 2.0*k3[j] + k4[j]);
                }
                for (int j = 0; j < ns; ++j) {
                    if (!std::isfinite(x[j])) { ok = false; break; }
                    x_guess(ns*i + j, c+1) = x[j];
                }
            }
        }
        if (ok) return true;
    }
    return false;
}

// Build the augmented problem for one scenario list. This is what spec.setup would have
// been, and it is generated instead.
static void robust_model_setup(Prob& problem, Alg& algorithm, const RobustModel& model,
                               const std::vector<RowVectorXd>& scenarios,
                               const RobustDesign& previous,
                               RobustAugmentation& aug, bool verbose)
{
    const int M   = (int) scenarios.size();
    const int ns  = model.nstates, nc = model.ncontrols;
    const int ne  = model.nevents, np = model.npath;

    aug.model     = &model;
    aug.scenarios = scenarios;

    problem.nphases   = 1;
    problem.nlinkages = 0;
    psopt_level1_setup(problem);

    problem.phases(1).nstates     = ns*M;
    problem.phases(1).ncontrols   = nc;
    problem.phases(1).nevents     = ne*M;
    problem.phases(1).npath       = np*M;
    problem.phases(1).nparameters = model.nparameters;
    problem.phases(1).nodes       = model.nodes;
    psopt_level2_setup(problem, algorithm);
    // The one moment at which the caller's options can be set, because the call above has
    // just reset every one of them to its default. See RobustModel::configure.
    if (model.configure) model.configure(algorithm, model.configure_data);

    if (problem.name.empty())        problem.name        = "robust design";
    if (problem.outfilename.empty()) problem.outfilename = "robust_design.txt";

    problem.robust_data = (void*) &aug;

    for (int i = 0; i < M; ++i)
        for (int j = 0; j < ns; ++j) {
            problem.phases(1).bounds.lower.states(ns*i + j) = model.states_lower(j);
            problem.phases(1).bounds.upper.states(ns*i + j) = model.states_upper(j);
        }
    for (int j = 0; j < nc; ++j) {
        problem.phases(1).bounds.lower.controls(j) = model.controls_lower(j);
        problem.phases(1).bounds.upper.controls(j) = model.controls_upper(j);
    }
    for (int j = 0; j < model.nparameters; ++j) {
        problem.phases(1).bounds.lower.parameters(j) = model.parameters_lower(j);
        problem.phases(1).bounds.upper.parameters(j) = model.parameters_upper(j);
    }

    // The tightening. Applied to the events and the path constraints, per scenario, and
    // reported once rather than per scenario.
    RowVectorXd emargin, pmargin;
    int e_one_sided = 0, p_one_sided = 0;
    (void) e_one_sided; (void) p_one_sided;   // reported once, by the driver
    if (ne > 0)
        robust_margins(model.events_lower, model.events_upper, model.tighten,
                       emargin, e_one_sided);
    if (np > 0)
        robust_margins(model.path_lower, model.path_upper, model.tighten,
                       pmargin, p_one_sided);

    for (int i = 0; i < M; ++i) {
        for (int j = 0; j < ne; ++j) {
            const double lo = model.events_lower(j), up = model.events_upper(j);
            const double mg = (up > lo) ? emargin(j) : 0.0;
            problem.phases(1).bounds.lower.events(ne*i + j) = lo + mg;
            problem.phases(1).bounds.upper.events(ne*i + j) = up - mg;
        }
        for (int j = 0; j < np; ++j) {
            const double lo = model.path_lower(j), up = model.path_upper(j);
            const double mg = (up > lo) ? pmargin(j) : 0.0;
            problem.phases(1).bounds.lower.path(np*i + j) = lo + mg;
            problem.phases(1).bounds.upper.path(np*i + j) = up - mg;
        }
    }

    problem.phases(1).bounds.lower.StartTime = model.t0_lower;
    problem.phases(1).bounds.upper.StartTime = model.t0_upper;
    problem.phases(1).bounds.lower.EndTime   = model.tf_lower;
    problem.phases(1).bounds.upper.EndTime   = model.tf_upper;

    problem.endpoint_cost  = model.endpoint_cost;
    problem.integrand_cost = model.integrand_cost;
    problem.dae            = &robust_augmented_dae;
    problem.events         = (ne > 0) ? &robust_augmented_events : 0;
    problem.linkages       = &robust_no_linkages;

    // The guess. The warm start when there is a previous design and it integrated
    // cleanly; the user's own nominal guess, replicated, otherwise.
    MatrixXd x_guess;
    bool warm = previous.valid
                && robust_warm_states(model, scenarios, previous, algorithm, x_guess);
    if (warm) {
        problem.phases(1).guess.controls = previous.controls;
        problem.phases(1).guess.time     = previous.time;
    } else {
        const int N = (int) model.guess_time.cols();
        x_guess = zeros(ns*M, N);
        for (int i = 0; i < M; ++i)
            for (int j = 0; j < ns; ++j)
                x_guess.row(ns*i + j) = model.guess_states.row(j);
        problem.phases(1).guess.controls = model.guess_controls;
        problem.phases(1).guess.time     = model.guess_time;
        if (verbose && previous.valid)
            printf("  the warm start did not stay finite at %d substeps or below; "
                   "using the nominal guess\n", robust_warm_substeps_max);
    }
    problem.phases(1).guess.states = x_guess;
    if (model.nparameters > 0 && (int) model.guess_parameters.size() == model.nparameters)
        problem.phases(1).guess.parameters = model.guess_parameters;
}

// What the caller is owed before the loop starts: whether the model is complete, and
// whether the tightening will actually reach every constraint it is meant to.
static bool robust_model_is_usable(const RobustModel& model, bool verbose)
{
    if (model.nstates <= 0 || !model.dae) {
        error_message("psopt_solve_robust: spec.model needs at least nstates and dae");
        return false;
    }
    if ((int) model.initial_state.size() != model.nstates) {
        error_message("psopt_solve_robust: spec.model.initial_state must have nstates "
                      "entries. The warm start integrates every scenario out of it, and "
                      "it is not inferred from the event bounds");
        return false;
    }
    if ((int) model.states_lower.size() != model.nstates ||
        (int) model.states_upper.size() != model.nstates) {
        error_message("psopt_solve_robust: spec.model needs nominal state bounds");
        return false;
    }
    if (model.nevents > 0 && !model.events) {
        error_message("psopt_solve_robust: spec.model declares nevents but no events "
                      "function");
        return false;
    }
    if (model.nodes.size() < 1 || (int) model.guess_time.cols() < 2) {
        error_message("psopt_solve_robust: spec.model needs nodes and a nominal guess");
        return false;
    }
    if (verbose) {
        RowVectorXd mg;
        int e1 = 0, p1 = 0;
        if (model.nevents > 0)
            robust_margins(model.events_lower, model.events_upper, model.tighten, mg, e1);
        if (model.npath > 0)
            robust_margins(model.path_lower, model.path_upper, model.tighten, mg, p1);
        if (e1 + p1 > 0)
            printf("\npsopt_solve_robust: %d event and %d path bound(s) are one-sided, "
                   "so tighten gives them\n  no margin. A design left sitting on such a "
                   "bound is where the between-scenario\n  overshoot appears; widen the "
                   "bound or state it two-sided if the loop will not converge.\n",
                   e1, p1);
    }
    return true;
}

// Build the problem for one scenario list, by whichever route the caller chose. Every
// call site in the driver goes through here, so the two paths cannot drift apart.
static void robust_do_setup(Prob& problem, Alg& algorithm, RobustSpec& spec,
                            const std::vector<RowVectorXd>& scenarios,
                            const RobustDesign& previous, RobustAugmentation& aug)
{
    if (spec.model)
        robust_model_setup(problem, algorithm, *spec.model, scenarios, previous, aug,
                           spec.verbose);
    else
        spec.setup(problem, algorithm, scenarios, previous, spec.user_data);
}

bool robust_prefer_cold(double cold_objective, double cold_worst,
                        double warm_objective, double warm_worst, double slack)
{
    const bool cert_cold = (cold_worst <= slack);
    const bool cert_warm = (warm_worst <= slack);
    if (cert_cold && cert_warm) return cold_objective < warm_objective;
    if (cert_cold)              return true;
    if (cert_warm)              return false;
    return cold_worst < warm_worst;
}

[[nodiscard]] int psopt_solve_robust(Sol& solution, RobustSpec& spec,
                                     Prob& problem, Alg& algorithm)
{
    if (!spec.violation) {
        error_message("psopt_solve_robust: spec.violation must be supplied. The "
                      "verification has to be independent of the transcription, and only "
                      "the caller knows what that means for their problem");
        return -1;
    }
    if (!spec.setup && !spec.model) {
        error_message("psopt_solve_robust: supply either spec.setup, which builds the "
                      "augmented problem, or spec.model, which is the nominal problem "
                      "the driver augments itself");
        return -1;
    }
    if (spec.model && !robust_model_is_usable(*spec.model, spec.verbose)) return -1;

    // Lives as long as this call, because the augmented dae and events read it through
    // problem.robust_data on every evaluation.
    RobustAugmentation aug;

    // Free the defect rows the transcription cannot fill. They are equality rows of
    // zeros, IPOPT counts equality rows against variables, and there are nstates of them
    // per phase -- which on an augmented problem is nstates PER SCENARIO, so the count
    // refuses a design with plenty of freedom left. On the two-link arm at 25 nodes it
    // refused above twelve scenarios where the problem has 51 degrees of freedom at
    // every count; freed, the same problem carries eighty. The option is off by default
    // in PSOPT because it moves the iterate path on problems that are already
    // rank-deficient; here it is the difference between a budget of tens and one of
    // hundreds. A caller who has set it explicitly keeps their choice.
    if (!algorithm.free_padded_defect_rows && !spec.keep_padded_defect_rows) {
        algorithm.free_padded_defect_rows = true;
        if (spec.verbose)
            printf("\npsopt_solve_robust: freeing the padded defect rows, which are "
                   "nstates per scenario\n  of counted equalities holding no dynamics. "
                   "Set spec.keep_padded_defect_rows to decline.\n");
    }

    // The starting scenario set: the caller's, if it filled one, otherwise the
    // unscented rule.
    std::vector<RowVectorXd> scenarios = spec.scenarios;
    if (scenarios.empty()) {
        RowVectorXd w;
        if (!robust_sigma_points(spec.uncertainty, scenarios, w)) {
            error_message("psopt_solve_robust: the unscented rule needs a negative "
                          "central weight in this dimension, which is admissible "
                          "in a quadrature and not in a scenario set. Fill "
                          "spec.scenarios with a low-discrepancy set instead.");
            return -1;
        }
    }

    spec.evaluations      = 0;
    spec.n_solves         = 0;
    spec.converged        = false;
    spec.budget_exhausted = false;

    RobustDesign best;                 // the last design that solved cleanly
    int          rc = 0;

    if (spec.verbose) {
        printf("\nRobust design: %d uncertain parameter(s), slack %.4g\n",
               robust_dimension(spec.uncertainty), spec.slack);
        printf("  %3s %4s %13s %13s\n", "it", "M", "objective", "worst found");
    }

    for (int it = 0; it < spec.max_iterations; ++it) {
        robust_do_setup(problem, algorithm, spec, scenarios, best, aug);

        // Before handing it to psopt(): is the guess even a number? A guess that
        // diverged while being built reaches the solver as NaN, and what comes back
        // is PSOPT's constraint-coverage guard reporting a defect in PSOPT. Saying
        // so here, in one sentence naming the array, saves the caller that hunt.
        const char* offender = "";
        if (!robust_guess_is_finite(problem, 1, &offender)) {
            printf("\npsopt_solve_robust: the initial guess for `%s` is not finite "
                   "at %d scenarios.\n  A guess built by integrating a previous "
                   "design can diverge, because the horizon it integrates over is a "
                   "decision variable and the robust final time grows as scenarios "
                   "are added: a step size that was ample at the start need not be "
                   "by the end.\n", offender, (int) scenarios.size());
            if (!best.valid) { solution.error_flag = 1; return 1; }
            scenarios.pop_back();
            robust_do_setup(problem, algorithm, spec, scenarios, best, aug);
            rc = psopt(solution, problem, algorithm);
            ++spec.n_solves;
            break;
        }

        rc = psopt(solution, problem, algorithm);
        ++spec.n_solves;

        if (!robust_trial_succeeded(solution)) {
            const int  dof     = robust_degrees_of_freedom(problem, 1, algorithm);
            const bool starved = (dof <= 0);
            char       buf[1400];
            if (starved)
                robust_dof_message(problem, 1, (int) scenarios.size(), algorithm,
                                   buf, sizeof buf);
            if (!best.valid) {
                if (starved)
                    printf("\npsopt_solve_robust: the first solve failed and the "
                           "transcription had no freedom left.\n  %s\n", buf);
                return rc;
            }
            if (spec.verbose) {
                if (starved)
                    printf("  %3d %4d  cannot carry another scenario.\n  %s\n",
                           it, (int) scenarios.size(), buf);
                else
                    printf("  %3d %4d  adding that scenario made the problem "
                           "unsolvable (NLP return code %d, error flag %d); "
                           "keeping the previous design\n",
                           it, (int) scenarios.size(), solution.nlp_return_code,
                           solution.error_flag);
            }
            // The failed solve has overwritten the caller's Sol, so the previous
            // design has to be put back. Re-solving the previous scenario list
            // warm-started from its own trajectory returns to it: that design is a
            // fixed point of the warm start, which a cold solve of the same set
            // would not be.
            scenarios.pop_back();
            robust_do_setup(problem, algorithm, spec, scenarios, best, aug);
            rc = psopt(solution, problem, algorithm);
            ++spec.n_solves;
            spec.budget_exhausted = starved;
            break;
        }

        best = design_of(solution, problem);

        OracleContext ctx;
        ctx.spec   = &spec;
        ctx.design = &best;
        ctx.calls  = 0;
        RowVectorXd  at;
        const double worst = robust_worst_case(spec.uncertainty, &oracle_bridge,
                                               &ctx, spec.n_seed, spec.n_refine,
                                               spec.seed + (unsigned) it, at);
        spec.evaluations += ctx.calls;

        if (spec.verbose) {
            printf("  %3d %4d %13.6f %13.3e  at", it, (int) scenarios.size(),
                   best.objective, worst);
            for (int i = 0; i < at.size(); ++i) printf(" %.4f", at(i));
            printf("\n");
        }

        spec.certificate    = worst;
        spec.certificate_at = at;

        if (worst <= spec.slack) { spec.converged = true; break; }
        if (it == spec.max_iterations - 1) break;

        scenarios.push_back(at);
    }

    // ---- the polish step -------------------------------------------------------
    //
    // Solve the final scenario set once more from the caller's own guess. Passing an
    // invalid RobustDesign as `previous` is how setup() is told to build its cold
    // guess, which is the same path it takes on the first iteration.
    //
    // The rule is the one the Python driver uses: prefer a design the caller's own
    // violation function certifies, and among certified designs the cheaper. A cold
    // answer that does not certify is not a design and is discarded whatever its
    // objective. See the note on RobustSpec::polish for what this is worth.
    if (spec.polish && best.valid) {
        const RobustDesign warm       = best;
        const double       warm_worst = spec.certificate;
        const RowVectorXd  warm_at    = spec.certificate_at;

        robust_do_setup(problem, algorithm, spec, scenarios, RobustDesign(), aug);
        const int rc_cold = psopt(solution, problem, algorithm);
        ++spec.n_solves;

        bool took_cold = false;
        if (robust_trial_succeeded(solution)) {
            const RobustDesign cold = design_of(solution, problem);

            OracleContext ctx;
            ctx.spec   = &spec;
            ctx.design = &cold;
            ctx.calls  = 0;
            RowVectorXd  cold_at;
            const double cold_worst =
                robust_worst_case(spec.uncertainty, &oracle_bridge, &ctx,
                                  spec.n_seed, spec.n_refine,
                                  spec.seed + 4242u, cold_at);
            spec.evaluations += ctx.calls;

            took_cold = robust_prefer_cold(cold.objective, cold_worst,
                                           warm.objective, warm_worst, spec.slack);

            if (spec.verbose) {
                printf("  polish: a cold solve of the same %d scenarios gives "
                       "%.6f (worst %.3e)\n", (int) scenarios.size(),
                       cold.objective, cold_worst);
                printf("          the warm chain gave %.6f (worst %.3e); keeping "
                       "the %s one\n", warm.objective, warm_worst,
                       took_cold ? "cold" : "warm");
            }

            if (took_cold) {
                best                = cold;
                spec.certificate    = cold_worst;
                spec.certificate_at = cold_at;
                spec.converged      = (cold_worst <= spec.slack);
                rc                  = rc_cold;
            }
        } else if (spec.verbose) {
            printf("  polish: the cold solve of the final %d scenarios failed "
                   "(NLP return code %d); keeping the warm one\n",
                   (int) scenarios.size(), solution.nlp_return_code);
        }

        if (!took_cold) {
            // The cold solve has overwritten the caller's Sol, so the warm design has
            // to be put back. Re-solving its own scenario list warm-started from its
            // own trajectory returns to it, that design being a fixed point of the
            // warm start, which is the same argument the failure path above uses.
            robust_do_setup(problem, algorithm, spec, scenarios, warm, aug);
            rc = psopt(solution, problem, algorithm);
            ++spec.n_solves;
            // spec.design stays the design that was certified, whatever the re-solve
            // returned. Reporting a certificate for one trajectory and handing back
            // another would be the one mistake this whole driver exists to avoid, and
            // the fixed-point argument is an argument and not a guarantee, so the
            // re-solve is checked against it instead of being trusted.
            best                = warm;
            spec.certificate    = warm_worst;
            spec.certificate_at = warm_at;
            if (spec.verbose && robust_trial_succeeded(solution)) {
                const double back = design_of(solution, problem).objective;
                if (fabs(back - warm.objective) > 1.0e-6*(1.0 + fabs(warm.objective)))
                    printf("  polish: putting the warm design back gave %.6f against "
                           "%.6f.\n          spec.design holds the certified one; "
                           "`solution` holds the re-solve.\n", back, warm.objective);
            }
        }
    }

    spec.scenarios = scenarios;
    spec.design    = best;

    if (spec.verbose) {
        printf("\n  scenarios in the final design : %d\n", (int) scenarios.size());
        printf("  objective                     : %.6f\n", best.objective);
        printf("  worst violation over the set  : %.3e\n", spec.certificate);
        printf("  found over                    : %ld evaluations of the caller's\n",
               spec.evaluations);
        printf("                                  independent integrator; nothing\n");
        printf("                                  worse was found, which is not a\n");
        printf("                                  proof that nothing worse exists\n");
        if (spec.budget_exhausted)
            printf("  STOPPED EARLY                 : the transcription could not\n"
                   "                                  carry another scenario\n");
        printf("  calls to psopt()              : %d\n", spec.n_solves);
    }
    return rc;
}

[[nodiscard]] int psopt_solve_robust(Sol& solution, RobustSpec& spec,
                                     RobustModel& model,
                                     Prob& problem, Alg& algorithm)
{
    spec.model = &model;
    return psopt_solve_robust(solution, spec, problem, algorithm);
}
