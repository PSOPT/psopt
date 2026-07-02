// include/integer_parameters.h
//
// Optimal control with integer (discrete-valued) static parameters, v1.
//
// A phase may declare one or more static parameters as integer (discrete-valued)
// by giving, for each, an admissible finite set of values. The dedicated driver
// psopt_solve_integer (src/integer_parameters.cxx) then enumerates the Cartesian
// product of all declared sets, solves the fixed-parameter OCP for each
// combination with the unmodified psopt(), and returns the best by cost. This is
// the exact solver for the small-integer-space regime; it is also the certified
// oracle against which the v2 convexification method is validated.
//
// This header holds the user-facing API (declare/reconstruct) and the pure,
// unit-testable building blocks (Cartesian-product enumerator, combination count,
// nearest-admissible snap). The driver itself is declared here and defined in
// src/integer_parameters.cxx.
//
// Copyright (c) Victor M. Becerra, 2026. Part of the PSOPT library (LGPL).

#ifndef PSOPT_INTEGER_PARAMETERS_H
#define PSOPT_INTEGER_PARAMETERS_H

#include "psopt.h"

#include <vector>
#include <cmath>
#include <limits>

//////////////////////////////////////////////////////////////////////////
///////////////////  Declaration  ////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////

// Declare a static parameter of a phase as integer (discrete-valued), giving the
// admissible finite set. A phase may declare several; each records the parameter's
// index in the phase parameter layout and its admissible values. Dormant until
// psopt_solve_integer is used - psopt() itself ignores these records.
inline void declare_integer_parameter(Prob& problem, int iphase,
                                       int parameter_index, const RowVectorXd& values)
{
    IntegerParameter ip;
    ip.parameter_index = parameter_index;
    ip.values          = values;
    problem.phases(iphase).integer_parameters.push_back(ip);
}

//////////////////////////////////////////////////////////////////////////
///////////////////  Building blocks (pure, unit-testable)  ///////////////
//////////////////////////////////////////////////////////////////////////

// Index of the admissible value closest to `value` (ties resolve to the lower
// index). Returns -1 if the value set is empty.
inline int nearest_admissible_index(double value, const RowVectorXd& values)
{
    const int M = static_cast<int>(values.size());
    if (M <= 0) return -1;
    int    best  = 0;
    double bestd = std::fabs(values(0) - value);
    for (int i = 1; i < M; ++i) {
        double d = std::fabs(values(i) - value);
        if (d < bestd) { bestd = d; best = i; }
    }
    return best;
}

// Product of the sizes of all declared integer-parameter value sets, across all
// phases, saturating at PSOPT_INTEGER_COMBINATIONS_SATURATE to avoid overflow so
// callers can safely compare against algorithm.max_integer_combinations. Empty
// value sets are skipped defensively. Returns 1 when nothing is declared (a single
// empty combination).
static const long PSOPT_INTEGER_COMBINATIONS_SATURATE = 1L << 40;   // ~1.1e12

inline long integer_parameter_combination_count(Prob& problem)
{
    long count = 1;
    for (int p = 1; p <= problem.nphases; ++p) {
        const std::vector<IntegerParameter>& ip = problem.phases(p).integer_parameters;
        for (size_t k = 0; k < ip.size(); ++k) {
            long m = static_cast<long>(ip[k].values.size());
            if (m <= 0) continue;
            if (count > PSOPT_INTEGER_COMBINATIONS_SATURATE / m)
                return PSOPT_INTEGER_COMBINATIONS_SATURATE;
            count *= m;
        }
    }
    return count;
}

// Advance a mixed-radix index vector `idx` in place given the per-digit radices
// `sizes` (idx[d] in [0, sizes[d])). The least-significant digit is idx[0].
// Returns true if there is a next combination, and false once the counter wraps
// back to all zeros (enumeration complete). A zero-length input returns false,
// i.e. a single (empty) combination.
inline bool integer_cartesian_advance(std::vector<int>& idx, const std::vector<int>& sizes)
{
    const int D = static_cast<int>(sizes.size());
    for (int d = 0; d < D; ++d) {
        if (++idx[d] < sizes[d]) return true;
        idx[d] = 0;
    }
    return false;
}

//////////////////////////////////////////////////////////////////////////
///////////////////  Reconstruction  /////////////////////////////////////
//////////////////////////////////////////////////////////////////////////

// Selected value of a declared integer parameter after a solve.
struct IntegerParameterReconstruction {
    double value;   // selected admissible value (NaN if not declared / unavailable)
    int    index;   // index into the admissible set (-1 if not declared)
};

// Read the selected value of a declared integer parameter from a solved problem.
// psopt_solve_integer fixes the parameter to one admissible value, so the solved
// solution carries it exactly; this reads it back and snaps to the nearest
// admissible value as a safety net. Returns {NaN, -1} if the phase declares no
// integer parameter at parameter_index.
inline IntegerParameterReconstruction
reconstruct_integer_parameter(Sol& solution, Prob& problem, int iphase, int parameter_index)
{
    IntegerParameterReconstruction r;
    r.value = std::numeric_limits<double>::quiet_NaN();
    r.index = -1;

    const std::vector<IntegerParameter>& ip = problem.phases(iphase).integer_parameters;
    const RowVectorXd* values = 0;
    for (size_t k = 0; k < ip.size(); ++k)
        if (ip[k].parameter_index == parameter_index && ip[k].values.size() > 0) {
            values = &ip[k].values; break;
        }
    if (!values) return r;

    MatrixXd psol   = solution.get_parameters_in_phase(iphase);
    double   solved = psol(parameter_index);
    r.index = nearest_admissible_index(solved, *values);
    r.value = (r.index >= 0) ? (*values)(r.index) : solved;
    return r;
}

//////////////////////////////////////////////////////////////////////////
///////////////////  Driver (defined in src/integer_parameters.cxx)  //////
//////////////////////////////////////////////////////////////////////////

// Dedicated driver for optimal control with integer (discrete-valued) static
// parameters (v1). Enumerates the Cartesian product of all declared admissible
// sets, solves the fixed-parameter OCP for each combination with the unmodified
// psopt(), and returns the best by cost in `solution`. Aborts (via error_message)
// if the number of combinations exceeds algorithm.max_integer_combinations.
// Returns the psopt() return code of the best solve (or of the last attempt if
// none succeeded).
[[nodiscard]] int psopt_solve_integer(Sol& solution, Prob& problem, Alg& algorithm);

#endif // PSOPT_INTEGER_PARAMETERS_H
