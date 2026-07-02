// src/integer_parameters.cxx
//
// Driver for optimal control with integer (discrete-valued) static parameters, v1.
//
// psopt_solve_integer enumerates the Cartesian product of the admissible sets of
// all declared integer static parameters, pins each combination (lower = upper =
// value), solves the resulting fixed-parameter OCP with the unmodified psopt(),
// and returns the best by cost. This is the exact solver for the small-integer
// regime and the certified oracle for the v2 convexification method.
//
// psopt() itself is never modified: it treats a pinned integer parameter as an
// ordinary fixed parameter and ignores the integer_parameters records entirely.
//
// Copyright (c) Victor M. Becerra, 2026. Part of the PSOPT library (LGPL).

#include "psopt.h"
#include "integer_parameters.h"

#include <vector>
#include <limits>
#include <cstdio>

// A trial solve is accepted only if psopt() reported no PSOPT-level error and the
// NLP itself converged. error_flag == 0 is not sufficient on its own: IPOPT may
// catch an internal exception and return a failure status while error_flag stays
// 0, so the NLP return code must also indicate success. Ipopt::Solve_Succeeded is
// 0 and Solved_To_Acceptable_Level is 1; both are genuine solutions.
static inline bool trial_succeeded(const Sol& s)
{
    return s.error_flag == 0 && (s.nlp_return_code == 0 || s.nlp_return_code == 1);
}

[[nodiscard]] int psopt_solve_integer(Sol& solution, Prob& problem, Alg& algorithm)
{
    solution.error_flag = 0;

    try {
        // 1. Flatten all declared integer parameters across phases. The pointers
        //    reference the phase records, which are not modified during the loop.
        struct Entry { int iphase; int index; const RowVectorXd* values; };
        std::vector<Entry> entries;
        for (int p = 1; p <= problem.nphases; ++p) {
            std::vector<IntegerParameter>& ip = problem.phases(p).integer_parameters;
            for (size_t k = 0; k < ip.size(); ++k) {
                if (ip[k].values.size() > 0) {
                    Entry e;
                    e.iphase = p;
                    e.index  = ip[k].parameter_index;
                    e.values = &ip[k].values;
                    entries.push_back(e);
                }
            }
        }

        // No integer parameters declared: behave as an ordinary single solve.
        if (entries.empty())
            return psopt(solution, problem, algorithm);

        // 2. Validate parameter indices against the phase parameter count.
        for (size_t d = 0; d < entries.size(); ++d) {
            const int np = problem.phases(entries[d].iphase).nparameters;
            if (entries[d].index < 0 || entries[d].index >= np) {
                char msg[256];
                snprintf(msg, sizeof(msg),
                    "psopt_solve_integer: integer parameter_index %d out of range for phase %d (nparameters = %d)",
                    entries[d].index, entries[d].iphase, np);
                error_message(msg);
            }
        }

        // 3. Combination-count guard (the small-integer-regime boundary).
        const long ncomb = integer_parameter_combination_count(problem);
        if (ncomb > (long) algorithm.max_integer_combinations) {
            char msg[320];
            snprintf(msg, sizeof(msg),
                "psopt_solve_integer: %ld integer-parameter combinations exceed algorithm.max_integer_combinations (%d); raise the limit or reduce the admissible sets",
                ncomb, algorithm.max_integer_combinations);
            error_message(msg);
        }

        // 4. Save the original bounds of the affected parameters.
        std::vector<double> save_lo(entries.size()), save_hi(entries.size());
        for (size_t d = 0; d < entries.size(); ++d) {
            save_lo[d] = problem.phases(entries[d].iphase).bounds.lower.parameters(entries[d].index);
            save_hi[d] = problem.phases(entries[d].iphase).bounds.upper.parameters(entries[d].index);
        }

        // 5. Enumerate the Cartesian product and record the best cleanly-solved
        //    combination by cost. No Sol is copied: sol_str owns raw arrays and
        //    frees them in its destructor but defines no copy constructor or copy
        //    assignment (a rule-of-three gap), so copying a Sol is unsafe. Instead
        //    we keep only the winning index vector and re-solve it once at the end.
        std::vector<int> sizes(entries.size());
        for (size_t d = 0; d < entries.size(); ++d)
            sizes[d] = static_cast<int>(entries[d].values->size());
        std::vector<int> idx(entries.size(), 0);
        std::vector<int> best_idx(entries.size(), 0);

        bool   have_best = false;
        double best_cost = std::numeric_limits<double>::infinity();

        do {
            for (size_t d = 0; d < entries.size(); ++d) {
                const double v = (*entries[d].values)(idx[d]);
                problem.phases(entries[d].iphase).bounds.lower.parameters(entries[d].index) = v;
                problem.phases(entries[d].iphase).bounds.upper.parameters(entries[d].index) = v;
            }

            Sol trial;
            (void) psopt(trial, problem, algorithm);

            if (trial_succeeded(trial) && (!have_best || trial.cost < best_cost)) {
                best_cost = trial.cost;
                best_idx  = idx;
                have_best = true;
            }
        } while (integer_cartesian_advance(idx, sizes));

        // 6. Re-solve the winning combination (or the first, if none solved) into
        //    the caller's solution, then restore the original bounds. The re-solve
        //    is deterministic, so it reproduces the recorded best cost exactly.
        const std::vector<int>& chosen = have_best ? best_idx : idx;   // idx is all-zero here
        for (size_t d = 0; d < entries.size(); ++d) {
            const double v = (*entries[d].values)(chosen[d]);
            problem.phases(entries[d].iphase).bounds.lower.parameters(entries[d].index) = v;
            problem.phases(entries[d].iphase).bounds.upper.parameters(entries[d].index) = v;
        }
        int ret_rc = psopt(solution, problem, algorithm);

        for (size_t d = 0; d < entries.size(); ++d) {
            problem.phases(entries[d].iphase).bounds.lower.parameters(entries[d].index) = save_lo[d];
            problem.phases(entries[d].iphase).bounds.upper.parameters(entries[d].index) = save_hi[d];
        }

        return ret_rc;
    }
    catch (ErrorHandler& /*handler*/) {
        // Mirror psopt(): report set-up/guard failures through solution.error_flag
        // rather than propagating out to the caller's main().
        solution.error_flag = 1;
        return solution.error_flag;
    }
}
