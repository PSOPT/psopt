// ---------------------------------------------------------------------------
// src/integer_controls.cxx
//
// Engine-side support for the ergonomic integer-control declaration API. When a
// phase declares a single integer control (via declare_integer_control), the
// IntegerControlExpansionGuard, constructed at the top of psopt(), performs the
// outer convexification transparently:
//
//   * it stashes the user's dae / integrand_cost callbacks and installs internal
//     wrappers in their place;
//   * it expands the phase control count to (nu-1)+M (the integer control slot
//     is dropped and M weight-controls are appended), adds one SOS1 equality
//     path constraint, and rewrites the control/path bounds, scale and guess to
//     match; and
//   * on destruction it restores the user layout (the solution is left in the
//     weights layout).
//
// The wrappers reconstruct the user control vector (continuous controls pass
// through; the integer slot is set to each admissible value in turn), call the
// user callback once per mode, and assemble the convex combination using the
// convex_combine / sos1_residual helpers.
//
// Scope: any number of integer controls per phase. With K of them, of sizes
// M_1..M_K, the convexification is over the Cartesian product of the admissible
// sets: the K integer slots are dropped from the control vector and P = prod(M_k)
// weight-controls are appended, with a single SOS1 equality over the P weights.
// The product, rather than a weight vector per control, is what keeps the
// relaxation convex in the weights and therefore tight. User path constraints are
// assumed independent of the integer controls (the path is taken from a
// representative mode); the integrand is convex-combined, so an integer-dependent
// running cost is handled correctly. See include/integer_controls.h.
// ---------------------------------------------------------------------------

#include "integer_controls.h"
#include <vector>
#include <cmath>
#include <cstdio>

using std::vector;

// --- layout helpers (during the solve, phase counts are the expanded ones) ---
namespace {

// The dae wrapper installed as problem.dae while an integer control is active.
void ic_dae_wrapper(adouble* derivatives, adouble* path, adouble* states,
                    adouble* controls, adouble* parameters, adouble& time,
                    adouble* xad, int iphase, Workspace* workspace)
{
    Prob&   problem = *workspace->problem;
    Phases& ph      = problem.phases(iphase);
    const std::vector<IntegerControl>& ic = ph.integer_controls;
    const int K = (int) ic.size();

    if (K <= 0) {   // this phase has no integer control: pass straight through
        problem.user_dae(derivatives, path, states, controls, parameters, time,
                         xad, iphase, workspace);
        return;
    }

    const int P          = integer_control_product_size(ic);
    const int nstates    = ph.nstates;
    const int ncont_int  = ph.ncontrols;           // expanded
    const int nu         = ncont_int - P + K;      // user control count
    const int npath_user = ph.npath - 1;           // user path count (SOS1 appended)

    const adouble* w = controls + (ncont_int - P); // weights = trailing P controls

    // Reconstruct the user control vector; the continuous controls pass through in
    // order, and the K integer slots are filled per product mode below.
    vector<bool> is_int(nu, false);
    for (int k = 0; k < K; ++k) is_int[ ic[k].control_index ] = true;
    vector<adouble> u_user(nu);
    {
        int src = 0;
        for (int j = 0; j < nu; ++j)
            if (!is_int[j]) u_user[j] = controls[src++];
    }

    // Evaluate the user dae at each product mode and convex-combine.
    vector<vector<adouble> > md(P, vector<adouble>(nstates));
    vector<const adouble*>   mdptr(P);
    vector<adouble>          path_i(npath_user > 0 ? npath_user : 1);
    vector<int>              modes;

    for (int p = 0; p < P; ++p) {
        integer_control_decode(ic, p, modes);
        for (int k = 0; k < K; ++k)
            u_user[ ic[k].control_index ] = ic[k].values( modes[k] );
        problem.user_dae(md[p].data(), path_i.data(), states, u_user.data(),
                         parameters, time, xad, iphase, workspace);
        mdptr[p] = md[p].data();
        if (p == 0)
            for (int q = 0; q < npath_user; ++q) path[q] = path_i[q];
    }

    convex_combine<adouble>(derivatives, nstates, P, w, mdptr.data());
    path[npath_user] = sos1_residual<adouble>(w, P);   // SOS1 equality (= 1)
}

// The integrand_cost wrapper installed as problem.integrand_cost while active.
adouble ic_integrand_wrapper(adouble* states, adouble* controls, adouble* parameters,
                             adouble& time, adouble* xad, int iphase, Workspace* workspace)
{
    Prob&   problem = *workspace->problem;
    Phases& ph      = problem.phases(iphase);
    const std::vector<IntegerControl>& ic = ph.integer_controls;
    const int K = (int) ic.size();

    if (K <= 0)
        return problem.user_integrand_cost(states, controls, parameters, time,
                                           xad, iphase, workspace);

    const int P         = integer_control_product_size(ic);
    const int ncont_int = ph.ncontrols;
    const int nu        = ncont_int - P + K;
    const adouble* w    = controls + (ncont_int - P);

    vector<bool> is_int(nu, false);
    for (int k = 0; k < K; ++k) is_int[ ic[k].control_index ] = true;
    vector<adouble> u_user(nu);
    {
        int src = 0;
        for (int j = 0; j < nu; ++j)
            if (!is_int[j]) u_user[j] = controls[src++];
    }

    adouble total = 0.0;
    vector<int> modes;
    for (int p = 0; p < P; ++p) {
        integer_control_decode(ic, p, modes);
        for (int k = 0; k < K; ++k)
            u_user[ ic[k].control_index ] = ic[k].values( modes[k] );
        total += w[p] * problem.user_integrand_cost(states, u_user.data(),
                                                    parameters, time, xad, iphase, workspace);
    }
    return total;
}

} // namespace

// ---------------------------------------------------------------------------
// Guard: expand on construction, restore on destruction.
// ---------------------------------------------------------------------------

IntegerControlExpansionGuard::IntegerControlExpansionGuard(Prob& problem)
    : problem_(problem), active_(false)
{
    const int nphases = problem.nphases;
    phase_active_.assign(nphases, false);

    bool any = false;
    for (int i = 1; i <= nphases; ++i)
        if (!problem.phases(i).integer_controls.empty()) any = true;
    if (!any) return;                       // complete no-op

    active_ = true;

    // Stash the user callbacks.
    problem.user_dae            = problem.dae;
    problem.user_integrand_cost = problem.integrand_cost;

    // Saved-state vectors are 1-indexed by phase for convenience.
    saved_ncontrols_.assign(nphases + 1, 0);
    saved_npath_.assign(nphases + 1, 0);
    saved_lower_controls_.resize(nphases + 1);
    saved_upper_controls_.resize(nphases + 1);
    saved_scale_controls_.resize(nphases + 1);
    saved_lower_path_.resize(nphases + 1);
    saved_upper_path_.resize(nphases + 1);
    saved_scale_path_.resize(nphases + 1);
    saved_guess_controls_.resize(nphases + 1);

    for (int i = 1; i <= nphases; ++i) {
        Phases& ph = problem.phases(i);
        const std::vector<IntegerControl>& ic = ph.integer_controls;
        const int K = (int) ic.size();
        if (K <= 0) { phase_active_[i - 1] = false; continue; }
        phase_active_[i - 1] = true;

        const int nu         = ph.ncontrols;      // user count before expansion
        const int npath_user = ph.npath;

        for (int k = 0; k < K; ++k) {
            if (ic[k].control_index < 0 || ic[k].control_index >= nu)
                error_message("declare_integer_control: control_index out of range");
            if (ic[k].values.size() < 1)
                error_message("declare_integer_control: an empty admissible set was declared");
            for (int j = 0; j < k; ++j)
                if (ic[j].control_index == ic[k].control_index)
                    error_message("declare_integer_control: the same control index was declared twice");
        }
        if (K > nu)
            error_message("declare_integer_control: more integer controls declared than the phase has controls");

        const int P = integer_control_product_size(ic);

        // Save the user layout for restoration.
        saved_ncontrols_[i]      = nu;
        saved_npath_[i]          = npath_user;
        saved_lower_controls_[i] = ph.bounds.lower.controls;
        saved_upper_controls_[i] = ph.bounds.upper.controls;
        saved_scale_controls_[i] = ph.scale.controls;
        saved_lower_path_[i]     = ph.bounds.lower.path;
        saved_upper_path_[i]     = ph.bounds.upper.path;
        saved_scale_path_[i]     = ph.scale.path;
        saved_guess_controls_[i] = ph.guess.controls;

        const int ncont_int = (nu - K) + P;
        const int npath_int = npath_user + 1;

        std::vector<bool> is_int(nu, false);
        for (int k = 0; k < K; ++k) is_int[ ic[k].control_index ] = true;

        if (K == 1) {
            printf("PSOPT: phase %d has an integer control (index %d, %d values). Its "
                   "control bounds are ignored; admissible values come from "
                   "declare_integer_control and %d weights are bounded to [0,1].\n",
                   i, ic[0].control_index, (int) ic[0].values.size(), P);
        } else {
            printf("PSOPT: phase %d has %d integer controls; the convexification is over "
                   "the product of their admissible sets, giving %d weights bounded to "
                   "[0,1] and one SOS1 constraint. Their control bounds are ignored.\n",
                   i, K, P);
        }

        // Control bounds/scale: drop the K integer slots, append P weight entries [0,1].
        MatrixXd lc(ncont_int, 1), uc(ncont_int, 1), sc(ncont_int, 1);
        int idx = 0;
        for (int k = 0; k < nu; ++k) {
            if (is_int[k]) continue;
            lc(idx, 0) = ph.bounds.lower.controls(k, 0);
            uc(idx, 0) = ph.bounds.upper.controls(k, 0);
            sc(idx, 0) = ph.scale.controls(k, 0);
            ++idx;
        }
        for (int j = 0; j < P; ++j) {
            lc(idx, 0) = 0.0; uc(idx, 0) = 1.0; sc(idx, 0) = 1.0; ++idx;
        }
        ph.bounds.lower.controls = lc;
        ph.bounds.upper.controls = uc;
        ph.scale.controls        = sc;

        // Path bounds/scale: append SOS1 equality [1,1] (scale 1).
        MatrixXd lp(npath_int, 1), up(npath_int, 1), sp(npath_int, 1);
        for (int p = 0; p < npath_user; ++p) {
            lp(p, 0) = ph.bounds.lower.path(p, 0);
            up(p, 0) = ph.bounds.upper.path(p, 0);
            sp(p, 0) = ph.scale.path(p, 0);
        }
        lp(npath_user, 0) = 1.0; up(npath_user, 0) = 1.0; sp(npath_user, 0) = 1.0;
        ph.bounds.lower.path = lp;
        ph.bounds.upper.path = up;
        ph.scale.path        = sp;

        // Guess: map each integer row to its nearest admissible mode, combine the K mode
        // indices into a product index, and one-hot expand over the P product modes.
        if (ph.guess.controls.size() > 0) {
            const int N = (int) ph.guess.controls.cols();
            RowVectorXi modeseq(N);
            for (int n = 0; n < N; ++n) {
                int prod = 0;
                for (int k = 0; k < K; ++k) {
                    const int Mk = (int) ic[k].values.size();
                    const double v = ph.guess.controls(ic[k].control_index, n);
                    int    b  = 0;
                    double bd = std::fabs(v - ic[k].values(0));
                    for (int j = 1; j < Mk; ++j) {
                        double d = std::fabs(v - ic[k].values(j));
                        if (d < bd) { bd = d; b = j; }
                    }
                    prod = prod * Mk + b;      // first declared control is the leading digit
                }
                modeseq(n) = prod;
            }
            MatrixXd Wg;
            integer_guess_to_weights(P, modeseq, Wg);   // P x N one-hot

            MatrixXd g(ncont_int, N);
            int r = 0;
            for (int k = 0; k < nu; ++k) {
                if (is_int[k]) continue;
                g.row(r) = ph.guess.controls.row(k);
                ++r;
            }
            for (int j = 0; j < P; ++j) { g.row(r) = Wg.row(j); ++r; }
            ph.guess.controls = g;
        }

        ph.ncontrols = ncont_int;
        ph.npath     = npath_int;
    }

    // Install the wrappers (integrand only if the user actually supplied one).
    problem.dae = &ic_dae_wrapper;
    if (problem.user_integrand_cost != nullptr)
        problem.integrand_cost = &ic_integrand_wrapper;
}

IntegerControlExpansionGuard::~IntegerControlExpansionGuard()
{
    if (!active_) return;

    const int nphases = problem_.nphases;
    for (int i = 1; i <= nphases; ++i) {
        if (!phase_active_[i - 1]) continue;
        Phases& ph = problem_.phases(i);
        ph.ncontrols             = saved_ncontrols_[i];
        ph.npath                 = saved_npath_[i];
        ph.bounds.lower.controls = saved_lower_controls_[i];
        ph.bounds.upper.controls = saved_upper_controls_[i];
        ph.scale.controls        = saved_scale_controls_[i];
        ph.bounds.lower.path     = saved_lower_path_[i];
        ph.bounds.upper.path     = saved_upper_path_[i];
        ph.scale.path            = saved_scale_path_[i];
        ph.guess.controls        = saved_guess_controls_[i];
    }

    problem_.dae                 = problem_.user_dae;
    problem_.integrand_cost      = problem_.user_integrand_cost;
    problem_.user_dae            = nullptr;
    problem_.user_integrand_cost = nullptr;
}
