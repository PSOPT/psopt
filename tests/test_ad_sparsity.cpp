//////////////////////////////////////////////////////////////////////////////
// test_ad_sparsity.cpp
//
// The cached sparsity pattern held by an AD handle, and the one condition under
// which it may outlive the tape it was computed from.
//
// psopt_ad::ad_record() replaces the tape in a handle. The sparsity pattern that
// goes with it is deliberately not discarded, because one caller depends on its
// survival: the IPOPT interface declares the Hessian's structure once, at
// get_nlp_info(), and then re-records the lambda-weighted tape at every single
// evaluation of the same problem. Throwing the pattern away there would mean a
// fresh sparsity sweep per evaluation.
//
// The case that must not be treated the same way is a change of shape. Every mesh
// refinement iteration gives the NLP a new number of variables and constraints, and
// a pattern computed for the previous mesh is then not a stale pattern for the same
// matrix, it is the pattern of a different matrix. CppAD is asked for the values at
// the positions the old pattern names, obtains them from the new tape, and returns
// them without complaint: a Jacobian of the wrong shape, wrong everywhere, and
// carrying no sign at all that anything is amiss. The symptom is a solver that
// converges on the first mesh and then behaves on the second as though the problem
// were badly conditioned, which is a diagnosis that leads a long way from the cause.
//
// So the rule is that the cache survives a re-record of the same shape and not a
// re-record of a different one, and it is the rule rather than any one caller's use
// of it that is worth pinning here.
//////////////////////////////////////////////////////////////////////////////

#include "gtest/gtest.h"
#include <psopt.h>
#include <cmath>
#include <vector>

#if PSOPT_AD_BACKEND == PSOPT_AD_CPPAD

namespace {

// f: R^3 -> R^2.  Its Jacobian has three non-zero entries.
void fun_a(const adouble* x, adouble* y)
{
    y[0] = x[0]*x[1];
    y[1] = x[2]*x[2];
}

// g: R^4 -> R^3, deliberately of a different shape and a different pattern, with a
// column (the last) that the pattern of f does not reach at all.
void fun_b(const adouble* x, adouble* y)
{
    y[0] = x[0] + 2.0*x[3];
    y[1] = x[1]*x[3];
    y[2] = x[2]*x[2] + x[3];
}

// The triplet as a dense matrix, so that a missing entry reads as a zero rather than
// as a shorter list.
std::vector<double> dense(const psopt_ad::SparseTriplet& T, int m, int n)
{
    std::vector<double> A((size_t)m*n, 0.0);
    for (int k = 0; k < T.nnz(); k++) A[(size_t)T.row[k]*n + T.col[k]] = T.val[k];
    return A;
}

} // namespace


TEST(AdDriver, SparsityCacheSurvivesARecordOfTheSameShape)
{
    psopt_ad::ADHandle h;
    const double xa[3] = {2.0, 3.0, 5.0};

    psopt_ad::ad_record(h, 3, 2, xa, fun_a);
    ASSERT_FALSE(h.jac_struct) << "a fresh handle cannot have a pattern yet";

    psopt_ad::ad_sparse_jacobian(h, xa, /*reuse=*/false);
    ASSERT_TRUE(h.jac_struct);

    // The same shape again: this is the IPOPT case, and the pattern must be kept.
    psopt_ad::ad_record(h, 3, 2, xa, fun_a);
    EXPECT_TRUE(h.jac_struct);
}

TEST(AdDriver, SparsityCacheIsDroppedWhenTheShapeChanges)
{
    psopt_ad::ADHandle h;

    const double xa[3] = {2.0, 3.0, 5.0};
    psopt_ad::ad_record(h, 3, 2, xa, fun_a);
    psopt_ad::ad_sparse_jacobian(h, xa, /*reuse=*/false);
    ASSERT_TRUE(h.jac_struct);

    const double xb[4] = {1.0, 2.0, 3.0, 4.0};
    psopt_ad::ad_record(h, 4, 3, xb, fun_b);
    EXPECT_FALSE(h.jac_struct) << "the pattern of the previous shape must not survive";

    // A caller that asks to reuse the pattern must still be given the right Jacobian:
    // there is nothing left to reuse, so the structure is recomputed.
    psopt_ad::SparseTriplet T = psopt_ad::ad_sparse_jacobian(h, xb, /*reuse=*/true);
    std::vector<double> J = dense(T, 3, 4);

    const double expect[12] = { 1.0, 0.0, 0.0, 2.0,
                                0.0, 4.0, 0.0, 2.0,
                                0.0, 0.0, 6.0, 1.0 };
    for (int k = 0; k < 12; k++)
        EXPECT_NEAR(J[k], expect[k], 1.0e-12) << "entry " << k/4 << "," << k%4;
}

// The same thing stated as the property that actually matters, without reference to
// the cache: two handles that have seen different histories must agree.
TEST(AdDriver, JacobianAfterAReshapeMatchesOneTakenFromAFreshHandle)
{
    const double xb[4] = {1.0, 2.0, 3.0, 4.0};

    psopt_ad::ADHandle fresh;
    psopt_ad::ad_record(fresh, 4, 3, xb, fun_b);
    std::vector<double> Jf = dense(psopt_ad::ad_sparse_jacobian(fresh, xb, false), 3, 4);

    psopt_ad::ADHandle reused;
    const double xa[3] = {2.0, 3.0, 5.0};
    psopt_ad::ad_record(reused, 3, 2, xa, fun_a);
    psopt_ad::ad_sparse_jacobian(reused, xa, false);
    psopt_ad::ad_record(reused, 4, 3, xb, fun_b);
    std::vector<double> Jr = dense(psopt_ad::ad_sparse_jacobian(reused, xb, true), 3, 4);

    for (int k = 0; k < 12; k++) EXPECT_NEAR(Jr[k], Jf[k], 1.0e-12) << "entry " << k;
}

#endif // PSOPT_AD_BACKEND == PSOPT_AD_CPPAD
