/*********************************************************************************************

This file is part of the PSOPT library, a software tool for computational optimal control

Copyright (C) 2009-2026 Victor M. Becerra

This library is free software; you can redistribute it and/or
modify it under the terms of the GNU Lesser General Public
License as published by the Free Software Foundation; either
version 2.1 of the License, or (at your option) any later version.

This library is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public
License along with this library; if not, write to the Free Software
Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA,
or visit http://www.gnu.org/licenses/

Author:    Professor Victor M. Becerra
Address:   University of Portsmouth
           School of Electrical and Mechanical Engineering
           Portsmouth PO1 3DJ
           United Kingdom
e-mail:    v.m.becerra@ieee.org

**********************************************************************************************/

// ============================================================================
//  psopt_ad_backend.hpp  —  pluggable AD backend (compile-time switch)
//  -DPSOPT_AD_BACKEND=PSOPT_AD_{CPPAD|AUTODIFF|XAD}; default CppAD.
//
//  CppAD (in-memory tapes) is the production backend. It replaced ADOL-C, whose
//  multi-gigabyte on-disk tapes corrupt/exhaust the disk at high node counts.
//  AUTODIFF and XAD are forward-mode cross-check backends (diagnostic only).
// ============================================================================
#ifndef PSOPT_AD_BACKEND_HPP
#define PSOPT_AD_BACKEND_HPP
#include <cmath>
#include <complex>

#define PSOPT_AD_CPPAD    2
#define PSOPT_AD_AUTODIFF 3
#define PSOPT_AD_XAD      4
#ifndef PSOPT_AD_BACKEND
#  define PSOPT_AD_BACKEND PSOPT_AD_CPPAD
#endif

#if   PSOPT_AD_BACKEND == PSOPT_AD_CPPAD
#  include <cppad/cppad.hpp>
   // adouble subclasses AD<double> so PSOPT's `x.value()` member calls keep compiling.
   // Layout-identical to AD<double> (single public base, no data members, no virtuals),
   // which lets ad_record() hand a reinterpret_cast<adouble*> view of an AD<double> taping
   // buffer (CppAD's vector type checks reject a subclass element type).
   class adouble : public CppAD::AD<double> {
   public:
     typedef CppAD::AD<double> ADBase;
     adouble() : ADBase() {}
     adouble(const ADBase& x) : ADBase(x) {}
     adouble(double x) : ADBase(x) {}
     adouble(int x) : ADBase(double(x)) {}
     double value() const { return CppAD::Value(CppAD::Var2Par(static_cast<const ADBase&>(*this))); }
   };
#  define PSOPT_AD_TAPED 1
#  define PSOPT_AD_NAME "CppAD"
#elif PSOPT_AD_BACKEND == PSOPT_AD_AUTODIFF
#  include <autodiff/forward/dual.hpp>
   typedef autodiff::dual adouble;
#  define PSOPT_AD_TAPED 0
#  define PSOPT_AD_NAME "autodiff(fwd)"
#elif PSOPT_AD_BACKEND == PSOPT_AD_XAD
#  include <XAD/XAD.hpp>
   // Forward-mode XAD. adouble subclasses xad::FReal<double> to add the `.value()` member
   // (so PSOPT's existing call sites keep compiling unchanged). XAD's expression templates
   // accept the derived type, so operators/elementary functions resolve normally.
   class adouble : public xad::FReal<double> {
   public:
     typedef xad::FReal<double> Base;
     using Base::FReal;        // inherit constructors (incl. from XAD expressions)
     using Base::operator=;    // inherit assignment from XAD expressions
     adouble() : Base() {}
     adouble(const Base& x) : Base(x) {}
     adouble(double x) : Base(x) {}
     adouble(int x) : Base(double(x)) {}
     double value() const { return xad::value(static_cast<const Base&>(*this)); }
   };
#  define PSOPT_AD_TAPED 0
#  define PSOPT_AD_NAME "XAD(fwd)"
#else
#  error "Unknown PSOPT_AD_BACKEND (expected CPPAD, AUTODIFF, or XAD)"
#endif

#if PSOPT_AD_BACKEND == PSOPT_AD_CPPAD
// adouble is a subclass of AD<double>; the derived->base conversion otherwise ties CppAD's
// pow overloads. Exact (adouble, scalar) matches win and resolve the ambiguity.
inline adouble pow(const adouble& x, double y){ return CppAD::pow(static_cast<const CppAD::AD<double>&>(x), CppAD::AD<double>(y)); }
inline adouble pow(const adouble& x, int    y){ return CppAD::pow(static_cast<const CppAD::AD<double>&>(x), CppAD::AD<double>((double)y)); }
#endif

#if PSOPT_AD_BACKEND == PSOPT_AD_AUTODIFF
// autodiff exposes abs() but not the C-style fabs() that PSOPT calls on adouble.
inline adouble fabs(const adouble& x){ using std::abs; return abs(x); }   // ADL finds autodiff::detail::abs
#endif

// ---- Guarded sqrt --------------------------------------------------------------------------
// ADOL-C defined d/dx sqrt(x)|_0 = 0 (a finite subgradient). CppAD / XAD / autodiff instead
// evaluate u'/(2*sqrt(u)) literally, so at u=0 they produce 0/0 = NaN. This bites whenever an
// initial guess sits on a coordinate singularity, e.g. X = sqrt(h*h+k*k) at an equatorial guess
// (h=k=0), where the single NaN partial poisons the constraint Jacobian and IPOPT fails.
//
// Adding a tiny epsilon to the radicand removes the singularity from the *derivative*: the inner
// sqrt is never evaluated at exactly 0, so no Inf/NaN Taylor coefficient is ever recorded, and
// the partial becomes 0 wherever the inner derivative is 0 (the correct subgradient). Because
// 1e-30 is far below the ULP of any argument >~1e-16, both the value and the derivative are
// bit-identical to sqrt(x) for every non-singular input, so results are otherwise unchanged.
// (A CondExp-style guard does NOT work: CppAD still tapes sqrt(0) and its NaN forward
//  coefficient trips CppAD's nan check even when the branch is not selected.)
#if   PSOPT_AD_BACKEND == PSOPT_AD_CPPAD
// adouble arithmetic returns the CppAD::AD<double> base, so user radicands such as (h*h+k*k)
// have type AD<double>, not adouble. Intercept at the base type: a non-template overload is
// preferred over CppAD's template sqrt, so this captures every recorded sqrt without ambiguity.
// The NaN must die *inside* sqrt (so a singular term contributes 0 to a larger sum) --
// sanitising the assembled Jacobian afterwards would wrongly zero the whole entry.
inline ::CppAD::AD<double> sqrt(const ::CppAD::AD<double>& x){
    return ::CppAD::sqrt( x + ::CppAD::AD<double>(1.0e-30) );
}
#elif PSOPT_AD_BACKEND == PSOPT_AD_XAD
// Forward-mode ET backends: this intercepts radicands that are a concrete adouble (e.g.
// sqrt(dot(v,v))). Inline-expression radicands (sqrt(h*h+k*k)) bind to the library's own
// sqrt() and would need model-level regularisation; these backends are diagnostic only.
inline adouble sqrt(const adouble& x){
    using std::sqrt;   // ADL also reaches xad::sqrt for the FReal expression
    return adouble( sqrt( static_cast<const adouble::Base&>(x) + adouble::Base(1.0e-30) ) );
}
#elif PSOPT_AD_BACKEND == PSOPT_AD_AUTODIFF
inline adouble sqrt(const adouble& x){
    using std::sqrt;   // ADL finds autodiff::detail::sqrt for the (dual + double) expression
    return sqrt( x + 1.0e-30 );
}
#endif

namespace psopt_ad {
#if   PSOPT_AD_BACKEND == PSOPT_AD_CPPAD
inline double value(const adouble& x){ return CppAD::Value(CppAD::Var2Par(x)); }
#elif PSOPT_AD_BACKEND == PSOPT_AD_AUTODIFF
inline double value(const adouble& x){ return autodiff::val(x); }
#elif PSOPT_AD_BACKEND == PSOPT_AD_XAD
inline double value(const adouble& x){ return xad::value(x); }
#endif
inline double value(double x){ return x; }
inline double value(const std::complex<double>& x){ return x.real(); }

#if !PSOPT_AD_TAPED
// Forward-mode directional-derivative seeding/extraction (one input direction at a time).
#  if   PSOPT_AD_BACKEND == PSOPT_AD_XAD
inline void   seed_deriv(adouble& x, double s){ xad::derivative(static_cast<adouble::Base&>(x)) = s; }
inline double get_deriv (const adouble& y){ return xad::derivative(static_cast<const adouble::Base&>(y)); }
#  elif PSOPT_AD_BACKEND == PSOPT_AD_AUTODIFF
inline void   seed_deriv(adouble& x, double s){ x.grad = s; }
inline double get_deriv (const adouble& y){ return y.grad; }   // seeded first-order derivative
#  endif
#endif

template<typename T> inline T smax0(const T& x){ using std::sqrt; return T(0.5)*(x + sqrt(x*x + T(1.0e-8))); }

template<typename T> inline void condassign(T& r, const T& cond, const T& a, const T& b){ r = (value(cond) > 0.0) ? a : b; }
#if PSOPT_AD_TAPED   // CppAD: a genuine taped conditional (selects a branch by the recorded value).
inline void condassign(adouble& r, const adouble& cond, const adouble& a, const adouble& b){
    r = CppAD::CondExpGt(cond, adouble(0.0), a, b);
}
#endif
template<typename T> inline T fmax(const T& a, const T& b){ T r, d=a-b; condassign(r,d,a,b); return r; }
template<typename T> inline T fmin(const T& a, const T& b){ T r, d=b-a; condassign(r,d,a,b); return r; }
} // namespace psopt_ad
#endif
