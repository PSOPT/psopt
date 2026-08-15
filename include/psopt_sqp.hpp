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

#ifndef PSOPT_SQP_HPP
#define PSOPT_SQP_HPP

//  psopt_sqp.hpp -- a sequential quadratic programming solver for PSOPT's NLP.
//
//  This is the dense stage of a two-stage plan. The algorithm here is the classical
//  one: a quadratic model of the Lagrangian built by damped BFGS updates, a quadratic
//  programming subproblem solved by qpOASES, and an l1 merit function with a
//  backtracking line search to globalise it. It descends from a dense SQP written by
//  the author in 2008 and is restated here against PSOPT's own function, gradient and
//  Jacobian machinery rather than the ADOL-C layer of the original.
//
//  The model of the problem is sparse throughout. The constraint Jacobian is kept in
//  compressed columns, on a pattern taken once from the recorded tape, and with
//  algorithm.hessian = "exact" the curvature is the exact Hessian of the Lagrangian,
//  also sparse, evaluated at every iterate through the same automatic differentiation
//  the IPOPT interface uses. With algorithm.hessian = "limited-memory", which remains
//  the default, the curvature is a dense n-by-n quasi-Newton model built by damped
//  BFGS updates: that model is what is dense, at O(n^2) storage, and it is the reason
//  the exact Hessian is worth having on a large mesh quite apart from its faster
//  convergence.
//
//  Two things are still dense, and are what a further stage would change:
//
//    * the quasi-Newton model, as above, when it is the one in use;
//    * qpOASES's own factorisations, when it is the QP solver in use. The
//      subproblem's matrices are handed over in sparse form, but the null-space
//      method inside qpOASES holds a dense n-by-n orthogonal factor whatever it is
//      given, so its memory is quadratic and its work per subproblem cubic in n.
//      algorithm.qp_solver replaces it with one of two proximal augmented-Lagrangian
//      methods, each of which factorises the KKT system sparsely and tolerates an
//      indefinite Hessian: "ProxQP", from INRIA's ProxSuite, header-only C++ over
//      Eigen under BSD-2-Clause (-DWITH_PROXQP=ON), or "QPALM", from KU Leuven, C
//      over its own LDL factorisation, under LGPL-3 (-DWITH_QPALM=ON). The licences
//      differ in consequence: linking QPALM makes the result distributable under
//      LGPL-3 rather than under PSOPT's LGPL-2.1.
//
//  Simple bounds are passed to the QP as bounds. The 2008 code expanded them into
//  2n general inequality rows with an identity block, which is harmless on a
//  three-variable problem and ruinous on a collocation mesh.
//
//  The globalisation depends on which model is in use, and has to. A line search on
//  the l1 merit function is enough for the quasi-Newton model, which is built to be
//  positive definite; the exact Hessian is indefinite wherever the problem is, and on
//  a minimum-time problem at a first iterate, where the multipliers are still zero, it
//  is exactly zero. It is therefore paired with a trust region and with a shift of the
//  diagonal that starts at one and halves at every iteration that does not need it.
//
//  Reference for the algorithm: Nocedal and Wright, Numerical Optimization, 2nd ed.,
//  chapter 18; the damped BFGS update is Powell (1978).

// This header is a companion to psopt.h and must be included after it: it uses the
// Eigen and PSOPT names that psopt.h brings into scope.

// Solve PSOPT's NLP by SQP. The signature matches NLP_interface's own so that the
// dispatch in NLP_interface.cxx is a single branch.
//
//   min f(x)  subject to  g_l <= g(x) <= g_u,  xlb <= x <= xub
//
// with the constraint bounds obtained from get_constraint_bounds(), so an equality is
// a constraint whose bounds coincide -- the same convention the IPOPT interface uses.
// On return x0 holds the solution and lambda the constraint multipliers in the sign
// convention of the IPOPT interface, that is, grad f + J' lambda = 0.
//
// Returns 0 on success and non-zero otherwise; the reason is placed in
// workspace->solution->status where the caller can report it.
int SQP_interface(Alg&         algorithm,
                  MatrixXd* x0,
                  double     (*f)(MatrixXd&, Workspace*),
                  void       (*g)(MatrixXd&, MatrixXd*, Workspace*),
                  int          nlp_ncons,
                  int          nlp_neq,
                  MatrixXd* xlb,
                  MatrixXd* xub,
                  MatrixXd* lambda,
                  int          hotflag,
                  int          iprint,
                  Workspace*   workspace,
                  void*        user_data);

#endif // PSOPT_SQP_HPP
