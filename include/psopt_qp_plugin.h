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

#ifndef PSOPT_QP_PLUGIN_H
#define PSOPT_QP_PLUGIN_H

/*  psopt_qp_plugin.h -- the interface between PSOPT's SQP solver and a QP backend.
 *
 *  Every QP library worth having comes with a sparse factorisation, and every sparse
 *  factorisation comes with an ordering, and the orderings are all the same handful of
 *  algorithms carrying the same handful of C symbol names. QPALM factorises through
 *  LADEL, which links SuiteSparse's AMD and exports amd_order and amd_l_order built for
 *  64-bit indices; ProxQP carries its own; GALAHAD brings SSIDS and another again;
 *  IPOPT has already brought whatever MUMPS was built against. Linked into one program
 *  they are one symbol with several definitions, the linker picks one, and a caller
 *  compiled against a different index width or a different control-array layout then
 *  writes past the end of something. The symptom is a segmentation fault a long way
 *  from the cause: with QPALM and ProxQP both linked, it was PSOPT's ProxQP test
 *  crashing inside ProxQP, which had not changed. qpOASES's BLAS replacements capturing
 *  dgemm_ for MUMPS was the same failure wearing a different coat.
 *
 *  Rather than keep discovering these one at a time, each backend is built as a plugin:
 *  a shared object that links its own QP library and its own linear algebra, compiles
 *  with hidden visibility so that none of it is exported, and offers the three
 *  functions below and nothing else. PSOPT loads plugins with RTLD_LOCAL, so what is
 *  inside one is invisible to the next, and two backends whose orderings would once
 *  have fought now do not know about each other.
 *
 *  The interface is C, and deliberately narrow: a subproblem in, a step and its
 *  multipliers out. It carries no C++ types, no Eigen, and no PSOPT types, so a plugin
 *  can be built against a different compiler or standard library than the one that
 *  built PSOPT, and so that adding a backend means writing one file.
 */

#ifdef __cplusplus
extern "C" {
#endif

/*  Bumped when the structures below change in a way an existing plugin would not
 *  survive. PSOPT refuses to load a plugin that reports a different version. */
#define PSOPT_QP_ABI_VERSION 1

/*  Status codes returned in psopt_qp_solution.status. An approximate solution is not a
 *  failure: the SQP's line search is the judge of whether a step is any good, and
 *  refusing a step because a QP hit its iteration limit sends the solver to a
 *  restoration it does not need. */
#define PSOPT_QP_SOLVED       0   /* solved to the requested tolerance          */
#define PSOPT_QP_APPROXIMATE  1   /* iteration limit reached; step is usable    */
#define PSOPT_QP_FAILED       2   /* no usable step                             */

/*  The subproblem
 *
 *      min  1/2 d' H d + g' d     s.t.  lbA <= A d <= ubA,   lb <= d <= ub
 *
 *  with H symmetric, given either sparsely, in compressed columns with both triangles
 *  stored, or densely; exactly one of the two is provided. A is m by n in compressed
 *  columns. An absent bound is +/- PSOPT_QP_INFINITY or beyond.
 *
 *  Index arrays are 64-bit signed, which is what every sparse library here uses at the
 *  sizes PSOPT reaches, and avoids the plugin having to guess.
 */
#define PSOPT_QP_INFINITY 1.0e20

typedef struct {
    int             abi_version;   /* PSOPT_QP_ABI_VERSION, checked by the plugin */

    int             n;             /* variables   */
    int             m;             /* general constraints; may be zero            */

    /* Hessian, sparse: column pointers (n+1), row indices, values. NULL if dense. */
    const long long* H_p;
    const long long* H_i;
    const double*    H_x;
    /* Hessian, dense, column-major n by n. NULL if sparse. */
    const double*    H_dense;

    const double*    g;            /* n */

    /* Constraint matrix, sparse, m by n. NULL when m == 0. */
    const long long* A_p;
    const long long* A_i;
    const double*    A_x;

    const double*    lbA;          /* m */
    const double*    ubA;          /* m */
    const double*    lb;           /* n */
    const double*    ub;           /* n */

    double           tolerance;    /* the NLP tolerance the SQP is working to */
    int              max_iter;     /* per-subproblem iteration budget         */
    int              nonconvex;    /* non-zero if H may be indefinite         */
} psopt_qp_problem;

/*  The answer, in PSOPT's sign convention: grad f + A' lambda - z = 0. Backends whose
 *  own convention differs -- qpOASES's duals carry the opposite sign, ProxQP's and
 *  QPALM's do not -- convert here, and each plugin's conversion is pinned by a test
 *  against a QP whose multiplier is known.
 *
 *  The arrays are allocated by the caller: d and z of length n, lambda of length m.
 */
typedef struct {
    double* d;
    double* lambda;
    double* z;
    int     iterations;
    int     status;
} psopt_qp_solution;

/*  The three functions a plugin exports, and the only symbols it may export. */

/*  The ABI the plugin was built against. */
int         psopt_qp_abi_version(void);

/*  A name for messages, e.g. "ProxQP". Owned by the plugin; not freed by the caller. */
const char* psopt_qp_name(void);

/*  Solve. Returns the same value placed in solution->status. Must not throw, must not
 *  abort, and must not call exit(): a subproblem a backend cannot take is an ordinary
 *  event and is reported, not escalated. */
int         psopt_qp_solve(const psopt_qp_problem* problem, psopt_qp_solution* solution);

#ifdef __cplusplus
}
#endif

#endif /* PSOPT_QP_PLUGIN_H */
