# Cam design — Example 15 of Chapter 2 of the book

This program does not use PSOPT. It is the cam shape optimisation problem worked
through in Chapter 2 of *Computational Optimal Control and Estimation*, solved by
calling IPOPT directly through its `TNLP` interface, with sparse first and second
derivatives supplied by CppAD. It is shipped with PSOPT so that a reader of the
book can build and run it; nothing in the library depends on it.

An earlier version used ADOL-C with graph colouring by ColPack. PSOPT now uses
CppAD throughout, so this example does too, and a reader needs one automatic
differentiation tool rather than two. The sparsity patterns it produces are
identical, entry for entry, to those the ADOL-C version produced.

## Building and running

It is built with the rest of the examples, that is with `-DBUILD_EXAMPLES=ON`,
provided the AD backend is CppAD (the default). Then

    ./ipopt_cam            # n = 1000, the case reported in the book
    ./ipopt_cam 200        # any other size

The program writes `ipopt.out` (IPOPT's own log, which is where the
statistics in the book's table for this example come from), `ipopt.m` (the solution, read by `plot_cam.m`), and
`jacobian.dat` and `Hessian.dat` (the sparsity patterns, read by `plotspy.m`).
The two MATLAB scripts draw the two figures for this example in the book.

## What to expect for n = 1000

64 iterations to `f(x*) = -4.2790674`, with no constraint violation. The
constraint Jacobian has 5,000 nonzeros out of 2,003,000 and the Hessian of the
Lagrangian 1,998 out of 1,000,000. The evaluation counts and the timing depend a
little on the version of IPOPT.
