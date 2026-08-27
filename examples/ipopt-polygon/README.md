# Largest small polygon — Example 14 of Chapter 2 of the book

This program does not use PSOPT. It is the largest small polygon problem worked
through in Chapter 2 of *Computational Optimal Control and Estimation*, solved by
calling IPOPT directly through its `TNLP` interface, with sparse first and second
derivatives supplied by CppAD. It is shipped with PSOPT so that a reader of the
book can build and run it; nothing in the library depends on it.

An earlier version used ADOL-C with graph colouring by ColPack. PSOPT now uses
CppAD throughout, so this example does too, and a reader needs one automatic
differentiation tool rather than two.

## Building and running

It is built with the rest of the examples, that is with `-DBUILD_EXAMPLES=ON`,
provided the AD backend is CppAD (the default). Then

    ./ipopt_polygon 16
    ./ipopt_polygon 50     # also the default

The program writes `ipopt.out` (IPOPT's own log, which is where the
statistics in the book's table for this example come from), `ipopt<nv>.m` (the solution), and `jacobian.dat` and
`Hessian.dat` (the sparsity patterns). Run both sizes and then `plotpolygon.m`
draws the figure of the two optimal polygons; `plotspy.m` draws the sparsity
figure.

## What to expect

Eleven iterations in both cases, to areas of 0.7718613 for nv = 16 and 0.7840773
for nv = 50. For nv = 50 the constraint Jacobian has 4,802 nonzeros out of
122,600 and the Hessian of the Lagrangian 4,857.

One option is not a detail here. The problem has many local maxima — a polygon
can be locally optimal without being the largest — and IPOPT's barrier-parameter
strategy decides which one is reached. With the `adaptive` strategy the program
sets, it converges in 11 iterations to the known optimum. With the `monotone`
default it converges instead, in 54 iterations, to a genuine but inferior
stationary point of area 0.7268687. That is worth knowing before reading too much
into a single run of any problem of this kind.
