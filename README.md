
PSOPT
=====

Copyright (C) 2009-2020 Victor M. Becerra


Introduction
------------

This is the PSOPT library, a software tool for computational [optimal control](http://www.scholarpedia.org/article/Optimal_control)

PSOPT is an open source optimal control package written in C++ that uses direct collocation methods. These methods solve optimal control problems by approximating the time-dependent variables using global or local polynomials. This allows to discretize the differential equations and continuous constraints over a grid of nodes, and to compute any integrals associated with the problem using well known quadrature formulas. [Nonlinear programming](https://en.wikipedia.org/wiki/Nonlinear_programming) then is used to find local optimal solutions. PSOPT is able to deal with problems with the following characteristics:

-  Single or multiphase problems
-  Continuous time nonlinear dynamics
-  General endpoint constraints
-  Nonlinear path constraints (equalities or inequalities) on states and/or control variables
- Integral constraints
-  Interior point constraints
-  Bounds on controls and state variables
-  General cost function with Lagrange and Mayer terms.
-  Free or fixed initial and final conditions
- Linear or nonlinear linkages between phases
-  Fixed or free initial time
-  Fixed or free final time
- Optimisation of static parameters
- Parameter estimation problems with sampled measurements • Differential equations with delayed variables.

The implementation has the following features:

- Automatic scaling
- Automatic first and second derivatives using the ADOL-C library
- Numerical differentiation by using sparse finite differences
- Automatic mesh refinement
- Automatic identification of the Jacobian and Hessian sparsity.
- DAE formulation, so that differential and algebraic constraints can be implemented in the same C++ function.

The PSOPT interface uses both Eigen3 (a linear algebra template library) and ADOL-C (an automatic differentiation library).

The PSOPT project webpage is <http://psopt.org>. The first release of PSOPT was published in 2009. This is release 5 of PSOPT.


License
----------


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

           School of Energy and Electronic Engineering
           
           Portsmouth PO1 3DJ
           
           United Kingdom

e-mail:    v.m.becerra@ieee.org


Getting Started
---------------

Please consult the detailed installation instructions in the [PSOPT PDF documentation](https://github.com/PSOPT/psopt/blob/master/doc/PSOPT_Manual_R5.pdf). In the following, we only summarize some main points.

### Dependencies

PSOPT requires the following libraries:

1. [IPOPT](https://github.com/coin-or/Ipopt )
2. [ADOL-C](https://github.com/coin-or/ADOL-C)
3. [EIGEN3](http://eigen.tuxfamily.org/)

Optionally, the user may wish to employ the following software
1. [SNOPT](http://www.sbsi-sol-optimize.com/manuals/SNOPT-Manual.pdf)
2. [GNUplot](http://www.gnuplot.info)

Moreover, PSOPT is built using [CMake](https://cmake.org/download/). CMake is an open-source, cross-platform family of tools designed to build, test and package software. The build process also requires the [pkg-config](https://www.freedesktop.org/wiki/Software/pkg-config/) tool.  

### Build

After installation of dependencies, a typical PSOPT build and installation on a Unix-like operating system follows these steps (please see the PDF documentation for futher details):

1. Download and extract the installation archive from the GItHub project page. Alternatively, the source code can be cloned using git, by issuing the following command: 
         `git clone https://github.com/PSOPT/psopt.git`
2.  `mkdir build; cd build`
3. `cmake -DBUILD_EXAMPLES=ON .. $ make`
4.  `sudo make install`


Getting Help
------------

* **[PSOPT Documentation](https://github.com/PSOPT/psopt/blob/master/doc/PSOPT_Manual_R5.pdf)** with installation instructions, background theory, examples and much more
 * **[Issue tracking system](https://github.com/PSOPT/psopt/issues/)**: If you believe you found a **bug** in the code, please use the issue tracking system.
   Please include as much information as possible, and if possible some (ideally simple) example code so that we can reproduce the error.
 * **[Mailing list](http://groups.google.com/group/psopt-users-group)**: subscribe to receive notifications about updates and to post questions and comments about PSOPT.


Please Cite Us
--------------

This software is provided for free in the hope that it may be useful to others, and we would very much like to hear about your experience with it. If you find PSOPT helpful for your work or research, please emai the author at v.m.becerra@ieee.org  to incorporate a feature on the PSOPT web page.

Given that a great deal of time and effort has gone into PSOPT's development, **please cite the following publication if you are using PSOPT for your own research**:

* Becerra, V.M. (2010). **Solving complex optimal control problems at no cost with PSOPT**. Proc. IEEE Multi-conference on Systems and Control, Yokohama, Japan, September 7-10, 2010, pp. 1391-1396.

