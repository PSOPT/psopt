
PSOPT
=====

Copyright (C) 2009-2025 Victor M. Becerra


Introduction
------------

This is the PSOPT library, a software tool for computational [optimal control](http://www.scholarpedia.org/article/Optimal_control)

PSOPT is an open source optimal control package written in C++ that uses [direct collocation methods](https://epubs.siam.org/doi/pdf/10.1137/16M1062569). These methods solve optimal control problems by approximating the time-dependent variables using global or local polynomials. This allows to discretize the differential equations and continuous constraints over a grid of nodes, and to compute any integrals associated with the problem using well known quadrature formulas. [Nonlinear programming](https://en.wikipedia.org/wiki/Nonlinear_programming) then is used to find local optimal solutions. PSOPT is able to deal with problems with the following characteristics:

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

- Choice between Legendre, Chebyshev, trapezoidal, or Hermite-Simpson based collocation
- Automatic scaling
- Automatic first and second derivatives using the ADOL-C library
- Numerical differentiation by using sparse finite differences
- Automatic mesh refinement
- Automatic identification of the Jacobian and Hessian sparsity.
- DAE formulation, so that differential and algebraic constraints can be implemented in the same C++ function.

The PSOPT interface uses both Eigen3 (a linear algebra template library) and ADOL-C (an automatic differentiation library).

The first release of PSOPT was published in 2009. 

The PSOPT website is [http://www.psopt.net](http://www.psopt.net).


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

Address:   
            
           University of Portsmouth

           School of Electrical and Mechanical Engineering
           
           Portsmouth PO1 3DJ
           
           United Kingdom

e-mail:    v.m.becerra@ieee.org


Rolling Release
---------------

From March 2025 PSOPT features a rolling release mode. Rolling release  is a concept in software development of frequently delivering updates to applications. This is in contrast to a standard or point release development model which uses software versions which replace the previous version. Users can download the latest source code from the GitHub repository. The documentation will also be updated  on a rolling release basis.

PSOPT documentation
-------------------

Please consult the [PSOPT User Manual (in PDF format)](https://github.com/PSOPT/psopt/blob/master/doc/PSOPT_Manual_RR.pdf) for further details on the software functionality and how to use it. 

There is also a [PSOPT Application Examples Document (in PDF format)](https://github.com/PSOPT/psopt/blob/master/doc/PSOPT_Application_Examples_Document_RR.pdf), which contains several application examples in various engineering/scientific domains, including their C++ code and results. 

Installation instructions are given below.


Installing PSOPT
----------------

Please consult the [PSOPT User Manual](https://github.com/PSOPT/psopt/blob/master/doc/PSOPT_Manual_RR.pdf) for further details on the software functionality and how to use it. 

PSOPT relies on three main software packages to perform a number of tasks: IPOPT, ADOL-C and EIGEN3. Some of these packages have their own dependencies.



**IPOPT**

IPOPT is an open-source C++ package for large-scale nonlinear optimization, which uses an interior point method. It is the default nonlinear programming algorithm used by PSOPT. IPOPT can be easily installed using a package manager in some, but not all, Linux distributions.

​	•	IPOPT repository:

​	https://github.com/coin-or/Ipopt

​	•	Version 3.12.12 is tested, but other versions may work.

​	https://www.coin-or.org/download/source/Ipopt/

​	•	Installation guide:

​	https://coin-or.github.io/Ipopt/INSTALL.html



**ADOL-C**

ADOL-C is a library for automatic differentiation of C++ code. It computes gradients and sparse Jacobians required by **PSOPT**.

A suitable version of ADOL-C can be easily installed using a package manager in some, but not all, Linux distributions. In some platforms, it may be necessary to manually install Adol-c and ColPack. The following commands should allow to perform a manual installation on various platforms:

```
$ wget --continue archive.ubuntu.com/ubuntu/pool/universe/a/adolc/adolc_2.7.2.orig.tar.xz 
$ tar -xf adolc_2.7.2.orig.tar.xz
$ cd ADOL-C-2.7.2
$ mkdir ./ThirdParty
$ cd ./ThirdParty
$ wget --continue http://archive.ubuntu.com/ubuntu/pool/universe/c/colpack/colpack_1.0.10.orig.tar.gz
$ tar zxvf colpack_1.0.10.orig.tar.gz
$ mv ColPack-1.0.10 ColPack
$ cd ColPack
$ autoreconf -fi
$ ./configure --prefix=/usr/local
$ make
$ sudo make install
$ cd ../..
$ autoreconf -fi
$ ./configure --prefix=/usr/local --enable-sparse --with-colpack=/usr/local
$ make
$ sudo make install
```

**EIGEN3**

[Eigen](http://eigen.tuxfamily.org/) is a lightweight, powerful linear algebra package for C++. Eigen is available on most major Linux distributions.



If necessary, Eigen can also be installed using CMake:

```
$ wget --continue https://gitlab.com/libeigen/eigen/-/archive/3.3.7/eigen-3.3.7.tar.gz
$ tar zxvf eigen-3.3.7.tar.gz
$ cd eigen-3.3.7
$ mkdir build
$ cd build
$ cmake ..
$ sudo make install
```



The following optional libraries can be employed for additional functionality.

**SNOPT**



[SNOPT](http://www.sbsi-sol-optimize.com/manuals/SNOPT-Manual.pdf) is an optimization algorithm for large-scale nonlinearly constrained problems based on sequential quadratic programming.



**GNUplot**



[GNUplot](http://www.gnuplot.info) is a portable, interactive data and function plotting utility. GNU plot is available on most Linux distributions. PSOPT includes a number of functions that allow to easily plot results using GNUplot.



**Building PSOPT**



PSOPT relies on [CMake](https://cmake.org/download/)  and '[pkg-config](https://en.wikipedia.org/wiki/Pkg-config)' for configuring builds, and on 'make' for managing compilation and linking.

CMake is an open-source tool for managing software builds. PSOPT requires CMake 3.12 or later. 

pkg-config is a helper tool used to provide the necessary details for compiling and linking a program to a library. It ensures that PSOPT’s dependencies are found correctly. pkg-config is available on most major Linux distributions. In particular, the build process expects to see pkg-config configuration files for IPOPT, ColPack and ADOL-C. These configuration files are usually installed under /usr/local/lib/pkgconfig or /usr/lib/pkgconfig. If these configuration files are not created during the build process for the above libraries, they can be created manually and be placed at the correct folder. If the pkg-config configuration files are being created manually, the contents of these files on the authors' computer are provided below as examples. Please note that the paths that are given in these files depend on the actual location where the different libraries have been installed.

For IPOPT (filename: ipopt.pc):

	prefix=/usr/local
	exec_prefix=${prefix}
	libdir=${exec_prefix}/lib
	includedir=${prefix}/include/coin-or
	Name: IPOPT
	Description: Interior Point Optimizer
	URL: https://github.com/coin-or/Ipopt
	Version: 3.13.2
	Cflags: -I${includedir}
	Libs: -L${libdir} -lipopt
	Requires.private: coinhsl coinmumps 



For ColPack (filename: ColPack.pc):

```
prefix=/usr/local
exec_prefix=${prefix}
libdir=${exec_prefix}/lib
includedir=${prefix}/include/ColPack

Name: ColPack
Version: 1.0.10 
Description: Graph Coloring Library 
Requires: 
Libs: -L${libdir} -lColPack -Wl,-rpath,${libdir} -Wl,-rpath,${libdir} 
Cflags: -I${includedir}
```

For ADOL-C (filename: adolc.pc):



```
prefix=/usr/local
exec_prefix=${prefix}
libdir=${exec_prefix}/lib
includedir=${prefix}/include

Name: adolc
Version: 2.6.3
Description: Algorithmic Differentiation Library for C/C++
Requires: 
Libs: -L${libdir} -ladolc -Wl,-rpath,${libdir} -lColPack -Wl,-rpath,${libdir} 
Cflags: -I${includedir}
\end{verbatim}
```

​	

For EIGEN3 (filename: eigen3.pc):

```
prefix=/usr/local
exec_prefix=${prefix}
libdir=${exec_prefix}/lib
includedir=${prefix}/include

Name: Eigen3
Version: 3.3.77
Description: Numerical linear algebra library for C++
Requires: 
Libs:  -Wl,-rpath,${libdir} -L$${libdir}  
Cflags: -I${includedir} -std=c++11
```

For SNOPT (filename: snopt7.pc):
	

```
prefix=/usr/local
exec_prefix=${prefix}
libdir=${exec_prefix}/lib
includedir=${prefix}/include/snopt7

Name: SNOPT7
Version: 7
Description: SNOPT NONLINEAR PROGRAMMING LIBRARY 
Requires: 
Libs: -L${libdir} -lsnopt7_cpp -Wl,-rpath,${libdir} -Wl,-rpath,${libdir} 
Cflags: -I${includedir}
```



**Tested Platforms**

PSOPT has been successfully compiled on:

​	•	Ubuntu Linux 24.04 LTS

​	•	OpenSUSE Linux 15.5 Leap and Tumbleweed

​	•	Arch Linux (latest versions as of 2025)

​	•	Manjaro Linux  (latest versions as of 2025)

​	•	MacOS Sequoia (MacPorts on Intel CPU)



**Installing Dependencies**



For **Ubuntu 24.04**:

```
$ sudo apt-get install git cmake gfortran g++ libboost-dev libboost-system-dev \
  coinor-libipopt-dev gnuplot libeigen3-dev libblas-dev liblapack-dev
```



Note that Adol-c and ColPack needs to be manually installed on the latest version of Ubuntu (version 24.04). See the instructions above.



For **Debian 12.9.0**:

```
$ su
$ apt-get install git cmake gfortran g++ libboost-dev libboost-system-dev \
  coinor-libipopt-dev gnuplot libeigen3-dev libblas-dev liblapack-dev
```

Note that Adol-c and ColPack needs to be manually installed on the latest version of Debian (version 11.9). See the instructions above.



For **OpenSUSE Leap 15.5 and Tumbleweed**:

```
$ sudo zypper install git gnuplot libboost_system1_66_0-devel eigen3-devel ColPack-devel \
  adolc-devel blas-devel lapack-devel Ipopt-devel cmake gcc-c++
```

For **Arch Linux / Manjaro**:

```
$ sudo pacman -Syu
$ sudo pacman -S git base-devel cmake gnuplot eigen boost blas lapack yay
$ yay -S coin-or-ipopt colpack adol-c
```



The use of the tool **yay** requires AUR support to be enabled on the package manager. On ARM64, it may be necessary to install [Anaconda]([https://www.anaconda.com/download), which provides gklib and IPOPT, as the installation script for IPOPT provided by AUR currently fails to build using yay.



For **MacOS**

It is also possible to build PSOPT on some versions of MacOS. This can be done by using MacPorts, which is an is an open-source system for compiling, installing, and upgrading open-source software on MacOS. 

To use MacPorts, download and install MacPorts from:

https://www.macports.org/install.php

To install the PSOPT dependencies, issue the following commands from a terminal window:

```
sudo port install cmake
sudo port install eigen3
sudo port install git
sudo port install ipopt
sudo port install ADOL-C
sudo port install gnuplot
```

The above is working with MacOS Sequoia with an Intel processor (date: 16 Feb 2025). Currently, there are issues with the ARM64 (M1 and above) processors that have MacOS Sequoia. The installation might work with previous versions of MacOS on ARM64. 



**Building and Installing PSOPT**

Once all dependencies are installed, PSOPT can be downloaded from GitHub, and built using CMake using the following commands.

```
$ git clone https://github.com/PSOPT/psopt.git
$ cd psopt
$ mkdir build
$ cd build
$ cmake -DBUILD_EXAMPLES=ON ..
$ make
$ sudo make install
```

If using SNOPT:

```
$ cmake -DBUILD_EXAMPLES=ON -DWITH_SNOPT_INTERFACE=ON ..
```

For debugging:

```
$ cmake -DBUILD_EXAMPLES=ON -DCMAKE_BUILD_TYPE=Debug ..
```

After installation, run at least one example to check that the build is working correctly:

```
$ cd build/examples/launch
$ ./launch
```


Running PSOPT witin a docker container
----------------

Docker containers are relatively small, standalone, executable software packages that include everything needed to run an application, such as code, runtime, libraries, and system tools. Containers are a form of operating system virtualisation. To use dockers containers, you need to install suitable software.  For instance, you can install Docker Desktop for [Windows 11](https://docs.docker.com/desktop/setup/install/windows-install/), [MacOS](https://docs.docker.com/desktop/setup/install/mac-install/), and various distributions of [Linux](https://docs.docker.com/desktop/setup/install/linux/).

The current distribution of PSOPT provides a Docker container file (Dockerfile). This provides an alternative way of installing and running PSOPT. 

The following are opportunities provided by the use of docker containers with PSOPT.

**Reproducible Environments:** A Docker container ensures everyone—whether on Linux, macOS, or Windows—runs PSOPT with the same OS libraries, compiler, and dependencies, eliminating configuration mismatches.

**Easier Setup:** Users avoid manually installing IPOPT, ADOL-C, COLPACK, EIGEN3, and other dependencies. A single docker pull or docker build command spins up a ready-to-run PSOPT environment.

**Continuous Integration (CI) Testing:** Automated pipelines (GitHub Actions, GitLab CI, etc.) can pull and test PSOPT in a Docker image, allowing fast and consistent builds without manually setting up each build agent. 

**Cloud or HPC Deployment:** Clusters often support container-based workloads. Docker images simplify running large-scale optimal control problems in cloud services or high-performance computing environments.

As it is not easy to get a docker to display graphical output (such as GNUplot plots), it is best to run PSOPT in headless mode (no graphical output) within the docker container, and visualise any graphical output from the host operating system (e.g. by opening any PDF files that PSOPT may have produced).

The steps to create a docker container and run PSOPT on the container are as follows:

1. Download [Dockerfile](https://github.com/PSOPT/psopt/blob/master/Dockerfile) from the PSOPT distribution, and place it in a folder. This Dockerfile uses [archlinux](https://hub.docker.com/_/archlinux/) as the base. This file clones the latest source code for PSOPT available from GitHub. If you have created your own version (for instance, to include your own examples or cases), you can modify the Dockerfile to copy your own source tree.

2. In your terminal, cd to the same folder where the Dockerfile is. The command to build the docker container (including PSOPT) is as follows: 
```
$ docker build --no-cache -t psopt-archlinux:latest .
```

3. Issue the following command to run the docker container interactively:
```
$ docker run -it psopt-archlinux:latest 
```
This will land you in the main 'psopt' folder. From there cd to 'build/examples' to run particular examples, etc.

4. Alternatively, you can use the following command to run the docker container interactively with a data connection to the host 

```
$ docker run -it --rm -v "$HOME/data:/data" psopt-archlinux:latest 
```
Here, the shared folder is "$HOME/data" as seen from the host, and "/data" as seen from the container.

From within the container, cd to 'build/examples' to run particular examples, etc.
Any output files must be manually copied to the folder /data from within the container. The copied files (e.g. PDFs or .txt files) appear within the corresponding directory of the host ($HOME/data). The host can send files to the container via the same folder.

Getting help
------------

* **[PSOPT Documentation](https://github.com/PSOPT/psopt/blob/master/doc/)** with information about the functionality and use of the software, background theory, examples, and more.
 * **[Issue tracking system](https://github.com/PSOPT/psopt/issues/)**: If you believe you found a **bug** in the code, please use the issue tracking system.
   Please include as much information as possible, and if possible some example code so that we can reproduce the error.
 * **[Mailing list](http://groups.google.com/group/psopt-users-group)**: subscribe to receive notifications about updates and to post questions and comments about PSOPT.


Please acknowledge this work
----------------------------

This software is provided for free in the hope that it may be useful to others, and we would very much like to hear about your experience with it. If you find PSOPT helpful for your work or research, please email the author at v.m.becerra@ieee.org  to incorporate a feature on the PSOPT web page.

Given that a great deal of time and effort has gone into PSOPT's development, **please cite the following publication if you are using PSOPT for your own research**:

* Becerra, V.M. (2010). [**Solving complex optimal control problems at no cost with PSOPT**](https://ieeexplore.ieee.org/document/5612676). Proc. IEEE Multi-conference on Systems and Control, Yokohama, Japan, September 7-10, 2010, pp. 1391-1396.

**BibTex entry:**

            @INPROCEEDINGS{5612676,  
            author={V. M. Becerra},  
            booktitle={2010 IEEE International Symposium on Computer-Aided Control System Design},          
            title={Solving complex optimal control problems at no cost with PSOPT},   
            year={2010},    
            pages={1391-1396},  
            doi={10.1109/CACSD.2010.5612676}}


