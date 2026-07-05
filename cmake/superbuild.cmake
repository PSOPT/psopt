# superbuild.cmake — download and build every dependency from source, then build
# PSOPT against them. Activated by -DPSOPT_SUPERBUILD=ON. Everything installs into
# ${CMAKE_BINARY_DIR}/_deps/install so nothing needs to be preinstalled.
#
# ABI note: C/C++ deps are built with the same compilers as PSOPT (libc++ on
# macOS) so ADOL-C/CasADi/PSOPT match; Fortran (MUMPS/HSL) uses gfortran.
include(ExternalProject)

set(SB_PREFIX   ${CMAKE_BINARY_DIR}/_deps)
set(SB_INSTALL  ${SB_PREFIX}/install)
file(MAKE_DIRECTORY ${SB_INSTALL})

set(SB_CC  ${CMAKE_C_COMPILER})
set(SB_CXX ${CMAKE_CXX_COMPILER})
find_program(SB_FC NAMES gfortran-mp-15 gfortran-15 gfortran REQUIRED)
set(SB_PKGCFG "${SB_INSTALL}/lib/pkgconfig")
# OpenBLAS for MUMPS/IPOPT LAPACK; prefer a found one, else let ThirdParty pick.
find_library(SB_OPENBLAS NAMES openblas HINTS /opt/local/lib /usr/local/lib /usr/lib)
if(SB_OPENBLAS)
  get_filename_component(_obdir ${SB_OPENBLAS} DIRECTORY)
  set(SB_BLAS_LFLAGS "-L${_obdir} -lopenblas")
else()
  set(SB_BLAS_LFLAGS "")   # ThirdParty scripts fall back to a reference BLAS
endif()

message(STATUS "PSOPT superbuild -> ${SB_INSTALL}")
message(STATUS "  CC=${SB_CC} CXX=${SB_CXX} FC=${SB_FC} BLAS='${SB_BLAS_LFLAGS}'")

set(_sb_common
  CMAKE_ARGS -DCMAKE_INSTALL_PREFIX=${SB_INSTALL} -DCMAKE_PREFIX_PATH=${SB_INSTALL}
             -DCMAKE_C_COMPILER=${SB_CC} -DCMAKE_CXX_COMPILER=${SB_CXX}
             -DBUILD_TESTING=OFF)

# ---- Eigen (CMake, header-only install) -------------------------------------
ExternalProject_Add(ep_eigen
  GIT_REPOSITORY https://gitlab.com/libeigen/eigen.git GIT_TAG 3.4.0 GIT_SHALLOW TRUE
  ${_sb_common} CMAKE_ARGS -DCMAKE_INSTALL_PREFIX=${SB_INSTALL} -DEIGEN_BUILD_DOC=OFF -DBUILD_TESTING=OFF)

# ---- ColPack + ADOL-C (autotools) -------------------------------------------
ExternalProject_Add(ep_colpack
  GIT_REPOSITORY https://github.com/CSCsw/ColPack.git GIT_TAG master GIT_SHALLOW TRUE
  CONFIGURE_COMMAND cd <SOURCE_DIR> && ./autoconf.sh && <SOURCE_DIR>/configure --prefix=${SB_INSTALL} CC=${SB_CC} CXX=${SB_CXX}
  BUILD_COMMAND cd <SOURCE_DIR> && make -j
  INSTALL_COMMAND cd <SOURCE_DIR> && make install
  BUILD_IN_SOURCE 1)
ExternalProject_Add(ep_adolc
  DEPENDS ep_colpack
  GIT_REPOSITORY https://github.com/coin-or/ADOL-C.git GIT_TAG master GIT_SHALLOW TRUE
  CONFIGURE_COMMAND cd <SOURCE_DIR> && autoreconf -fi && <SOURCE_DIR>/configure --prefix=${SB_INSTALL} --with-colpack=${SB_INSTALL} --enable-sparse CC=${SB_CC} CXX=${SB_CXX}
  BUILD_COMMAND cd <SOURCE_DIR> && make -j
  INSTALL_COMMAND cd <SOURCE_DIR> && make install
  BUILD_IN_SOURCE 1)

# ---- MUMPS (sequential, via COIN-OR ThirdParty) -----------------------------
ExternalProject_Add(ep_mumps
  GIT_REPOSITORY https://github.com/coin-or-tools/ThirdParty-Mumps.git GIT_TAG stable/3.0 GIT_SHALLOW TRUE
  CONFIGURE_COMMAND cd <SOURCE_DIR> && ./get.Mumps && <SOURCE_DIR>/configure --prefix=${SB_INSTALL} CC=${SB_CC} FC=${SB_FC} --with-lapack-lflags=${SB_BLAS_LFLAGS} --with-blas-lflags=${SB_BLAS_LFLAGS}
  BUILD_COMMAND cd <SOURCE_DIR> && make -j
  INSTALL_COMMAND cd <SOURCE_DIR> && make install
  BUILD_IN_SOURCE 1)

# ---- HSL (opt-in, license-gated; needs user-provided Coin-HSL source) -------
set(_ipopt_hsl_dep "")
set(_ipopt_hsl_flag "--without-hsl")
if(PSOPT_WITH_HSL)
  if(NOT COINHSL_SOURCE_DIR OR NOT EXISTS "${COINHSL_SOURCE_DIR}")
    message(FATAL_ERROR
      "PSOPT_WITH_HSL=ON requires -DCOINHSL_SOURCE_DIR=<path to coinhsl-x.y.z>.\n"
      "Register for a (free academic) licence and download Coin-HSL at\n"
      "  https://licences.stfc.ac.uk/product/coin-hsl")
  endif()
  ExternalProject_Add(ep_hsl
    GIT_REPOSITORY https://github.com/coin-or-tools/ThirdParty-HSL.git GIT_TAG stable/2.2 GIT_SHALLOW TRUE
    CONFIGURE_COMMAND ${CMAKE_COMMAND} -E copy_directory ${COINHSL_SOURCE_DIR} <SOURCE_DIR>/coinhsl
            COMMAND cd <SOURCE_DIR> && <SOURCE_DIR>/configure --prefix=${SB_INSTALL} CC=${SB_CC} FC=${SB_FC} --with-lapack-lflags=${SB_BLAS_LFLAGS}
    BUILD_COMMAND cd <SOURCE_DIR> && make -j
    INSTALL_COMMAND cd <SOURCE_DIR> && make install
    BUILD_IN_SOURCE 1)
  set(_ipopt_hsl_dep ep_hsl)
  set(_ipopt_hsl_flag "--with-hsl")
endif()

# ---- IPOPT (autotools; uses coinmumps + optional coinhsl via pkg-config) ----
ExternalProject_Add(ep_ipopt
  DEPENDS ep_mumps ${_ipopt_hsl_dep}
  GIT_REPOSITORY https://github.com/coin-or/Ipopt.git GIT_TAG stable/3.14 GIT_SHALLOW TRUE
  CONFIGURE_COMMAND ${CMAKE_COMMAND} -E env PKG_CONFIG_PATH=${SB_PKGCFG}
      <SOURCE_DIR>/configure --prefix=${SB_INSTALL}
      CC=${SB_CC} CXX=${SB_CXX} FC=${SB_FC}
      --with-mumps ${_ipopt_hsl_flag}
      --with-lapack-lflags=${SB_BLAS_LFLAGS} --with-blas-lflags=${SB_BLAS_LFLAGS}
      --disable-java --without-asl
  BUILD_COMMAND make -j
  INSTALL_COMMAND make install)

# ---- CasADi (CMake; with the IPOPT plugin, optional HSL) --------------------
set(_casadi_dep "")
if(PSOPT_WITH_CASADI)
  ExternalProject_Add(ep_casadi
    DEPENDS ep_ipopt
    GIT_REPOSITORY https://github.com/casadi/casadi.git GIT_TAG main GIT_SHALLOW TRUE
    ${_sb_common}
    CMAKE_ARGS -DCMAKE_INSTALL_PREFIX=${SB_INSTALL} -DCMAKE_PREFIX_PATH=${SB_INSTALL}
               -DWITH_IPOPT=ON -DWITH_HSL=${PSOPT_WITH_HSL} -DWITH_PYTHON=OFF -DBUILD_TESTING=OFF)
  set(_casadi_dep ep_casadi)
endif()

# ---- PSOPT itself (re-invoke this CMake with superbuild off, deps resolved) --
ExternalProject_Add(ep_psopt
  DEPENDS ep_eigen ep_adolc ep_ipopt ${_casadi_dep}
  SOURCE_DIR ${CMAKE_SOURCE_DIR}
  CMAKE_ARGS
    -DPSOPT_SUPERBUILD=OFF
    -DCMAKE_INSTALL_PREFIX=${CMAKE_INSTALL_PREFIX}
    -DCMAKE_PREFIX_PATH=${SB_INSTALL}
    -DIPOPT_ROOT=${SB_INSTALL}
    -DCMAKE_C_COMPILER=${SB_CC} -DCMAKE_CXX_COMPILER=${SB_CXX}
    -DPSOPT_WITH_MUMPS=${PSOPT_WITH_MUMPS} -DPSOPT_WITH_HSL=${PSOPT_WITH_HSL}
    -DPSOPT_WITH_CASADI=${PSOPT_WITH_CASADI} -DPSOPT_WITH_SNOPT=${PSOPT_WITH_SNOPT}
    -DPSOPT_DEFAULT_LINEAR_SOLVER=${PSOPT_DEFAULT_LINEAR_SOLVER}
    -DBUILD_EXAMPLES=${BUILD_EXAMPLES} -DBUILD_TESTS=${BUILD_TESTS} -DHEADLESS=${HEADLESS}
  INSTALL_COMMAND "")

message(STATUS "Superbuild configured. Run the build to download+compile all dependencies and PSOPT.")
