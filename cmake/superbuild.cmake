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
# Legacy C deps (METIS 4.0, older MUMPS/HSL C) predate C99 strictness; clang 21
# makes implicit declarations/int a hard error. Relax it for those builds.
set(SB_LEGACY_CFLAGS "-Wno-implicit-function-declaration -Wno-implicit-int -Wno-error=implicit-function-declaration")
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

# Apple clang needs explicit libomp flags for deps that require OpenMP (ColPack).
set(SB_OMP_ARGS "")
if(APPLE)
  find_path(SB_OMP_INC omp.h HINTS /opt/local/include/libomp /opt/local/include /usr/local/include /opt/homebrew/include)
  find_library(SB_OMP_LIB NAMES omp libomp HINTS /opt/local/lib/libomp /opt/local/lib /usr/local/lib /opt/homebrew/lib)
  if(SB_OMP_INC AND SB_OMP_LIB)
    set(SB_OMP_ARGS
      "-DOpenMP_C_FLAGS=-Xpreprocessor -fopenmp -I${SB_OMP_INC}"   "-DOpenMP_C_LIB_NAMES=omp"
      "-DOpenMP_CXX_FLAGS=-Xpreprocessor -fopenmp -I${SB_OMP_INC}" "-DOpenMP_CXX_LIB_NAMES=omp"
      "-DOpenMP_omp_LIBRARY=${SB_OMP_LIB}")
  endif()
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
  # Use ColPack's CMake build (build/cmake) — builds just the library, avoiding
  # the automake path's broken Examples targets.
  SOURCE_SUBDIR build/cmake
  CMAKE_ARGS -DCMAKE_INSTALL_PREFIX=${SB_INSTALL} -DCMAKE_PREFIX_PATH=${SB_INSTALL}
             -DCMAKE_C_COMPILER=${SB_CC} -DCMAKE_CXX_COMPILER=${SB_CXX}
             -DBUILD_SHARED_LIBS=ON -DBUILD_TESTING=OFF ${SB_OMP_ARGS})
ExternalProject_Add(ep_adolc
  DEPENDS ep_colpack
  GIT_REPOSITORY https://github.com/coin-or/ADOL-C.git GIT_TAG master GIT_SHALLOW TRUE
  CONFIGURE_COMMAND cd <SOURCE_DIR> && autoreconf -fi && <SOURCE_DIR>/configure --prefix=${SB_INSTALL} --with-colpack=${SB_INSTALL} --enable-sparse CC=${SB_CC} CXX=${SB_CXX}
  BUILD_COMMAND cd <SOURCE_DIR> && make -j
  INSTALL_COMMAND cd <SOURCE_DIR> && make install
  BUILD_IN_SOURCE 1)

# ---- METIS (fill-reducing ordering; speeds MUMPS and the HSL solvers) --------
# ThirdParty-Metis installs libcoinmetis + coinmetis.pc; MUMPS and HSL pick it
# up automatically via pkg-config (PKG_CONFIG_PATH below).
ExternalProject_Add(ep_metis
  GIT_REPOSITORY https://github.com/coin-or-tools/ThirdParty-Metis.git GIT_TAG stable/2.0 GIT_SHALLOW TRUE
  CONFIGURE_COMMAND cd <SOURCE_DIR> && ./get.Metis && <SOURCE_DIR>/configure --prefix=${SB_INSTALL} CC=${SB_CC} "CFLAGS=${SB_LEGACY_CFLAGS}"
  BUILD_COMMAND cd <SOURCE_DIR> && make -j
  INSTALL_COMMAND cd <SOURCE_DIR> && make install
  BUILD_IN_SOURCE 1)

# ---- MUMPS (sequential, via COIN-OR ThirdParty; uses METIS ordering) --------
ExternalProject_Add(ep_mumps
  DEPENDS ep_metis
  GIT_REPOSITORY https://github.com/coin-or-tools/ThirdParty-Mumps.git GIT_TAG stable/3.0 GIT_SHALLOW TRUE
  CONFIGURE_COMMAND cd <SOURCE_DIR> && ./get.Mumps
          COMMAND ${CMAKE_COMMAND} -E env PKG_CONFIG_PATH=${SB_PKGCFG} <SOURCE_DIR>/configure --prefix=${SB_INSTALL} CC=${SB_CC} FC=${SB_FC} --with-metis --with-lapack-lflags=${SB_BLAS_LFLAGS} --with-blas-lflags=${SB_BLAS_LFLAGS}
  BUILD_COMMAND cd <SOURCE_DIR> && make -j
  INSTALL_COMMAND cd <SOURCE_DIR> && make install
  BUILD_IN_SOURCE 1)

# ---- HSL (opt-in, NON-FREE licence; notify then auto-download+build) --------
# Coin-HSL (MA27/MA57/MA86/MA97) is NOT free software. We do not vendor it; the
# superbuild fetches it automatically ONLY after the user acknowledges the
# licence (PSOPT_HSL_ACCEPT_LICENSE=ON) or supplies a local tree
# (COINHSL_SOURCE_DIR). The download URL is configurable (PSOPT_HSL_URL) so a
# personalised/licensed link can be used.
set(_ipopt_hsl_dep "")
set(_ipopt_hsl_flag "--without-hsl")
set(_hsl_extra_dep "")
if(PSOPT_WITH_HSL)
  message(WARNING
    "\n============================ HSL LICENCE NOTICE ============================\n"
    "PSOPT_WITH_HSL=ON. The HSL / Coin-HSL linear solvers (MA27/MA57/MA86/MA97)\n"
    "are NON-FREE, licensed software from STFC. They are NOT redistributed by\n"
    "PSOPT. A free licence for academic/personal use is available at:\n"
    "    https://licences.stfc.ac.uk/product/coin-hsl\n"
    "By setting PSOPT_HSL_ACCEPT_LICENSE=ON you confirm you have obtained a\n"
    "licence and accept its terms; the source is then downloaded automatically.\n"
    "===========================================================================\n")

  if(COINHSL_SOURCE_DIR AND EXISTS "${COINHSL_SOURCE_DIR}")
    # Use a local Coin-HSL tree the user already downloaded.
    set(_hsl_populate ${CMAKE_COMMAND} -E copy_directory ${COINHSL_SOURCE_DIR} <SOURCE_DIR>/coinhsl)
    message(STATUS "HSL: using local source ${COINHSL_SOURCE_DIR}")
  elseif(PSOPT_HSL_ACCEPT_LICENSE)
    # Guard against the common mistake of pointing at the ThirdParty-HSL build
    # wrapper (which contains NO HSL source) instead of the Coin-HSL tarball.
    if(PSOPT_HSL_URL MATCHES "ThirdParty-HSL" OR PSOPT_HSL_URL MATCHES "\\.git$")
      message(FATAL_ERROR
        "PSOPT_HSL_URL='${PSOPT_HSL_URL}' points at the ThirdParty-HSL build\n"
        "wrapper, which does NOT contain HSL solver source (no MA27/57/86/97).\n"
        "That wrapper is already used internally. Set PSOPT_HSL_URL to a Coin-HSL\n"
        "SOURCE tarball (coinhsl-x.y.z.tar.gz), e.g. your download link from:\n"
        "  https://licences.stfc.ac.uk/product/coin-hsl            (full HSL), or\n"
        "  https://licences.stfc.ac.uk/product/coin-hsl-archive    (free MA27/MA28/MC19)")
    endif()
    # Auto-download the Coin-HSL archive, unpack it into ThirdParty-HSL/coinhsl.
    message(STATUS "HSL: licence accepted — will auto-download ${PSOPT_HSL_URL}")
    ExternalProject_Add(ep_coinhsl_src
      URL ${PSOPT_HSL_URL}
      DOWNLOAD_NO_EXTRACT FALSE
      CONFIGURE_COMMAND "" BUILD_COMMAND "" INSTALL_COMMAND ""
      SOURCE_DIR ${SB_PREFIX}/coinhsl-src)
    set(_hsl_populate ${CMAKE_COMMAND} -E copy_directory ${SB_PREFIX}/coinhsl-src <SOURCE_DIR>/coinhsl)
    set(_hsl_extra_dep ep_coinhsl_src)
  else()
    message(FATAL_ERROR
      "PSOPT_WITH_HSL=ON but the non-free HSL licence has not been acknowledged.\n"
      "Either:\n"
      "  * accept the licence and auto-download:  -DPSOPT_HSL_ACCEPT_LICENSE=ON "
      "(optionally -DPSOPT_HSL_URL=<your licensed link>), or\n"
      "  * point at a local Coin-HSL tree:        -DCOINHSL_SOURCE_DIR=<path>\n"
      "Obtain a (free academic) licence at https://licences.stfc.ac.uk/product/coin-hsl")
  endif()

  ExternalProject_Add(ep_hsl
    DEPENDS ${_hsl_extra_dep} ep_metis
    GIT_REPOSITORY https://github.com/coin-or-tools/ThirdParty-HSL.git GIT_TAG stable/2.2 GIT_SHALLOW TRUE
    CONFIGURE_COMMAND ${_hsl_populate}
            COMMAND ${CMAKE_COMMAND} -E env PKG_CONFIG_PATH=${SB_PKGCFG} <SOURCE_DIR>/configure --prefix=${SB_INSTALL} CC=${SB_CC} FC=${SB_FC} --with-metis --with-lapack-lflags=${SB_BLAS_LFLAGS}
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
# With HSL built, default to the OpenMP-parallel ma97 solver.
set(_psopt_default_solver ${PSOPT_DEFAULT_LINEAR_SOLVER})
if(PSOPT_WITH_HSL AND _psopt_default_solver STREQUAL "mumps")
  set(_psopt_default_solver "ma97")
endif()
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
    -DPSOPT_WITH_OPENMP=${PSOPT_WITH_OPENMP} -DPSOPT_WITH_MPI=${PSOPT_WITH_MPI}
    -DPSOPT_DEFAULT_LINEAR_SOLVER=${_psopt_default_solver}
    -DBUILD_EXAMPLES=${BUILD_EXAMPLES} -DBUILD_TESTS=${BUILD_TESTS} -DHEADLESS=${HEADLESS}
  INSTALL_COMMAND "")

message(STATUS "Superbuild configured. Run the build to download+compile all dependencies and PSOPT.")
