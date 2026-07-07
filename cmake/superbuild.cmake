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

# Legacy C deps (METIS 4.0, older MUMPS/HSL C) predate C99 strictness; modern
# gcc also makes implicit declarations a hard error. Relax it for those builds.
set(SB_LEGACY_CFLAGS "-Wno-implicit-function-declaration -Wno-implicit-int -Wno-error=implicit-function-declaration")
# Fortran is OPTIONAL. It is needed only for the Fortran solver stack (MUMPS/IPOPT,
# and the opt-in HSL/SPRAL/SCIP) and — transitively — for PSOPT's sparse ADOL-C
# derivatives (sparse_jac), which the IPOPT path uses. Without a Fortran compiler the
# superbuild builds the C++-only stack (Eigen + ADOL-C WITHOUT sparse + CasADi), which
# PSOPT drives via finite-difference derivatives and a CasADi C++ solver. ColPack (an
# ADOL-C sparse dependency) is therefore only required when Fortran is found; it is
# always git-downloadable, but is only built/linked on the Fortran path.
if(MSVC)
  # MSVC has no gcc/gfortran. Use the MSVC C/C++ compiler for the CMake-based deps
  # (Eigen/CasADi); ADOL-C 2.7.2 has no CMake, so it is built via its Visual Studio
  # solution (msbuild). This is a C++-only, Fortran-free superbuild — no ColPack/
  # METIS/MUMPS/IPOPT. EXPERIMENTAL/UNTESTED; expect to iterate on paths/flags.
  set(SB_CC  "${CMAKE_C_COMPILER}")
  set(SB_CXX "${CMAKE_CXX_COMPILER}")
  set(SB_HAVE_FORTRAN FALSE)
  set(SB_LEGACY_CFLAGS "")
  message(STATUS "MSVC superbuild — C++-only (Eigen + ADOL-C via msbuild + CasADi); "
                 "no Fortran/IPOPT/ColPack. EXPERIMENTAL/UNTESTED.")
else()
  # MacPorts/MinGW GNU compilers: gcc is far more tolerant of the old scientific C/C++
  # (ColPack/METIS/ADOL-C/MUMPS), has native OpenMP, and keeps the libstdc++ ABI
  # self-consistent (the superbuild compiles every dependency AND PSOPT itself).
  find_program(SB_CC  NAMES gcc-mp-15 gcc-mp-14 gcc-mp-13 gcc REQUIRED)
  find_program(SB_CXX NAMES g++-mp-15 g++-mp-14 g++-mp-13 g++ REQUIRED)
  find_program(SB_FC NAMES gfortran-mp-15 gfortran-15 gfortran)
  if(SB_FC)
    set(SB_HAVE_FORTRAN TRUE)
  else()
    set(SB_HAVE_FORTRAN FALSE)
    message(STATUS "No Fortran compiler found — building the C++-only superbuild stack "
                   "(no ColPack/METIS/MUMPS/IPOPT); PSOPT will use the CasADi backend.")
  endif()
endif()
set(SB_PKGCFG "${SB_INSTALL}/lib/pkgconfig")
# OpenBLAS for MUMPS/IPOPT LAPACK; prefer a found one, else let ThirdParty pick.
# For PSOPT_WITH_INT64 (ILP64) prefer a 64-bit-integer BLAS (openblas64 /
# libopenblas with INTERFACE64) and build MUMPS with 8-byte Fortran integers.
set(SB_INT64_MUMPS_FLAGS "")
if(PSOPT_WITH_INT64)
  find_library(SB_OPENBLAS NAMES openblas64 openblas_ilp64 openblas HINTS /mingw64/lib /opt/local/lib /usr/local/lib /usr/lib)
  # MUMPS with 64-bit ordinals: 8-byte default Fortran integer + INTSIZE64.
  set(SB_INT64_MUMPS_FLAGS "ADD_FCFLAGS=-fdefault-integer-8" "ADD_CFLAGS=-DINTSIZE64")
  message(STATUS "PSOPT_WITH_INT64=ON — building MUMPS/BLAS with 64-bit integers "
                 "(full ILP64 IPOPT also needs a 64-bit-Index IPOPT; see docs)")
else()
  find_library(SB_OPENBLAS NAMES openblas HINTS /mingw64/lib /opt/local/lib /usr/local/lib /usr/lib)
endif()
if(SB_OPENBLAS)
  get_filename_component(_obdir ${SB_OPENBLAS} DIRECTORY)
  get_filename_component(_obname ${SB_OPENBLAS} NAME_WE)
  string(REGEX REPLACE "^lib" "" _obname "${_obname}")
  set(SB_BLAS_LFLAGS "-L${_obdir} -l${_obname}")
else()
  set(SB_BLAS_LFLAGS "")   # ThirdParty scripts fall back to a reference BLAS
endif()

# GNU gcc supports OpenMP natively (-fopenmp), so no libomp flag gymnastics.
set(SB_OMP_ARGS "")

message(STATUS "PSOPT superbuild -> ${SB_INSTALL}")
message(STATUS "  CC=${SB_CC} CXX=${SB_CXX} FC=${SB_FC} BLAS='${SB_BLAS_LFLAGS}'")

set(_sb_common
  CMAKE_ARGS -DCMAKE_INSTALL_PREFIX=${SB_INSTALL} -DCMAKE_PREFIX_PATH=${SB_INSTALL}
             -DCMAKE_C_COMPILER=${SB_CC} -DCMAKE_CXX_COMPILER=${SB_CXX}
             -DBUILD_TESTING=OFF)

if(CMAKE_TOOLCHAIN_FILE)
  list(APPEND _sb_common "-DCMAKE_TOOLCHAIN_FILE=${CMAKE_TOOLCHAIN_FILE}")
endif()
if(VCPKG_TARGET_TRIPLET)
  list(APPEND _sb_common "-DVCPKG_TARGET_TRIPLET=${VCPKG_TARGET_TRIPLET}")
endif()

# ---- Eigen (CMake, header-only install) -------------------------------------
ExternalProject_Add(ep_eigen
  GIT_REPOSITORY https://gitlab.com/libeigen/eigen.git GIT_TAG 3.4.0 GIT_SHALLOW TRUE
  ${_sb_common} CMAKE_ARGS -DCMAKE_INSTALL_PREFIX=${SB_INSTALL} -DEIGEN_BUILD_DOC=OFF -DBUILD_TESTING=OFF)

# ---- ColPack + ADOL-C (autotools) -------------------------------------------
# ColPack supplies ADOL-C's sparse Jacobian/Hessian drivers (sparse_jac/sparse_hess),
# which PSOPT uses on the IPOPT/Fortran path. It is only required when Fortran is
# present; on the C++-only (CasADi/finite-difference) path ADOL-C is built WITHOUT
# sparse and ColPack is skipped. (ColPack itself is always git-downloadable.)
if(SB_HAVE_FORTRAN)
  # Pin ColPack to v1.0.10 — the release whose GenerateSeedHessian/Jacobian API
  # matches ADOL-C 2.7.2 (master ColPack changed those signatures). Top-level
  # autotools; builds cleanly under GNU gcc.
  ExternalProject_Add(ep_colpack
    GIT_REPOSITORY https://github.com/CSCsw/ColPack.git GIT_TAG v1.0.10 GIT_SHALLOW TRUE
    # --disable-dependency-tracking: newer automake (e.g. Ubuntu 26.04) fails to
    # bootstrap ColPack's dep-tracking makefile fragments; disabling is safe here.
    CONFIGURE_COMMAND cd <SOURCE_DIR> && autoreconf -fi && ./configure --prefix=${SB_INSTALL} --disable-dependency-tracking --disable-shared CC=${SB_CC} CXX=${SB_CXX}
    BUILD_COMMAND cd <SOURCE_DIR> && make -j
    INSTALL_COMMAND cd <SOURCE_DIR> && make install
    BUILD_IN_SOURCE 1)
  # Pin ADOL-C to release 2.7.2: it still ships the classic <adolc/adouble.h> that
  # PSOPT includes (ADOL-C master removed it), has a committed ./configure, and
  # installs adolc.pc. Autotools build with sparse drivers against our ColPack.
  # ADOL-C 2.7.2 looks for libColPack under <colpack>/lib64; our ColPack installs
  # to lib, so provide a lib64 -> lib symlink first (idempotent).
  # macOS rejects the undefined ColPack symbols ADOL-C 2.7.2 leaves in
  # libadolc.dylib (its Makefile doesn't add -lColPack to the shared lib, which
  # only Linux tolerates). Link ColPack explicitly via LDFLAGS/LIBS.
  # --without-boost: Boost 1.69+ made Boost.System header-only, so ADOL-C 2.7.2's
  # configure fails linking -lboost_system on newer distros (e.g. Ubuntu 26.04).
  # Boost is only an optional pool-allocator optimisation for ADOL-C; drop it.
  # NB: do NOT put linker force-load flags (e.g. --whole-archive) in LIBS — configure
  # link-tests use $LIBS, and force-loading all of ColPack makes those tests fail, giving
  # false negatives for ColPack AND trunc/fmax (the latter then makes common.h define a
  # trunc macro that clashes with MinGW's math.h). Plain -lColPack keeps configure happy.
  #
  # On MinGW, build ADOL-C STATIC (--disable-shared): a shared libadolc.dll must fully
  # resolve ColPack's (mutually-referencing) sparse-driver symbols at build time, and
  # single-pass -lColPack drops them ("undefined reference to ColPack::..."). A static
  # libadolc.a defers that to PSOPT's final executable link, where ColPack (added via
  # Findadolc) is pulled in and its internal deps resolve. macOS/Linux keep shared.
  if(WIN32 AND NOT MSVC)
    set(_adolc_libtype --disable-shared --enable-static)
    # MinGW declares trunc as an inline in <math.h>, which conflicts with autoconf's fake
    # 'char trunc()' prototype, so AC_CHECK_FUNCS(trunc) fails and HAVE_TRUNC stays unset ->
    # common.h then #defines a trunc() macro that clashes with <math.h> ("expected ')' before
    # 'double'"). Force the correct result via the configure cache. (Static build too, above.)
    set(_adolc_win_cache ac_cv_func_trunc=yes)
  else()
    set(_adolc_libtype "")
    set(_adolc_win_cache "")
  endif()
  ExternalProject_Add(ep_adolc
    DEPENDS ep_colpack
    GIT_REPOSITORY https://github.com/coin-or/ADOL-C.git GIT_TAG releases/2.7.2 GIT_SHALLOW TRUE
    CONFIGURE_COMMAND ln -sfn lib ${SB_INSTALL}/lib64 && cd <SOURCE_DIR> && ./configure --prefix=${SB_INSTALL} --with-colpack=${SB_INSTALL} --enable-sparse ${_adolc_libtype} --disable-dependency-tracking --without-boost CC=${SB_CC} CXX=${SB_CXX} ${_adolc_win_cache} "LDFLAGS=-L${SB_INSTALL}/lib -Wl,-rpath,${SB_INSTALL}/lib" "LIBS=-lColPack"
    BUILD_COMMAND cd <SOURCE_DIR> && make -j
    INSTALL_COMMAND cd <SOURCE_DIR> && make install
    BUILD_IN_SOURCE 1)
elseif(MSVC)
  # No Fortran, MSVC: ADOL-C 2.7.2 has no CMake. Build its Visual Studio solution with
  # msbuild and copy the artifacts into the install prefix. ColPack/sparse are not used
  # (CasADi path is finite-difference). EXPERIMENTAL/UNTESTED — the .sln output paths and
  # ADOL-C's MSVC config.h may need adjustment for your VS version.
  if(MSVC_VERSION GREATER_EQUAL 1930)
    set(_adolc_toolset "v143")
  elseif(MSVC_VERSION GREATER_EQUAL 1920)
    set(_adolc_toolset "v142")
  else()
    set(_adolc_toolset "v140")
  endif()
  ExternalProject_Add(ep_adolc
    GIT_REPOSITORY https://github.com/coin-or/ADOL-C.git GIT_TAG releases/2.7.2 GIT_SHALLOW TRUE
    CONFIGURE_COMMAND ${CMAKE_COMMAND} -DFILE=<SOURCE_DIR>/MSVisualStudio/v14/x64/nosparse/config.h -P ${CMAKE_SOURCE_DIR}/cmake/strip_boost.cmake
            COMMAND ${CMAKE_COMMAND} -DFILE=<SOURCE_DIR>/MSVisualStudio/v14/nosparse/config.h -P ${CMAKE_SOURCE_DIR}/cmake/strip_boost.cmake
            COMMAND ${CMAKE_COMMAND} -DFILE=<SOURCE_DIR>/ADOL-C/include/adolc/internal/adolc_settings.h -P ${CMAKE_SOURCE_DIR}/cmake/strip_boost.cmake
    BUILD_COMMAND msbuild <SOURCE_DIR>/MSVisualStudio/v14/adolc.vcxproj /p:Configuration=nosparse /p:Platform=x64 /p:WindowsTargetPlatformVersion=10.0 /p:PlatformToolset=${_adolc_toolset} /p:IntDir=x64/nosparse/
    INSTALL_COMMAND ${CMAKE_COMMAND} -E copy_directory <SOURCE_DIR>/ADOL-C/include ${SB_INSTALL}/include
            COMMAND ${CMAKE_COMMAND} -E make_directory ${SB_INSTALL}/lib ${SB_INSTALL}/bin
            COMMAND ${CMAKE_COMMAND} -E copy <SOURCE_DIR>/MSVisualStudio/v14/x64/nosparse/adolc.lib ${SB_INSTALL}/lib/adolc.lib
            COMMAND ${CMAKE_COMMAND} -E copy_if_different <SOURCE_DIR>/MSVisualStudio/v14/x64/nosparse/adolc.dll ${SB_INSTALL}/bin/adolc.dll
    BUILD_IN_SOURCE 1)
else()
  # No Fortran, MinGW/MSYS2/Unix: autotools build of ADOL-C WITHOUT sparse (no ColPack).
  ExternalProject_Add(ep_adolc
    GIT_REPOSITORY https://github.com/coin-or/ADOL-C.git GIT_TAG releases/2.7.2 GIT_SHALLOW TRUE
    CONFIGURE_COMMAND cd <SOURCE_DIR> && ./configure --prefix=${SB_INSTALL} --disable-sparse --disable-dependency-tracking --without-boost CC=${SB_CC} CXX=${SB_CXX}
    BUILD_COMMAND cd <SOURCE_DIR> && make -j
    INSTALL_COMMAND cd <SOURCE_DIR> && make install
    BUILD_IN_SOURCE 1)
endif()

# ==== Fortran solver stack (METIS/MUMPS/IPOPT + opt-in HSL/SPRAL/PARDISO) ======
# Everything below needs a Fortran compiler; skipped entirely on the C++-only path.
# MUMPS/HSL/IPOPT configure via `cmake -E env ... <SOURCE_DIR>/configure`, which execs the
# script directly. Unix runs it via the shebang; Windows cannot exec a shell script
# ("inappropriate file type or format"), so run it through sh explicitly there.
set(_sb_mumps_deps ep_metis)
if(WIN32 AND NOT MSVC)
  # ThirdParty-Mumps depends on COIN-OR BuildTools m4 macros (coin.m4,
  # coin_chk_lapack.m4, coin_chk_libhdr.m4, coin_fortran.m4). They are not part
  # of the ThirdParty-Mumps checkout, so fetch the BuildTools repo once and feed
  # its m4 directory to autoreconf.
  ExternalProject_Add(ep_buildtools
    GIT_REPOSITORY https://github.com/coin-or-tools/BuildTools.git
    GIT_TAG 743a4f662032b246661042902da1a2e374b956c5
    GIT_SHALLOW TRUE
    CONFIGURE_COMMAND ""
    BUILD_COMMAND ""
    INSTALL_COMMAND ""
    BUILD_IN_SOURCE 1)
  ExternalProject_Get_Property(ep_buildtools source_dir)
  set(_sb_buildtools_m4 "${source_dir}")
  set(_sb_sh sh)
  # Clean up mock symlinks so autoreconf/libtoolize/automake copy actual system versions.
  set(_sb_reconf rm -f config.guess config.sub install-sh compile missing ar-lib ltmain.sh depcomp && autoreconf -I ${_sb_buildtools_m4} -fi)
  list(APPEND _sb_mumps_deps ep_buildtools)
else()
  set(_sb_sh "")
  set(_sb_reconf ${CMAKE_COMMAND} -E true)
endif()
set(_psopt_ipopt_dep "")
set(_casadi_ipopt_dep "")
set(_casadi_with_ipopt OFF)
set(_psopt_scip_dep "")
if(SB_HAVE_FORTRAN)

# ---- METIS (fill-reducing ordering; speeds MUMPS and the HSL solvers) --------
# ThirdParty-Metis installs libcoinmetis + coinmetis.pc; MUMPS and HSL pick it
# up automatically via pkg-config (PKG_CONFIG_PATH below).
ExternalProject_Add(ep_metis
  GIT_REPOSITORY https://github.com/coin-or-tools/ThirdParty-Metis.git GIT_TAG stable/2.0 GIT_SHALLOW TRUE
  CONFIGURE_COMMAND cd <SOURCE_DIR> && ./get.Metis && ./configure --prefix=${SB_INSTALL} CC=${SB_CC} "CFLAGS=${SB_LEGACY_CFLAGS}"
  BUILD_COMMAND cd <SOURCE_DIR> && make -j
  INSTALL_COMMAND cd <SOURCE_DIR> && make install
  BUILD_IN_SOURCE 1)

# ---- MUMPS (sequential, via COIN-OR ThirdParty; uses METIS ordering) --------
ExternalProject_Add(ep_mumps
  DEPENDS ${_sb_mumps_deps}
  GIT_REPOSITORY https://github.com/coin-or-tools/ThirdParty-Mumps.git GIT_TAG stable/3.0 GIT_SHALLOW TRUE
  CONFIGURE_COMMAND ./get.Mumps && ${_sb_reconf}
          COMMAND ${CMAKE_COMMAND} -E env PKG_CONFIG_PATH=${SB_PKGCFG} ${_sb_sh} ./configure --prefix=${SB_INSTALL} CC=${SB_CC} FC=${SB_FC} --with-metis --with-lapack-lflags=${SB_BLAS_LFLAGS} --with-blas-lflags=${SB_BLAS_LFLAGS} ${SB_INT64_MUMPS_FLAGS}
  # MUMPS' Makefile does not order its Fortran module dependencies, so a parallel
  # `make -j` intermittently fails ("Cannot open module file 'mumps_lr_common.mod'").
  # Try parallel first, then fall back to a serial `make` which resolves the modules
  # in order and finishes whatever the parallel pass missed.
  BUILD_COMMAND cd <SOURCE_DIR> && sh -c "make -j || make -j1"
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
            COMMAND ${CMAKE_COMMAND} -E env PKG_CONFIG_PATH=${SB_PKGCFG} ${_sb_sh} <SOURCE_DIR>/configure --prefix=${SB_INSTALL} CC=${SB_CC} FC=${SB_FC} --with-metis --with-lapack-lflags=${SB_BLAS_LFLAGS}
    BUILD_COMMAND cd <SOURCE_DIR> && make -j
    INSTALL_COMMAND cd <SOURCE_DIR> && make install
    BUILD_IN_SOURCE 1)
  set(_ipopt_hsl_dep ep_hsl)
  set(_ipopt_hsl_flag "--with-hsl")
endif()

# PARDISO: build IPOPT with the (commercial Panua/Intel-MKL) PARDISO solver when
# the user supplies its link flags in PSOPT_PARDISO_LFLAGS. PARDISO/MKL are NOT
# redistributed and are typically unavailable on arm64 macOS; without them
# IPOPT is MUMPS(+optional HSL)-only. Once built, select at runtime with
# algorithm.ipopt_linear_solver="pardiso" (Panua) or "pardisomkl" (MKL).
set(_ipopt_pardiso_flag "--without-pardiso")
if(PSOPT_PARDISO_LFLAGS)
  set(_ipopt_pardiso_flag "--with-pardiso=${PSOPT_PARDISO_LFLAGS}")
  message(STATUS "IPOPT: building with PARDISO (${PSOPT_PARDISO_LFLAGS})")
endif()

# SPRAL (ssids): open-source OpenMP-threaded sparse symmetric-indefinite solver
# from STFC, a free alternative to PARDISO/ma97. Built with Meson; needs METIS 5,
# BLAS/LAPACK (OpenBLAS), and hwloc (MacPorts /opt/local by default). Select at
# runtime with algorithm.ipopt_linear_solver="spral". Verified: bryson_denham ->
# 3.999539, obstacle -> 4.571044 (same optima as MUMPS). NOTE: SSIDS REQUIRES an
# OpenMP environment or it aborts (info.flag=-53) -- run with, e.g.:
#   OMP_CANCELLATION=TRUE OMP_NESTED=TRUE OMP_PROC_BIND=TRUE OMP_STACKSIZE=64M
# Its threaded factorization pays off on LARGE problems; on small ones it is
# comparable to (slightly slower than) serial MUMPS.
set(_ipopt_spral_flag "--without-spral")
set(_ipopt_spral_cflag "")
set(_ipopt_spral_dep "")
set(PSOPT_SPRAL_DEP_PREFIX "/opt/local" CACHE PATH "Prefix providing METIS 5 / OpenBLAS / hwloc for SPRAL")
if(PSOPT_WITH_SPRAL)
  ExternalProject_Add(ep_spral
    GIT_REPOSITORY https://github.com/ralna/spral.git GIT_TAG master GIT_SHALLOW TRUE
    CONFIGURE_COMMAND ${CMAKE_COMMAND} -E env CC=${SB_CC} CXX=${SB_CXX} FC=${SB_FC}
        meson setup <BINARY_DIR> <SOURCE_DIR> --prefix=${SB_INSTALL} -Ddefault_library=shared
        -Dlibblas=openblas -Dliblapack=openblas -Dlibmetis=metis -Dlibhwloc=hwloc
        -Dlibmetis_path=${PSOPT_SPRAL_DEP_PREFIX}/lib -Dlibblas_path=${PSOPT_SPRAL_DEP_PREFIX}/lib -Dlibhwloc_path=${PSOPT_SPRAL_DEP_PREFIX}/lib
        -Dopenmp=true -Dexamples=false -Dtests=false -Dbinaries=false
    BUILD_COMMAND meson compile -C <BINARY_DIR>
    INSTALL_COMMAND meson install -C <BINARY_DIR>)
  set(_ipopt_spral_dep ep_spral)
  set(_ipopt_spral_flag "--with-spral-lflags=-L${SB_INSTALL}/lib -lspral -Wl,-rpath,${SB_INSTALL}/lib -L${PSOPT_SPRAL_DEP_PREFIX}/lib -lopenblas -lmetis -lhwloc")
  set(_ipopt_spral_cflag "--with-spral-cflags=-I${SB_INSTALL}/include")
endif()

# ---- IPOPT (autotools; uses coinmumps + optional coinhsl via pkg-config) ----
ExternalProject_Add(ep_ipopt
  DEPENDS ep_mumps ${_ipopt_hsl_dep} ${_ipopt_spral_dep}
  GIT_REPOSITORY https://github.com/coin-or/Ipopt.git GIT_TAG stable/3.14 GIT_SHALLOW TRUE
  CONFIGURE_COMMAND ${CMAKE_COMMAND} -E env PKG_CONFIG_PATH=${SB_PKGCFG}
      ${_sb_sh} ./configure --prefix=${SB_INSTALL}
      CC=${SB_CC} CXX=${SB_CXX} FC=${SB_FC}
      --with-mumps ${_ipopt_hsl_flag} ${_ipopt_pardiso_flag} ${_ipopt_spral_flag} ${_ipopt_spral_cflag}
      --with-lapack-lflags=${SB_BLAS_LFLAGS} --with-blas-lflags=${SB_BLAS_LFLAGS}
      --disable-java --without-asl
  BUILD_COMMAND make -j
  INSTALL_COMMAND make install
  BUILD_IN_SOURCE 1)
set(_psopt_ipopt_dep ep_ipopt)
set(_casadi_ipopt_dep ep_ipopt)
set(_casadi_with_ipopt ON)

# ---- SoPlex + SCIP (opt-in; free MILP/MINLP for mixed-integer optimal control) ----
# SoPlex is SCIP's LP backend; SCIP uses the superbuild IPOPT as its NLP relaxation
# solver for MINLP. GMP (both) and readline (SCIP shell) come from
# PSOPT_SCIP_DEP_PREFIX (MacPorts /opt/local by default; a Homebrew/apt prefix works too).
set(_psopt_scip_dep "")
if(PSOPT_WITH_SCIP)
  ExternalProject_Add(ep_soplex
    GIT_REPOSITORY https://github.com/scipopt/soplex.git GIT_TAG master GIT_SHALLOW TRUE
    CMAKE_ARGS -DCMAKE_INSTALL_PREFIX=${SB_INSTALL}
               -DCMAKE_PREFIX_PATH=${SB_INSTALL}\;${PSOPT_SCIP_DEP_PREFIX}
               -DCMAKE_C_COMPILER=${SB_CC} -DCMAKE_CXX_COMPILER=${SB_CXX}
               -DCMAKE_BUILD_TYPE=Release -DGMP=on -DBOOST=off -DBUILD_TESTING=OFF)
  ExternalProject_Add(ep_scip
    DEPENDS ep_soplex ep_ipopt
    GIT_REPOSITORY https://github.com/scipopt/scip.git GIT_TAG master GIT_SHALLOW TRUE
    CMAKE_ARGS -DCMAKE_INSTALL_PREFIX=${SB_INSTALL}
               -DCMAKE_PREFIX_PATH=${SB_INSTALL}\;${PSOPT_SCIP_DEP_PREFIX}
               -DCMAKE_C_COMPILER=${SB_CC} -DCMAKE_CXX_COMPILER=${SB_CXX}
               -DCMAKE_BUILD_TYPE=Release
               -DSOPLEX_DIR=${SB_INSTALL} -DGMP=on -DREADLINE=on
               -DIPOPT=on -DIPOPT_DIR=${SB_INSTALL}
               # AMPL=off: SCIP's bundled AMPL .nl reader (amplmp/src/os.cpp) pulls in
               # Apple <mach/message.h>, which MacPorts gcc-15 cannot parse. Not needed
               # for the transcription-based MIOC use case.
               -DAMPL=off -DZIMPL=off -DPAPILO=off -DSYM=none -DBUILD_TESTING=OFF)
  set(_psopt_scip_dep ep_scip)
endif()

endif()  # SB_HAVE_FORTRAN — end of the Fortran solver stack (METIS/MUMPS/IPOPT/…)

# ---- CasADi (CMake; with the IPOPT plugin, optional HSL) --------------------
set(_casadi_dep "")
if(PSOPT_WITH_CASADI)
  set(_casadi_build_lapack OFF)
  if(MSVC AND NOT CMAKE_TOOLCHAIN_FILE)
    set(_casadi_build_lapack ON)
  endif()
  ExternalProject_Add(ep_casadi
    DEPENDS ${_casadi_ipopt_dep}
    GIT_REPOSITORY https://github.com/casadi/casadi.git GIT_TAG main GIT_SHALLOW TRUE
    ${_sb_common}
    # WITH_IPOPT only when the Fortran path built IPOPT; on the C++-only path CasADi
    # drives the NLP with its own solvers (sqpmethod+qpOASES / blockSQP).
    CMAKE_ARGS -DCMAKE_INSTALL_PREFIX=${SB_INSTALL} -DCMAKE_PREFIX_PATH=${SB_INSTALL}
               -DBLA_VENDOR=Generic
               -DWITH_IPOPT=${_casadi_with_ipopt} -DWITH_HSL=${PSOPT_WITH_HSL} -DWITH_PYTHON=OFF -DBUILD_TESTING=OFF
               -DWITH_BUILD_LAPACK=${_casadi_build_lapack}
               # Alternative NLP/QP solver plugins so algorithm.casadi_solver can be
               # sqpmethod (needs a QP solver), blocksqp, fatrop, etc. — not just ipopt.
               -DWITH_LAPACK=ON                            # required by qpOASES
               -DWITH_QPOASES=ON -DWITH_BUILD_QPOASES=ON   # qpOASES QP -> enables sqpmethod
               -DWITH_BLOCKSQP=ON)                         # blockSQP NLP
               # (OSQP/fatrop omitted: CasADi's bundled external_projects builds
               #  of libosqp/libblasfeo fail on macOS 26. qpOASES builds inline
               #  and enables sqpmethod; blockSQP is built in.)
  set(_casadi_dep ep_casadi)
endif()

# ---- PSOPT itself (re-invoke this CMake with superbuild off, deps resolved) --
# With HSL built, default to the OpenMP-parallel ma97 solver.
set(_psopt_default_solver ${PSOPT_DEFAULT_LINEAR_SOLVER})
if(PSOPT_WITH_HSL AND _psopt_default_solver STREQUAL "mumps")
  set(_psopt_default_solver "ma97")
endif()
# On MSVC, PSOPT's own remaining deps (Boost/OpenBLAS) come from vcpkg — forward its
# toolchain so ep_psopt can find them.
set(_sb_psopt_extra "")
if(MSVC AND DEFINED ENV{VCPKG_ROOT})
  set(_sb_psopt_extra "-DCMAKE_TOOLCHAIN_FILE=$ENV{VCPKG_ROOT}/scripts/buildsystems/vcpkg.cmake")
endif()
ExternalProject_Add(ep_psopt
  DEPENDS ep_eigen ep_adolc ${_psopt_ipopt_dep} ${_casadi_dep} ${_psopt_scip_dep}
  SOURCE_DIR ${CMAKE_SOURCE_DIR}
  CMAKE_ARGS
    -DPSOPT_SUPERBUILD=OFF
    -DCMAKE_INSTALL_PREFIX=${CMAKE_INSTALL_PREFIX}
    -DCMAKE_PREFIX_PATH=${SB_INSTALL}
    -DIPOPT_ROOT=${SB_INSTALL}
    # IPOPT is built only on the Fortran path, so tell PSOPT whether to expect it.
    -DPSOPT_WITH_IPOPT=${SB_HAVE_FORTRAN}
    ${_sb_psopt_extra}
    -DCMAKE_C_COMPILER=${SB_CC} -DCMAKE_CXX_COMPILER=${SB_CXX}
    -DPSOPT_WITH_MUMPS=${PSOPT_WITH_MUMPS} -DPSOPT_WITH_HSL=${PSOPT_WITH_HSL}
    -DPSOPT_WITH_CASADI=${PSOPT_WITH_CASADI} -DPSOPT_WITH_SNOPT=${PSOPT_WITH_SNOPT}
    -DPSOPT_WITH_OPENMP=${PSOPT_WITH_OPENMP} -DPSOPT_WITH_MPI=${PSOPT_WITH_MPI}
    -DPSOPT_WITH_INT64=${PSOPT_WITH_INT64}
    -DPSOPT_DEFAULT_LINEAR_SOLVER=${_psopt_default_solver}
    -DBUILD_EXAMPLES=${BUILD_EXAMPLES} -DBUILD_TESTS=${BUILD_TESTS} -DHEADLESS=${HEADLESS}
  INSTALL_COMMAND "")

message(STATUS "Superbuild configured. Run the build to download+compile all dependencies and PSOPT.")
