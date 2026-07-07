cmake_minimum_required(VERSION 3.16)

# ────────────────────────────────────────────────────────────────
# VERSIONING (unchanged behaviour: derive version from git tag/date)
# ────────────────────────────────────────────────────────────────
execute_process(COMMAND git describe --tags --abbrev=0
  WORKING_DIRECTORY ${CMAKE_SOURCE_DIR} OUTPUT_VARIABLE PSOPT_TAG_RAW
  ERROR_QUIET OUTPUT_STRIP_TRAILING_WHITESPACE)
if(PSOPT_TAG_RAW STREQUAL "")
  string(TIMESTAMP PSOPT_VERSION "20%y.%m.%d")
else()
  string(REPLACE "-" "." PSOPT_VERSION "${PSOPT_TAG_RAW}")
endif()
execute_process(COMMAND git describe --tags --always
  WORKING_DIRECTORY ${CMAKE_SOURCE_DIR} OUTPUT_VARIABLE PSOPT_RELEASE_STRING
  ERROR_QUIET OUTPUT_STRIP_TRAILING_WHITESPACE)
string(TIMESTAMP PSOPT_BUILD_DATE "%Y-%m-%d %H:%M:%S")

# ────────────────────────────────────────────────────────────────
# ROBUST COMPILER SELECTION  (must precede project())
#   Preference:  GNU  ›  Clang  ›  MSVC.
#   Searches package managers and OSes: MacPorts, Homebrew, Fink, Nix,
#   MSYS2 / native Windows, the BSDs (Free/Open/Net/DragonFly), and Linux.
#   An explicit -DCMAKE_CXX_COMPILER=... (or $CXX / a toolchain file) always
#   wins; set -DPSOPT_AUTO_COMPILER=OFF to disable this heuristic entirely.
#   NB on macOS: GNU is favoured as requested, but Apple's SDK headers can clash
#   with GCC (e.g. GoogleTest's <mach/*.h>); force clang there with
#   -DCMAKE_CXX_COMPILER=clang++ if a GCC build of a given target fails.
# ────────────────────────────────────────────────────────────────
option(PSOPT_AUTO_COMPILER "Auto-pick a working GNU/Clang/MSVC toolchain if none is specified" ON)
if(PSOPT_AUTO_COMPILER AND NOT CMAKE_CXX_COMPILER AND NOT DEFINED ENV{CXX} AND NOT DEFINED ENV{CC}
   AND NOT CMAKE_TOOLCHAIN_FILE)
  set(_pm_bindirs
      /opt/local/bin                                   # MacPorts   (favoured on macOS)
      /opt/homebrew/bin /usr/local/bin                 # Homebrew (arm64 / x86_64)
      /sw/bin                                          # Fink
      "$ENV{HOME}/.nix-profile/bin"                    # Nix (per-user profile)
      /run/current-system/sw/bin                       # Nix (NixOS system)
      /nix/var/nix/profiles/default/bin                # Nix (default profile)
      C:/msys64/ucrt64/bin C:/msys64/mingw64/bin       # MSYS2 (Windows, out-of-shell)
      /ucrt64/bin /mingw64/bin                         # MSYS2 (in MSYS shell)
      /usr/pkg/bin                                     # NetBSD pkgsrc
      /usr/local/bin /usr/bin)                         # *BSD / Linux / generic
  # (C++|C) candidate pairs, most-preferred first.  GNU is favoured EXCEPT on
  # macOS, where Apple's SDK headers clash with GCC (e.g. GoogleTest's
  # <mach/*.h>) -- there clang is preferred so `cmake -B build && ctest` works
  # out of the box (GCC still reachable via -DCMAKE_CXX_COMPILER=g++-mp-15).
  # MSVC is the implicit fallback (leave unset).
  if(CMAKE_HOST_APPLE)
    set(_cxx_cc_pairs
        "clang++|clang"                                                  # Apple / MacPorts clang
        "clang++-20|clang-20" "clang++-19|clang-19" "clang++-18|clang-18"
        "g++-mp-15|gcc-mp-15" "g++-mp-14|gcc-mp-14"                       # GCC fallback (MacPorts)
        "g++-15|gcc-15" "g++-14|gcc-14" "g++|gcc")
  else()
    set(_cxx_cc_pairs
        "g++-mp-15|gcc-mp-15" "g++-mp-14|gcc-mp-14" "g++-mp-13|gcc-mp-13" # MacPorts (non-Apple hosts too)
        "g++-15|gcc-15" "g++-14|gcc-14" "g++-13|gcc-13" "g++-12|gcc-12" "g++-11|gcc-11"
        "g++|gcc" "eg++|egcc"                                             # generic GNU / OpenBSD
        "clang++|clang" "clang++-20|clang-20" "clang++-19|clang-19"
        "clang++-18|clang-18" "clang++-17|clang-17")
  endif()
  foreach(_pair IN LISTS _cxx_cc_pairs)
    string(REPLACE "|" ";" _pc "${_pair}")
    list(GET _pc 0 _cxx_name)
    list(GET _pc 1 _cc_name)
    find_program(_PSOPT_CXX NAMES ${_cxx_name} HINTS ${_pm_bindirs})
    find_program(_PSOPT_CC  NAMES ${_cc_name}  HINTS ${_pm_bindirs})
    if(_PSOPT_CXX AND _PSOPT_CC)
      set(CMAKE_CXX_COMPILER "${_PSOPT_CXX}" CACHE FILEPATH "C++ compiler auto-selected by PSOPT")
      set(CMAKE_C_COMPILER   "${_PSOPT_CC}"  CACHE FILEPATH "C compiler auto-selected by PSOPT")
      message(STATUS "PSOPT auto-selected toolchain: ${_PSOPT_CXX} | ${_PSOPT_CC}")
      break()
    endif()
    unset(_PSOPT_CXX CACHE)
    unset(_PSOPT_CC CACHE)
  endforeach()
  if(NOT CMAKE_CXX_COMPILER)
    message(STATUS "PSOPT: no GNU/Clang toolchain found in known prefixes; "
                   "using CMake default (e.g. MSVC on Windows, cc/c++ elsewhere)")
  endif()
endif()

project(PSOPT VERSION ${PSOPT_VERSION} LANGUAGES C CXX)

list(APPEND CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/cmake)
include(GNUInstallDirs)

# ────────────────────────────────────────────────────────────────
# OPTIONS — choose the solver stack; build any missing dep ourselves
# ────────────────────────────────────────────────────────────────
option(PSOPT_WITH_IPOPT   "IPOPT NLP backend (headers + link). OFF -> CasADi-only, no IPOPT/Fortran" ON)
option(PSOPT_WITH_MUMPS   "IPOPT + MUMPS linear solver"                        ON)
option(PSOPT_WITH_HSL     "IPOPT + HSL (MA27/57/97); needs Coin-HSL source"    OFF)
option(PSOPT_WITH_CASADI  "Full CasADi NLP backend (CasADi -> IPOPT[+HSL])"    OFF)
option(PSOPT_WITH_SNOPT   "SNOPT NLP backend"                                  OFF)
option(PSOPT_WITH_OPENMP  "Enable OpenMP (threaded BLAS/HSL/ColPack, parallel MUMPS)" OFF)
option(PSOPT_WITH_MPI     "Enable MPI (distributed-memory MUMPS; see macOS-26 caveat)" OFF)
option(PSOPT_SUPERBUILD   "Download and build missing dependencies from source" OFF)
option(PSOPT_WITH_INT64   "64-bit integer addressing (ILP64) for very large NLPs" OFF)
option(PSOPT_WITH_SPRAL   "Build the open-source SPRAL (ssids) threaded sparse solver into IPOPT" OFF)
option(PSOPT_WITH_SCIP    "Build SoPlex + SCIP (free MILP/MINLP solver) for mixed-integer optimal control" OFF)
set(PSOPT_SCIP_DEP_PREFIX "/opt/local" CACHE PATH
    "Prefix providing GMP and readline for the SCIP/SoPlex build (MacPorts /opt/local, or Homebrew/apt prefix)")
option(BUILD_EXAMPLES     "Build examples from the example subdirectory."      OFF)
option(BUILD_TESTS        "Build tests for PSOPT."                             OFF)
option(HEADLESS           "Gnuplot writes PDFs instead of spawning windows."   OFF)
set(PSOPT_DEFAULT_LINEAR_SOLVER "mumps" CACHE STRING
    "Default IPOPT linear_solver (mumps, ma27, ma57, ma86, ma97, pardiso)")
# PSOPT_WITH_IPOPT consistency: MUMPS is IPOPT's solver, and SCIP's MINLP driver
# resolves NLPs through IPOPT — both are meaningless without it. A no-IPOPT build must
# therefore also drop MUMPS/SCIP and provide an alternative NLP backend (CasADi).
if(NOT PSOPT_WITH_IPOPT)
  if(PSOPT_WITH_MUMPS)
    message(STATUS "PSOPT_WITH_IPOPT=OFF -> forcing PSOPT_WITH_MUMPS=OFF (MUMPS is IPOPT's solver)")
    set(PSOPT_WITH_MUMPS OFF CACHE BOOL "" FORCE)
  endif()
  if(PSOPT_WITH_SCIP)
    message(FATAL_ERROR "PSOPT_WITH_SCIP=ON requires PSOPT_WITH_IPOPT=ON (SCIP resolves NLP subproblems through IPOPT).")
  endif()
  if(NOT PSOPT_WITH_CASADI)
    message(FATAL_ERROR "PSOPT_WITH_IPOPT=OFF needs an alternative NLP backend: also set -DPSOPT_WITH_CASADI=ON.")
  endif()
endif()

set(COINHSL_SOURCE_DIR "" CACHE PATH "Path to a local Coin-HSL source tree (skips download)")
set(PSOPT_PARDISO_LFLAGS "" CACHE STRING
    "Link flags for a PARDISO library (commercial Panua PARDISO or Intel MKL) to build IPOPT with; then select algorithm.ipopt_linear_solver=\"pardiso\" or \"pardisomkl\"")
option(PSOPT_HSL_ACCEPT_LICENSE "I accept the non-free HSL licence; auto-download Coin-HSL" OFF)
set(PSOPT_HSL_URL
    "https://www.hsl.rl.ac.uk/download/coinhsl/2023.11.17/coinhsl-2023.11.17.tar.gz"
    CACHE STRING "Coin-HSL source archive URL (override with your licensed/personalised link)")

# Backwards-compat alias for the old option name.
if(DEFINED WITH_SNOPT_INTERFACE)
  set(PSOPT_WITH_SNOPT ${WITH_SNOPT_INTERFACE})
endif()

# ────────────────────────────────────────────────────────────────
# SUPERBUILD DISPATCH: build all deps (+ PSOPT itself) via ExternalProject
# ────────────────────────────────────────────────────────────────
if(PSOPT_SUPERBUILD)
  include(superbuild)   # cmake/superbuild.cmake — sets up ExternalProjects and returns
  return()
endif()

# ────────────────────────────────────────────────────────────────
# DEPENDENCIES (find, or light fetch; heavy deps come from the superbuild)
# ────────────────────────────────────────────────────────────────
include(dependencies)   # defines Eigen3::Eigen, adolc, IPOPT::ipopt, (casadi)

# ────────────────────────────────────────────────────────────────
# OPTIONAL BACKEND: SNOPT (existing object-library pattern)
# ────────────────────────────────────────────────────────────────
if(PSOPT_WITH_SNOPT)
  add_subdirectory(snopt-interface)
endif()

# ────────────────────────────────────────────────────────────────
# PSOPT LIBRARY TARGET
# ────────────────────────────────────────────────────────────────
file(GLOB SRC    CONFIGURE_DEPENDS src/*.cxx)
file(GLOB HEADER CONFIGURE_DEPENDS include/*.h)
if(NOT PSOPT_WITH_IPOPT)
  # Exclude the IPOPT backend TU; the IPOPT branch in NLP_interface.cxx is #ifdef'd out.
  list(REMOVE_ITEM SRC ${CMAKE_CURRENT_SOURCE_DIR}/src/IPOPT_interface.cxx)
endif()
if(PSOPT_WITH_CASADI)
  # CASADI_interface.cxx is compiled in; the CASADI branch is guarded by USE_CASADI
  list(APPEND SRC ${CMAKE_CURRENT_SOURCE_DIR}/src/CASADI_interface.cxx)
  list(REMOVE_DUPLICATES SRC)
endif()

add_library(psopt
  $<$<BOOL:$<TARGET_EXISTS:PSOPT_SNOPT_interface>>:$<TARGET_OBJECTS:PSOPT_SNOPT_interface>>
  ${SRC} ${HEADER})
add_library(PSOPT::psopt ALIAS psopt)
set_target_properties(psopt PROPERTIES OUTPUT_NAME PSOPT)

target_include_directories(psopt PUBLIC
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
  $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}/PSOPT>)

target_compile_features(psopt PUBLIC cxx_std_17)

target_compile_definitions(psopt PUBLIC
  PSOPT_RELEASE_STRING="${PSOPT_RELEASE_STRING}"
  PSOPT_BUILD_DATE="${PSOPT_BUILD_DATE}"
  HAVE_CSTDDEF
  PSOPT_DEFAULT_LINEAR_SOLVER="${PSOPT_DEFAULT_LINEAR_SOLVER}"
  $<$<BOOL:${PSOPT_WITH_IPOPT}>:USE_IPOPT>
  $<$<BOOL:${PSOPT_WITH_SNOPT}>:USE_SNOPT>
  $<$<BOOL:${PSOPT_WITH_CASADI}>:USE_CASADI>
  $<$<BOOL:${PSOPT_WITH_SCIP}>:USE_SCIP>
  $<$<BOOL:${PSOPT_WITH_HSL}>:USE_HSL>
  $<$<BOOL:${PSOPT_WITH_INT64}>:PSOPT_INT64>
  $<$<BOOL:${HEADLESS}>:HEADLESS>)

# All three appear in the public header psopt.h (Eigen/Dense, adolc/*, IpTNLP),
# so they are PUBLIC usage requirements of consumers that include it.
# ADOL-C: link its namespaced CMake-config target directly when available
# (superbuild), else the plain "adolc" from the pkg-config Find module (MacPorts).
if(TARGET adolc::adolc)
  target_link_libraries(psopt PUBLIC adolc::adolc)
else()
  target_link_libraries(psopt PUBLIC adolc)
endif()
target_link_libraries(psopt
  PUBLIC  Eigen3::Eigen
  $<$<BOOL:${PSOPT_WITH_IPOPT}>:IPOPT::ipopt>
  $<$<BOOL:${PSOPT_WITH_CASADI}>:casadi::casadi>
  $<$<BOOL:${PSOPT_WITH_SCIP}>:libscip>
  $<$<BOOL:${PSOPT_WITH_OPENMP}>:OpenMP::OpenMP_CXX>
  $<$<BOOL:${PSOPT_WITH_MPI}>:MPI::MPI_CXX>)
if(BLAS_FOUND)
  target_link_libraries(psopt PUBLIC ${BLAS_LIBRARIES})
endif()
if(LAPACK_FOUND)
  target_link_libraries(psopt PUBLIC ${LAPACK_LIBRARIES})
endif()

# ────────────────────────────────────────────────────────────────
# EXAMPLES / TESTS
# ────────────────────────────────────────────────────────────────
if(BUILD_TESTS OR BUILD_EXAMPLES)
  enable_testing()   # must precede any add_test() in the subdirectories below
endif()
if(BUILD_EXAMPLES)
  add_subdirectory(examples)
endif()
if(BUILD_TESTS)
  add_subdirectory(tests)
endif()

# ────────────────────────────────────────────────────────────────
# INSTALL & PACKAGE (config now re-finds deps via find_dependency)
# ────────────────────────────────────────────────────────────────
include(CMakePackageConfigHelpers)
install(TARGETS psopt EXPORT PSOPTTargets
  LIBRARY  DESTINATION ${CMAKE_INSTALL_LIBDIR}
  ARCHIVE  DESTINATION ${CMAKE_INSTALL_LIBDIR}
  RUNTIME  DESTINATION ${CMAKE_INSTALL_BINDIR})
install(DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/include/" DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/PSOPT)

set(PSOPT_CMAKE_DIR ${CMAKE_INSTALL_LIBDIR}/cmake/PSOPT)
install(EXPORT PSOPTTargets NAMESPACE PSOPT:: DESTINATION ${PSOPT_CMAKE_DIR})
# bundle our Find modules so find_dependency(IPOPT/adolc) works downstream
install(FILES cmake/FindIPOPT.cmake cmake/Findadolc.cmake cmake/FindColPack.cmake
        DESTINATION ${PSOPT_CMAKE_DIR})

configure_package_config_file(cmake/PSOPTConfig.cmake.in
  "${CMAKE_CURRENT_BINARY_DIR}/PSOPTConfig.cmake"
  INSTALL_DESTINATION ${PSOPT_CMAKE_DIR})
write_basic_package_version_file(
  "${CMAKE_CURRENT_BINARY_DIR}/PSOPTConfigVersion.cmake"
  VERSION ${PSOPT_VERSION} COMPATIBILITY SameMajorVersion)
install(FILES
  "${CMAKE_CURRENT_BINARY_DIR}/PSOPTConfig.cmake"
  "${CMAKE_CURRENT_BINARY_DIR}/PSOPTConfigVersion.cmake"
  DESTINATION ${PSOPT_CMAKE_DIR})

message(STATUS "PSOPT ${PSOPT_VERSION}  solver stack: "
  "MUMPS=${PSOPT_WITH_MUMPS} HSL=${PSOPT_WITH_HSL} CASADI=${PSOPT_WITH_CASADI} SNOPT=${PSOPT_WITH_SNOPT} "
  "default_linear_solver=${PSOPT_DEFAULT_LINEAR_SOLVER}")
message(STATUS "PSOPT parallelism: OPENMP=${PSOPT_WITH_OPENMP} MPI=${PSOPT_WITH_MPI}  "
  "INT64=${PSOPT_WITH_INT64}  BLAS=${BLAS_FOUND} LAPACK=${LAPACK_FOUND}")
