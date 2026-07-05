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

project(PSOPT VERSION ${PSOPT_VERSION} LANGUAGES C CXX)

list(APPEND CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/cmake)
include(GNUInstallDirs)

# ────────────────────────────────────────────────────────────────
# OPTIONS — choose the solver stack; build any missing dep ourselves
# ────────────────────────────────────────────────────────────────
option(PSOPT_WITH_MUMPS   "IPOPT + MUMPS linear solver"                        ON)
option(PSOPT_WITH_HSL     "IPOPT + HSL (MA27/57/97); needs Coin-HSL source"    OFF)
option(PSOPT_WITH_CASADI  "Full CasADi NLP backend (CasADi -> IPOPT[+HSL])"    OFF)
option(PSOPT_WITH_SNOPT   "SNOPT NLP backend"                                  OFF)
option(PSOPT_WITH_OPENMP  "Enable OpenMP (threaded BLAS/HSL/ColPack, parallel MUMPS)" OFF)
option(PSOPT_WITH_MPI     "Enable MPI (distributed-memory MUMPS; see macOS-26 caveat)" OFF)
option(PSOPT_SUPERBUILD   "Download and build missing dependencies from source" OFF)
option(PSOPT_WITH_INT64   "64-bit integer addressing (ILP64) for very large NLPs" OFF)
option(BUILD_EXAMPLES     "Build examples from the example subdirectory."      ON)
option(BUILD_TESTS        "Build tests for PSOPT."                             ON)
option(HEADLESS           "Gnuplot writes PDFs instead of spawning windows."   OFF)
set(PSOPT_DEFAULT_LINEAR_SOLVER "mumps" CACHE STRING
    "Default IPOPT linear_solver (mumps, ma27, ma57, ma86, ma97, pardiso)")
set(COINHSL_SOURCE_DIR "" CACHE PATH "Path to a local Coin-HSL source tree (skips download)")
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
  $<$<BOOL:${PSOPT_WITH_SNOPT}>:USE_SNOPT>
  $<$<BOOL:${PSOPT_WITH_CASADI}>:USE_CASADI>
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
  PUBLIC  Eigen3::Eigen IPOPT::ipopt
  $<$<BOOL:${PSOPT_WITH_CASADI}>:casadi::casadi>
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
