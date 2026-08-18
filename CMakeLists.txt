cmake_minimum_required(VERSION 3.12)

# ────────────────────────────────────────────────────────────────
# VERSIONING
# ────────────────────────────────────────────────────────────────

# Get the latest tag for use as version (e.g. 2025-04-21)   
execute_process(
  COMMAND git describe --tags --abbrev=0
  WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
  OUTPUT_VARIABLE PSOPT_TAG_RAW
  ERROR_QUIET
  OUTPUT_STRIP_TRAILING_WHITESPACE
)

# Convert to CMake-friendly version format: 2025.04.21
if(PSOPT_TAG_RAW STREQUAL "")
  string(TIMESTAMP PSOPT_VERSION "20%y.%m.%d")
else()
  string(REPLACE "-" "." PSOPT_VERSION "${PSOPT_TAG_RAW}")
endif()

# Full Git description: tag + commits ahead + short hash
execute_process(
  COMMAND git describe --tags --always
  WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
  OUTPUT_VARIABLE PSOPT_RELEASE_STRING
  OUTPUT_STRIP_TRAILING_WHITESPACE
)

# Get build date
string(TIMESTAMP PSOPT_BUILD_DATE "%Y-%m-%d %H:%M:%S")

# Debug prints
message(STATUS "PSOPT Git Version String: ${PSOPT_RELEASE_STRING}")
message(STATUS "PSOPT Version from Tag or Date: ${PSOPT_VERSION}")
message(STATUS "PSOPT Build Date: ${PSOPT_BUILD_DATE}")

# Define macros for use in C++ code
add_definitions(
  -DPSOPT_RELEASE_STRING="${PSOPT_RELEASE_STRING}"
  -DPSOPT_BUILD_DATE="${PSOPT_BUILD_DATE}"
)

# Set project metadata
project(PSOPT VERSION ${PSOPT_VERSION} LANGUAGES CXX)

# ────────────────────────────────────────────────────────────────
# OPTIONS & DEPENDENCIES
# ────────────────────────────────────────────────────────────────

option(BUILD_EXAMPLES "Build examples from the example subdirectory." OFF)
option(BUILD_TESTS "Build tests for PSOPT." OFF)
option(HEADLESS "Headless option (Gnuplot does not spawn windows with plots, but generates PDFs)" OFF)
option(PSOPT_STRICT_WARNINGS "Build the PSOPT library with -Wall -Wextra -Werror (for CI). Default OFF so normal builds are unaffected." OFF)
option(PSOPT_PYTHON "Build the Python binding (python/). Builds libPSOPT position-independent, which the extension requires." OFF)
option(PSOPT_ALLOW_ENV_OVERRIDES "Let algorithm settings be overridden from the environment (PSOPT_NLP_METHOD and friends). For benchmarking one binary many ways; every override is announced on stdout. OFF by default." OFF)
include(GNUInstallDirs)
option(WITH_SQP "Build PSOPT's own SQP solver (algorithm.nlp_method = \"SQP\"). Requires at least one QP backend, e.g. WITH_GALAHAD." OFF)
option(WITH_PROXQP "Add ProxQP as a QP backend for the SQP solver (algorithm.qp_solver = \"ProxQP\"). Requires proxsuite; implies WITH_SQP." OFF)
option(WITH_QPALM "Add QPALM as a QP backend for the SQP solver (algorithm.qp_solver = \"QPALM\"). Requires QPALM; implies WITH_SQP." OFF)
option(WITH_OSQP "Add OSQP as a QP backend for the SQP solver (algorithm.qp_solver = \"OSQP\"). Requires OSQP; implies WITH_SQP." OFF)
option(WITH_GALAHAD "Add GALAHAD's QPA as a QP backend for the SQP solver (algorithm.qp_solver = \"GALAHAD\"). Requires GALAHAD; implies WITH_SQP." OFF)

# Automatic-differentiation backend. CPPAD (in-memory tapes) is the production backend; it
# replaced ADOL-C, whose multi-GB on-disk tapes corrupt/exhaust disk at high node counts.
# AUTODIFF and XAD are forward-mode cross-check backends (diagnostic only).
set(PSOPT_AD_BACKEND "CPPAD" CACHE STRING "AD backend: CPPAD (default), AUTODIFF, or XAD")
set_property(CACHE PSOPT_AD_BACKEND PROPERTY STRINGS CPPAD AUTODIFF XAD)

list(APPEND CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/cmake)
find_package(Eigen3 REQUIRED NO_MODULE)
find_package(PkgConfig REQUIRED)
pkg_check_modules(ipopt REQUIRED IMPORTED_TARGET ipopt)

if(PSOPT_AD_BACKEND STREQUAL "CPPAD")
  find_path(CPPAD_INCLUDE_DIR cppad/cppad.hpp HINTS ENV CPPAD_DIR PATH_SUFFIXES include)
  find_library(CPPAD_LIBRARY NAMES cppad_lib HINTS ENV CPPAD_DIR PATH_SUFFIXES lib lib64 cppad_lib build/cppad_lib)
  if(NOT CPPAD_INCLUDE_DIR OR NOT CPPAD_LIBRARY)
    message(FATAL_ERROR "PSOPT_AD_BACKEND=CPPAD but CppAD not found. "
      "Set -DCPPAD_INCLUDE_DIR=<dir with cppad/cppad.hpp> and -DCPPAD_LIBRARY=<libcppad_lib>, or CPPAD_DIR.")
  endif()
  message(STATUS "PSOPT AD backend: CppAD  (${CPPAD_LIBRARY})")
elseif(PSOPT_AD_BACKEND STREQUAL "AUTODIFF")
  find_path(AUTODIFF_INCLUDE_DIR autodiff/forward/dual.hpp HINTS ENV AUTODIFF_DIR PATH_SUFFIXES include)
  if(NOT AUTODIFF_INCLUDE_DIR)
    message(FATAL_ERROR "PSOPT_AD_BACKEND=AUTODIFF but autodiff not found. "
      "Set -DAUTODIFF_INCLUDE_DIR=<dir with autodiff/forward/dual.hpp>, or AUTODIFF_DIR.")
  endif()
  message(STATUS "PSOPT AD backend: autodiff (forward)  (${AUTODIFF_INCLUDE_DIR})")
elseif(PSOPT_AD_BACKEND STREQUAL "XAD")
  find_path(XAD_INCLUDE_DIR XAD/XAD.hpp HINTS ENV XAD_DIR PATH_SUFFIXES include)
  find_library(XAD_LIBRARY NAMES xad HINTS ENV XAD_DIR PATH_SUFFIXES lib lib64)
  if(NOT XAD_INCLUDE_DIR OR NOT XAD_LIBRARY)
    message(FATAL_ERROR "PSOPT_AD_BACKEND=XAD but XAD not found. "
      "Set -DXAD_INCLUDE_DIR=<dir with XAD/XAD.hpp> and -DXAD_LIBRARY=<libxad>, or XAD_DIR.")
  endif()
  message(STATUS "PSOPT AD backend: XAD (forward)  (${XAD_LIBRARY})")
else()
  message(FATAL_ERROR "Unknown PSOPT_AD_BACKEND='${PSOPT_AD_BACKEND}' (expected CPPAD, AUTODIFF, or XAD).")
endif()

# ────────────────────────────────────────────────────────────────
# BUILD TARGETS
# ────────────────────────────────────────────────────────────────

file(GLOB SRC src/*.cxx)
file(GLOB HEADER include/*.h)

add_library(${PROJECT_NAME}
  ${SRC} ${HEADER}
)

target_include_directories(${PROJECT_NAME}
  PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
    $<INSTALL_INTERFACE:include/${PROJECT_NAME}>
)

target_link_libraries(${PROJECT_NAME}
  PRIVATE Eigen3::Eigen
  PUBLIC PkgConfig::ipopt
)

# Select the AD backend at compile time. PUBLIC so examples/user code linking PSOPT inherit the
# define, CppAD headers, and the in-memory-tape support library.
if(PSOPT_AD_BACKEND STREQUAL "CPPAD")
  target_compile_definitions(${PROJECT_NAME} PUBLIC PSOPT_AD_BACKEND=2)
  target_include_directories(${PROJECT_NAME} PUBLIC ${CPPAD_INCLUDE_DIR})
  target_link_libraries(${PROJECT_NAME} PUBLIC ${CPPAD_LIBRARY})
elseif(PSOPT_AD_BACKEND STREQUAL "AUTODIFF")
  target_compile_definitions(${PROJECT_NAME} PUBLIC PSOPT_AD_BACKEND=3)
  target_include_directories(${PROJECT_NAME} PUBLIC ${AUTODIFF_INCLUDE_DIR})
elseif(PSOPT_AD_BACKEND STREQUAL "XAD")
  target_compile_definitions(${PROJECT_NAME} PUBLIC PSOPT_AD_BACKEND=4)
  target_include_directories(${PROJECT_NAME} PUBLIC ${XAD_INCLUDE_DIR})
  target_link_libraries(${PROJECT_NAME} PUBLIC ${XAD_LIBRARY})
endif()

if(PSOPT_ALLOW_ENV_OVERRIDES)
  target_compile_definitions(${PROJECT_NAME} PUBLIC PSOPT_ALLOW_ENV_OVERRIDES)
  message(STATUS "PSOPT: algorithm settings may be overridden from the environment "
                 "(benchmarking build; every override is announced)")
endif()

# PSOPT's own SQP solver. The quadratic programming subproblem goes to one of the
# backends below, each of which is sparse. The solver is opt-in: with WITH_SQP=OFF the
# file compiles to a stub that reports the option was not built, and PSOPT gains no new
# dependency.
if(WITH_SQP)
  if(NOT (WITH_PROXQP OR WITH_QPALM OR WITH_OSQP OR WITH_GALAHAD))
    message(FATAL_ERROR
      "WITH_SQP=ON but no QP backend was selected. The SQP sends every subproblem to a "
      "backend and has none of its own, so at least one of WITH_GALAHAD, WITH_PROXQP, "
      "WITH_QPALM or WITH_OSQP must also be ON. WITH_GALAHAD is the one to choose: it is "
      "sparse, BSD-3 licensed, and the configuration the solver has been measured "
      "against. scripts/build_galahad.sh will install it.")
  endif()

  # Where to look for the SQP's dependencies beyond the compiler's defaults. MacPorts
  # installs under /opt/local and Homebrew under /opt/homebrew on Apple Silicon; neither
  # is searched by find_library out of the box, so a Mac that has the dependency still
  # reports it missing. The per-dependency *_DIR variables continue to take precedence.
  set(PSOPT_SEARCH_PREFIXES /usr/local /usr)
  if(APPLE)
    list(INSERT PSOPT_SEARCH_PREFIXES 0 /opt/local /opt/homebrew)
  endif()
  set(PSOPT_EXTRA_PREFIXES_INCLUDE "")
  set(PSOPT_EXTRA_PREFIXES_LIB "")
  foreach(_p ${PSOPT_SEARCH_PREFIXES})
    list(APPEND PSOPT_EXTRA_PREFIXES_INCLUDE ${_p}/include)
    list(APPEND PSOPT_EXTRA_PREFIXES_LIB     ${_p}/lib)
  endforeach()

  # The QP backend plugins link their own linear algebra, and GALAHAD's needs these.
  find_package(LAPACK REQUIRED)
  find_package(BLAS REQUIRED)

  find_program(PSOPT_NM_EXECUTABLE nm)

  # MUMPS supplies the inertia of the KKT matrix, which the SQP needs to decide whether
  # the Hessian has to be modified -- Betts (2.41). IPOPT's default linear solver is the
  # same sequential MUMPS, so nothing new is being installed here; both the header and
  # the library have to be found all the same.
  #
  # The header is found here rather than left to the compiler's default include path,
  # which is where it happens to sit on Debian and Ubuntu but not on macOS or on any
  # build where IPOPT and MUMPS were installed under their own prefix.
  # PATH_SUFFIXES matters as much as the hints. coinbrew -- which is what README.md tells
  # macOS users to build IPOPT and MUMPS with -- installs the header as
  # <prefix>/include/coin-or/mumps/dmumps_c.h, so a bare <prefix>/include search misses
  # it. Listing the suffixes also means CMAKE_PREFIX_PATH alone is enough, which is
  # already set to that prefix in the macOS build instructions.
  find_path(MUMPS_INCLUDE_DIR dmumps_c.h
            HINTS ${MUMPS_DIR}/include ${MUMPS_DIR}
                  ${IPOPT_DIR}/include
                  ${ipopt_INCLUDE_DIRS}
                  ${PSOPT_EXTRA_PREFIXES_INCLUDE}
            PATH_SUFFIXES coin-or/mumps coin/mumps mumps coin-or coin)
  if(NOT MUMPS_INCLUDE_DIR)
    message(FATAL_ERROR
      "WITH_SQP=ON but dmumps_c.h was not found. PSOPT's SQP reads the inertia of the "
      "KKT matrix from MUMPS, which IPOPT already provides as its default linear solver, "
      "so this is almost always a matter of pointing at an existing installation rather "
      "than installing anything new. Set MUMPS_DIR to the prefix containing "
      "include/dmumps_c.h (on macOS, the prefix you built MUMPS into when following the "
      "IPOPT instructions in README.md).")
  endif()
  message(STATUS "PSOPT SQP: MUMPS headers at ${MUMPS_INCLUDE_DIR}")

  # And the library, which is a separate question from the header and was for a while
  # answered by accident. On Debian and Ubuntu `pkg-config --libs ipopt` lists
  # -ldmumps_seq itself, so linking PkgConfig::ipopt resolves dmumps_c and nothing here
  # had to ask for it. An IPOPT built with coinbrew does not list it -- the dependency is
  # recorded inside libipopt instead -- and macOS's two-level namespace will not resolve a
  # symbol through an indirect dylib the way ELF will. The result is a link that succeeds
  # everywhere it was tested and fails on a Mac with three undefined _dmumps_c.
  #
  # Name order matters. coinmumps is what coinbrew's ThirdParty/Mumps produces;
  # dmumps_seq is Debian's sequential build and is preferred over plain dmumps, which on
  # the same system is the MPI build and would drag in an MPI runtime that IPOPT is not
  # using.
  find_library(MUMPS_LIBRARY NAMES coinmumps dmumps_seq dmumps
               HINTS ${MUMPS_DIR}/lib ${MUMPS_DIR}
                     ${IPOPT_DIR}/lib
                     ${ipopt_LIBRARY_DIRS}
                     ${PSOPT_EXTRA_PREFIXES_LIB})
  if(MUMPS_LIBRARY)
    message(STATUS "PSOPT SQP: MUMPS library at ${MUMPS_LIBRARY}")
  else()
    # Not fatal: if IPOPT's own link line carries MUMPS, as it does on Debian and Ubuntu,
    # the symbol resolves without any help from us.
    string(FIND "${ipopt_LIBRARIES}" "mumps" _psopt_ipopt_has_mumps)
    if(_psopt_ipopt_has_mumps GREATER -1)
      message(STATUS "PSOPT SQP: MUMPS comes from IPOPT's own link line")
    else()
      message(WARNING
        "WITH_SQP=ON but no MUMPS library was found, and IPOPT's pkg-config link line "
        "does not mention one either. If the build fails with an undefined dmumps_c, set "
        "MUMPS_LIBRARY to the library holding it -- libcoinmumps if you built IPOPT with "
        "coinbrew, libdmumps_seq on Debian and Ubuntu.")
    endif()
  endif()

  target_compile_definitions(${PROJECT_NAME} PUBLIC USE_SQP)
  target_include_directories(${PROJECT_NAME} PUBLIC ${MUMPS_INCLUDE_DIR})
  target_link_libraries(${PROJECT_NAME} PUBLIC ${LAPACK_LIBRARIES} ${BLAS_LIBRARIES})
  if(MUMPS_LIBRARY)
    target_link_libraries(${PROJECT_NAME} PUBLIC ${MUMPS_LIBRARY})
  endif()

  # ---- QP backend plugins ----------------------------------------------------
  # Each backend is a shared object that links its own QP library and its own linear
  # algebra, exports the three entry points of psopt_qp_plugin.h and nothing else, and
  # is loaded at run time with RTLD_LOCAL. The reason is in that header: every one of
  # these libraries carries an ordering, the orderings share C symbol names, and linked
  # together they corrupt each other. Kept in separate plugins they coexist, which is
  # what lets a build offer more than one.
  #
  # -fvisibility=hidden hides this project's own symbols; --exclude-libs,ALL hides the
  # ones that arrive from the static libraries the plugin links. Both are needed, and
  # the CTest check below verifies the result rather than trusting it.
  #
  # Neither flag exists on macOS: Apple's ld64 has no --exclude-libs and does not
  # understand -z. It has something better for this purpose -- -exported_symbols_list
  # takes a file naming everything the module may export, which is a whitelist rather
  # than two overlapping suppressions, and so cannot be defeated by a static library
  # arriving through a path nobody thought about. The symbol names in that file carry
  # Mach-O's leading underscore.
  if(APPLE)
    set(PSOPT_QP_PLUGIN_EXPORTS "${CMAKE_BINARY_DIR}/psopt_qp_plugin.exp")
    file(WRITE "${PSOPT_QP_PLUGIN_EXPORTS}"
         "_psopt_qp_abi_version\n_psopt_qp_name\n_psopt_qp_solve\n_psopt_qp_environment_ok\n")
  endif()
  set(PSOPT_QP_PLUGIN_DIR "${CMAKE_INSTALL_PREFIX}/${CMAKE_INSTALL_LIBDIR}/psopt")
  target_compile_definitions(${PROJECT_NAME}
      PRIVATE PSOPT_QP_PLUGIN_INSTALL_DIR="${PSOPT_QP_PLUGIN_DIR}"
              PSOPT_QP_PLUGIN_BUILD_DIR="${CMAKE_BINARY_DIR}/qp_plugins")
  target_link_libraries(${PROJECT_NAME} PUBLIC ${CMAKE_DL_LIBS})

  function(psopt_add_qp_plugin name source)
    add_library(psopt_qp_${name} MODULE ${source})
    set_target_properties(psopt_qp_${name} PROPERTIES
        PREFIX "lib"
        CXX_VISIBILITY_PRESET hidden
        VISIBILITY_INLINES_HIDDEN ON
        LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/qp_plugins)
    target_include_directories(psopt_qp_${name} PRIVATE ${CMAKE_SOURCE_DIR}/include)
    target_include_directories(psopt_qp_${name} PRIVATE ${MUMPS_INCLUDE_DIR})
    target_link_libraries(psopt_qp_${name} PRIVATE Eigen3::Eigen ${ARGN})
    if(APPLE)
      target_link_options(psopt_qp_${name} PRIVATE
          "LINKER:-exported_symbols_list,${PSOPT_QP_PLUGIN_EXPORTS}"
          "LINKER:-undefined,error")
      set_property(TARGET psopt_qp_${name} APPEND PROPERTY
          LINK_DEPENDS "${PSOPT_QP_PLUGIN_EXPORTS}")
    else()
      target_link_options(psopt_qp_${name} PRIVATE
          "LINKER:--exclude-libs,ALL" "LINKER:-z,defs")
    endif()
    install(TARGETS psopt_qp_${name} LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}/psopt)
  endfunction()

  if(WITH_PROXQP)
    find_path(PROXSUITE_INCLUDE_DIR proxsuite/proxqp/sparse/sparse.hpp
              HINTS ${PROXSUITE_DIR}/include ${PSOPT_EXTRA_PREFIXES_INCLUDE})
    if(NOT PROXSUITE_INCLUDE_DIR)
      message(FATAL_ERROR "WITH_PROXQP=ON but proxsuite was not found. Set PROXSUITE_DIR to its install prefix.")
    endif()
    message(STATUS "PSOPT SQP: ProxQP backend plugin enabled (proxsuite at ${PROXSUITE_INCLUDE_DIR})")
    psopt_add_qp_plugin(proxqp ${CMAKE_SOURCE_DIR}/src/qp_backends/proxqp_plugin.cxx)
    target_compile_definitions(${PROJECT_NAME} PUBLIC USE_PROXQP)
    target_include_directories(psopt_qp_proxqp PRIVATE ${PROXSUITE_INCLUDE_DIR})
  endif()

  if(WITH_QPALM)
    find_package(QPALM QUIET)
    if(NOT QPALM_FOUND)
      message(FATAL_ERROR "WITH_QPALM=ON but QPALM was not found. Set QPALM_DIR to the directory containing QPALMConfig.cmake.")
    endif()
    # QPALM's headers choose their index width from LADEL_64BIT_INDICES and it must
    # match the library's own build, or the C++ interface fails to link with an
    # undefined reference to eigen_to_ladel_copy taking a matrix of the other width.
    find_library(LADEL_LIBRARY NAMES ladel HINTS ${QPALM_DIR}/../.. ${PSOPT_EXTRA_PREFIXES_LIB})
    find_library(LADEL_AMD_LIBRARY NAMES ladel_amd HINTS ${QPALM_DIR}/../.. ${PSOPT_EXTRA_PREFIXES_LIB})
    if(NOT LADEL_LIBRARY)
      message(FATAL_ERROR "WITH_QPALM=ON but LADEL, which QPALM factorises with, was not found.")
    endif()
    message(STATUS "PSOPT SQP: QPALM backend plugin enabled")
    psopt_add_qp_plugin(qpalm ${CMAKE_SOURCE_DIR}/src/qp_backends/qpalm_plugin.cxx
                        qpalm_cxx qpalm ${LADEL_LIBRARY} ${LADEL_AMD_LIBRARY})
    target_compile_definitions(${PROJECT_NAME} PUBLIC USE_QPALM)
    target_compile_definitions(psopt_qp_qpalm PRIVATE LADEL_64BIT_INDICES)
  endif()

  if(WITH_OSQP)
    find_package(osqp QUIET)
    if(NOT osqp_FOUND)
      message(FATAL_ERROR "WITH_OSQP=ON but OSQP was not found. Set osqp_DIR to the directory containing osqpConfig.cmake.")
    endif()
    message(STATUS "PSOPT SQP: OSQP backend plugin enabled")
    psopt_add_qp_plugin(osqp ${CMAKE_SOURCE_DIR}/src/qp_backends/osqp_plugin.cxx osqp::osqpstatic)
    target_compile_definitions(${PROJECT_NAME} PUBLIC USE_OSQP)
  endif()

  if(WITH_GALAHAD)
    find_path(GALAHAD_INCLUDE_DIR galahad_qpa.h
              HINTS ${GALAHAD_DIR}/include ${PSOPT_EXTRA_PREFIXES_INCLUDE})
    find_library(GALAHAD_LIBRARY NAMES galahad_double galahad
              HINTS ${GALAHAD_DIR}/lib ${PSOPT_EXTRA_PREFIXES_LIB})
    if(NOT GALAHAD_INCLUDE_DIR OR NOT GALAHAD_LIBRARY)
      message(FATAL_ERROR "WITH_GALAHAD=ON but GALAHAD was not found. Set GALAHAD_DIR to its install prefix.")
    endif()

    # GALAHAD is Fortran, so its plugin needs the Fortran runtime. Naming "gfortran" and
    # hoping the linker finds it works on a Linux distribution that installs libgfortran
    # into a default search path and nowhere else: on macOS it lives beside whichever
    # gcc MacPorts or Homebrew put there, under a versioned directory the C++ linker has
    # never heard of. Asking CMake for the Fortran compiler's own implicit link line
    # gets the right library and the right -L on both.
    enable_language(Fortran)
    set(PSOPT_FORTRAN_RUNTIME "")
    foreach(_d ${CMAKE_Fortran_IMPLICIT_LINK_DIRECTORIES})
      list(APPEND PSOPT_FORTRAN_RUNTIME "-L${_d}")
    endforeach()
    if(CMAKE_Fortran_IMPLICIT_LINK_LIBRARIES)
      foreach(_l ${CMAKE_Fortran_IMPLICIT_LINK_LIBRARIES})
        if(IS_ABSOLUTE "${_l}")
          list(APPEND PSOPT_FORTRAN_RUNTIME "${_l}")
        else()
          list(APPEND PSOPT_FORTRAN_RUNTIME "-l${_l}")
        endif()
      endforeach()
    else()
      list(APPEND PSOPT_FORTRAN_RUNTIME "-lgfortran")
    endif()
    message(STATUS "PSOPT SQP: Fortran runtime for the GALAHAD plugin: ${PSOPT_FORTRAN_RUNTIME}")
    message(STATUS "PSOPT SQP: GALAHAD backend plugin enabled (${GALAHAD_LIBRARY})")
    enable_language(Fortran)
    psopt_add_qp_plugin(galahad ${CMAKE_SOURCE_DIR}/src/qp_backends/galahad_plugin.cxx
                        ${GALAHAD_LIBRARY} ${PSOPT_FORTRAN_RUNTIME}
                        ${LAPACK_LIBRARIES} ${BLAS_LIBRARIES})
    target_compile_definitions(${PROJECT_NAME} PUBLIC USE_GALAHAD)
    target_include_directories(psopt_qp_galahad PRIVATE ${GALAHAD_INCLUDE_DIR})
  endif()

  # The check that the isolation actually holds. A plugin that exports anything beyond
  # its three entry points has failed at the one job it was given, and the failure would
  # otherwise show up as a crash somewhere else entirely.
  if(WITH_PROXQP OR WITH_QPALM OR WITH_OSQP OR WITH_GALAHAD)
    find_program(PSOPT_NM_EXECUTABLE nm)
    if(PSOPT_NM_EXECUTABLE)
      if(APPLE)
        set(PSOPT_PLUGIN_FORMAT macho)
      else()
        set(PSOPT_PLUGIN_FORMAT elf)
      endif()
      add_test(NAME qp_plugins_export_nothing_else
               COMMAND ${CMAKE_COMMAND}
                       -DPLUGIN_DIR=${CMAKE_BINARY_DIR}/qp_plugins
                       -DNM=${PSOPT_NM_EXECUTABLE}
                       -DPLUGIN_FORMAT=${PSOPT_PLUGIN_FORMAT}
                       -P ${CMAKE_SOURCE_DIR}/cmake/CheckPluginSymbols.cmake)
    endif()
  endif()
elseif(WITH_PROXQP OR WITH_QPALM OR WITH_OSQP OR WITH_GALAHAD)
  message(FATAL_ERROR "The QP backend options require WITH_SQP=ON: they are backends of PSOPT's SQP solver.")
endif()

add_definitions(-DPSOPT_RELEASE_STRING="${PSOPT_RELEASE_STRING}" -DHAVE_CSTDDEF)
target_compile_features(${PROJECT_NAME} PUBLIC cxx_std_17)

# Opt-in strict warnings for CI. Applied only to the PSOPT library target so
# dependency headers and user/example code are not affected.
#  - unused-parameter is silenced: the fixed dae/cost/events/linkage callback
#    signatures force parameters that individual problems may not use.
#  - maybe-uninitialized is kept visible but non-fatal: it is prone to
#    optimizer/compiler-version false positives (e.g. inside Eigen).
if(PSOPT_STRICT_WARNINGS)
  if(CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
    target_compile_options(${PROJECT_NAME} PRIVATE
      -Wall -Wextra -Werror
      -Wno-unused-parameter
      -Wno-error=maybe-uninitialized
    )
  endif()
endif()

if(${HEADLESS})
    add_definitions(-DPSOPT_RELEASE_STRING="${PSOPT_RELEASE_STRING}" -DHAVE_CSTDDEF -DHEADLESS)
endif()

if(${BUILD_EXAMPLES})
    add_subdirectory(examples/)
endif()

if(${BUILD_TESTS})
    enable_testing()
    add_subdirectory(tests/)
endif()

# The Python extension is a shared module that whole-archives libPSOPT, so the
# library must be position-independent; a default static build is not, and the
# link fails with a relocation error against the thread-local print_level. This
# is the in-tree route documented in python/README.md.
if(${PSOPT_PYTHON})
    set_property(TARGET PSOPT PROPERTY POSITION_INDEPENDENT_CODE ON)
    add_subdirectory(python/)
endif()

# ────────────────────────────────────────────────────────────────
# INSTALL & PACKAGE
# ────────────────────────────────────────────────────────────────

set(PACKAGE_FIND_NAME "${PROJECT_NAME}")
install(TARGETS ${PROJECT_NAME} EXPORT ${PROJECT_NAME} DESTINATION "lib/${PROJECT_NAME}")
install(DIRECTORY "${CMAKE_SOURCE_DIR}/include/" DESTINATION "include/${PROJECT_NAME}")
install(FILES psopt-config.cmake DESTINATION "share/cmake/${PROJECT_NAME}/")
install(EXPORT ${PROJECT_NAME} DESTINATION "share/cmake/${PROJECT_NAME}")
