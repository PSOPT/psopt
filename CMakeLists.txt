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
option(WITH_SNOPT_INTERFACE "Build interface for the SNOPT optimizer." OFF)
option(BUILD_TESTS "Build tests for PSOPT." OFF)
option(HEADLESS "Headless option (Gnuplot does not spawn windows with plots, but generates PDFs)" OFF)
option(PSOPT_STRICT_WARNINGS "Build the PSOPT library with -Wall -Wextra -Werror (for CI). Default OFF so normal builds are unaffected." OFF)

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

if(${WITH_SNOPT_INTERFACE})
    add_subdirectory(snopt-interface/)
endif()

file(GLOB SRC src/*.cxx)
file(GLOB HEADER include/*.h)

add_library(${PROJECT_NAME}
  $<$<BOOL:$<TARGET_EXISTS:PSOPT_SNOPT_interface>>:$<TARGET_OBJECTS:PSOPT_SNOPT_interface>>
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

if(${WITH_SNOPT_INTERFACE})
    add_definitions(-DPSOPT_RELEASE_STRING="${PSOPT_RELEASE_STRING}" -DHAVE_CSTDDEF -DUSE_SNOPT)
endif()

if(${HEADLESS})
    add_definitions(-DPSOPT_RELEASE_STRING="${PSOPT_RELEASE_STRING}" -DHAVE_CSTDDEF -DHEADLESS)
endif()

if(${BUILD_EXAMPLES})
    add_subdirectory(examples/)
endif()

if(${BUILD_TESTS})
    add_subdirectory(tests/)
endif()

# ────────────────────────────────────────────────────────────────
# INSTALL & PACKAGE
# ────────────────────────────────────────────────────────────────

set(PACKAGE_FIND_NAME "${PROJECT_NAME}")
install(TARGETS ${PROJECT_NAME} EXPORT ${PROJECT_NAME} DESTINATION "lib/${PROJECT_NAME}")
install(DIRECTORY "${CMAKE_SOURCE_DIR}/include/" DESTINATION "include/${PROJECT_NAME}")
install(FILES psopt-config.cmake DESTINATION "share/cmake/${PROJECT_NAME}/")
install(EXPORT ${PROJECT_NAME} DESTINATION "share/cmake/${PROJECT_NAME}")
