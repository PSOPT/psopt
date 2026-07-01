# ──────────────────────────────────────────────────────────────────────────
# PSOPT Python binding
#
# Builds the `_psopt` pybind11 extension and generates psopt/_toolchain.py from
# the compiler and flags CMake used, so the JIT codegen path needs no
# hand-written, machine-specific configuration.
#
# Two configure modes:
#   * In-tree   : pulled in by the top-level build with -DPSOPT_PYTHON=ON, where
#                 the PSOPT target already exists and is built position-independent.
#   * Standalone: configured on its own against a prebuilt PIC libPSOPT, e.g.
#                   cmake -S python -B build/python \
#                         -DPSOPT_ROOT=<src root> -DPSOPT_LIB=<.../libPSOPT.a> \
#                         -Dpybind11_DIR=$(python -m pybind11 --cmakedir)
#
# The extension whole-archives libPSOPT (so auto_link is retained) and exports
# its dynamic symbols (-rdynamic); psopt/__init__.py imports it with RTLD_GLOBAL
# so each dlopened math .so resolves those symbols at load time.  The extension
# is built with default symbol visibility for the same reason.
# ──────────────────────────────────────────────────────────────────────────
cmake_minimum_required(VERSION 3.15)

if(NOT TARGET PSOPT)
    project(psopt_python LANGUAGES CXX)
    set(PSOPT_ROOT "" CACHE PATH     "PSOPT source root (standalone build)")
    set(PSOPT_LIB  "" CACHE FILEPATH "Prebuilt PIC libPSOPT.a (standalone build)")
    if(NOT EXISTS "${PSOPT_LIB}" OR NOT IS_DIRECTORY "${PSOPT_ROOT}")
        message(FATAL_ERROR
            "Standalone build needs -DPSOPT_ROOT=<src> -DPSOPT_LIB=<.../libPSOPT.a>")
    endif()
    set(_PSOPT_ARCHIVE "${PSOPT_LIB}")
    set(_PSOPT_INCLUDE "${PSOPT_ROOT}/include")
else()
    # In-tree: link the freshly built archive and use the in-tree headers.
    set(_PSOPT_ARCHIVE "$<TARGET_FILE:PSOPT>")
    set(_PSOPT_INCLUDE "${CMAKE_SOURCE_DIR}/include")
endif()

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

find_package(pybind11 CONFIG REQUIRED)
find_package(Eigen3 REQUIRED)
find_package(PkgConfig REQUIRED)
pkg_check_modules(ipopt REQUIRED IMPORTED_TARGET ipopt)
find_library(CPPAD_LIBRARY NAMES cppad_lib PATHS /usr/local/lib REQUIRED)

# Resolve the Eigen include directory for both the build and the generated toolchain.
if(DEFINED EIGEN3_INCLUDE_DIR AND EIGEN3_INCLUDE_DIR)
    set(_EIGEN_INCLUDE "${EIGEN3_INCLUDE_DIR}")
else()
    get_target_property(_EIGEN_INCLUDE Eigen3::Eigen INTERFACE_INCLUDE_DIRECTORIES)
endif()

# ── the extension ─────────────────────────────────────────────────────────
pybind11_add_module(_psopt MODULE src/_psopt_driver.cpp)

# Default visibility: pybind11 defaults to hidden, which would strip the
# whole-archived PSOPT symbols from the dynamic table even with -rdynamic.
set_target_properties(_psopt PROPERTIES
    CXX_VISIBILITY_PRESET default
    VISIBILITY_INLINES_HIDDEN OFF
    LIBRARY_OUTPUT_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/psopt")

target_compile_definitions(_psopt PRIVATE
    HAVE_CSTDDEF
    PSOPT_RELEASE_STRING="${PSOPT_VERSION}"
    PSOPT_BUILD_DATE="${PSOPT_BUILD_DATE}"
    PSOPT_AD_NAME="CppAD")

target_include_directories(_psopt PRIVATE
    "${_PSOPT_INCLUDE}" "${_EIGEN_INCLUDE}" ${ipopt_INCLUDE_DIRS} /usr/local/include)

# whole-archive libPSOPT (retain auto_link) + export dynamic symbols (-rdynamic).
target_link_options(_psopt PRIVATE
    "LINKER:--whole-archive" "${_PSOPT_ARCHIVE}" "LINKER:--no-whole-archive"
    "LINKER:--export-dynamic")

target_link_libraries(_psopt PRIVATE
    PkgConfig::ipopt Eigen3::Eigen "${CPPAD_LIBRARY}" ${CMAKE_DL_LIBS} pthread)

if(TARGET PSOPT)
    add_dependencies(_psopt PSOPT)
endif()

# ── generate psopt/_toolchain.py ──────────────────────────────────────────
# Same compiler/flags/includes the JIT math .so needs; it links only CppAD.
set(PSOPT_PY_CXX      "${CMAKE_CXX_COMPILER}")
set(PSOPT_PY_CXXFLAGS "\"-O2\", \"-std=c++17\"")
set(PSOPT_PY_INCLUDES "\"-I${_PSOPT_INCLUDE}\", \"-I${_EIGEN_INCLUDE}\"")
foreach(_d ${ipopt_INCLUDE_DIRS})
    string(APPEND PSOPT_PY_INCLUDES ", \"-I${_d}\"")
endforeach()
string(APPEND PSOPT_PY_INCLUDES ", \"-I/usr/local/include\", \"-DHAVE_CSTDDEF\"")
set(PSOPT_PY_LIBS     "\"${CPPAD_LIBRARY}\"")

configure_file(
    "${CMAKE_CURRENT_SOURCE_DIR}/psopt/_toolchain.py.in"
    "${CMAKE_CURRENT_SOURCE_DIR}/psopt/_toolchain.py"
    @ONLY)
