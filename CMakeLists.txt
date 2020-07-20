cmake_minimum_required (VERSION 3.10)
project(PSOPT_ROOT)

option(BUILD_EXAMPLES "Build examples for dmatrix and PSOPT" OFF)
option(BUILD_DOCS "Build documentation for dmatrix and PSOPT" OFF)

add_subdirectory(dmatrix/)
add_subdirectory(PSOPT/)
