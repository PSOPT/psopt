cmake_minimum_required (VERSION 3.10)
project(ipopt_cam LANGUAGES CXX)

# This example does not use PSOPT. It is the worked example of Chapter 2 of
# "Computational Optimal Control and Estimation", which calls IPOPT directly
# through its TNLP interface and gets sparse derivatives from CppAD, and it is
# shipped so that a reader of the book can build and run it. It therefore links
# IPOPT and CppAD but not PSOPT, and it is built only when CppAD is the AD
# backend, since that is the only configuration in which CppAD is looked for.

if(NOT PSOPT_AD_BACKEND STREQUAL "CPPAD")
    message(STATUS "ipopt_cam: skipped (needs PSOPT_AD_BACKEND=CPPAD)")
    return()
endif()

add_executable(${PROJECT_NAME} ipopt_cam.cxx)
target_include_directories(${PROJECT_NAME} PRIVATE ${CPPAD_INCLUDE_DIR})
target_link_libraries(${PROJECT_NAME} PRIVATE ${CPPAD_LIBRARY} PkgConfig::ipopt)
