cmake_minimum_required(VERSION 3.16)
project(PSOPT_Tests LANGUAGES CXX)

find_package(Threads REQUIRED)

# GoogleTest: use a system install if present, else fetch it.
find_package(GTest QUIET)
if(NOT GTest_FOUND)
  message(STATUS "GoogleTest not found — fetching release-1.14.0")
  include(FetchContent)
  FetchContent_Declare(googletest
    GIT_REPOSITORY https://github.com/google/googletest.git
    GIT_TAG release-1.14.0 GIT_SHALLOW TRUE)
  set(gtest_force_shared_crt ON CACHE BOOL "" FORCE)
  set(INSTALL_GTEST OFF CACHE BOOL "" FORCE)
  FetchContent_MakeAvailable(googletest)
endif()

file(GLOB SRC "*.cpp")
add_executable(${PROJECT_NAME} ${SRC})
add_dependencies(${PROJECT_NAME} psopt)

# PSOPT::psopt carries Eigen/ADOL-C/IPOPT as PUBLIC usage requirements.
target_link_libraries(${PROJECT_NAME} PRIVATE
    GTest::gtest_main
    Threads::Threads
    PSOPT::psopt)

enable_testing()
add_test(NAME psopt_tests COMMAND ${PROJECT_NAME})
