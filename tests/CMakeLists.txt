cmake_minimum_required(VERSION 3.16)
project(PSOPT_Tests LANGUAGES CXX)

find_package(Threads REQUIRED)

# Always build GoogleTest from source with THIS project's compiler: a system
# GoogleTest may use a different C++ ABI (e.g. libc++/Abseil) than a gcc/libstdc++
# PSOPT build, causing link failures. FetchContent guarantees an ABI match.
include(FetchContent)
FetchContent_Declare(googletest
  GIT_REPOSITORY https://github.com/google/googletest.git
  GIT_TAG v1.15.2 GIT_SHALLOW TRUE)
set(gtest_force_shared_crt ON CACHE BOOL "" FORCE)
set(INSTALL_GTEST OFF CACHE BOOL "" FORCE)
FetchContent_MakeAvailable(googletest)

file(GLOB SRC CONFIGURE_DEPENDS "*.cpp")
add_executable(${PROJECT_NAME} ${SRC})
add_dependencies(${PROJECT_NAME} psopt)

# PSOPT::psopt carries Eigen/ADOL-C/IPOPT as PUBLIC usage requirements.
target_link_libraries(${PROJECT_NAME} PRIVATE
    GTest::gtest_main
    Threads::Threads
    PSOPT::psopt)

enable_testing()
add_test(NAME psopt_tests COMMAND ${PROJECT_NAME})
