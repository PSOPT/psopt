cmake_minimum_required (VERSION 3.10)
project(breakwell LANGUAGES CXX)

add_executable(${PROJECT_NAME} breakwell.cxx)
add_dependencies(${PROJECT_NAME} PSOPT)
target_include_directories(${PROJECT_NAME} PRIVATE ${adolc_INCLUDE_DIRS})
target_link_libraries(${PROJECT_NAME} PUBLIC PSOPT PRIVATE PkgConfig::ipopt Eigen3::Eigen)
