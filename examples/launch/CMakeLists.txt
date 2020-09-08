cmake_minimum_required (VERSION 3.10)
project(launch LANGUAGES CXX)

add_executable(${PROJECT_NAME} launch.cxx)
add_dependencies(${PROJECT_NAME} PSOPT)
target_link_libraries(${PROJECT_NAME} PUBLIC PSOPT PRIVATE PkgConfig::ipopt adolc Eigen3::Eigen)
