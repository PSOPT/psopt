cmake_minimum_required (VERSION 3.10)
project(mixed_ic LANGUAGES CXX)

add_executable(${PROJECT_NAME} mixed_ic.cxx)
add_dependencies(${PROJECT_NAME} PSOPT)
target_link_libraries(${PROJECT_NAME} PUBLIC PSOPT PRIVATE PkgConfig::ipopt Eigen3::Eigen)
