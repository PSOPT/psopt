cmake_minimum_required (VERSION 3.10)
project(mineng_target_set LANGUAGES CXX)

add_executable(${PROJECT_NAME} mineng_target_set.cxx)
add_dependencies(${PROJECT_NAME} PSOPT)
target_link_libraries(${PROJECT_NAME} PUBLIC PSOPT PRIVATE PkgConfig::ipopt Eigen3::Eigen)
