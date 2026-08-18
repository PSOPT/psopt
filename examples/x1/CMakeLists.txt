cmake_minimum_required (VERSION 3.10)
project(x1_example LANGUAGES CXX)

add_executable(${PROJECT_NAME} x1_example.cxx)
add_dependencies(${PROJECT_NAME} PSOPT)
target_link_libraries(${PROJECT_NAME} PUBLIC PSOPT PRIVATE PkgConfig::ipopt Eigen3::Eigen)
