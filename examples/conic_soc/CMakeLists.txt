cmake_minimum_required (VERSION 3.10)
project(conic_soc LANGUAGES CXX)

add_executable(${PROJECT_NAME} conic_soc.cxx)
add_dependencies(${PROJECT_NAME} PSOPT)
target_link_libraries(${PROJECT_NAME} PUBLIC PSOPT PRIVATE PkgConfig::ipopt Eigen3::Eigen)
