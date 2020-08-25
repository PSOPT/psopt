cmake_minimum_required (VERSION 3.10)
project(delay1 LANGUAGES CXX)

find_package(Eigen3 REQUIRED NO_MODULE)
find_package(PkgConfig REQUIRED)
find_package(adolc REQUIRED)
pkg_check_modules(ipopt REQUIRED ipopt)
find_package(ColPack REQUIRED)

add_executable(${PROJECT_NAME} delay1.cxx)

target_include_directories(${PROJECT_NAME} PRIVATE ${ipopt_INCLUDE_DIRS})
target_link_libraries(${PROJECT_NAME} PUBLIC PSOPT PRIVATE ${ipopt_LIBRARIES} adolc ColPack Eigen3::Eigen)