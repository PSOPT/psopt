cmake_minimum_required (VERSION 3.10)
project(alpine LANGUAGES CXX)

find_package(Eigen3 REQUIRED NO_MODULE)
find_package(adolc REQUIRED)
find_package(ColPack REQUIRED)
find_package(PkgConfig REQUIRED)
pkg_check_modules(ipopt REQUIRED ipopt)

add_executable(${PROJECT_NAME} alpine.cxx)
target_include_directories(${PROJECT_NAME}
PRIVATE 
    ${ipopt_INCLUDE_DIRS})

set(CMAKE_CXX_FLAGS "-std=c++11 -g")

target_link_libraries(${PROJECT_NAME} PUBLIC PSOPT PRIVATE -L${ipopt_LIBRARY_DIRS} ipopt adolc ColPack Eigen3::Eigen)
