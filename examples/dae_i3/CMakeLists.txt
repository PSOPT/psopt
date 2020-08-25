cmake_minimum_required (VERSION 3.10)
project(dae_i3 LANGUAGES CXX)

add_executable(${PROJECT_NAME} dae_i3.cxx)
find_package(Eigen3 REQUIRED NO_MODULE)
find_package(PkgConfig REQUIRED)
pkg_check_modules(adolc adolc)
pkg_check_modules(ipopt REQUIRED ipopt)
pkg_check_modules(ColPack REQUIRED ColPack)


target_include_directories(${PROJECT_NAME}
PUBLIC 
    $<BUILD_INTERFACE:${CMAKE_SOURCE_DIR}/include/>
    $<INSTALL_INTERFACE:include/${PROJECT_NAME}>
PRIVATE 
    ${ipopt_INCLUDE_DIRS}
    ${adolc_INCLUDE_DIRS}
    ${ColPack_INCLUDE_DIRS})

set(CMAKE_CXX_FLAGS "-std=c++11 -g")

target_link_libraries(${PROJECT_NAME} PUBLIC PSOPT PRIVATE -L${ipopt_LIBRARY_DIRS} ipopt adolc ColPack Eigen3::Eigen)
