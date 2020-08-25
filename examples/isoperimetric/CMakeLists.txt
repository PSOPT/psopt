cmake_minimum_required (VERSION 3.10)
project(isop LANGUAGES CXX)

add_executable(${PROJECT_NAME} isoperimetric.cxx)
find_package(Eigen3 REQUIRED NO_MODULE)
find_package(PkgConfig REQUIRED)
find_package(adolc REQUIRED)
pkg_check_modules(ipopt REQUIRED ipopt)
find_package(ColPack REQUIRED)


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
