cmake_minimum_required (VERSION 3.10)
project(PSOPT_SNOPT_interface LANGUAGES CXX)

find_package(snopt7 REQUIRED)

file(GLOB SRC src/*.cpp)
file(GLOB HEADER include/*.h)

add_library(${PROJECT_NAME} ${SRC} ${HEADER})
target_include_directories(${PROJECT_NAME}
PUBLIC 
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include/>
    $<INSTALL_INTERFACE:include/${PROJECT_NAME}>)
target_link_directories(${PROJECT_NAME} PRIVATE snopt7)

