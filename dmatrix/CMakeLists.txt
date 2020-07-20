cmake_minimum_required (VERSION 3.10)
project(dmatrix LANGUAGES CXX)

file(GLOB SRC src/*.cxx)
file(GLOB HEADER include/*.h)

add_library(${PROJECT_NAME} ${SRC} ${HEADER})
target_include_directories(${PROJECT_NAME} PUBLIC include/)
set(${PROJECT_NAME} ${PROJECT_NAME} CACHE INTERNAL "${PROJECT_NAME} library target.")

if(${BUILD_EXAMPLES})
	add_subdirectory(examples/)
endif()