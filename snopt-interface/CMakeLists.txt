cmake_minimum_required (VERSION 3.10)
project(PSOPT_SNOPT_interface LANGUAGES Fortran C CXX)

file(GLOB SRC src/*.cpp src/*.c src/*.f90)
file(GLOB HEADER include/*.h include/*.hpp)

add_library(${PROJECT_NAME} OBJECT ${SRC} ${HEADER})
target_include_directories(${PROJECT_NAME}
PUBLIC 
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include/>
    $<INSTALL_INTERFACE:include/${PROJECT_NAME}>)
target_link_libraries(${PROJECT_NAME} PRIVATE gfortran ${snopt7_LIBRARIES})

if(${BUILD_EXAMPLES})
	add_subdirectory(cppexamples/)
	add_subdirectory(cexamples/)
endif()
