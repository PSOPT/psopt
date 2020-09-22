cmake_minimum_required (VERSION 3.10)
project(PSOPT_SNOPT_interface LANGUAGES Fortran C CXX)

find_package(snopt7 REQUIRED)
message(STATUS "SNOPT interface will be build.")

file(GLOB SRC src/*.cpp src/*.c src/*.f90)
file(GLOB HEADER include/*.h include/*.hpp)

add_library(${PROJECT_NAME} OBJECT ${SRC} ${HEADER})
target_include_directories(${PROJECT_NAME}
PUBLIC 
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include/>
    $<INSTALL_INTERFACE:include/PSOPT/snopt>)
target_link_libraries(${PROJECT_NAME} PRIVATE gfortran adolc ${snopt7_LIBRARIES})
target_compile_features(${PROJECT_NAME} PUBLIC cxx_delegating_constructors)

if(${BUILD_EXAMPLES})
    add_subdirectory(cppexamples/)
    add_subdirectory(cexamples/)
endif()

install(DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/include/" DESTINATION "include/PSOPT/snopt")
