cmake_minimum_required (VERSION 3.10)
project(PSOPT_SNOPT_interface LANGUAGES CXX)

find_package(snopt7 REQUIRED)

file(GLOB SRC src/*.cpp)
file(GLOB HEADER include/*.h include/*.hpp)

add_library(${PROJECT_NAME} STATIC ${SRC} ${HEADER})
target_include_directories(${PROJECT_NAME}
PUBLIC 
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include/>
    $<INSTALL_INTERFACE:include/${PROJECT_NAME}>)
target_link_libraries(${PROJECT_NAME} PRIVATE gfortran ${snopt7_LIBRARIES})

if(${BUILD_EXAMPLES})
	add_subdirectory(cppexamples/)
endif()
