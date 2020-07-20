cmake_minimum_required (VERSION 3.10)
project(PSOPT VERSION 4.0.1 LANGUAGES CXX)

find_package(PkgConfig REQUIRED)
pkg_check_modules(adolc adolc)

if(NOT TARGET ${adolc})
	message(STATUS "ADOL-C will be downloaded and added to the project...")
	# Download and unpack adolc at configure time
	configure_file(CMakeLists-adolc.txt.in adolc-download/CMakeLists.txt)
	execute_process(COMMAND ${CMAKE_COMMAND} -G "${CMAKE_GENERATOR}" .
  	RESULT_VARIABLE result
  	WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/adolc-download )
	if(result)
	  message(FATAL_ERROR "CMake step for adolc failed: ${result}")
	endif()
	execute_process(COMMAND ${CMAKE_COMMAND} --build .
	  RESULT_VARIABLE result
	  WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/adolc-download )
	if(result)
	  message(FATAL_ERROR "Build step for adolc failed: ${result}")
	endif()

	add_library(adolc SHARED IMPORTED)
    target_include_directories(adolc INTERFACE ${CMAKE_CURRENT_BINARY_DIR}/adolc-build/include/)
	set_target_properties(adolc PROPERTIES IMPORTED_LOCATION ${CMAKE_CURRENT_BINARY_DIR}/adolc-build/lib64)
 
    set(CMAKE_CXX_FLAGS "-std=c++11")
endif()

file(GLOB SRC src/*.cxx)
file(GLOB HEADER "src/*.h")

add_library(${PROJECT_NAME} ${SRC} ${HEADER})
add_dependencies(${PROJECT_NAME} dmatrix)
target_include_directories(${PROJECT_NAME} PUBLIC src/)
target_link_libraries(${PROJECT_NAME} PUBLIC adolc PRIVATE dmatrix)
