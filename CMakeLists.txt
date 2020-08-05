cmake_minimum_required (VERSION 3.10)
project(PSOPT VERSION 4.0.1 LANGUAGES CXX)

option(BUILD_EXAMPLES "Build examples from the example subdirectory." OFF)

find_package(Eigen3 REQUIRED NO_MODULE)
find_package(PkgConfig REQUIRED)
pkg_check_modules(adolc adolc)
pkg_check_modules(ipopt REQUIRED ipopt)

if(NOT TARGET ${adolc})
    list(APPEND CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake/")
    find_package(ColPack REQUIRED)
	find_package(Python2 REQUIRED COMPONENTS Development)

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
target_include_directories(${PROJECT_NAME} PUBLIC src/ PRIVATE ${ipopt_INCLUDE_DIRS})
target_link_libraries(${PROJECT_NAME} PUBLIC adolc ${ipopt_LIBRARIES} PRIVATE Eigen3::Eigen)
add_definitions(-DPSOPT_RELEASE_STRING="${PSOPT_VERSION}")

if(${BUILD_EXAMPLES})
	add_subdirectory(examples/)
endif()
