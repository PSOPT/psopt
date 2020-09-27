cmake_minimum_required(VERSION 3.10)
project(PSOPT_Tests)

find_package(GTest)
find_package(Threads REQUIRED)

if(NOT ${GTest_FOUND})
	message(STATUS "Google Test not found. It will be downloaded.")
	# Download and unpack googletest at configure time
	configure_file(CMakeLists.txt.in googletest-download/CMakeLists.txt)
	execute_process(COMMAND ${CMAKE_COMMAND} -G "${CMAKE_GENERATOR}" .
	  RESULT_VARIABLE result
	  WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/googletest-download )
	if(result)
	  message(FATAL_ERROR "CMake step for googletest failed: ${result}")
	endif()
	execute_process(COMMAND ${CMAKE_COMMAND} --build .
	  RESULT_VARIABLE result
	  WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/googletest-download )
	if(result)
	  message(FATAL_ERROR "Build step for googletest failed: ${result}")
	endif()
	
	add_subdirectory(${CMAKE_CURRENT_BINARY_DIR}/googletest-src
	  ${CMAKE_CURRENT_BINARY_DIR}/googletest-build
	  EXCLUDE_FROM_ALL)
endif()
		 
file(GLOB SRC "*.cpp")
add_executable(${PROJECT_NAME} ${SRC})
add_dependencies(${PROJECT_NAME} PSOPT gtest)

target_link_libraries(${PROJECT_NAME} PRIVATE
    gtest_main
    Threads::Threads
    Eigen3::Eigen
    adolc
    PkgConfig::ipopt 
    PSOPT
    )

if(${WITH_SNOPT_INTERFACE})
    enable_language(Fortran)
    target_link_libraries(${PROJECT_NAME} PRIVATE gfortran ${snopt7_LIBRARIES} PSOPT_SNOPT_interface)
endif()

enable_testing()
add_test(test ${PROJECT_NAME})
