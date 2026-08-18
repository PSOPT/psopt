cmake_minimum_required(VERSION 3.10)
project(PSOPT_Tests)

find_package(GTest)
find_package(Threads REQUIRED)

if(NOT GTest_FOUND)
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
	set(PSOPT_GTEST_MAIN gtest_main)
else()
	set(PSOPT_GTEST_MAIN GTest::gtest_main)
endif()
		 
file(GLOB SRC "*.cpp")
add_executable(${PROJECT_NAME} ${SRC})
add_dependencies(${PROJECT_NAME} PSOPT)

target_link_libraries(${PROJECT_NAME} PRIVATE
    ${PSOPT_GTEST_MAIN}
    Threads::Threads
    Eigen3::Eigen
    PkgConfig::ipopt 
    PSOPT
    )

enable_testing()
add_test(test ${PROJECT_NAME})

# A second pass with the environment GALAHAD's linear solver requires. It cannot be set
# from inside the process -- the OpenMP runtime reads these once, when it initialises --
# so the tests that need it skip themselves in the ordinary pass and run here.
add_test(NAME test_galahad_environment COMMAND ${PROJECT_NAME} --gtest_filter=SQPSolver.*)
set_tests_properties(test_galahad_environment PROPERTIES
    ENVIRONMENT "OMP_CANCELLATION=TRUE;OMP_PROC_BIND=TRUE")
