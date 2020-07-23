cmake_minimum_required (VERSION 3.10)
project(dmatrix LANGUAGES CXX)

file(GLOB SRC src/*.cxx)
file(GLOB HEADER include/*.h)

# Download and unpack SuiteSparse at configure time

configure_file(CMakeLists-SuiteSparse.txt.in SuiteSparse-download/CMakeLists.txt)
execute_process(COMMAND ${CMAKE_COMMAND} -G "${CMAKE_GENERATOR}" .
  RESULT_VARIABLE result
  WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/SuiteSparse-download )

if(result)
	message(FATAL_ERROR "CMake step for SuiteSparse failed: ${result}")
endif()

execute_process(COMMAND ${CMAKE_COMMAND} --build .
	RESULT_VARIABLE result
	WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/SuiteSparse-download )

if(result)
	message(FATAL_ERROR "Build step for SuiteSparse failed: ${result}")
endif()

execute_process(COMMAND make -j4
	RESULT_VARIABLE result
	WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/SuiteSparse-src/)

if(result)
	message(FATAL_ERROR "Install step for SuiteSparse failed: ${result}")
endif()


add_library(SuiteSparse STATIC IMPORTED)
target_include_directories(SuiteSparse INTERFACE ${CMAKE_CURRENT_BINARY_DIR}/SuiteSparse-src/include/ ${CMAKE_CURRENT_BINARY_DIR}/SuiteSparse-src/include/)
set_target_properties(SuiteSparse PROPERTIES IMPORTED_LOCATION ${CMAKE_CURRENT_BINARY_DIR}/SuiteSparse-src/lib)

# --- download and unpack lusol
configure_file(CMakeLists-lusol.txt.in lusol-download/CMakeLists.txt)
execute_process(COMMAND ${CMAKE_COMMAND} -G "${CMAKE_GENERATOR}" .
  RESULT_VARIABLE result
  WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/lusol-download )

if(result)
	message(FATAL_ERROR "CMake step for lusol failed: ${result}")
endif()

execute_process(COMMAND ${CMAKE_COMMAND} --build .
	RESULT_VARIABLE result
	WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/lusol-download )

if(result)
	message(FATAL_ERROR "Build step for lusol failed: ${result}")
endif()

file(GLOB LUSOL_SRC ${CMAKE_CURRENT_BINARY_DIR}/lusol-src/csrc/*.c)
list(REMOVE_ITEM LUSOL_SRC ${CMAKE_CURRENT_BINARY_DIR}/lusol-src/csrc/LU1LOT.c)
file(GLOB LUSOL_HEADERS ${CMAKE_CURRENT_BINARY_DIR}/lusol-src/csrc/*.h)
add_library(lusol ${LUSOL_SRC} ${LUSOL_HEADERS})
target_include_directories(lusol PUBLIC ${CMAKE_CURRENT_BINARY_DIR}/lusol-src/csrc/)
# ---

add_library(${PROJECT_NAME} ${SRC} ${HEADER})
target_include_directories(${PROJECT_NAME} PUBLIC include/ PRIVATE ${SuiteSparse_INCLUDE_DIRS})
target_link_libraries(${PROJECT_NAME} PRIVATE SuiteSparse lusol)

set(${PROJECT_NAME} ${PROJECT_NAME} CACHE INTERNAL "${PROJECT_NAME} library target.")

if(${BUILD_EXAMPLES})
	add_subdirectory(examples/)
endif()
