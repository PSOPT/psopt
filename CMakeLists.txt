cmake_minimum_required (VERSION 3.12)
string(TIMESTAMP BUILD_DATE "%Y.%m.%d")

execute_process(
  COMMAND git describe --tags --abbrev=0
  WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
  OUTPUT_VARIABLE PSOPT_RELEASE_DATE
  OUTPUT_STRIP_TRAILING_WHITESPACE
) 

message(STATUS "PSOPT Rolling Release Date: ${PSOPT_RELEASE_DATE}")

# Define it as a version macro
add_definitions(-DPSOPT_VERSION__="${PSOPT_RELEASE_DATE}")


project(PSOPT VERSION ${PSOPT_VERSION__} LANGUAGES CXX)

option(BUILD_EXAMPLES "Build examples from the example subdirectory." OFF)
option(WITH_SNOPT_INTERFACE "Build interface for the SNOPT optimizer." OFF)
option(BUILD_TESTS "Build tests for PSOPT." OFF)
option(HEADLESS "Headless option (Gnuplot does not spawn windows with plots, but generates PDFs)" OFF)

list(APPEND CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/cmake)
find_package(Eigen3 REQUIRED NO_MODULE)
find_package(adolc REQUIRED)
find_package(PkgConfig REQUIRED)
pkg_check_modules(ipopt REQUIRED IMPORTED_TARGET ipopt)

if(${WITH_SNOPT_INTERFACE})
    add_subdirectory(snopt-interface/)
endif()

file(GLOB SRC src/*.cxx)
file(GLOB HEADER include/*.h)

add_library(${PROJECT_NAME} $<$<BOOL:$<TARGET_EXISTS:PSOPT_SNOPT_interface>>:$<TARGET_OBJECTS:PSOPT_SNOPT_interface>> ${SRC} ${HEADER})
target_include_directories(${PROJECT_NAME}
PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
    $<INSTALL_INTERFACE:include/${PROJECT_NAME}>
)

target_link_libraries(${PROJECT_NAME} PRIVATE adolc Eigen3::Eigen PUBLIC PkgConfig::ipopt)
add_definitions(-DPSOPT_RELEASE_STRING="${PSOPT_VERSION}" -DHAVE_CSTDDEF)

target_compile_features(${PROJECT_NAME} PUBLIC cxx_delegating_constructors)

if(${WITH_SNOPT_INTERFACE})
    add_definitions(-DPSOPT_RELEASE_STRING="${PSOPT_VERSION}" -DHAVE_CSTDDEF -DUSE_SNOPT)
endif()

if(${HEADLESS})
    add_definitions(-DPSOPT_RELEASE_STRING="${PSOPT_VERSION}" -DHAVE_CSTDDEF -DHEADLESS)
endif()

if(${BUILD_EXAMPLES})
	add_subdirectory(examples/)
endif()

if(${BUILD_TESTS})
	add_subdirectory(tests/)
endif()


# INSTALLING and PAKAGING
set(PACKAGE_FIND_NAME "${PROJECT_NAME}")
install(TARGETS ${PROJECT_NAME} EXPORT ${PROJECT_NAME} DESTINATION "lib/${PROJECT_NAME}")
install(DIRECTORY "${CMAKE_SOURCE_DIR}/include/" DESTINATION "include/${PROJECT_NAME}")
install(FILES psopt-config.cmake DESTINATION "share/cmake/${PROJECT_NAME}/")
install(EXPORT ${PROJECT_NAME} DESTINATION "share/cmake/${PROJECT_NAME}")
