cmake_minimum_required (VERSION 3.10)
project(PSOPT VERSION 5.00 LANGUAGES CXX)

option(BUILD_EXAMPLES "Build examples from the example subdirectory." OFF)
option(WITH_SNOPT_INTERFACE "Build interface for the SNOPT optimizer." OFF)

list(APPEND CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/cmake)
find_package(Eigen3 REQUIRED NO_MODULE)
find_package(adolc REQUIRED)
find_package(PkgConfig REQUIRED)
pkg_check_modules(ipopt REQUIRED IMPORTED_TARGET ipopt)

if(${WITH_SNOPT_INTERFACE})
    find_package(snopt7 REQUIRED)
    message(STATUS "SNOPT7 found. SNOPT interface will be build.")
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

target_link_libraries(${PROJECT_NAME} PRIVATE adolc PkgConfig::ipopt Eigen3::Eigen)
add_definitions(-DPSOPT_RELEASE_STRING="${PSOPT_VERSION}" -DHAVE_CSTDDEF)

if(${BUILD_EXAMPLES})
	add_subdirectory(examples/)
    if(${WITH_SNOPT_INTERFACE})
        add_subdirectory(snopt-interface/cppexamples/)
        add_subdirectory(snopt-interface/cexamples/)
    endif()
endif()

# INSTALLING and PAKAGING
set(PACKAGE_FIND_NAME "${PROJECT_NAME}")
install(TARGETS ${PROJECT_NAME} EXPORT ${PROJECT_NAME} DESTINATION "lib/${PROJECT_NAME}")
install(DIRECTORY "${CMAKE_SOURCE_DIR}/include/" DESTINATION "include/${PROJECT_NAME}")
install(FILES psopt-config.cmake DESTINATION "share/cmake/${PROJECT_NAME}/")
install(EXPORT ${PROJECT_NAME} DESTINATION "share/cmake/${PROJECT_NAME}")
