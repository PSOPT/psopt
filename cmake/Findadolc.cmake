
# - Try to find adolc
# Once done, this will define
#
#  adolc_FOUND - system has adolc
#  adolc_INCLUDE_DIRS - the adolc include directories
#  adolc_LIBRARIES - link these to use adolc

find_package(PkgConfig)
include(FindPackageHandleStandardArgs)

pkg_check_modules(adolc QUIET IMPORTED_TARGET adolc)

if(${adolc_FOUND}) # if Adolc could be found by pkgconfig
    add_library(adolc ALIAS PkgConfig::adolc)
else()  # is it already installed locally by this file?
    # sometimes, AdolC will be downloaded each time the user calls cmake. prevent this by searching compiled files in the build dir
    find_path(adolc_INCLUDE_DIR adolc.h
            HINTS ${CMAKE_BINARY_DIR}/adolc-build/include/
            NO_CMAKE_PATH)

    find_library(adolc_LIBRARY adolc
                HINTS ${CMAKE_BINARY_DIR}/adolc-build/lib64
                NO_CMAKE_PATH)
                
    if(NOT EXISTS ${adolc_INCLUDE_DIR}) # if everything fails, download it
        message(STATUS "AdolC has not been installed on this system and will be automatically added to this project.")
        
        find_package(ColPack REQUIRED)
        find_package(Python2 REQUIRED COMPONENTS Development)

        # Download and unpack adolc at configure time
        configure_file(cmake/CMakeLists-adolc.txt.in ${CMAKE_BINARY_DIR}/adolc-download/CMakeLists.txt)
        
        execute_process(COMMAND ${CMAKE_COMMAND} -G "${CMAKE_GENERATOR}" .
        RESULT_VARIABLE result
        WORKING_DIRECTORY ${CMAKE_BINARY_DIR}/adolc-download/ )
        if(result)
            message(FATAL_ERROR "CMake step for adolc failed: ${result}")
        endif()
        
        execute_process(COMMAND ${CMAKE_COMMAND} --build .
        RESULT_VARIABLE result
        WORKING_DIRECTORY ${CMAKE_BINARY_DIR}/adolc-download )
        if(result)
            message(FATAL_ERROR "Build step for adolc failed: ${result}")
        endif()

        add_library(adolc SHARED IMPORTED)
        target_include_directories(adolc INTERFACE ${CMAKE_BINARY_DIR}/adolc-build/include)
        set_target_properties(adolc PROPERTIES IMPORTED_LOCATION ${CMAKE_BINARY_DIR}/adolc-build/lib64)
        
        find_package_handle_standard_args(adolc DEFAULT_MSG
                                        ${adolc_LIBRARIES} ${adolc_INCLUDE_DIRS})
    else()
        # handle the QUIETLY and REQUIRED arguments and set adolc_FOUND to TRUE
        # if all listed variables are TRUE
        find_package_handle_standard_args(adolc  DEFAULT_MSG
                                        adolc_LIBRARY adolc_INCLUDE_DIR)

        mark_as_advanced(adolc_INCLUDE_DIR adolc_LIBRARY )

        add_library(adolc SHARED IMPORTED)
        target_include_directories(adolc INTERFACE ${adolc_INCLUDE_DIR})
        set_target_properties(adolc PROPERTIES IMPORTED_LOCATION ${CMAKE_BINARY_DIR}/adolc-build/lib64)
    endif()
endif()

