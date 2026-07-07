
# - Try to find adolc
# Once done, this will define
#
#  adolc_FOUND - system has adolc
#  adolc_INCLUDE_DIRS - the adolc include directories
#  adolc_LIBRARIES - link these to use adolc

find_package(PkgConfig)
include(FindPackageHandleStandardArgs)

pkg_check_modules(adolc QUIET IMPORTED_TARGET GLOBAL adolc)

if(${adolc_FOUND}) # if Adolc could be found by pkgconfig
    add_library(adolc ALIAS PkgConfig::adolc)
    # ADOL-C's sparse drivers use ColPack. A statically-linked libadolc (e.g. the MinGW
    # superbuild) leaves those symbols for the final executable link, so ColPack must be
    # on the link line. Pull it into the adolc link interface when present (harmless where
    # a shared libadolc already embeds it).
    find_package(ColPack QUIET)
    if(ColPack_FOUND)
        set_property(TARGET PkgConfig::adolc APPEND PROPERTY INTERFACE_LINK_LIBRARIES ${ColPack_LIBRARIES})
    endif()
else()  # no pkg-config (e.g. MSVC/Windows): search a user-supplied ADOL-C, then build dir
    # Honor ADOLC_ROOT / adolc_ROOT (e.g. a Visual Studio ADOL-C build) as well as the
    # in-tree build dir. NB: ColPack is only needed to *build* ADOL-C's sparse drivers;
    # linking a prebuilt ADOL-C (sparse or not — e.g. the Fortran-free MSVC/CasADi path)
    # does NOT require finding ColPack here. It is therefore only pulled in below in the
    # genuine download-and-build fallback.
    find_path(adolc_INCLUDE_DIR adolc/adouble.h
            HINTS ${ADOLC_ROOT} ${adolc_ROOT} $ENV{ADOLC_ROOT} ${CMAKE_BINARY_DIR}/adolc-build
            PATH_SUFFIXES include)

    find_library(adolc_LIBRARY NAMES adolc
            HINTS ${ADOLC_ROOT} ${adolc_ROOT} $ENV{ADOLC_ROOT} ${CMAKE_BINARY_DIR}/adolc-build
            PATH_SUFFIXES lib lib64)

    if(adolc_INCLUDE_DIR AND adolc_LIBRARY)  # prebuilt ADOL-C found (ADOLC_ROOT or build dir)
        find_package_handle_standard_args(adolc DEFAULT_MSG adolc_LIBRARY adolc_INCLUDE_DIR)
        mark_as_advanced(adolc_INCLUDE_DIR adolc_LIBRARY)
        add_library(adolc UNKNOWN IMPORTED)
        set_target_properties(adolc PROPERTIES
            IMPORTED_LOCATION ${adolc_LIBRARY}
            INTERFACE_INCLUDE_DIRECTORIES ${adolc_INCLUDE_DIR})
        find_package(ColPack QUIET)
        if(ColPack_FOUND)
            set_property(TARGET adolc APPEND PROPERTY INTERFACE_LINK_LIBRARIES ${ColPack_LIBRARIES})
        endif()
    elseif(MSVC OR CMAKE_HOST_WIN32)
        # The autotools download-build fallback below cannot run on MSVC. Fail with a
        # clear instruction instead of a confusing ColPack error.
        message(FATAL_ERROR
          "ADOL-C not found. On Windows/MSVC, build ADOL-C from its Visual Studio solution "
          "(ADOL-C/MSVisualStudio/v14/adolc.sln, x64/Release) and set:\n"
          "    -DADOLC_ROOT=<prefix>\n"
          "so this module finds <prefix>/include/adolc/adouble.h and <prefix>/lib/adolc.lib "
          "(or point -Dadolc_INCLUDE_DIR=<...> and -Dadolc_LIBRARY=<...> directly). "
          "ColPack is NOT needed for a Fortran-free MSVC/CasADi build (ADOL-C built without "
          "--enable-sparse). The autotools download-build fallback does not work on MSVC.")
    else()  # nothing found: download + build ADOL-C (Unix/autotools; needs ColPack + Python2)
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

        add_library(adolc UNKNOWN IMPORTED)
        target_include_directories(adolc INTERFACE ${CMAKE_BINARY_DIR}/adolc-build/include)
        set_target_properties(adolc PROPERTIES IMPORTED_LOCATION ${CMAKE_BINARY_DIR}/adolc-build/lib64)
    endif()
endif()

