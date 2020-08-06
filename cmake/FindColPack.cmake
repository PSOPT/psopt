# - Try to find ColPack
# Once done, this will define
#
#  ColPack_FOUND - system has ColPack
#  ColPack_INCLUDE_DIRS - the ColPack include directories
#  ColPack_LIBRARIES - link these to use ColPack

find_package(PkgConfig)
pkg_check_modules(ColPack QUIET ColPack)

set(ColPack_DEFINITIONS ${ColPack_CFLAGS_OTHER})

find_path(ColPack_INCLUDE_DIR ColPack/ColPackHeaders.h
          HINTS ${ColPack_INCLUDEDIR} ${ColPack_INCLUDE_DIRS}
          PATH_SUFFIXES libxml2 )

find_library(ColPack_LIBRARY NAMES ColPack
             HINTS ${ColPack_LIBDIR} ${ColPack_LIBRARY_DIRS} )

include(FindPackageHandleStandardArgs)
# handle the QUIETLY and REQUIRED arguments and set ColPack_FOUND to TRUE
# if all listed variables are TRUE
find_package_handle_standard_args(ColPack  DEFAULT_MSG
                                  ColPack_LIBRARY ColPack_INCLUDE_DIR)

mark_as_advanced(ColPack_INCLUDE_DIR ColPack_LIBRARY )

set(ColPack_LIBRARIES ${ColPackLIBRARY} )
set(ColPack_INCLUDE_DIRS ${ColPack_INCLUDE_DIR} )
