# - Try to find snopt7
# Once done, this will define
#
#  snopt7_FOUND - system has snopt7
#  snopt7_INCLUDE_DIRS - the snopt7 include directories
#  snopt7_LIBRARIES - link these to use snopt7

find_package(PkgConfig)
pkg_check_modules(snopt7 QUIET snopt7)

set(snopt7_DEFINITIONS ${snopt7_CFLAGS_OTHER})

find_library(snopt7_LIBRARY NAMES snopt7 snopt7_cpp
             HINTS ${snopt7_LIBDIR} ${snopt7_LIBRARY_DIRS} )

include(FindPackageHandleStandardArgs)
# handle the QUIETLY and REQUIRED arguments and set snopt7_FOUND to TRUE
# if all listed variables are TRUE
find_package_handle_standard_args(snopt7  DEFAULT_MSG
                                  snopt7_LIBRARY)

mark_as_advanced(snopt7_LIBRARY)

set(snopt7_LIBRARIES ${snopt7_LIBRARY})
