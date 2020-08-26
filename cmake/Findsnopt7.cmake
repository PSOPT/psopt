# - Try to find snopt7
# Once done, this will define
#
#  snopt7_FOUND - system has snopt7
#  snopt7_LIBRARIES - link these to use snopt7

include(FindPackageHandleStandardArgs)

find_library(snopt7_LIBRARY NAMES snopt7_cpp)

# handle the QUIETLY and REQUIRED arguments and set snopt7_FOUND to TRUE
# if all listed variables are TRUE
find_package_handle_standard_args(snopt7  DEFAULT_MSG
                                  snopt7_LIBRARY)
mark_as_advanced(snopt7_LIBRARY)
set(snopt7_LIBRARIES ${snopt7_LIBRARY})

message(STATUS "Found snopt7_cpp: ${snopt7_LIBRARY}")
