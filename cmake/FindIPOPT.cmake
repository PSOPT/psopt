# FindIPOPT.cmake - locate a COIN-OR IPOPT installation and expose IPOPT::ipopt.
#
# Search order:
#   1. IPOPT_ROOT / IPOPT_DIR (CMake var or env) prefix.
#   2. pkg-config "ipopt" (honours PKG_CONFIG_PATH, e.g. /opt/ipopt/lib/pkgconfig).
#   3. Standard system prefixes.
#
# Sets IPOPT_FOUND and an imported target IPOPT::ipopt with include dirs + libs.

if(TARGET IPOPT::ipopt)
  set(IPOPT_FOUND TRUE)
  return()
endif()

set(_ipopt_hints "")
foreach(_h IPOPT_ROOT IPOPT_DIR)
  if(DEFINED ${_h})
    list(APPEND _ipopt_hints "${${_h}}")
  endif()
  if(DEFINED ENV{${_h}})
    list(APPEND _ipopt_hints "$ENV{${_h}}")
  endif()
endforeach()

# 1) pkg-config (adds its result as a fallback for include/lib discovery)
find_package(PkgConfig QUIET)
if(PkgConfig_FOUND)
  foreach(_h ${_ipopt_hints})
    if(EXISTS "${_h}/lib/pkgconfig")
      set(ENV{PKG_CONFIG_PATH} "${_h}/lib/pkgconfig:$ENV{PKG_CONFIG_PATH}")
    endif()
  endforeach()
  pkg_check_modules(_PC_IPOPT QUIET ipopt)
endif()

find_path(IPOPT_INCLUDE_DIR
  NAMES IpIpoptApplication.hpp coin-or/IpIpoptApplication.hpp
  HINTS ${_ipopt_hints} ${_PC_IPOPT_INCLUDE_DIRS}
  PATH_SUFFIXES include include/coin-or coin-or)

find_library(IPOPT_LIBRARY
  NAMES ipopt
  HINTS ${_ipopt_hints} ${_PC_IPOPT_LIBRARY_DIRS}
  PATH_SUFFIXES lib lib64)

set(IPOPT_VERSION "${_PC_IPOPT_VERSION}")

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(IPOPT
  REQUIRED_VARS IPOPT_LIBRARY IPOPT_INCLUDE_DIR
  VERSION_VAR IPOPT_VERSION)

if(IPOPT_FOUND)
  add_library(IPOPT::ipopt UNKNOWN IMPORTED)
  set_target_properties(IPOPT::ipopt PROPERTIES
    IMPORTED_LOCATION "${IPOPT_LIBRARY}"
    INTERFACE_INCLUDE_DIRECTORIES "${IPOPT_INCLUDE_DIR}")
  # propagate any extra pkg-config link libraries (mumps, blas, gfortran, ...)
  if(_PC_IPOPT_LIBRARIES)
    set_target_properties(IPOPT::ipopt PROPERTIES
      INTERFACE_LINK_LIBRARIES "${_PC_IPOPT_LDFLAGS}")
  endif()
endif()

mark_as_advanced(IPOPT_INCLUDE_DIR IPOPT_LIBRARY)
