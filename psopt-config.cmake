# Config file for the PSOPT package.
#
# Usage from a consumer project:
#
#     find_package(PSOPT REQUIRED)
#     target_link_libraries(myprog PRIVATE PSOPT)
#
# The exported target's link interface names IMPORTED targets that belong to
# PSOPT's own dependencies -- Eigen3::Eigen and PkgConfig::ipopt.  An IMPORTED
# target does not survive into the consumer's scope on its own, so those
# dependencies have to be found again here, BEFORE the exported target file is
# included.  Without this the consumer's find_package(PSOPT) fails at configure
# time with "the link interface of target PSOPT contains Eigen3::Eigen but the
# target was not found", which is what it did until September 2026: nothing had
# ever consumed the installed package, every example being built inside the
# source tree.
#
# It also defines, for projects that prefer variables to targets:
#   PSOPT_INCLUDE_DIRS   include directories for PSOPT
#   PSOPT_LIBRARIES      the target to link against

include(CMakeFindDependencyMacro)

find_dependency(Eigen3 NO_MODULE)

find_dependency(PkgConfig)
pkg_check_modules(ipopt REQUIRED IMPORTED_TARGET ipopt)

include("${CMAKE_CURRENT_LIST_DIR}/PSOPT.cmake")

get_target_property(PSOPT_INCLUDE_DIRS PSOPT INTERFACE_INCLUDE_DIRECTORIES)
set(PSOPT_LIBRARIES PSOPT)
