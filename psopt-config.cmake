# - Config file for the <APP> package
# It defines the following variables
#  PSOPT_INCLUDE_DIRS - include directories for <APP>
#  PSOPT_LIBRARIES    - libraries to link against
#  PSOPT_EXECUTABLE   - the bar executable
 
# Compute paths
get_filename_component(SELF_DIR "${CMAKE_CURRENT_LIST_FILE}" PATH)
include(${SELF_DIR}/PSOPT.cmake)
