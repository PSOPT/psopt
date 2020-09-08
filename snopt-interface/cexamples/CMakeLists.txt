cmake_minimum_required (VERSION 3.10)
project(PSOPT_SNOPT_interface_C_examples LANGUAGES C)

add_executable(toyB toyB.c)
add_dependencies(toyB PSOPT_SNOPT_interface)
target_link_libraries(toyB PRIVATE gfortran ${snopt7_LIBRARIES} PSOPT)

add_executable(toyC toyC.c)
add_dependencies(toyC PSOPT_SNOPT_interface)
target_link_libraries(toyC PRIVATE gfortran ${snopt7_LIBRARIES} PSOPT)
