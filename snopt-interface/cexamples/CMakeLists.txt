cmake_minimum_required (VERSION 3.10)
project(PSOPT_SNOPT_interface_C_examples LANGUAGES C Fortran)

add_executable(toyA toyA.c)
add_dependencies(toyA PSOPT)
target_link_libraries(toyA PRIVATE gfortran ${snopt7_LIBRARIES} PSOPT_SNOPT_interface)

add_executable(toyB toyB.c)
add_dependencies(toyB PSOPT)
target_link_libraries(toyB PRIVATE gfortran ${snopt7_LIBRARIES} PSOPT_SNOPT_interface)

add_executable(toyC toyC.c)
add_dependencies(toyC PSOPT)
target_link_libraries(toyC PRIVATE gfortran ${snopt7_LIBRARIES} PSOPT_SNOPT_interface)

add_executable(toyC_w toyC_w.c)
add_dependencies(toyC_w PSOPT)
target_link_libraries(toyC_w PRIVATE gfortran ${snopt7_LIBRARIES} PSOPT_SNOPT_interface)
