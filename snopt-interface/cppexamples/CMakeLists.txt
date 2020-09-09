cmake_minimum_required (VERSION 3.10)
project(PSOPT_SNOPT_interface_CPP_examples LANGUAGES Fortran CXX)

add_executable(sntoya sntoya.cpp sntoyA.hpp)
add_dependencies(sntoya PSOPT_SNOPT_interface)
target_link_libraries(sntoya PRIVATE gfortran ${snopt7_LIBRARIES} PSOPT)

add_executable(catmixa catmixa.cpp)
add_dependencies(catmixa PSOPT_SNOPT_interface)
target_link_libraries(catmixa PRIVATE gfortran ${snopt7_LIBRARIES} PSOPT)

add_executable(sntoya-jac sntoya-jac.cpp sntoyA.hpp)
add_dependencies(sntoya-jac PSOPT_SNOPT_interface)
target_link_libraries(sntoya-jac PRIVATE gfortran ${snopt7_LIBRARIES} PSOPT)

add_executable(sntoya2 sntoya2.cpp sntoyA.hpp)
add_dependencies(sntoya2 PSOPT_SNOPT_interface)
target_link_libraries(sntoya2 PRIVATE gfortran ${snopt7_LIBRARIES} PSOPT)

add_executable(sntoyB sntoyB.cpp sntoyB.hpp)
add_dependencies(sntoyB PSOPT_SNOPT_interface)
target_link_libraries(sntoyB PRIVATE gfortran ${snopt7_LIBRARIES} PSOPT)

add_executable(sntoyC sntoyC.cpp sntoyC.hpp)
add_dependencies(sntoyC PSOPT_SNOPT_interface)
target_link_libraries(sntoyC PRIVATE gfortran ${snopt7_LIBRARIES} PSOPT)

