cmake_minimum_required (VERSION 3.10)
project(PSOPT_SNOPT_interface_CPP_examples LANGUAGES Fortran CXX)

add_executable(sntoya sntoya.cpp)
add_dependencies(sntoya PSOPT)
target_link_libraries(sntoya PRIVATE  gfortran ${snopt7_LIBRARIES} PSOPT_SNOPT_interface)

add_executable(catmixa catmixa.cpp)
add_dependencies(catmixa PSOPT)
target_link_libraries(catmixa PRIVATE  gfortran ${snopt7_LIBRARIES} PSOPT_SNOPT_interface)

add_executable(catmixa2 catmixa2.cpp)
add_dependencies(catmixa2 PSOPT)
target_link_libraries(catmixa2 PRIVATE  gfortran ${snopt7_LIBRARIES} PSOPT_SNOPT_interface)

add_executable(catmixa_stop catmixa_stop.cpp)
add_dependencies(catmixa_stop PSOPT)
target_link_libraries(catmixa_stop PRIVATE  gfortran ${snopt7_LIBRARIES} PSOPT_SNOPT_interface)

add_executable(sntoyB sntoyB.cpp)
add_dependencies(sntoyB PSOPT)
target_link_libraries(sntoyB PRIVATE  gfortran ${snopt7_LIBRARIES} PSOPT_SNOPT_interface)

add_executable(sntoyC sntoyC.cpp)
add_dependencies(sntoyC PSOPT)
target_link_libraries(sntoyC PRIVATE  gfortran ${snopt7_LIBRARIES} PSOPT_SNOPT_interface)

add_executable(hs118 hs118.cpp)
add_dependencies(hs118 PSOPT)
target_link_libraries(hs118 PRIVATE  gfortran ${snopt7_LIBRARIES} PSOPT_SNOPT_interface)

