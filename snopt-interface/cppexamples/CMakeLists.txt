cmake_minimum_required (VERSION 3.10)
project(PSOPT_SNOPT_interface_examples LANGUAGES CXX)

add_executable(sntoya-jac sntoya-jac.cpp)
add_dependencies(sntoya-jac PSOPT_SNOPT_interface)
target_link_directories(sntoya-jac PUBLIC ${snopt7_LINK_DIRECTORIES})
target_link_libraries(sntoya-jac PRIVATE PSOPT_SNOPT_interface gfortran)

