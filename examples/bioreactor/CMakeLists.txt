cmake_minimum_required (VERSION 3.10)
project(bioreactor LANGUAGES CXX)

add_executable(${PROJECT_NAME} bioreactor.cxx)
target_link_libraries(bioreactor PUBLIC PSOPT  )
