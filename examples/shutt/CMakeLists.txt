cmake_minimum_required (VERSION 3.10)
project(shutt LANGUAGES CXX)

add_executable(${PROJECT_NAME} shuttle_reentry1.cxx)
target_link_libraries(shutt PUBLIC PSOPT  )
