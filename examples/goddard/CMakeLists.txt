cmake_minimum_required (VERSION 3.10)
project(goddard LANGUAGES CXX)

add_executable(${PROJECT_NAME} goddard.cxx)
target_link_libraries(goddard PUBLIC PSOPT  )
