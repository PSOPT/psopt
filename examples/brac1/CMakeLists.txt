cmake_minimum_required (VERSION 3.10)
project(brac1 LANGUAGES CXX)

add_executable(${PROJECT_NAME} brac1.cxx)
target_link_libraries(brac1 PUBLIC PSOPT  )
