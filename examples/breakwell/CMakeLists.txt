cmake_minimum_required (VERSION 3.10)
project(breakwell LANGUAGES CXX)

add_executable(${PROJECT_NAME} breakwell.cxx)
target_link_libraries(breakwell PUBLIC PSOPT  )
