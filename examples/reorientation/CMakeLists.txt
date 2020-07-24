cmake_minimum_required (VERSION 3.10)
project(reorientation LANGUAGES CXX)

add_executable(${PROJECT_NAME} reorientation.cxx)
target_link_libraries(reorientation PUBLIC PSOPT  )
