cmake_minimum_required (VERSION 3.10)
project(isop LANGUAGES CXX)

add_executable(${PROJECT_NAME} isoperimetric.cxx)
target_link_libraries(isop PUBLIC PSOPT  )
