cmake_minimum_required (VERSION 3.10)
project(delay1 LANGUAGES CXX)

add_executable(${PROJECT_NAME} delay1.cxx)
target_link_libraries(delay1 PUBLIC PSOPT  )
