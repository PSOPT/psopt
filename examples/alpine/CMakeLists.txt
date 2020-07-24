cmake_minimum_required (VERSION 3.10)
project(alpine LANGUAGES CXX)

add_executable(${PROJECT_NAME} alpine.cxx)
target_link_libraries(${PROJECT_NAME} PUBLIC PSOPT  )
