cmake_minimum_required (VERSION 3.10)
project(crane LANGUAGES CXX)

add_executable(${PROJECT_NAME} crane.cxx)
target_link_libraries(crane PUBLIC PSOPT  )
