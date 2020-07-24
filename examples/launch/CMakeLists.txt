cmake_minimum_required (VERSION 3.10)
project(launch LANGUAGES CXX)

add_executable(${PROJECT_NAME} launch.cxx)
target_link_libraries(launch PUBLIC PSOPT  )
