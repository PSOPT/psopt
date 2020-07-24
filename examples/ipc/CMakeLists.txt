cmake_minimum_required (VERSION 3.10)
project(ipc LANGUAGES CXX)

add_executable(${PROJECT_NAME} interior_point.cxx)
target_link_libraries(ipc PUBLIC PSOPT  )
