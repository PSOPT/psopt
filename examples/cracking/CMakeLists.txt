cmake_minimum_required (VERSION 3.10)
project(cracking LANGUAGES CXX)

add_executable(${PROJECT_NAME} cracking.cxx)
target_link_libraries(cracking PUBLIC PSOPT  )
