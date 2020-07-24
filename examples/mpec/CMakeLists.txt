cmake_minimum_required (VERSION 3.10)
project(mpec LANGUAGES CXX)

add_executable(${PROJECT_NAME} mpec.cxx)
target_link_libraries(mpec PUBLIC PSOPT  )
