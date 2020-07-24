cmake_minimum_required (VERSION 3.10)
project(notorious LANGUAGES CXX)

add_executable(${PROJECT_NAME} notorious.cxx)
target_link_libraries(notorious PUBLIC PSOPT  )
