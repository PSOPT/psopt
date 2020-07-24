cmake_minimum_required (VERSION 3.10)
project(manutec LANGUAGES CXX)

add_executable(${PROJECT_NAME} manutec.cxx)
target_link_libraries(manutec PUBLIC PSOPT  )
