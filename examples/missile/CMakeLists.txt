cmake_minimum_required (VERSION 3.10)
project(missile LANGUAGES CXX)

add_executable(${PROJECT_NAME} missile.cxx)
target_link_libraries(missile PUBLIC PSOPT  )
