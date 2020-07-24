cmake_minimum_required (VERSION 3.10)
project(moon LANGUAGES CXX)

add_executable(${PROJECT_NAME} moon.cxx)
target_link_libraries(moon PUBLIC PSOPT  )
