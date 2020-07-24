cmake_minimum_required (VERSION 3.10)
project(predator LANGUAGES CXX)

add_executable(${PROJECT_NAME} predator.cxx)
target_link_libraries(predator PUBLIC PSOPT  )
