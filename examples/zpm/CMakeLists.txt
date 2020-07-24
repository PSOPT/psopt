cmake_minimum_required (VERSION 3.10)
project(zpm LANGUAGES CXX)

add_executable(${PROJECT_NAME} zpm.cxx)
target_link_libraries(zpm PUBLIC PSOPT  )
