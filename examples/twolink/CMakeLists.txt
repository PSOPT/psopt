cmake_minimum_required (VERSION 3.10)
project(twolink LANGUAGES CXX)

add_executable(${PROJECT_NAME} twolinkarm.cxx)
target_link_libraries(twolink PUBLIC PSOPT  )
