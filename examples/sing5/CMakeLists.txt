cmake_minimum_required (VERSION 3.10)
project(sing5 LANGUAGES CXX)

add_executable(${PROJECT_NAME} singular5.cxx)
target_link_libraries(sing5 PUBLIC PSOPT  )
