cmake_minimum_required (VERSION 3.10)
project(param2 LANGUAGES CXX)

add_executable(${PROJECT_NAME} param2.cxx)
target_link_libraries(param2 PUBLIC PSOPT  )
