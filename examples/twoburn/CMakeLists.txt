cmake_minimum_required (VERSION 3.10)
project(twoburn LANGUAGES CXX)

add_executable(${PROJECT_NAME} twoburn.cxx)
target_link_libraries(twoburn PUBLIC PSOPT  )
