cmake_minimum_required (VERSION 3.10)
project(lambert LANGUAGES CXX)

add_executable(${PROJECT_NAME} lambert.cxx)
target_link_libraries(lambert PUBLIC PSOPT  )
