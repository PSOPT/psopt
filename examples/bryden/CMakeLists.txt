cmake_minimum_required (VERSION 3.10)
project(bryden LANGUAGES CXX)

add_executable(${PROJECT_NAME} bryson_denham.cxx)
target_link_libraries(bryden PUBLIC PSOPT  )
