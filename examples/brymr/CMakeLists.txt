cmake_minimum_required (VERSION 3.10)
project(brymr LANGUAGES CXX)

add_executable(${PROJECT_NAME} bryson_max_range.cxx)
target_link_libraries(brymr PUBLIC PSOPT  )
