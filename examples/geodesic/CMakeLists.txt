cmake_minimum_required (VERSION 3.10)
project(geodesic LANGUAGES CXX)

add_executable(${PROJECT_NAME} geodesic.cxx)
target_link_libraries(geodesic PUBLIC PSOPT  )
