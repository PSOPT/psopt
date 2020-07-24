cmake_minimum_required (VERSION 3.10)
project(hyper LANGUAGES CXX)

add_executable(${PROJECT_NAME} hypersensitive.cxx)
target_link_libraries(hyper PUBLIC PSOPT  )
