cmake_minimum_required (VERSION 3.10)
project(twophsc LANGUAGES CXX)

add_executable(${PROJECT_NAME} twophase_schwartz.cxx)
target_link_libraries(twophsc PUBLIC PSOPT  )
