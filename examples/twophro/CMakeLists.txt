cmake_minimum_required (VERSION 3.10)
project(twophro LANGUAGES CXX)

add_executable(${PROJECT_NAME} twophase_robot.cxx)
target_link_libraries(twophro PUBLIC PSOPT  )
