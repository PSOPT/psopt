cmake_minimum_required (VERSION 3.10)
project(lowthr LANGUAGES CXX)

add_executable(${PROJECT_NAME} low_thrust.cxx)
target_link_libraries(lowthr PUBLIC PSOPT  )
