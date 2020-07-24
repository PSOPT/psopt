cmake_minimum_required (VERSION 3.10)
project(user LANGUAGES CXX)

add_executable(${PROJECT_NAME} user.cxx)
target_link_libraries(user PUBLIC PSOPT)
