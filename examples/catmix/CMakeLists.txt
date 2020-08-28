cmake_minimum_required (VERSION 3.10)
project(catmix LANGUAGES CXX)

add_executable(${PROJECT_NAME} catmix.cxx)
add_dependencies(${PROJECT_NAME} PSOPT)
target_include_directories(${PROJECT_NAME} PRIVATE ${ipopt_INCLUDE_DIRS})
target_link_libraries(${PROJECT_NAME} PUBLIC PSOPT PRIVATE ${ipopt_LIBRARIES} adolc ColPack Eigen3::Eigen)

