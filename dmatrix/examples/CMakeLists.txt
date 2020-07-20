cmake_minimum_required (VERSION 3.10)
project(dmatrix_examples LANGUAGES CXX)

add_executable(colon "colon.cxx")
target_link_libraries(colon PUBLIC dmatrix)

add_executable(expm "expm.cxx")
target_link_libraries(expm PUBLIC dmatrix)

add_executable(rand "rand.cxx")
target_link_libraries(rand PUBLIC dmatrix)

add_executable(ricc "ricc.cxx")
target_link_libraries(ricc PUBLIC dmatrix)

add_executable(sparse "sparse.cxx")
target_link_libraries(sparse PUBLIC dmatrix)

add_executable(template "template.cxx")
target_link_libraries(template PUBLIC dmatrix)
