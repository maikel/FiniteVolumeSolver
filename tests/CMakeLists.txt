cmake_minimum_required(VERSION 3.10)

project(HyperbolicSolver_tests CXX)

add_executable(test.core.span test.core.span.cpp)
target_link_libraries(test.core.span Catch2)
target_include_directories(test.core.span PRIVATE ../include)
target_compile_features(test.core.span PUBLIC cxx_std_14)

add_executable(test.ParallelBuffer test.ParallelBuffer.cpp)
target_link_libraries(test.ParallelBuffer PRIVATE HyperbolicSolver Catch2)
