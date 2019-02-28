cmake_minimum_required(VERSION 3.10)

project(HyperbolicSolver_tests CXX)

add_executable(test.core.span test.core.span.cpp)
target_link_libraries(test.core.span HyperbolicSolver.Core Catch2)

add_executable(test.core.mdspan test.core.mdspan.cpp)
target_link_libraries(test.core.mdspan HyperbolicSolver.Core Catch2)

# add_executable(test.VariableDescription test.VariableDescription.cpp)
# target_link_libraries(test.VariableDescription HyperbolicSolver.Core Eigen3::Eigen)

# add_executable(test.Eigen test.Eigen.cpp)
# target_link_libraries(test.Eigen HyperbolicSolver.Core Eigen3::Eigen Catch2)