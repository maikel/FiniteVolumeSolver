cmake_minimum_required(VERSION 3.10)

project(HyperbolicSolver_examples CXX)

find_package(benchmark)
if (TARGET benchmark::benchmark)
  add_executable(benchmark.TChem benchmark.TChem.cpp)
  target_link_libraries(benchmark.TChem HyperbolicSolver.IdealGas benchmark::benchmark)
endif()

add_executable(FlameMasterReactor_Advance FlameMasterReactor_Advance.cpp)
target_link_libraries(FlameMasterReactor_Advance PRIVATE HyperbolicSolver.IdealGas)

add_executable(FlameMasterSolver FlameMasterSolver.cpp)
target_link_libraries(FlameMasterSolver PRIVATE HyperbolicSolver.IdealGas)

add_executable(TChemSolver TChemSolver.cpp)
target_link_libraries(TChemSolver PRIVATE HyperbolicSolver.IdealGas)
