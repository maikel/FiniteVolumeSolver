cmake_minimum_required(VERSION 3.10)

project(ReactiveEulerSolver_examples CXX)

find_package(benchmark QUIET)
if (FUB_WITH_SUNDIALS AND SAMRAI_FOUND AND TARGET benchmark::benchmark)
  add_executable(benchmark.TChem.Zhao benchmark.TChem.Zhao.cpp)
  target_link_libraries(benchmark.TChem.Zhao ReactiveEulerSolver.SAMRAI.IdealGas benchmark::benchmark)

  add_executable(benchmark.TChem.Aramco benchmark.TChem.Aramco.cpp)
  target_link_libraries(benchmark.TChem.Aramco ReactiveEulerSolver.SAMRAI.IdealGas benchmark::benchmark)
endif()

add_executable(FlameMasterReactor_Advance FlameMasterReactor_Advance.cpp)
target_link_libraries(FlameMasterReactor_Advance PRIVATE ReactiveEulerSolver.IdealGas)

add_executable(Jacobian Jacobian.cpp)
target_link_libraries(Jacobian PRIVATE ReactiveEulerSolver.IdealGas)

if (SAMRAI_FOUND)
  add_executable(KineticDriver KineticDriver.cpp)
  target_link_libraries(KineticDriver PRIVATE ReactiveEulerSolver.SAMRAI.IdealGas)

  add_executable(FlameMasterSolver FlameMasterSolver.cpp)
  target_link_libraries(FlameMasterSolver PRIVATE ReactiveEulerSolver.SAMRAI.IdealGas)

  add_executable(TChemSolver TChemSolver.cpp)
  target_link_libraries(TChemSolver PRIVATE ReactiveEulerSolver.SAMRAI.IdealGas)
endif()

if (AMReX_FOUND AND TARGET AMReX::amrex) 
  add_executable(amrex.hello_world AMReX/hello_world.cpp)
  target_link_libraries(amrex.hello_world AMReX::amrex)
endif()