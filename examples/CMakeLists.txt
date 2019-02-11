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
  add_executable(SAMRAI.Advection SAMRAI/Advection.cpp)
  target_link_libraries(SAMRAI.Advection ReactiveEulerSolver.SAMRAI)

  add_executable(KineticDriver KineticDriver.cpp)
  target_link_libraries(KineticDriver PRIVATE ReactiveEulerSolver.SAMRAI.IdealGas)

  add_executable(FlameMasterSolver FlameMasterSolver.cpp)
  target_link_libraries(FlameMasterSolver PRIVATE ReactiveEulerSolver.SAMRAI.IdealGas)

  add_executable(TChemSolver TChemSolver.cpp)
  target_link_libraries(TChemSolver PRIVATE ReactiveEulerSolver.SAMRAI.IdealGas)
endif()

if (AMReX_FOUND AND TARGET AMReX::amrex) 
  add_executable(AMReX.Advection AMReX/Advection.cpp)
  target_link_libraries(AMReX.Advection ReactiveEulerSolver.AMReX)
endif()

add_executable(eigen3.tensor_concatenate Eigen/TensorConcatenate.cpp)
target_link_libraries(eigen3.tensor_concatenate Eigen3::Eigen)
target_compile_features(eigen3.tensor_concatenate PRIVATE cxx_std_14)

add_executable(equations.advection Equations/Advection.cpp)
target_link_libraries(equations.advection PUBLIC ReactiveEulerSolver.Equations benchmark::benchmark)

add_executable(equations.burgers Equations/Burgers.cpp)
target_link_libraries(equations.burgers PUBLIC ReactiveEulerSolver.Equations benchmark::benchmark)

add_executable(equations.shallow_water Equations/ShallowWater.cpp)
target_link_libraries(equations.shallow_water PUBLIC ReactiveEulerSolver.Equations benchmark::benchmark)

add_executable(equations.plot_shallow_water Equations/ShallowWater_Driver.cpp)
target_link_libraries(equations.plot_shallow_water PUBLIC ReactiveEulerSolver.Equations benchmark::benchmark)