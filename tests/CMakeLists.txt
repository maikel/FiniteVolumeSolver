cmake_minimum_required(VERSION 3.10)

project(ReactiveEulerSolver_tests CXX)

add_executable(test.core.span test.core.span.cpp)
target_link_libraries(test.core.span Catch2)
target_include_directories(test.core.span PRIVATE ../include)
target_compile_features(test.core.span PRIVATE cxx_std_14)

add_executable(test.core.mdspan test.core.mdspan.cpp)
target_link_libraries(test.core.mdspan Catch2)
target_include_directories(test.core.mdspan PRIVATE ../include)
target_compile_features(test.core.mdspan PRIVATE cxx_std_14)

add_executable(test.SimpleGrid test.SimpleGrid.cpp)
target_link_libraries(test.SimpleGrid Catch2 Eigen3::Eigen Boost::boost)
target_include_directories(test.SimpleGrid PRIVATE ../include)
target_compile_features(test.SimpleGrid PRIVATE cxx_std_17)

if (SAMRAI_FOUND)
  add_executable(test.ParallelBuffer test.ParallelBuffer.cpp)
  target_link_libraries(test.ParallelBuffer 
    PRIVATE 
      ReactiveEulerSolver.SAMRAI 
      Catch2
  )
endif()

add_executable(test.RadauSolver test.ode_solver.Radau.cpp)
target_link_libraries(test.RadauSolver PRIVATE ReactiveEulerSolver.IdealGas)

add_executable(test.TChemReactor test.TChemReactor.cpp)
target_link_libraries(test.TChemReactor PRIVATE ReactiveEulerSolver.IdealGas)

add_executable(test.ideal_gas.IdealGasEquation 
  test.ideal_gas.IdealGasEquation.cpp)
target_link_libraries(test.ideal_gas.IdealGasEquation 
  PRIVATE
    ReactiveEulerSolver.IdealGas
    Eigen3::Eigen
    Catch2)


find_package(benchmark QUIET)
if (TARGET benchmark::benchmark)
  add_executable(benchmark.SAMRAI.HlleGodunovMethod benchmark.HlleGodunovMethod.cpp)
  target_link_libraries(benchmark.SAMRAI.HlleGodunovMethod
    PRIVATE
    ReactiveEulerSolver.SAMRAI.IdealGas
    benchmark::benchmark
  )
endif()