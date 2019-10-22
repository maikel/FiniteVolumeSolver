cmake_minimum_required(VERSION 3.10)

project(FiniteVolumeSolver_tests CXX)

add_executable(test.geometry.Polygon test.Geometry.Polygon)
target_link_libraries(test.geometry.Polygon FiniteVolumeSolver.Solver)

add_executable(test.AMReX.average_down test.AMReX.average_down.cpp)
target_link_libraries(test.AMReX.average_down AMReX::amrex)

if (TARGET Catch2)
  add_executable(test.FlameMasterReactor test.FlameMasterReactor.cpp)
  target_link_libraries(test.FlameMasterReactor FiniteVolumeSolver.Solver Catch2)

  add_executable(test.core.span test.core.span.cpp)
  target_link_libraries(test.core.span FiniteVolumeSolver.Core Catch2)

  add_executable(test.core.mdspan test.core.mdspan.cpp)
  target_link_libraries(test.core.mdspan FiniteVolumeSolver.Core Catch2)

  add_executable(test.PerfectGas.ExactRiemannProblem test.PerfectGas.ExactRiemannProblem.cpp)
  target_link_libraries(test.PerfectGas.ExactRiemannProblem FiniteVolumeSolver.Solver)

  add_executable(test.HyperbolicSplitPatchIntegrator test.HyperbolicSplitPatchIntegrator.cpp)
  target_link_libraries(test.HyperbolicSplitPatchIntegrator FiniteVolumeSolver.Solver Catch2)

  add_executable(test.State test.State.cpp)
  target_link_libraries(test.State FiniteVolumeSolver.Solver Catch2)

  add_executable(test.simd_FluxMethod test.simd_FluxMethod.cpp)
  target_link_libraries(test.simd_FluxMethod Catch2 FiniteVolumeSolver.AMReX)
endif()

if (AMReX_FOUND)
  add_executable(test.AMReX.n_proper test.AMReX.n_proper.cpp)
  target_link_libraries(test.AMReX.n_proper FiniteVolumeSolver.AMReX)
endif()

if (AMReX_FOUND AND AMREX_DIM EQUAL 2 AND AMREX_ENABLE_EB)
  add_executable(test.AMReX.MultiCutFab test.AMReX.MultiCutFab.cpp)
  target_link_libraries(test.AMReX.MultiCutFab FiniteVolumeSolver.AMReX)

  add_executable(test.AMReX.EB.Wedge_nan test.AMReX.EB.plane_nan.cpp)
  target_link_libraries(test.AMReX.EB.Wedge_nan FiniteVolumeSolver.AMReX)
  target_compile_options(test.AMReX.EB.Wedge_nan
      PRIVATE 
      $<$<OR:$<CXX_COMPILER_ID:AppleClang>,$<CXX_COMPILER_ID:Clang>,$<CXX_COMPILER_ID:GNU>>:-Wpedantic>
      $<$<OR:$<CXX_COMPILER_ID:AppleClang>,$<CXX_COMPILER_ID:Clang>,$<CXX_COMPILER_ID:GNU>>:-Wall>
      $<$<OR:$<CXX_COMPILER_ID:AppleClang>,$<CXX_COMPILER_ID:Clang>,$<CXX_COMPILER_ID:GNU>>:-Wextra>
      $<$<OR:$<CXX_COMPILER_ID:AppleClang>,$<CXX_COMPILER_ID:Clang>,$<CXX_COMPILER_ID:GNU>>:-Wconversion>)
endif()

if (AMReX_FOUND AND AMREX_DIM EQUAL 3 AND AMREX_ENABLE_EB)
  add_executable(test.AMReX.EB.Cylinder test.AMReX.EB.Cylinder.cpp)
  target_link_libraries(test.AMReX.EB.Cylinder FiniteVolumeSolver.AMReX)
  target_compile_options(test.AMReX.EB.Cylinder
      PRIVATE 
      $<$<OR:$<CXX_COMPILER_ID:AppleClang>,$<CXX_COMPILER_ID:Clang>,$<CXX_COMPILER_ID:GNU>>:-Wpedantic>
      $<$<OR:$<CXX_COMPILER_ID:AppleClang>,$<CXX_COMPILER_ID:Clang>,$<CXX_COMPILER_ID:GNU>>:-Wall>
      $<$<OR:$<CXX_COMPILER_ID:AppleClang>,$<CXX_COMPILER_ID:Clang>,$<CXX_COMPILER_ID:GNU>>:-Wextra>)
endif()

find_package(benchmark QUIET)
if (TARGET benchmark::benchmark)
  add_executable(benchmark.HyperbolicSplitPatchIntegrator benchmark/HyperbolicSplitPatchIntegrator.cpp)
  target_link_libraries(benchmark.HyperbolicSplitPatchIntegrator FiniteVolumeSolver.Solver benchmark::benchmark)

  add_executable(benchmark.Advection benchmark/Advection.cpp)
  target_link_libraries(benchmark.Advection FiniteVolumeSolver.Solver benchmark::benchmark)

  add_executable(benchmark.ShallowWater benchmark/ShallowWater.cpp)
  target_link_libraries(benchmark.ShallowWater FiniteVolumeSolver.Solver benchmark::benchmark)

  add_executable(benchmark.PerfectGas benchmark/PerfectGas.cpp)
  target_link_libraries(benchmark.PerfectGas FiniteVolumeSolver.Solver benchmark::benchmark)
endif()

# endif()
# add_executable(test.VariableDescription test.VariableDescription.cpp)
# target_link_libraries(test.VariableDescription FiniteVolumeSolver.Core Eigen3::Eigen)

# add_executable(test.Eigen test.Eigen.cpp)
# target_link_libraries(test.Eigen FiniteVolumeSolver.Core Eigen3::Eigen Catch2)
