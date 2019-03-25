cmake_minimum_required(VERSION 3.10)

project(HyperbolicSolver_examples CXX)

if (AMReX_FOUND AND AMREX_DIM EQUAL 2) 
  add_executable(AMReX.Advection2d_Godunov AMReX/Advection2d_Godunov.cpp)
  target_link_libraries(AMReX.Advection2d_Godunov HyperbolicSolver.AMReX)

  add_executable(AMReX.Advection2d_MUSCLHancock AMReX/Advection2d_MUSCLHancock.cpp)
  target_link_libraries(AMReX.Advection2d_MUSCLHancock HyperbolicSolver.AMReX)

  add_executable(AMReX.PerfectGas2d AMReX/PerfectGas2d.cpp)
  target_link_libraries(AMReX.PerfectGas2d HyperbolicSolver.AMReX)

  add_executable(AMReX.ShallowWater2d AMReX/ShallowWater2d.cpp)
  target_link_libraries(AMReX.ShallowWater2d HyperbolicSolver.AMReX)

  if (AMREX_ENABLE_EB) 
    add_executable(AMReX.EB.LinearShock2d AMReX/EB/LinearShock2d.cpp)
    target_link_libraries(AMReX.EB.LinearShock2d HyperbolicSolver.AMReX)
    target_compile_options(AMReX.EB.LinearShock2d
      PRIVATE 
      $<$<OR:$<CXX_COMPILER_ID:AppleClang>,$<CXX_COMPILER_ID:Clang>,$<CXX_COMPILER_ID:GNU>>:-Wpedantic>
      $<$<OR:$<CXX_COMPILER_ID:AppleClang>,$<CXX_COMPILER_ID:Clang>,$<CXX_COMPILER_ID:GNU>>:-Wall>
      $<$<OR:$<CXX_COMPILER_ID:AppleClang>,$<CXX_COMPILER_ID:Clang>,$<CXX_COMPILER_ID:GNU>>:-Wextra>
      $<$<OR:$<CXX_COMPILER_ID:AppleClang>,$<CXX_COMPILER_ID:Clang>,$<CXX_COMPILER_ID:GNU>>:-Wconversion>)
  endif()

  add_executable(NodalPoisson AMReX/NodalPoisson.cpp)
  target_link_libraries(NodalPoisson HyperbolicSolver.AMReX)
endif()

if (AMReX_FOUND AND AMREX_DIM EQUAL 3) 
  add_executable(AMReX.PerfectGas3d AMReX/PerfectGas3d.cpp)
  target_link_libraries(AMReX.PerfectGas3d HyperbolicSolver.AMReX)

  if (AMREX_ENABLE_EB) 
    add_executable(AMReX.EB.LinearShock3d AMReX/EB/LinearShock3d.cpp)
    target_link_libraries(AMReX.EB.LinearShock3d HyperbolicSolver.AMReX)
    target_compile_options(AMReX.EB.LinearShock3d
      PRIVATE 
      $<$<OR:$<CXX_COMPILER_ID:AppleClang>,$<CXX_COMPILER_ID:Clang>,$<CXX_COMPILER_ID:GNU>>:-Wpedantic>
      $<$<OR:$<CXX_COMPILER_ID:AppleClang>,$<CXX_COMPILER_ID:Clang>,$<CXX_COMPILER_ID:GNU>>:-Wall>
      $<$<OR:$<CXX_COMPILER_ID:AppleClang>,$<CXX_COMPILER_ID:Clang>,$<CXX_COMPILER_ID:GNU>>:-Wextra>
      $<$<OR:$<CXX_COMPILER_ID:AppleClang>,$<CXX_COMPILER_ID:Clang>,$<CXX_COMPILER_ID:GNU>>:-Wconversion>
      $<$<OR:$<CXX_COMPILER_ID:AppleClang>,$<CXX_COMPILER_ID:Clang>,$<CXX_COMPILER_ID:GNU>>:-Wsign-conversion>
      $<$<CXX_COMPILER_ID:AppleClang>:-Wno-missing-braces>)
  endif()
endif()

if (SAMRAI_FOUND)
  add_executable(SAMRAI.Advection2d SAMRAI/Advection2d.cpp)
  target_link_libraries(SAMRAI.Advection2d HyperbolicSolver.SAMRAI)

  add_executable(SAMRAI.Burgers1d SAMRAI/Burgers1d.cpp)
  target_link_libraries(SAMRAI.Burgers1d HyperbolicSolver.SAMRAI)

  add_executable(SAMRAI.ShallowWater2d SAMRAI/ShallowWater2d.cpp)
  target_link_libraries(SAMRAI.ShallowWater2d HyperbolicSolver.SAMRAI)

  add_executable(SAMRAI.PerfectGas2d SAMRAI/PerfectGas2d.cpp)
  target_link_libraries(SAMRAI.PerfectGas2d HyperbolicSolver.SAMRAI)

  add_executable(SAMRAI.PerfectGas3d SAMRAI/PerfectGas3d.cpp)
  target_link_libraries(SAMRAI.PerfectGas3d HyperbolicSolver.SAMRAI)
endif()

add_executable(small_eps small_eps.cpp)
target_link_libraries(small_eps fmt::fmt-header-only)
target_compile_features(small_eps PUBLIC cxx_std_17)