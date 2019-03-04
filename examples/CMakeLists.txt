cmake_minimum_required(VERSION 3.10)

project(HyperbolicSolver_examples CXX)

if (AMReX_FOUND) 
  add_executable(AMReX.Advection2d_Godunov AMReX/Advection2d_Godunov.cpp)
  target_link_libraries(AMReX.Advection2d_Godunov HyperbolicSolver.AMReX)

  add_executable(AMReX.Advection2d_MUSCLHancock AMReX/Advection2d_MUSCLHancock.cpp)
  target_link_libraries(AMReX.Advection2d_MUSCLHancock HyperbolicSolver.AMReX)

  add_executable(AMReX.PerfectGas2d AMReX/PerfectGas2d.cpp)
  target_link_libraries(AMReX.PerfectGas2d HyperbolicSolver.AMReX)

  add_executable(AMReX.ShallowWater2d AMReX/ShallowWater2d.cpp)
  target_link_libraries(AMReX.ShallowWater2d HyperbolicSolver.AMReX)

  add_executable(NodalPoisson AMReX/NodalPoisson.cpp)
  target_link_libraries(NodalPoisson HyperbolicSolver.AMReX)
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