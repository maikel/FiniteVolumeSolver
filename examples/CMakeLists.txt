cmake_minimum_required(VERSION 3.10)

project(FiniteVolumeSolver_examples CXX)

function(add_executable_with_warnings TARGET)
  add_executable(${ARGV})
  target_compile_options("${TARGET}"
    PRIVATE
      $<$<OR:$<CXX_COMPILER_ID:AppleClang>,$<CXX_COMPILER_ID:Clang>,$<CXX_COMPILER_ID:GNU>>:-Wpedantic>
      $<$<OR:$<CXX_COMPILER_ID:AppleClang>,$<CXX_COMPILER_ID:Clang>,$<CXX_COMPILER_ID:GNU>>:-Wall>
      $<$<OR:$<CXX_COMPILER_ID:AppleClang>,$<CXX_COMPILER_ID:Clang>,$<CXX_COMPILER_ID:GNU>>:-Wextra>)
endfunction()

if (AMReX_FOUND)
  add_executable_with_warnings(AMReX.PerfectGas1d AMReX/PerfectGas1d.cpp)
  target_link_libraries(AMReX.PerfectGas1d FiniteVolumeSolver.AMReX)

  add_executable_with_warnings(AMReX.IdealGasMix AMReX/IdealGasMix.cpp)
  target_link_libraries(AMReX.IdealGasMix FiniteVolumeSolver.AMReX Boost::program_options)

  add_executable_with_warnings(AMReX.AxialSourceTerm AMReX/AxialVariation.cpp)
  target_link_libraries(AMReX.AxialSourceTerm FiniteVolumeSolver.AMReX Boost::program_options)

  add_executable_with_warnings(AMReX.Advection2d_Godunov AMReX/Advection2d_Godunov.cpp)
  target_link_libraries(AMReX.Advection2d_Godunov FiniteVolumeSolver.AMReX)

  add_executable_with_warnings(AMReX.Advection2d_MUSCLHancock AMReX/Advection2d_MUSCLHancock.cpp)
  target_link_libraries(AMReX.Advection2d_MUSCLHancock FiniteVolumeSolver.AMReX)

  add_executable_with_warnings(AMReX.GreshoVortex AMReX/GreshoVortex.cpp)
  target_link_libraries(AMReX.GreshoVortex FiniteVolumeSolver.AMReX)

  add_executable_with_warnings(AMReX.ShallowWater2d AMReX/ShallowWater2d.cpp)
  target_link_libraries(AMReX.ShallowWater2d FiniteVolumeSolver.AMReX)
endif()

if (AMReX_FOUND AND AMReX_2D_FOUND)
  add_executable(NodalPoisson AMReX/NodalPoisson.cpp)
  target_link_libraries(NodalPoisson FiniteVolumeSolver.AMReX)

  if (AMReX_EB_FOUND)
    add_executable_with_warnings(AMReX.EB.Divider2D AMReX/EB/Divider2D.cpp)
    target_link_libraries(AMReX.EB.Divider2D FiniteVolumeSolver.AMReX.EB)

    add_executable_with_warnings(AMReX.EB.LinearShock2d AMReX/EB/LinearShock2d.cpp)
    target_link_libraries(AMReX.EB.LinearShock2d FiniteVolumeSolver.AMReX.EB)

    add_executable_with_warnings(AMReX.EB.Schardin AMReX/EB/Schardin.cpp)
    target_link_libraries(AMReX.EB.Schardin FiniteVolumeSolver.AMReX.EB)

    add_executable_with_warnings(AMReX.EB.Wedge AMReX/EB/Wedge.cpp)
    target_link_libraries(AMReX.EB.Wedge FiniteVolumeSolver.AMReX.EB)

    add_executable_with_warnings(AMReX.EB.Ramp AMReX/EB/Ramp.cpp)
    target_link_libraries(AMReX.EB.Ramp FiniteVolumeSolver.AMReX.EB)

    add_executable_with_warnings(AMReX.EB.Ramp2 AMReX/EB/Ramp2.cpp)
    target_link_libraries(AMReX.EB.Ramp2 FiniteVolumeSolver.AMReX.EB)

    add_executable_with_warnings(AMReX.EB.Tube AMReX/EB/Tube.cpp)
    target_link_libraries(AMReX.EB.Tube FiniteVolumeSolver.AMReX.EB)

    add_executable_with_warnings(AMReX.EB.CoupledLinearShock2d AMReX/EB/CoupledLinearShock2d.cpp)
    target_link_libraries(AMReX.EB.CoupledLinearShock2d FiniteVolumeSolver.AMReX.EB)
  endif()
endif()

if (AMReX_FOUND AND AMReX_3D_FOUND)
  if (AMReX_EB_FOUND)
    add_executable_with_warnings(AMReX.EB.LinearShock3d AMReX/EB/LinearShock3d.cpp)
    target_link_libraries(AMReX.EB.LinearShock3d FiniteVolumeSolver.AMReX.EB)

    add_executable_with_warnings(AMReX.EB.CoupledLinearShock3d AMReX/EB/CoupledLinearShock3d.cpp)
    target_link_libraries(AMReX.EB.CoupledLinearShock3d FiniteVolumeSolver.AMReX.EB)

    add_executable_with_warnings(AMReX.EB.SingleTubePlenum AMReX/EB/SingleTubePlenum.cpp)
    target_link_libraries(AMReX.EB.SingleTubePlenum FiniteVolumeSolver.AMReX.EB)

    add_executable_with_warnings(AMReX.EB.SingleTubeRing AMReX/EB/SingleTubeRing.cpp)
    target_link_libraries(AMReX.EB.SingleTubeRing FiniteVolumeSolver.AMReX.EB)

    add_executable_with_warnings(AMReX.EB.MultiTube AMReX/EB/MultiTube.cpp)
    target_link_libraries(AMReX.EB.MultiTube FiniteVolumeSolver.AMReX.EB)
    
    add_executable_with_warnings(AMReX.EB.MultiTube_Compressor AMReX/EB/MultiTube_Compressor.cpp)
    target_link_libraries(AMReX.EB.MultiTube_Compressor FiniteVolumeSolver.AMReX.EB)

    add_executable_with_warnings(AMReX.EB.Divider AMReX/EB/Divider.cpp)
    target_link_libraries(AMReX.EB.Divider FiniteVolumeSolver.AMReX.EB)
  endif()
endif()

if (SAMRAI_FOUND)
  add_executable_with_warnings(SAMRAI.Advection2d SAMRAI/Advection2d.cpp)
  target_link_libraries(SAMRAI.Advection2d FiniteVolumeSolver.SAMRAI)

  add_executable_with_warnings(SAMRAI.Burgers1d SAMRAI/Burgers1d.cpp)
  target_link_libraries(SAMRAI.Burgers1d FiniteVolumeSolver.SAMRAI)

  add_executable_with_warnings(SAMRAI.ShallowWater2d SAMRAI/ShallowWater2d.cpp)
  target_link_libraries(SAMRAI.ShallowWater2d FiniteVolumeSolver.SAMRAI)

  add_executable_with_warnings(SAMRAI.PerfectGas2d SAMRAI/PerfectGas2d.cpp)
  target_link_libraries(SAMRAI.PerfectGas2d FiniteVolumeSolver.SAMRAI)

  add_executable_with_warnings(SAMRAI.PerfectGas3d SAMRAI/PerfectGas3d.cpp)
  target_link_libraries(SAMRAI.PerfectGas3d FiniteVolumeSolver.SAMRAI)

  add_executable_with_warnings(SAMRAI.VariableDatabase SAMRAI/VariableDatabase.cpp)
  target_link_libraries(SAMRAI.VariableDatabase FiniteVolumeSolver.SAMRAI)

  add_executable_with_warnings(SAMRAI.PatchHierarchy SAMRAI/PatchHierarchyExample.cpp)
  target_link_libraries(SAMRAI.PatchHierarchy FiniteVolumeSolver.SAMRAI FiniteVolumeSolver.AMReX)
endif()

add_executable_with_warnings(small_eps small_eps.cpp)
target_link_libraries(small_eps fmt::fmt-header-only)
target_compile_features(small_eps PUBLIC cxx_std_17)
