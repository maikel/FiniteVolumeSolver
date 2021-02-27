cmake_minimum_required(VERSION 3.10)

project(FiniteVolumeSolver_examples CXX)

function(add_executable_with_warnings TARGET)
  add_executable(${ARGV})
  target_compile_options("${TARGET}"
    PRIVATE
      $<$<OR:$<CXX_COMPILER_ID:AppleClang>,$<CXX_COMPILER_ID:Clang>,$<CXX_COMPILER_ID:GNU>>:-Wpedantic>
      $<$<OR:$<CXX_COMPILER_ID:AppleClang>,$<CXX_COMPILER_ID:Clang>,$<CXX_COMPILER_ID:GNU>>:-Wall>
      $<$<OR:$<CXX_COMPILER_ID:AppleClang>,$<CXX_COMPILER_ID:Clang>,$<CXX_COMPILER_ID:GNU>>:-Wextra>)
  set_property(TARGET "${TARGET}" PROPERTY CXX_VISIBILITY_PRESET "hidden")
  set_property(TARGET "${TARGET}" PROPERTY VISIBILITY_INLINES_HIDDEN "hidden")
endfunction()

add_executable_with_warnings(ExactRiemannProblem ExactRiemannProblem.cpp)
target_link_libraries(ExactRiemannProblem FiniteVolumeSolver.Solver MPI::MPI_CXX MPI::MPI_C)

if (AMReX_FOUND)
  add_executable_with_warnings(AMReX.PerfectGas1d AMReX/PerfectGas1d.cpp)
  target_link_libraries(AMReX.PerfectGas1d FiniteVolumeSolver.AMReX)

  add_executable_with_warnings(AMReX.PerfectGas1d_Char AMReX/PerfectGas1d_Char.cpp)
  target_link_libraries(AMReX.PerfectGas1d_Char FiniteVolumeSolver.AMReX)

  add_executable_with_warnings(AMReX.SEC_ResonantDemo_1D AMReX/SEC_ResonantDemo_1D.cpp)
  target_link_libraries(AMReX.SEC_ResonantDemo_1D FiniteVolumeSolver.AMReX)

  add_executable_with_warnings(AMReX.SEC_Arrhenius AMReX/SEC_Arrhenius.cpp)
  target_link_libraries(AMReX.SEC_Arrhenius FiniteVolumeSolver.AMReX)

  add_executable_with_warnings(AMReX.SEC_Arrhenius_Deflagration AMReX/SEC_Arrhenius_Deflagration.cpp)
  target_link_libraries(AMReX.SEC_Arrhenius_Deflagration FiniteVolumeSolver.AMReX)

  add_executable_with_warnings(AMReX.PerfectGas1d_ThirdOrder AMReX/PerfectGas1d_ThirdOrder.cpp)
  target_link_libraries(AMReX.PerfectGas1d_ThirdOrder FiniteVolumeSolver.AMReX)

  add_executable_with_warnings(AMReX.IdealGasMix AMReX/IdealGasMix.cpp)
  target_link_libraries(AMReX.IdealGasMix FiniteVolumeSolver.AMReX)

  add_executable_with_warnings(AMReX.IdealGasMix_Advection AMReX/IdealGasMix_Advection.cpp)
  target_link_libraries(AMReX.IdealGasMix_Advection FiniteVolumeSolver.AMReX)

  add_executable_with_warnings(AMReX.IdealGasMix2 AMReX/IdealGasMix2.cpp)
  target_link_libraries(AMReX.IdealGasMix2 FiniteVolumeSolver.AMReX)

  add_executable_with_warnings(AMReX.AxialSourceTerm AMReX/AxialVariation.cpp)
  target_link_libraries(AMReX.AxialSourceTerm FiniteVolumeSolver.AMReX)

  add_executable_with_warnings(AMReX.Advection2d_Godunov AMReX/Advection2d_Godunov.cpp)
  target_link_libraries(AMReX.Advection2d_Godunov FiniteVolumeSolver.AMReX)

  add_executable_with_warnings(AMReX.Advection2d_MUSCLHancock AMReX/Advection2d_MUSCLHancock.cpp)
  target_link_libraries(AMReX.Advection2d_MUSCLHancock FiniteVolumeSolver.AMReX)

  add_executable_with_warnings(AMReX.GreshoVortex AMReX/GreshoVortex.cpp)
  target_link_libraries(AMReX.GreshoVortex FiniteVolumeSolver.AMReX)

  add_executable_with_warnings(AMReX.ShallowWater2d AMReX/ShallowWater2d.cpp)
  target_link_libraries(AMReX.ShallowWater2d FiniteVolumeSolver.AMReX)

  add_executable_with_warnings(AMReX.SodShockTube2d AMReX/SodShockTube2d.cpp)
  target_link_libraries(AMReX.SodShockTube2d FiniteVolumeSolver.AMReX)

  add_executable_with_warnings(AMReX.SodShockTube2d_IdealGas AMReX/SodShockTube2d_IdealGas.cpp)
  target_link_libraries(AMReX.SodShockTube2d_IdealGas FiniteVolumeSolver.AMReX)
endif()

if (AMReX_FOUND AND AMReX_2D_FOUND)

  add_executable_with_warnings(AMReX.InitialHierarchy AMReX/InitialHierarchy.cpp)
  target_link_libraries(AMReX.InitialHierarchy FiniteVolumeSolver.AMReX)

  add_executable_with_warnings(AMReX.NodalProjectionSingleLevel AMReX/NodalProjectionSingleLevel.cpp)
  target_link_libraries(AMReX.NodalProjectionSingleLevel FiniteVolumeSolver.AMReX)

  add_executable_with_warnings(AMReX.CompressibleAdvection AMReX/CompressibleAdvection.cpp)
  target_link_libraries(AMReX.CompressibleAdvection FiniteVolumeSolver.AMReX)

  add_executable_with_warnings(AMReX.BK19PseudoInc_TravellingVortex AMReX/BK19PseudoInc_TravellingVortex.cpp)
  target_link_libraries(AMReX.BK19PseudoInc_TravellingVortex FiniteVolumeSolver.AMReX)

  add_executable_with_warnings(AMReX.BK19PseudoInc_VortexCoriolis AMReX/BK19PseudoInc_VortexCoriolis.cpp)
  target_link_libraries(AMReX.BK19PseudoInc_VortexCoriolis FiniteVolumeSolver.AMReX)

  add_executable_with_warnings(AMReX.BK19Compressible_TravellingVortex AMReX/BK19Compressible_TravellingVortex.cpp)
  target_link_libraries(AMReX.BK19Compressible_TravellingVortex FiniteVolumeSolver.AMReX)

  add_executable_with_warnings(AMReX.BK19Compressible_AcousticWave AMReX/BK19Compressible_AcousticWave.cpp)
  target_link_libraries(AMReX.BK19Compressible_AcousticWave FiniteVolumeSolver.AMReX)

  add_executable_with_warnings(AMReX.BK19Compressible_VortexCoriolis AMReX/BK19Compressible_VortexCoriolis.cpp)
  target_link_libraries(AMReX.BK19Compressible_VortexCoriolis FiniteVolumeSolver.AMReX)

  if (AMReX_EB_FOUND)
    add_executable_with_warnings(AMReX.EB.ConvergentNozzle2D AMReX/EB/2D/ConvergentNozzle.cpp)
    target_link_libraries(AMReX.EB.ConvergentNozzle2D FiniteVolumeSolver.AMReX.EB)

    add_executable_with_warnings(AMReX.EB.SEC_Plenum AMReX/EB/2D/SEC_Plenum.cpp)
    target_link_libraries(AMReX.EB.SEC_Plenum FiniteVolumeSolver.AMReX.EB)

    add_executable_with_warnings(AMReX.EB.SEC_Plenum_Arrhenius AMReX/EB/2D/SEC_Plenum_Arrhenius.cpp)
    target_link_libraries(AMReX.EB.SEC_Plenum_Arrhenius FiniteVolumeSolver.AMReX.EB range-v3)

    add_executable_with_warnings(AMReX.EB.SEC_Plenum_Arrhenius_slanted AMReX/EB/2D/SEC_Plenum_Arrhenius_slanted.cpp)
    target_link_libraries(AMReX.EB.SEC_Plenum_Arrhenius_slanted FiniteVolumeSolver.AMReX.EB)

    add_executable_with_warnings(AMReX.EB.ConvergentNozzleAxiSymmetric AMReX/EB/2D/ConvergentNozzleAxiSymmetric.cpp)
    target_link_libraries(AMReX.EB.ConvergentNozzleAxiSymmetric FiniteVolumeSolver.AMReX.EB)

    add_executable_with_warnings(AMReX.EB.Divider2D AMReX/EB/2D/Divider2D.cpp)
    target_link_libraries(AMReX.EB.Divider2D FiniteVolumeSolver.AMReX.EB)

    add_executable_with_warnings(AMReX.EB.Divider2D_AxiSymmetric AMReX/EB/2D/Divider2D_AxiSymmetric.cpp)
    target_link_libraries(AMReX.EB.Divider2D_AxiSymmetric FiniteVolumeSolver.AMReX.EB)

    add_executable_with_warnings(AMReX.EB.ShockTube AMReX/EB/2D/ShockTube.cpp)
    target_link_libraries(AMReX.EB.ShockTube FiniteVolumeSolver.AMReX.EB)

    add_executable_with_warnings(AMReX.EB.LinearShock2d AMReX/EB/2D/LinearShock2d.cpp)
    target_link_libraries(AMReX.EB.LinearShock2d FiniteVolumeSolver.AMReX.EB)

    add_executable_with_warnings(AMReX.EB.LinearShock2d_PerfectGasMix AMReX/EB/2D/LinearShock2d_PerfectGasMix.cpp)
    target_link_libraries(AMReX.EB.LinearShock2d_PerfectGasMix FiniteVolumeSolver.AMReX.EB)

    add_executable_with_warnings(AMReX.EB.LinearShockAxiSymmetric AMReX/EB/2D/LinearShockAxiSymmetric.cpp)
    target_link_libraries(AMReX.EB.LinearShockAxiSymmetric FiniteVolumeSolver.AMReX.EB)

    add_executable_with_warnings(AMReX.EB.Schardin AMReX/EB/2D/Schardin.cpp)
    target_link_libraries(AMReX.EB.Schardin FiniteVolumeSolver.AMReX.EB)

    add_executable_with_warnings(AMReX.EB.Wedge AMReX/EB/2D/Wedge.cpp)
    target_link_libraries(AMReX.EB.Wedge FiniteVolumeSolver.AMReX.EB)

    add_executable_with_warnings(AMReX.EB.Ramp AMReX/EB/2D/Ramp.cpp)
    target_link_libraries(AMReX.EB.Ramp FiniteVolumeSolver.AMReX.EB)

    add_executable_with_warnings(AMReX.EB.Ramp2 AMReX/EB/2D/Ramp2.cpp)
    target_link_libraries(AMReX.EB.Ramp2 FiniteVolumeSolver.AMReX.EB)

    add_executable_with_warnings(AMReX.EB.SlopedChannel AMReX/EB/2D/SlopedChannel.cpp)
    target_link_libraries(AMReX.EB.SlopedChannel FiniteVolumeSolver.AMReX.EB)

    add_executable_with_warnings(AMReX.EB.SlopedChannel_new AMReX/EB/2D/SlopedChannel_new.cpp)
    target_link_libraries(AMReX.EB.SlopedChannel_new FiniteVolumeSolver.AMReX.EB)

    add_executable_with_warnings(AMReX.EB.Tube AMReX/EB/2D/Tube.cpp)
    target_link_libraries(AMReX.EB.Tube FiniteVolumeSolver.AMReX.EB)

    add_executable_with_warnings(AMReX.EB.CoupledLinearShock2d AMReX/EB/2D/CoupledLinearShock2d.cpp)
    target_link_libraries(AMReX.EB.CoupledLinearShock2d FiniteVolumeSolver.AMReX.EB)
  endif()
endif()

if (AMReX_FOUND AND AMReX_3D_FOUND)
  if (AMReX_EB_FOUND)
    add_executable_with_warnings(AMReX.EB.LinearShock3d AMReX/EB/3D/LinearShock3d.cpp)
    target_link_libraries(AMReX.EB.LinearShock3d FiniteVolumeSolver.AMReX.EB)

    add_executable_with_warnings(AMReX.EB.ConvergentNozzle AMReX/EB/3D/ConvergentNozzle.cpp)
    target_link_libraries(AMReX.EB.ConvergentNozzle FiniteVolumeSolver.AMReX.EB)

    add_executable_with_warnings(AMReX.EB.CoupledLinearShock3d AMReX/EB/3D/CoupledLinearShock3d.cpp)
    target_link_libraries(AMReX.EB.CoupledLinearShock3d FiniteVolumeSolver.AMReX.EB)

    add_executable_with_warnings(AMReX.EB.SingleTubePlenum AMReX/EB/3D/SingleTubePlenum.cpp)
    target_link_libraries(AMReX.EB.SingleTubePlenum FiniteVolumeSolver.AMReX.EB)

    add_executable_with_warnings(AMReX.EB.SingleTubeRing AMReX/EB/3D/SingleTubeRing.cpp)
    target_link_libraries(AMReX.EB.SingleTubeRing FiniteVolumeSolver.AMReX.EB)

    add_executable_with_warnings(AMReX.EB.MultiTube AMReX/EB/3D/MultiTube.cpp)
    target_link_libraries(AMReX.EB.MultiTube FiniteVolumeSolver.AMReX.EB)

    add_executable_with_warnings(AMReX.EB.MultiTube_blocking AMReX/EB/3D/MultiTube_blocking.cpp)
    target_link_libraries(AMReX.EB.MultiTube_blocking FiniteVolumeSolver.AMReX.EB)

    add_executable_with_warnings(AMReX.EB.MultiTube_PerfectGas AMReX/EB/3D/MultiTube_PerfectGas.cpp)
    target_link_libraries(AMReX.EB.MultiTube_PerfectGas FiniteVolumeSolver.AMReX.EB)

    #add_executable_with_warnings(AMReX.EB.MultiTube_Compressor AMReX/EB/3D/MultiTube_Compressor.cpp)
    #target_link_libraries(AMReX.EB.MultiTube_Compressor FiniteVolumeSolver.AMReX.EB)

    add_executable_with_warnings(AMReX.EB.Divider AMReX/EB/3D/Divider.cpp)
    target_link_libraries(AMReX.EB.Divider FiniteVolumeSolver.AMReX.EB)
  endif()
endif()

if (FUB_WITH_SAMRAI AND SAMRAI_FOUND)
  add_executable_with_warnings(SAMRAI.Advection2d_Godunov SAMRAI/Advection2d_Godunov.cpp)
  target_link_libraries(SAMRAI.Advection2d_Godunov FiniteVolumeSolver.SAMRAI)

#  add_executable_with_warnings(SAMRAI.Burgers1d SAMRAI/Burgers1d.cpp)
#  target_link_libraries(SAMRAI.Burgers1d FiniteVolumeSolver.SAMRAI)

#  add_executable_with_warnings(SAMRAI.ShallowWater2d SAMRAI/ShallowWater2d.cpp)
#  target_link_libraries(SAMRAI.ShallowWater2d FiniteVolumeSolver.SAMRAI)

#  add_executable_with_warnings(SAMRAI.PerfectGas2d SAMRAI/PerfectGas2d.cpp)
#  target_link_libraries(SAMRAI.PerfectGas2d FiniteVolumeSolver.SAMRAI)

#  add_executable_with_warnings(SAMRAI.PerfectGas3d SAMRAI/PerfectGas3d.cpp)
#  target_link_libraries(SAMRAI.PerfectGas3d FiniteVolumeSolver.SAMRAI)

  add_executable_with_warnings(SAMRAI.VariableDatabase SAMRAI/VariableDatabase.cpp)
  target_link_libraries(SAMRAI.VariableDatabase FiniteVolumeSolver.SAMRAI)

  add_executable_with_warnings(SAMRAI.PatchHierarchy SAMRAI/PatchHierarchyExample.cpp)
  target_link_libraries(SAMRAI.PatchHierarchy FiniteVolumeSolver.SAMRAI)
endif()

add_executable_with_warnings(small_eps small_eps.cpp)
target_link_libraries(small_eps fmt::fmt-header-only)
target_compile_features(small_eps PUBLIC cxx_std_17)
