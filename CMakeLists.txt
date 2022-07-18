cmake_minimum_required(VERSION 3.13)
list(APPEND CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake")
project(FiniteVolumeSolver C CXX Fortran)

# Set a default build type if none was specified
set(default_build_type "Release")
if(EXISTS "${CMAKE_SOURCE_DIR}/.git")
  set(default_build_type "Debug")
endif()

if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
  message(STATUS "Setting build type to '${default_build_type}' as none was specified.")
  set(CMAKE_BUILD_TYPE "${default_build_type}" CACHE STRING "Choose the type of build." FORCE)
  # Set the possible values of build type for cmake-gui
  set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS
    "Debug" "Release" "MinSizeRel" "RelWithDebInfo")
else()
  message(STATUS "Setting build type to '${CMAKE_BUILD_TYPE}'.")
endif()

if (NOT CMAKE_VISIBILITY_INLINES_HIDDEN)
  set(CMAKE_VISIBILITY_INLINES_HIDDEN "hidden" CACHE STRING "Choose visibility for inline symbols" FORCE)
endif()
message(STATUS "Setting the visbility for inline symbols to '${CMAKE_VISIBILITY_INLINES_HIDDEN}'.")

# Detect if we use conan to build this project and its dependencies.
# If yes, we include the neccessary build information, i.e. include paths and
# linker options.
if (EXISTS ${CMAKE_BINARY_DIR}/conanbuildinfo.cmake)
  include(${CMAKE_BINARY_DIR}/conanbuildinfo.cmake)
  conan_basic_setup(TARGETS NO_OUTPUT_DIRS)
  set(Boost_NO_BOOST_CMAKE ON)
  set(Boost_NO_SYSTEM_PATHS ON)
  set(BOOST_ROOT ${CONAN_BOOST_ROOT})
endif()

find_package(range-v3 REQUIRED)
if (NOT TARGET range-v3)
  message(FATAL "Could not find range-v3 target.")
endif()

# Try to find SAMRAI (optional)
find_package(Boost REQUIRED COMPONENTS filesystem container program_options serialization log)

find_package(SAMRAI QUIET)
if (SAMRAI_FOUND)
  message(STATUS "SAMRAI version: ${SAMRAI_VERSION}")
  option(FUB_WITH_SAMRAI "Compile with SAMRAI modules." ON)
else()
  message(STATUS "SAMRAI: Not Found.")
  option(FUB_WITH_SAMRAI "Compile with SAMRAI modules." OFF)
endif()

find_package(MPI REQUIRED COMPONENTS C CXX)

# Try to find AMReX (optional)
find_package(AMReX QUIET)


if (AMReX_FOUND)
  message(STATUS "AMReX version: ${AMReX_VERSION}")
  option(FUB_WITH_AMREX "Compile with AMReX modules." ON)
else()
  message(STATUS "AMReX: Not Found.")
  option(FUB_WITH_AMREX "Compile with AMReX modules." OFF)
endif()

if (TARGET CONAN_PKG::Eigen3)
  add_library(Eigen3::Eigen IMPORTED INTERFACE)
  target_link_libraries(Eigen3::Eigen INTERFACE CONAN_PKG::Eigen3)
elseif (TARGET CONAN_PKG::eigen)
  add_library(Eigen3::Eigen IMPORTED INTERFACE)
  target_link_libraries(Eigen3::Eigen INTERFACE CONAN_PKG::eigen)
else()
  find_package(Eigen3 3.3 REQUIRED NO_MODULE)
endif()

# Attempt to find Sundials (which contains an optional ODE solver)
find_package(Sundials QUIET)

if (TARGET Sundials::cvode)
  option(FUB_WITH_SUNDIALS "Compile with support for CVode" OFF)
else()
  option(FUB_WITH_SUNDIALS "Compile with support for CVode" OFF)
endif()

# Attempt to find Intel TBB (optional multi threading support)
find_package(TBB QUIET)

if (TARGET TBB::tbb)
  option(FUB_WITH_TBB "Use TBB for shared memory parallelism." ON)
else()
  option(FUB_WITH_TBB "Use TBB for shared memory parallelism." OFF)
endif()

# TChem's source code is included in our product.
add_subdirectory(third_party/TChem2)

# Check for OpenMP (optional multi threading support)
include(CheckFortranCompilerFlag)
check_fortran_compiler_flag("-fopenmp" FUB_WITH_FORTRAN_OPENMP_SUPPORT)
include(CheckCXXCompilerFlag)
check_cxx_compiler_flag("-fopenmp" FUB_WITH_CXX_OPENMP_SUPPORT)

if (FUB_WITH_CXX_OPENMP_SUPPORT)
  option(FUB_WITH_OPENMP "Fallback to OpenMP for shared memory parallelism." ON)
else()
  option(FUB_WITH_OPENMP "Fallback to OpenMP for shared memory parallelism." OFF)
endif()

set(FUB_DEFAULT_CHUNK_SIZE "8" CACHE STRING "The Array width for vectorization routines")

# In a multi-threaded environment we have to specify '-fopenmp' else we risk
# data races in the fortran implementation of RADAU.
add_library(RadauSolverFortran
  src/ode_solver/decsol.f
  src/ode_solver/dc_decsol.f
  src/ode_solver/radau.f)
if (FUB_WITH_FORTRAN_OPENMP_SUPPORT)
 target_compile_options(RadauSolverFortran PRIVATE "-fopenmp")
endif()

install(
  TARGETS RadauSolverFortran
  EXPORT FiniteVolumeSolver
  ARCHIVE DESTINATION ${CMAKE_INSTALL_PREFIX}/lib
  LIBRARY DESTINATION ${CMAKE_INSTALL_PREFIX}/lib
)

option(FUB_WARN_CONVERSION "Compile with warnings regarding implicit sign conversion." OFF)

# Write information about found libraries into fub/config.hpp
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/include/fub/config.hpp.in
               ${CMAKE_BINARY_DIR}/include/fub/config.hpp)

# Add a custom command that produces version.cpp, plus
# a dummy output that's not actually produced, in order
# to force version.cmake to always be re-run before the build
ADD_CUSTOM_COMMAND(
    OUTPUT ${CMAKE_CURRENT_BINARY_DIR}/version.cpp
           ${CMAKE_CURRENT_BINARY_DIR}/_version.cpp
    COMMAND ${CMAKE_COMMAND} -P
            ${CMAKE_CURRENT_SOURCE_DIR}/cmake/version.cmake)
################################################################################
# [FiniteVolumeSolver.Core]

add_library(FiniteVolumeSolver.Core INTERFACE)
target_include_directories(FiniteVolumeSolver.Core
  INTERFACE
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
    $<BUILD_INTERFACE:${CMAKE_BINARY_DIR}/include>
    $<INSTALL_INTERFACE:${CMAKE_INSTALL_PREFIX}/include>)
target_link_libraries(FiniteVolumeSolver.Core INTERFACE Boost::boost)
target_compile_features(FiniteVolumeSolver.Core INTERFACE cxx_std_17)

install(
  TARGETS FiniteVolumeSolver.Core
  EXPORT FiniteVolumeSolver
  ARCHIVE DESTINATION ${CMAKE_INSTALL_PREFIX}/lib
  LIBRARY DESTINATION ${CMAKE_INSTALL_PREFIX}/lib
)

################################################################################
# [FiniteVolumeSolver.OdeSolver]

add_library(FiniteVolumeSolver.OdeSolver
  include/fub/ode_solver/OdeSolver.hpp
  include/fub/ode_solver/RadauSolver.hpp
  include/fub/ode_solver/OdeSolverFactory.hpp
  src/ode_solver/RadauSolver.cpp
  src/ode_solver/OdeSolverFactory.cpp)

# If we found Sundials::cvode then we include its sources to the target.
if (Sundials_FOUND AND FUB_WITH_SUNDIALS)
  target_sources(FiniteVolumeSolver.OdeSolver
    PRIVATE
      src/ode_solver/CVodeSolver.cpp)
  target_link_libraries(FiniteVolumeSolver.OdeSolver PUBLIC Sundials::cvode)
endif()

target_link_libraries(FiniteVolumeSolver.OdeSolver
  PUBLIC
    FiniteVolumeSolver.Core
    RadauSolverFortran)
target_include_directories(FiniteVolumeSolver.OdeSolver
  PRIVATE
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>)
install(
  TARGETS FiniteVolumeSolver.OdeSolver
  EXPORT FiniteVolumeSolver
  ARCHIVE DESTINATION ${CMAKE_INSTALL_PREFIX}/lib
  LIBRARY DESTINATION ${CMAKE_INSTALL_PREFIX}/lib
)

find_package(fmt REQUIRED)
find_package(Vc REQUIRED)
find_package(HDF5 REQUIRED COMPONENTS C)
find_package(Python COMPONENTS Interpreter Development REQUIRED)
if (TARGET CONAN_PKG::pybind11)
  message(STATUS "Use pybind11 from conan repository")
  message(STATUS "pybind11 configs: ${CONAN_BUILD_MODULES_PATHS_PYBIND11}")
  find_package(pybind11 CONFIGS "pybind11Install.cmake" REQUIRED)
else()
  find_package(pybind11 REQUIRED)
endif()

add_library(FiniteVolumeSolver.Solver
  include/fub/boundary_condition/TransmissiveBoundary.hpp
  include/fub/boundary_condition/ReflectiveBoundary.hpp
  include/fub/equations/ideal_gas_mix/mechanism/Burke2012.hpp
  include/fub/equations/ideal_gas_mix/FlameMasterReactor.hpp
  include/fub/equations/ideal_gas_mix/KineticSourceTerm.hpp
  include/fub/equations/Advection.hpp
  include/fub/equations/Burgers.hpp
  include/fub/equations/IdealGasMix.hpp
  include/fub/equations/PerfectGas.hpp
  include/fub/equations/PerfectGasMix.hpp
  include/fub/equations/perfect_gas/ExactRiemannSolver.hpp
  include/fub/equations/perfect_gas/EinfeldtSignalVelocities.hpp
  include/fub/equations/perfect_gas/HllemMethod.hpp
  include/fub/equations/perfect_gas/MusclHancockPrimMethod.hpp
  include/fub/equations/perfect_gas/MusclHancockMethod.hpp
  include/fub/equations/perfect_gas/ThirdOrderRungeKuttaMethod.hpp
  include/fub/equations/perfect_gas/GodunovMethod.hpp
  include/fub/equations/ShallowWater.hpp
  include/fub/ext/Eigen.hpp
  include/fub/ext/outcome.hpp
  include/fub/ext/Log.hpp
  include/fub/ext/Version.hpp
  include/fub/ext/Mpi.hpp
  include/fub/PoolAllocator.hpp
  src/PoolAllocator.cpp
  include/fub/flux_method/FluxMethod.hpp
  include/fub/flux_method/GodunovMethod.hpp
  include/fub/flux_method/HllMethod.hpp
  include/fub/flux_method/MusclHancockMethod.hpp
  include/fub/cutcell_method/KbnStabilisation.hpp
  include/fub/cutcell_method/MyStabilisation.hpp
  include/fub/cutcell_method/HGridReconstruction.hpp
  include/fub/initial_data/Constant.hpp
  include/fub/initial_data/ShockMachnumber.hpp
  include/fub/tagging_method/GradientDetector.hpp
  include/fub/tagging_method/TagBuffer.hpp
  include/fub/tagging_method/TagCutCells.hpp
  include/fub/geometry/Geometry.hpp
  include/fub/geometry/Halfspace.hpp
  include/fub/geometry/ExpandTube.hpp
  include/fub/geometry/Polygon.hpp
  include/fub/geometry/Cone.hpp
  include/fub/output/BasicOutput.hpp
  include/fub/output/OutputAtFrequencyOrInterval.hpp
  include/fub/output/AsOutput.hpp
  include/fub/AnyInitialData.hpp
  include/fub/AnyBoundaryCondition.hpp
  include/fub/AnyTaggingMethod.hpp
  include/fub/Solver.hpp
  include/fub/CompleteFromCons.hpp
  include/fub/CutCellData.hpp
  include/fub/Direction.hpp
  include/fub/Duration.hpp
  include/fub/EinfeldtSignalVelocities.hpp
  include/fub/Execution.hpp
  include/fub/ForEach.hpp
  include/fub/split_method/GodunovSplitting.hpp
  include/fub/split_method/SplittingMethod.hpp
  include/fub/split_method/StrangSplitting.hpp
  include/fub/solver/DimensionalSplitLevelIntegrator.hpp
  include/fub/solver/NoSubcycleSolver.hpp
  include/fub/solver/SolverFacade.hpp
  include/fub/solver/SplitSystemSourceLevelIntegrator.hpp
  include/fub/solver/SubcycleCoarseFirstSolver.hpp
  include/fub/solver/SubcycleFineFirstSolver.hpp
  include/fub/HyperbolicPatchIntegrator.hpp
  include/fub/HyperbolicMethod.hpp
  include/fub/NewtonIteration.hpp
  include/fub/PatchDataView.hpp
  include/fub/RunSimulation.hpp
  include/fub/State.hpp
  include/fub/StateArray.hpp
  include/fub/StateRow.hpp
  include/fub/TimeStepError.hpp
  include/fub/split_method/SplittingMethod.hpp
  include/fub/split_method/GodunovSplitting.hpp
  include/fub/split_method/StrangSplitting.hpp
  include/fub/split_method/StrangSplittingLumped.hpp
  include/fub/ext/ProgramOptions.hpp
  include/fub/ext/CopyInputFile.hpp
  include/fub/counter/Counter.hpp
  include/fub/counter/CounterRegistry.hpp
  include/fub/counter/CounterResult.hpp
  include/fub/counter/Timer.hpp
  src/ext/Log.cpp
  src/ext/Version.cpp
  src/boundary_condition/ReflectiveBoundary.cpp
  src/equations/Advection.cpp
  src/equations/advection/GodunovMethod.cpp
  src/equations/advection/MusclHancockMethod.cpp
  src/equations/Burgers.cpp
  src/equations/burgers/GodunovMethod.cpp
  src/equations/burgers/MusclHancockMethod.cpp
  src/equations/ShallowWater.cpp
  src/equations/shallow_water/GodunovMethod.cpp
  src/equations/shallow_water/HllMethod.cpp
  src/equations/shallow_water/MusclHancockMethod.cpp
  src/equations/PerfectGas.cpp
  src/equations/perfect_gas/EinfeldtSignalVelocities_1d.cpp
  src/equations/perfect_gas/EinfeldtSignalVelocities_2d.cpp
  src/equations/perfect_gas/EinfeldtSignalVelocities_3d.cpp
  src/equations/perfect_gas/ExactRiemannSolver_1d.cpp
  src/equations/perfect_gas/ExactRiemannSolver_2d.cpp
  src/equations/perfect_gas/ExactRiemannSolver_3d.cpp
  src/equations/perfect_gas/GodunovMethod_1d.cpp
  src/equations/perfect_gas/GodunovMethod_2d.cpp
  src/equations/perfect_gas/GodunovMethod_3d.cpp
  src/equations/perfect_gas/HllemMethod_1d.cpp
  src/equations/perfect_gas/HllemMethod_2d.cpp
  src/equations/perfect_gas/HllemMethod_3d.cpp
  src/equations/perfect_gas/ThirdOrderRungeKuttaMethod_1d.cpp
  src/equations/perfect_gas/MusclHancockMethod_1d.cpp
  src/equations/perfect_gas/MusclHancockMethod_2d.cpp
  src/equations/perfect_gas/MusclHancockMethod_3d.cpp
  src/equations/perfect_gas/MusclHancockPrimMethod_1d.cpp
  src/equations/perfect_gas/MusclHancockPrimMethod_2d.cpp
  src/equations/perfect_gas/MusclHancockPrimMethod_3d.cpp
  src/equations/perfect_gas/MusclHancockCharacteristicMethod_1d.cpp
  src/equations/perfect_gas/MusclHancockCharacteristicMethod_2d.cpp
  src/equations/perfect_gas/MusclHancockCharacteristicMethod_3d.cpp
  src/equations/PerfectGasMix.cpp
  src/equations/perfect_gas_mix/HllemMethod_1d.cpp
  src/equations/perfect_gas_mix/HllemMethod_2d.cpp
  src/equations/perfect_gas_mix/EinfeldtSignalVelocities_1d.cpp
  src/equations/perfect_gas_mix/EinfeldtSignalVelocities_2d.cpp
  src/equations/perfect_gas_mix/EinfeldtSignalVelocities_3d.cpp
  src/equations/IdealGasMix.cpp
  src/equations/ideal_gas_mix/EinfeldtSignalVelocities_1d.cpp
  src/equations/ideal_gas_mix/EinfeldtSignalVelocities_2d.cpp
  src/equations/ideal_gas_mix/EinfeldtSignalVelocities_3d.cpp
  src/equations/ideal_gas_mix/HlleMethod_1d.cpp
  src/equations/ideal_gas_mix/HlleMethod_2d.cpp
  src/equations/ideal_gas_mix/HlleMethod_3d.cpp
  src/equations/ideal_gas_mix/MusclHancockCharacteristicMethod_1d.cpp
  src/equations/ideal_gas_mix/MusclHancockCharacteristicMethod_2d.cpp
  src/equations/ideal_gas_mix/MusclHancockCharacteristicMethod_3d.cpp
  src/equations/ideal_gas_mix/MusclHancockPrimMethod_1d.cpp
  src/equations/ideal_gas_mix/MusclHancockPrimMethod_2d.cpp
  src/equations/ideal_gas_mix/MusclHancockPrimMethod_3d.cpp
  src/equations/ideal_gas_mix/FlameMasterReactor.cpp
  src/geometry/Cone.cpp
  src/equations/ideal_gas_mix/mechanism/Burke2012.cpp
  src/cutcell_method/KbnStabilisation.cpp
  src/cutcell_method/MyStabilisation.cpp
  src/cutcell_method/HGridReconstruction.cpp
  src/tagging_method/TagBuffer.cpp
  src/tagging_method/TagCutCells.cpp
  src/ext/Mpi.cpp
  src/split_method/GodunovSplitting.cpp
  src/split_method/StrangSplitting.cpp
  src/split_method/StrangSplittingLumped.cpp
  src/geometry/Halfspace.cpp
  src/geometry/Polygon.cpp
  src/geometry/ExpandTube.cpp
  src/CutCellData.cpp
  src/HyperbolicPatchIntegrator.cpp
  src/RunSimulation.cpp
  src/ext/ProgramOptions.cpp
  src/ext/CopyInputFile.cpp
  src/counter/Counter.cpp
  src/counter/CounterRegistry.cpp
  src/counter/Timer.cpp
  ${CMAKE_CURRENT_BINARY_DIR}/version.cpp)
target_link_libraries(FiniteVolumeSolver.Solver
  PUBLIC
    FiniteVolumeSolver.Core
    FiniteVolumeSolver.OdeSolver
    MPI::MPI_C
    fmt::fmt-header-only
    Boost::container
    Boost::log
    Boost::program_options
    Boost::serialization
    Eigen3::Eigen
    Vc::Vc
    pybind11::embed
    range-v3
    CONAN_PKG::cgal)
target_include_directories(FiniteVolumeSolver.Solver
  PRIVATE
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>)
target_compile_options(FiniteVolumeSolver.Solver
  PRIVATE
  $<$<OR:$<CXX_COMPILER_ID:AppleClang>,$<CXX_COMPILER_ID:Clang>,$<CXX_COMPILER_ID:GNU>>:-Wpedantic>
  $<$<OR:$<CXX_COMPILER_ID:AppleClang>,$<CXX_COMPILER_ID:Clang>,$<CXX_COMPILER_ID:GNU>>:-Wall>
  $<$<OR:$<CXX_COMPILER_ID:AppleClang>,$<CXX_COMPILER_ID:Clang>,$<CXX_COMPILER_ID:GNU>>:-Wextra>
  $<$<AND:$<BOOL:${FUB_WARN_CONVERSION}>,$<OR:$<CXX_COMPILER_ID:AppleClang>,$<CXX_COMPILER_ID:Clang>,$<CXX_COMPILER_ID:GNU>>>:-Wconversion>
  $<$<AND:$<BOOL:${FUB_WARN_CONVERSION}>,$<OR:$<CXX_COMPILER_ID:AppleClang>,$<CXX_COMPILER_ID:Clang>,$<CXX_COMPILER_ID:GNU>>>:-Wsign-conversion>)
set_property(TARGET FiniteVolumeSolver.Solver PROPERTY CXX_VISIBILITY_PRESET "hidden")
set_property(TARGET FiniteVolumeSolver.Solver PROPERTY VISIBILITY_INLINES_HIDDEN "hidden")
install(
  TARGETS FiniteVolumeSolver.Solver
  EXPORT FiniteVolumeSolver
  ARCHIVE DESTINATION ${CMAKE_INSTALL_PREFIX}/lib
  LIBRARY DESTINATION ${CMAKE_INSTALL_PREFIX}/lib)

################################################################################
# [FiniteVolumeSolver.AMReX]

if (AMReX_FOUND AND TARGET AMReX::amrex)
  add_library(FiniteVolumeSolver.AMReX
    include/fub/AMReX.hpp
    include/fub/AMReX/boundary_condition/BoundarySet.hpp
    include/fub/AMReX/boundary_condition/TransmissiveBoundary.hpp
    include/fub/AMReX/boundary_condition/ReflectiveBoundary.hpp
    include/fub/AMReX/boundary_condition/IsentropicPressureBoundary.hpp
    include/fub/AMReX/boundary_condition/PressureValveBoundary.hpp
    include/fub/AMReX/boundary_condition/MassflowBoundary.hpp
    include/fub/AMReX/initial_data/ConstantData.hpp
    include/fub/AMReX/AxialSourceTerm.hpp
    include/fub/AMReX/DiffusionSourceTerm.hpp
    src/AMReX/DiffusionSourceTerm.cpp
    include/fub/AMReX/GriddingAlgorithm.hpp
    include/fub/AMReX/IgniteDetonation.hpp
    include/fub/AMReX/PatchHierarchy.hpp
    include/fub/AMReX/Reconstruction.hpp
    include/fub/AMReX/TimeIntegrator.hpp
    include/fub/AMReX/IntegratorContext.hpp
    include/fub/AMReX/FluxMethodAdapter.hpp
    include/fub/AMReX/ScopeGuard.hpp
    include/fub/AMReX/ViewFArrayBox.hpp
    include/fub/AMReX/MLMG/MLNodeHelmholtz.hpp
    include/fub/AMReX/MLMG/MLNodeHelmDualCstVel.hpp
    include/fub/AMReX/MLMG/MLNodeHelmDualLinVel.hpp
    include/fub/AMReX/tagging_method/TagBuffer.hpp
    include/fub/AMReX/tagging_method/ConstantRegion.hpp
    include/fub/geometry/Geometry.hpp
    include/fub/geometry/PolymorphicGeometry.hpp
    include/fub/geometry/Union.hpp
    include/fub/AMReX/output/WriteHdf5.hpp
    include/fub/AMReX/output/WritePlotfiles.hpp
    include/fub/equations/CompressibleAdvection.hpp
    include/fub/AMReX/output/DebugOutput.hpp
    src/AMReX/output/DebugOutput.cpp
    src/equations/CompressibleAdvection.cpp
    src/AMReX/output/WriteHdf5.cpp
    src/equations/ideal_gas_mix/KineticSourceTerm.cpp
    src/equations/perfect_gas_mix/IgnitionDelayKinetics.cpp
    src/equations/perfect_gas_mix/ArrheniusKinetics.cpp
    src/equations/perfect_gas_mix/PlenaControl.cpp
    src/equations/perfect_gas_mix/FluxMethodFactory.cpp
    src/AMReX/boundary_condition/BoundarySet.cpp
    src/AMReX/boundary_condition/TransmissiveBoundary.cpp
    src/AMReX/boundary_condition/IsentropicPressureBoundary.cpp
    src/AMReX/boundary_condition/PressureValveBoundary.cpp
    src/AMReX/boundary_condition/MassflowBoundary.cpp
    src/AMReX/boundary_condition/TurbinePlenumBoundaryCondition.cpp
    src/AMReX/AxialSourceTerm.cpp
    src/AMReX/IgniteDetonation.cpp
    src/AMReX/ScopeGuard.cpp
    src/AMReX/CartesianGridGeometry.cpp
    src/AMReX/ViewFArrayBox.cpp
    src/AMReX/MultiFabUtilities.cpp
    src/AMReX/MLMG/MLNodeHelmDualCstVel.cpp
    src/AMReX/MLMG/MLNodeHelmDualCstVel_K.cpp
    src/AMReX/MLMG/MLNodeHelmDualLinVel.cpp
    src/AMReX/MLMG/MLNodeHelmDualLinVel_K.cpp
    src/AMReX/PatchHierarchy.cpp
    src/AMReX/GriddingAlgorithm.cpp
    src/AMReX/TimeIntegrator.cpp
    src/AMReX/IntegratorContext.cpp
    src/AMReX/tagging_method/TagBuffer.cpp
    src/AMReX/tagging_method/ConstantRegion.cpp
    src/geometry/Union.cpp
    src/AMReX/output/WriteHdf5Impl.hpp
    src/AMReX/output/WriteHdf5Impl.cpp)
    # src/AMReX_MultiBlockBoundary.H
    # src/AMReX_MultiBlockBoundaryImpl.H
    # src/AMReX_MultiBlockBoundary.cpp)
  target_include_directories(FiniteVolumeSolver.AMReX
    PUBLIC
      $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
      $<INSTALL_INTERFACE:${CMAKE_INSTALL_PREFIX}/include>
    PRIVATE
      ${CMAKE_CURRENT_SOURCE_DIR}
  )
  target_compile_options(FiniteVolumeSolver.AMReX
    PRIVATE
      $<$<OR:$<CXX_COMPILER_ID:AppleClang>,$<CXX_COMPILER_ID:Clang>,$<CXX_COMPILER_ID:GNU>>:-Wpedantic>
      $<$<OR:$<CXX_COMPILER_ID:AppleClang>,$<CXX_COMPILER_ID:Clang>,$<CXX_COMPILER_ID:GNU>>:-Wall>
      $<$<OR:$<CXX_COMPILER_ID:AppleClang>,$<CXX_COMPILER_ID:Clang>,$<CXX_COMPILER_ID:GNU>>:-Wextra>
      $<$<AND:$<BOOL:${FUB_WARN_CONVERSION}>,$<OR:$<CXX_COMPILER_ID:AppleClang>,$<CXX_COMPILER_ID:Clang>,$<CXX_COMPILER_ID:GNU>>>:-Wconversion>
      $<$<AND:$<BOOL:${FUB_WARN_CONVERSION}>,$<OR:$<CXX_COMPILER_ID:AppleClang>,$<CXX_COMPILER_ID:Clang>,$<CXX_COMPILER_ID:GNU>>>:-Wsign-conversion>)
  set_property(TARGET FiniteVolumeSolver.AMReX PROPERTY CXX_VISIBILITY_PRESET "hidden")
  set_property(TARGET FiniteVolumeSolver.AMReX PROPERTY VISIBILITY_INLINES_HIDDEN "hidden")
  target_link_libraries(FiniteVolumeSolver.AMReX
    PUBLIC
      FiniteVolumeSolver.Solver
      AMReX::amrex
      Boost::boost
      Boost::filesystem
      Boost::serialization
      Eigen3::Eigen
      CONAN_PKG::HDF5)
  target_include_directories(FiniteVolumeSolver.AMReX PUBLIC ${HDF5_C_INCLUDE_DIR})

  if (AMReX_2D_FOUND)
    target_sources(FiniteVolumeSolver.AMReX PRIVATE
      include/fub/AMReX/solver/BK19LevelIntegrator.hpp
      include/fub/AMReX/CompressibleAdvectionIntegratorContext.hpp
      src/AMReX/solver/BK19LevelIntegrator.cpp
      src/AMReX/CompressibleAdvectionIntegratorContext.cpp)
  endif()

  if (AMReX_EB_FOUND)
    target_sources(FiniteVolumeSolver.AMReX PRIVATE
        include/fub/AMReX/cutcell/AllRegularIndexSpace.hpp
        src/AMReX/cutcell/AllRegularIndexSpace.cpp)

    add_library(FiniteVolumeSolver.AMReX.EB
        src/geometry/ExpandTube.cpp
        include/fub/AMReX_CutCell.hpp
        include/fub/AMReX/cutcell/FluxMethod.hpp
        include/fub/AMReX/cutcell/GriddingAlgorithm.hpp
        include/fub/AMReX/cutcell/IndexSpace.hpp
        include/fub/AMReX/cutcell/PatchHierarchy.hpp
        include/fub/AMReX/cutcell/Reconstruction.hpp
        include/fub/AMReX/cutcell/IntegratorContext.hpp
        include/fub/AMReX/cutcell/TimeIntegrator.hpp
        include/fub/AMReX/cutcell/AxiSymmetricSourceTerm.hpp
        include/fub/AMReX/cutcell/initial_data/RiemannProblem.hpp
        include/fub/AMReX/cutcell/initial_data/InterpolateFrom1d.hpp
        include/fub/AMReX/cutcell/tagging_method/ConstantRegion.hpp
        include/fub/AMReX/cutcell/tagging_method/GradientDetector.hpp
        include/fub/AMReX/cutcell/tagging_method/TagBuffer.hpp
        include/fub/AMReX/cutcell/tagging_method/TagCutCells.hpp
        include/fub/AMReX/cutcell/boundary_condition/BoundarySet.hpp
        include/fub/AMReX/cutcell/boundary_condition/TransmissiveBoundary.hpp
        include/fub/AMReX/cutcell/boundary_condition/IsentropicPressureBoundary.hpp
        include/fub/AMReX/cutcell/boundary_condition/MassflowBoundary.hpp
        include/fub/AMReX/cutcell/boundary_condition/ReflectiveBoundary.hpp
        include/fub/AMReX/cutcell/boundary_condition/TurbineMassflowBoundary.hpp
        include/fub/AMReX/cutcell/boundary_condition/PressureOutflowBoundary.hpp
        include/fub/AMReX/cutcell/output/WriteHdf5.hpp
        include/fub/AMReX/cutcell/output/WritePlotfiles.hpp
        include/fub/AMReX/cutcell/output/PerfectGasProbesOutput.hpp
        include/fub/AMReX/cutcell/output/DebugOutput.hpp
        src/equations/perfect_gas_mix/CutCellMethodFactory.cpp
        src/AMReX/cutcell/output/DebugOutput.cpp
        include/fub/AMReX/multi_block/MultiBlockBoundary.hpp
        include/fub/AMReX/multi_block/MultiBlockBoundary2.hpp
        include/fub/AMReX/multi_block/MultiBlockGriddingAlgorithm.hpp
        include/fub/AMReX/multi_block/MultiBlockGriddingAlgorithm2.hpp
        include/fub/AMReX/multi_block/MultiBlockIntegratorContext.hpp
        include/fub/AMReX/multi_block/MultiBlockKineticSourceTerm.hpp
        include/fub/AMReX/multi_block/MultiBlockIgniteDetonation.hpp
        include/fub/AMReX/multi_block/MultiWriteHdf5.hpp
        include/fub/AMReX/multi_block/MultiBlockPlotfileOutput.hpp
        include/fub/AMReX/cutcell/output/LogProbes.hpp
        src/AMReX/cutcell/output/LogProbes.cpp
        src/AMReX/multi_block/MultiBlockPlotfileOutput.cpp
        src/AMReX/multi_block/MultiWriteHdf5.cpp
        src/AMReX/cutcell/AxiSymmetricSourceTerm.cpp
        src/AMReX/cutcell/IntegratorContext.cpp
        src/AMReX/cutcell/PatchHierarchy.cpp
        src/AMReX/cutcell/FillCutCellData.cpp
        src/AMReX/cutcell/GriddingAlgorithm.cpp
        src/AMReX/cutcell/TimeIntegrator.cpp
        src/AMReX/cutcell/boundary_condition/BoundarySet.cpp
        src/AMReX/cutcell/boundary_condition/TransmissiveBoundary.cpp
        src/AMReX/cutcell/boundary_condition/IsentropicPressureBoundary.cpp
        src/AMReX/cutcell/boundary_condition/MassflowBoundary.cpp
        src/AMReX/cutcell/boundary_condition/MachnumberBoundary.cpp
        src/AMReX/cutcell/boundary_condition/TurbineMassflowBoundary.cpp
        src/AMReX/cutcell/boundary_condition/PressureOutflowBoundary.cpp
        src/AMReX/cutcell/boundary_condition/PressureOutflowBoundary2.cpp
        src/AMReX/cutcell/tagging_method/ConstantRegion.cpp
        src/AMReX/cutcell/tagging_method/TagBuffer.cpp
        src/AMReX/cutcell/tagging_method/TagCutCells.cpp
        src/AMReX/cutcell/initial_data/InterpolateFrom1d.cpp
        src/AMReX/cutcell/output/WriteHdf5.cpp
        src/AMReX/cutcell/output/PerfectGasProbesOutput.cpp
        src/AMReX/multi_block/MultiBlockBoundary.cpp
        src/AMReX/multi_block/MultiBlockBoundary2.cpp
        src/AMReX/multi_block/MultiBlockGriddingAlgorithm.cpp
        src/AMReX/multi_block/MultiBlockGriddingAlgorithm2.cpp
        src/AMReX/multi_block/MultiBlockIntegratorContext.cpp
        src/AMReX/multi_block/MultiBlockIntegratorContext2.cpp
        src/AMReX/multi_block/MultiBlockKineticSourceTerm.cpp
        src/AMReX/multi_block/MultiBlockIgniteDetonation.cpp
        src/AMReX/cutcell/output/WriteHdf5Impl.hpp
        src/AMReX/cutcell/output/WriteHdf5Impl.cpp)
    target_compile_options(FiniteVolumeSolver.AMReX.EB
      PRIVATE
        $<$<OR:$<CXX_COMPILER_ID:AppleClang>,$<CXX_COMPILER_ID:Clang>,$<CXX_COMPILER_ID:GNU>>:-Wpedantic>
        $<$<OR:$<CXX_COMPILER_ID:AppleClang>,$<CXX_COMPILER_ID:Clang>,$<CXX_COMPILER_ID:GNU>>:-Wall>
        $<$<OR:$<CXX_COMPILER_ID:AppleClang>,$<CXX_COMPILER_ID:Clang>,$<CXX_COMPILER_ID:GNU>>:-Wextra>
        $<$<AND:$<BOOL:${FUB_WARN_CONVERSION}>,$<OR:$<CXX_COMPILER_ID:AppleClang>,$<CXX_COMPILER_ID:Clang>,$<CXX_COMPILER_ID:GNU>>>:-Wconversion>
        $<$<AND:$<BOOL:${FUB_WARN_CONVERSION}>,$<OR:$<CXX_COMPILER_ID:AppleClang>,$<CXX_COMPILER_ID:Clang>,$<CXX_COMPILER_ID:GNU>>>:-Wsign-conversion>)
    set_property(TARGET FiniteVolumeSolver.AMReX.EB PROPERTY CXX_VISIBILITY_PRESET "hidden")
    set_property(TARGET FiniteVolumeSolver.AMReX.EB PROPERTY VISIBILITY_INLINES_HIDDEN "hidden")
    target_link_libraries(FiniteVolumeSolver.AMReX.EB PUBLIC FiniteVolumeSolver.AMReX)
    target_include_directories(FiniteVolumeSolver.AMReX.EB PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})

    if (AMReX_2D_FOUND)
      target_sources(FiniteVolumeSolver.AMReX.EB PRIVATE
        src/equations/perfect_gas_mix/PlenaControlOutput.cpp)
    endif()
  endif()

  if (FUB_WITH_OPENMP AND AMREX_ENABLE_OMP)
    target_compile_options(FiniteVolumeSolver.AMReX PUBLIC "-fopenmp")
  endif()
endif()

################################################################################
# [FiniteVolumeSolver.SAMRAI]

if (SAMRAI_FOUND)
  add_library(FiniteVolumeSolver.SAMRAI
    include/fub/SAMRAI/PatchHierarchy.hpp
    include/fub/SAMRAI/IntegratorContext.hpp
    include/fub/SAMRAI/RegisterVariables.hpp
    include/fub/SAMRAI/ScopeGuard.hpp
    include/fub/SAMRAI/GriddingAlgorithm.hpp
    include/fub/SAMRAI/HyperbolicTimeIntegrator.hpp
    include/fub/SAMRAI/FluxMethodAdapter.hpp
    include/fub/SAMRAI/CompleteFromCons.hpp
    include/fub/SAMRAI/ViewPatch.hpp
    include/fub/SAMRAI/tagging_method/GradientDetector.hpp
    include/fub/SAMRAI/tagging_method/ConstantBox.hpp
    src/SAMRAI/ScopeGuard.cpp
    src/SAMRAI/ViewPatch.cpp
    src/SAMRAI/RegisterVariables.cpp
    src/SAMRAI/CartesianPatchHierarchy.cpp
    src/SAMRAI/PatchHierarchy.cpp
    src/SAMRAI/GriddingAlgorithm.cpp
    src/SAMRAI/IntegratorContext.cpp
    src/SAMRAI/HyperbolicTimeIntegrator.cpp
    src/SAMRAI/tagging_method/ConstantBox.cpp)
  target_include_directories(FiniteVolumeSolver.SAMRAI
    PUBLIC
      $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
      $<INSTALL_INTERFACE:${CMAKE_INSTALL_PREFIX}/include>
    PRIVATE
      ${CMAKE_CURRENT_SOURCE_DIR})
  target_link_libraries(FiniteVolumeSolver.SAMRAI
    PUBLIC
      FiniteVolumeSolver.Solver
      SAMRAI_appu SAMRAI_xfer SAMRAI_mesh SAMRAI_geom
      SAMRAI_pdat SAMRAI_hier SAMRAI_tbox range-v3)
  target_compile_options(FiniteVolumeSolver.SAMRAI
    PRIVATE
      $<$<OR:$<CXX_COMPILER_ID:AppleClang>,$<CXX_COMPILER_ID:Clang>,$<CXX_COMPILER_ID:GNU>>:-Wpedantic>
      $<$<OR:$<CXX_COMPILER_ID:AppleClang>,$<CXX_COMPILER_ID:Clang>,$<CXX_COMPILER_ID:GNU>>:-Wall>
      $<$<OR:$<CXX_COMPILER_ID:AppleClang>,$<CXX_COMPILER_ID:Clang>,$<CXX_COMPILER_ID:GNU>>:-Wextra>
      $<$<AND:$<BOOL:${FUB_WARN_CONVERSION}>,$<OR:$<CXX_COMPILER_ID:AppleClang>,$<CXX_COMPILER_ID:Clang>,$<CXX_COMPILER_ID:GNU>>>:-Wconversion>
      $<$<AND:$<BOOL:${FUB_WARN_CONVERSION}>,$<OR:$<CXX_COMPILER_ID:AppleClang>,$<CXX_COMPILER_ID:Clang>,$<CXX_COMPILER_ID:GNU>>>:-Wsign-conversion>)
  install(
    TARGETS FiniteVolumeSolver.SAMRAI
    EXPORT FiniteVolumeSolver
    ARCHIVE DESTINATION ${CMAKE_INSTALL_PREFIX}/lib
    LIBRARY DESTINATION ${CMAKE_INSTALL_PREFIX}/lib)
  install(
    DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/include/fub/SAMRAI
    DESTINATION ${CMAKE_INSTALL_PREFIX}/include/fub/
    FILES_MATCHING PATTERN *.hpp)
endif()

################################################################################
# Doxygen target

find_package(Doxygen QUIET)
if (DOXYGEN_FOUND)
    # set input and output files
    set(DOXYGEN_IN ${PROJECT_SOURCE_DIR}/docs/Doxyfile.in)
    set(DOXYGEN_OUT ${CMAKE_CURRENT_BINARY_DIR}/Doxyfile)

    # request to configure the file
    set(FUB_SOURCE_DIR "${PROJECT_SOURCE_DIR}")
    configure_file(${DOXYGEN_IN} ${DOXYGEN_OUT} @ONLY)
    message("-- Doyxgen: Generated a Doxyfile at ${DOXYGEN_OUT}.")

    add_custom_target( doc 
        COMMAND ${DOXYGEN_EXECUTABLE} ${DOXYGEN_OUT}
        WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
        COMMENT "Generating API documentation with Doxygen"
        VERBATIM )
else (DOXYGEN_FOUND)
  message("-- Doyxgen: Doxygen needs to be installed to generate the doxygen documentation")
endif (DOXYGEN_FOUND)

################################################################################
# Examples

option(FUB_EXAMPLES "Build Examples." ON)
if (FUB_EXAMPLES)
  add_subdirectory(examples)
endif()

################################################################################
#  Tests

include(CTest)
if (EXISTS "${CMAKE_SOURCE_DIR}/third_party/Catch2/.git" AND BUILD_TESTING)
  # We use Catch2 for unit tests. We include that via a git submodule.
  add_library(Catch2 INTERFACE)
  target_include_directories(Catch2 INTERFACE third_party/Catch2/single_include)
  add_subdirectory(tests)
endif()

add_library(AddOpenMPFlag INTERFACE)
if (FUB_WITH_OPENMP)
  target_compile_options(AddOpenMPFlag INTERFACE "-fopenmp")
endif()

install(
  EXPORT FiniteVolumeSolver
  NAMESPACE FiniteVolumeSolver
  DESTINATION ${CMAKE_INSTALL_PREFIX}/cmake)
