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
  set(CMAKE_BUILD_TYPE "${default_build_type}" CACHE
      STRING "Choose the type of build." FORCE)
  # Set the possible values of build type for cmake-gui
  set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS
    "Debug" "Release" "MinSizeRel" "RelWithDebInfo")
else()
  message(STATUS "Setting build type to '${CMAKE_BUILD_TYPE}'.")
endif()

# Detect if we use conan to build this project and its dependencies.
# If yes, we include the neccessary build information, i.e. include paths and
# linker options.
if (EXISTS ${CMAKE_BINARY_DIR}/conanbuildinfo.cmake)
  include(${CMAKE_BINARY_DIR}/conanbuildinfo.cmake)
  conan_basic_setup()
endif()

# Try to find SAMRAI (optional)
find_package(Boost REQUIRED COMPONENTS filesystem container)

find_package(SAMRAI QUIET)
if (SAMRAI_FOUND)
  message(STATUS "SAMRAI version: ${SAMRAI_VERSION}")
  option(FUB_WITH_SAMRAI "Compile with SAMRAI modules." ON)
else()
  message(STATUS "SAMRAI: Not Found.")
  option(FUB_WITH_SAMRAI "Compile with SAMRAI modules." OFF)
endif()

# Try to find AMReX (optional)
find_package(AMReX QUIET)

if (AMReX_FOUND)
  message(STATUS "AMReX version: ${AMReX_VERSION}")
  option(FUB_WITH_AMREX "Compile with AMReX modules." ON)
else()
  message(STATUS "AMReX: Not Found.")
  option(FUB_WITH_AMREX "Compile with AMReX modules." OFF)
endif()

find_package(Eigen3 3.3 REQUIRED NO_MODULE)

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

# Write information about found libraries into fub/config.hpp
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/include/fub/config.hpp.in
               ${CMAKE_BINARY_DIR}/include/fub/config.hpp)

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
  include/fub/equations/ShallowWater.hpp
  include/fub/ext/Eigen.hpp
  include/fub/ext/outcome.hpp
  include/fub/flux_method/FluxMethod.hpp
  include/fub/flux_method/GodunovMethod.hpp
  include/fub/flux_method/HllMethod.hpp
  include/fub/flux_method/MusclHancockMethod.hpp
  include/fub/cutcell_method/KbnStabilisation.hpp
  include/fub/initial_data/Constant.hpp
  include/fub/initial_data/ShockMachnumber.hpp
  include/fub/tagging/GradientDetector.hpp
  include/fub/tagging/TagBuffer.hpp
  include/fub/tagging/TagCutCells.hpp
  include/fub/geometry/Geometry.hpp
  include/fub/geometry/Halfspace.hpp
  include/fub/geometry/ExpandTube.hpp
  include/fub/Solver.hpp
  include/fub/CompleteFromCons.hpp
  include/fub/Direction.hpp
  include/fub/Duration.hpp
  include/fub/EinfeldtSignalVelocities.hpp
  include/fub/Execution.hpp
  include/fub/ForEach.hpp
  include/fub/HyperbolicSplitLevelIntegrator.hpp
  include/fub/HyperbolicPatchIntegrator.hpp
  include/fub/HyperbolicSplitSystemSolver.hpp
  include/fub/HyperbolicMethod.hpp
  include/fub/SplitSystemSourceSolver.hpp
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
  src/equations/perfect_gas/GodunovMethod.cpp
  src/equations/perfect_gas/MusclHancockMethod.cpp
  src/equations/IdealGasMix.cpp
  src/equations/ideal_gas_mix/HllMethod.cpp
  src/equations/ideal_gas_mix/MusclHancockMethod.cpp
  src/equations/ideal_gas_mix/FlameMasterReactor.cpp
  src/equations/ideal_gas_mix/mechanism/Burke2012.cpp
  src/tagging/TagBuffer.cpp
  src/tagging/TagCutCells.cpp
  src/CartesianCoordinates.cpp
  src/split_method/GodunovSplitting.cpp
  src/split_method/StrangSplitting.cpp
  src/geometry/Halfspace.cpp
  src/geometry/Polygon.cpp
  src/geometry/ExpandTube.cpp
  src/HyperbolicPatchIntegrator.cpp
  src/RunSimulation.cpp)
target_link_libraries(FiniteVolumeSolver.Solver
  PUBLIC
    FiniteVolumeSolver.Core
    FiniteVolumeSolver.OdeSolver
    fmt::fmt-header-only
    Boost::container
    Eigen3::Eigen
    Vc::Vc)
target_compile_options(FiniteVolumeSolver.Solver 
  PRIVATE 
  $<$<OR:$<CXX_COMPILER_ID:AppleClang>,$<CXX_COMPILER_ID:Clang>,$<CXX_COMPILER_ID:GNU>>:-Wpedantic>
  $<$<OR:$<CXX_COMPILER_ID:AppleClang>,$<CXX_COMPILER_ID:Clang>,$<CXX_COMPILER_ID:GNU>>:-Wall>
  $<$<OR:$<CXX_COMPILER_ID:AppleClang>,$<CXX_COMPILER_ID:Clang>,$<CXX_COMPILER_ID:GNU>>:-Wextra>)
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
    include/fub/AMReX/boundary_condition/IsentropicBoundary.hpp
    include/fub/AMReX/GriddingAlgorithm.hpp
          include/fub/AMReX/PatchHierarchyExample.hpp
    include/fub/AMReX/Reconstruction.hpp
    include/fub/AMReX/RunSimulation.hpp
    include/fub/AMReX/BoundaryCondition.hpp
    include/fub/AMReX/Tagging.hpp
    include/fub/AMReX/TimeIntegrator.hpp
    include/fub/AMReX/IntegratorContext.hpp
    include/fub/AMReX/InitialData.hpp
    include/fub/AMReX/FluxMethod.hpp
    include/fub/AMReX/ScopeGuard.hpp
    include/fub/AMReX/ViewFArrayBox.hpp
    include/fub/AMReX/tagging/TagBuffer.hpp
    include/fub/AMReX/tagging/ConstantRegion.hpp
    src/equations/ideal_gas_mix/KineticSourceTerm.cpp
    src/AMReX/boundary_condition/BoundarySet.cpp
    src/AMReX/boundary_condition/TransmissiveBoundary.cpp
    src/AMReX/boundary_condition/IsentropicBoundary.cpp
    src/AMReX/FluxMethod.cpp
    src/AMReX/ScopeGuard.cpp
    src/AMReX/CartesianGridGeometry.cpp
    src/AMReX/ViewFArrayBox.cpp
    src/AMReX/PatchHierarchy.cpp
    src/AMReX/GriddingAlgorithm.cpp
    src/AMReX/BoundaryCondition.cpp
    src/AMReX/TimeIntegrator.cpp
    src/AMReX/IntegratorContext.cpp
    src/AMReX/tagging/TagBuffer.cpp
    src/AMReX//tagging/ConstantRegion.cpp)
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
      $<$<OR:$<CXX_COMPILER_ID:AppleClang>,$<CXX_COMPILER_ID:Clang>,$<CXX_COMPILER_ID:GNU>>:-Wextra>)
  target_link_libraries(FiniteVolumeSolver.AMReX
    PUBLIC
      FiniteVolumeSolver.Solver
      AMReX::amrex
      Boost::boost
      Eigen3::Eigen)

  if (AMReX_EB_FOUND) 
    add_library(FiniteVolumeSolver.AMReX.EB
        src/CutCellData.cpp
        src/geometry/ExpandTube.cpp
        include/fub/AMReX_CutCell.hpp
        include/fub/AMReX/cutcell/BoundaryCondition.hpp
        include/fub/AMReX/cutcell/FluxMethod.hpp
        include/fub/AMReX/cutcell/GriddingAlgorithm.hpp
        include/fub/AMReX/cutcell/IndexSpace.hpp
        include/fub/AMReX/cutcell/InitialData.hpp
        include/fub/AMReX/cutcell/PatchHierarchy.hpp
        include/fub/AMReX/cutcell/Reconstruction.hpp
        include/fub/AMReX/cutcell/IntegratorContext.hpp
        include/fub/AMReX/cutcell/Tagging.hpp
        include/fub/AMReX/cutcell/TimeIntegrator.hpp
        include/fub/AMReX/cutcell/initial_data/RiemannProblem.hpp
        include/fub/AMReX/cutcell/tagging/ConstantRegion.hpp
        include/fub/AMReX/cutcell/tagging/GradientDetector.hpp
        include/fub/AMReX/cutcell/tagging/TagBuffer.hpp
        include/fub/AMReX/cutcell/tagging/TagCutCells.hpp
        include/fub/AMReX/cutcell/boundary_condition/BoundarySet.hpp
        include/fub/AMReX/cutcell/boundary_condition/TransmissiveBoundary.hpp
        include/fub/AMReX/multi_block/MultiBlockBoundary.hpp
        include/fub/AMReX/multi_block/MultiBlockGriddingAlgorithm.hpp
        include/fub/AMReX/multi_block/MultiBlockIntegratorContext.hpp
        include/fub/AMReX/multi_block/MultiBlockKineticSourceTerm.hpp
        src/AMReX/cutcell/BoundaryCondition.cpp
        src/AMReX/cutcell/IntegratorContext.cpp
        src/AMReX/cutcell/PatchHierarchy.cpp
        src/AMReX/cutcell/FillCutCellData.cpp
        src/AMReX/cutcell/GriddingAlgorithm.cpp
        src/AMReX/cutcell/TimeIntegrator.cpp
        src/AMReX/cutcell/boundary_condition/BoundarySet.cpp
        src/AMReX/cutcell/boundary_condition/TransmissiveBoundary.cpp
        src/AMReX/cutcell/tagging/ConstantRegion.cpp
        src/AMReX/cutcell/tagging/TagBuffer.cpp
        src/AMReX/cutcell/tagging/TagCutCells.cpp
        src/AMReX/multi_block/MultiBlockBoundary.cpp
        src/AMReX/multi_block/MultiBlockGriddingAlgorithm.cpp
        src/AMReX/multi_block/MultiBlockIntegratorContext.cpp
        src/AMReX/multi_block/MultiBlockKineticSourceTerm.cpp)
    target_compile_options(FiniteVolumeSolver.AMReX.EB
      PRIVATE
        $<$<OR:$<CXX_COMPILER_ID:AppleClang>,$<CXX_COMPILER_ID:Clang>,$<CXX_COMPILER_ID:GNU>>:-Wpedantic>
        $<$<OR:$<CXX_COMPILER_ID:AppleClang>,$<CXX_COMPILER_ID:Clang>,$<CXX_COMPILER_ID:GNU>>:-Wall>
        $<$<OR:$<CXX_COMPILER_ID:AppleClang>,$<CXX_COMPILER_ID:Clang>,$<CXX_COMPILER_ID:GNU>>:-Wextra>)
    target_link_libraries(FiniteVolumeSolver.AMReX.EB PUBLIC FiniteVolumeSolver.AMReX)
  endif()

  if (FUB_WITH_OPENMP AND AMREX_ENABLE_OMP)
    target_compile_options(FiniteVolumeSolver.AMReX PUBLIC "-fopenmp")
  endif()
endif()

################################################################################
# [FiniteVolumeSolver.SAMRAI]

if (SAMRAI_FOUND)
  find_package(HDF5 REQUIRED)

  add_library(FiniteVolumeSolver.SAMRAI
    src/SAMRAI/PatchHierarchy.cpp
    src/SAMRAI/ScopeGuard.cpp
    src/SAMRAI/RegisterVariables.cpp
    src/SAMRAI/GriddingAlgorithm.cpp
    src/SAMRAI/ViewPatch.cpp
    src/SAMRAI/tagging/ConstantBox.cpp)
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
      SAMRAI_pdat SAMRAI_hier SAMRAI_tbox)
  target_compile_options(FiniteVolumeSolver.SAMRAI 
    PRIVATE 
      $<$<OR:$<CXX_COMPILER_ID:AppleClang>,$<CXX_COMPILER_ID:Clang>,$<CXX_COMPILER_ID:GNU>>:-Wpedantic>
      $<$<OR:$<CXX_COMPILER_ID:AppleClang>,$<CXX_COMPILER_ID:Clang>,$<CXX_COMPILER_ID:GNU>>:-Wall>
      $<$<OR:$<CXX_COMPILER_ID:AppleClang>,$<CXX_COMPILER_ID:Clang>,$<CXX_COMPILER_ID:GNU>>:-Wextra>)
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
# Examples

option(FUB_EXAMPLES "Build Examples." ON)
if (FUB_EXAMPLES)
  add_subdirectory(examples)
endif()

################################################################################
#  Tests

include(CTest)
if (BUILD_TESTING)
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
