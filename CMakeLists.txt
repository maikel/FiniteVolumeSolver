cmake_minimum_required(VERSION 3.13)
list(APPEND CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake")
project(ReactiveEulerSolver CXX Fortran)

# Detect if we use conan to build this project and its dependencies.
# If yes, we include the neccessary build information, i.e. include paths and
# linker options.
if (EXISTS ${CMAKE_BINARY_DIR}/conanbuildinfo.cmake)
  include(${CMAKE_BINARY_DIR}/conanbuildinfo.cmake)
  conan_basic_setup()
endif()

# Attempt to find SAMRAI (optinonal)
find_package(SAMRAI QUIET)

if (SAMRAI_FOUND)
  message(STATUS "SAMRAI version: ${SAMRAI_VERSION}")
  option(FUB_WITH_SAMRAI "Compile with SAMRAI modules." ON)
  find_package(Boost REQUIRED COMPONENTS filesystem)
else()
  message(STATUS "SAMRAI: Not Found.")
  option(FUB_WITH_SAMRAI "Compile with SAMRAI modules." OFF)
  find_package(Boost REQUIRED)
endif()

# Try to find AMReX
find_package(AMReX QUIET)

if (AMReX_FOUND)
  message(STATUS "AMReX version: ${AMReX_VERSION}")
  option(FUB_WITH_AMREX "Compile with AMReX modules." ON)
else()
  message(STATUS "AMReX: Not Found.")
  option(FUB_WITH_AMREX "Compile with AMReX modules." OFF)
endif()

# We currently use Eigen3 only in one place, to compute cell coordinates. 
# Note: Consider to drop that dependency, or use it more.
find_package(Eigen3 3.3 REQUIRED NO_MODULE)

# Attempt to find Sundials (which contains an optional ODE solver)
find_package(Sundials QUIET)

if (TARGET Sundials::cvode)
  option(FUB_WITH_SUNDIALS "Compile with support for CVode" ON)
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
  EXPORT ReactiveEulerSolver
  ARCHIVE DESTINATION ${CMAKE_INSTALL_PREFIX}/lib
  LIBRARY DESTINATION ${CMAKE_INSTALL_PREFIX}/lib
)

# Write information about found libraries into fub/config.hpp
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/include/fub/config.hpp.in 
               ${CMAKE_BINARY_DIR}/include/fub/config.hpp)

################################################################################
# [ReactiveEulerSolver.OdeSolver]

add_library(ReactiveEulerSolver.OdeSolver 
  src/ode_solver/RadauSolver.cpp
  src/ode_solver/OdeSolverFactory.cpp)

# If we found Sundials::cvode then we include its sources to the target.
if (Sundials_FOUND AND FUB_WITH_SUNDIALS)
  target_sources(ReactiveEulerSolver.OdeSolver 
    PRIVATE 
      src/ode_solver/CVodeSolver.cpp)
  target_link_libraries(ReactiveEulerSolver.OdeSolver PUBLIC Sundials::cvode)
endif()

target_link_libraries(ReactiveEulerSolver.OdeSolver
  PUBLIC 
    RadauSolverFortran)
target_include_directories(ReactiveEulerSolver.OdeSolver 
  PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
    $<BUILD_INTERFACE:${CMAKE_BINARY_DIR}/include>
    $<INSTALL_INTERFACE:${CMAKE_INSTALL_PREFIX}/include>
  PRIVATE
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>)
target_compile_features(ReactiveEulerSolver.OdeSolver PUBLIC cxx_std_14)
install(
  TARGETS ReactiveEulerSolver.OdeSolver
  EXPORT ReactiveEulerSolver
  ARCHIVE DESTINATION ${CMAKE_INSTALL_PREFIX}/lib
  LIBRARY DESTINATION ${CMAKE_INSTALL_PREFIX}/lib
)

################################################################################
# [ReactiveEulerSolver.IdealGas]

add_library(ReactiveEulerSolver.IdealGas
  src/ideal_gas/FlameMasterReactor.cpp
  src/ideal_gas/TChemReactor.cpp
  src/ideal_gas/mechanism/Burke2012.cpp
  src/ideal_gas/mechanism/Gri30.cpp
  src/ideal_gas/mechanism/Zhao2008Dme.cpp
  src/ideal_gas/mechanism/AramcoMech_DMEonly_74spec.cpp
  src/ideal_gas/mechanism/AramcoMech_woC3H4_91spec.cpp)
target_link_libraries(ReactiveEulerSolver.IdealGas 
  PUBLIC
    ReactiveEulerSolver.OdeSolver
    Boost::boost
    TChem)
target_include_directories(ReactiveEulerSolver.IdealGas
  PRIVATE 
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>)
target_compile_features(ReactiveEulerSolver.IdealGas PUBLIC cxx_std_14)
install(
  TARGETS ReactiveEulerSolver.IdealGas
  EXPORT ReactiveEulerSolver
  ARCHIVE DESTINATION ${CMAKE_INSTALL_PREFIX}/lib
  LIBRARY DESTINATION ${CMAKE_INSTALL_PREFIX}/lib
)
install(
  DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/include/fub/core
  DESTINATION ${CMAKE_INSTALL_PREFIX}/include/fub/
  FILES_MATCHING PATTERN *.hpp
)
install(
  DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/include/fub/ode_solver
  DESTINATION ${CMAKE_INSTALL_PREFIX}/include/fub/
  FILES_MATCHING PATTERN *.hpp
)
install(
  DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/include/fub/ideal_gas
  DESTINATION ${CMAKE_INSTALL_PREFIX}/include/fub/
  FILES_MATCHING PATTERN *.hpp
)

################################################################################
# [ReactiveEulerSolver.SAMRAI]

if (SAMRAI_FOUND)
  add_library(ReactiveEulerSolver.SAMRAI
    src/SAMRAI/ScopeGuard.cpp
    src/SAMRAI/utility.cpp
    src/SAMRAI/DimensionalSplitTimeIntegrator.cpp
    src/SAMRAI/DimensionalSplitSystemSolver.cpp
    src/SAMRAI/DimensionalSplitSystemSourceSolver.cpp
    src/SAMRAI/DimensionalSplitBoundaryCondition.cpp
    src/SAMRAI/SourceTermIntegrator.cpp
    src/SAMRAI/GodunovSplitting.cpp
    src/SAMRAI/geometry/Halfspace.cpp
    src/SAMRAI/output/GnuplotWriter.cpp
    src/SAMRAI/tagging/GradientDetector.cpp
    src/SAMRAI/initial_data/RiemannProblem.cpp)
  target_include_directories(ReactiveEulerSolver.SAMRAI 
    PUBLIC 
      $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
      $<INSTALL_INTERFACE:${CMAKE_INSTALL_PREFIX}/include>
    PRIVATE 
      ${CMAKE_CURRENT_SOURCE_DIR}
  )
  target_link_libraries(ReactiveEulerSolver.SAMRAI
    PUBLIC
      SAMRAI_appu SAMRAI_mesh SAMRAI_geom 
      Boost::boost
      Boost::filesystem
      Eigen3::Eigen)
  target_compile_features(ReactiveEulerSolver.SAMRAI PUBLIC cxx_std_14)
  install(
    TARGETS ReactiveEulerSolver.SAMRAI
    EXPORT ReactiveEulerSolver
    ARCHIVE DESTINATION ${CMAKE_INSTALL_PREFIX}/lib
    LIBRARY DESTINATION ${CMAKE_INSTALL_PREFIX}/lib
  )
  install(
    DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/include/fub/SAMRAI
    DESTINATION ${CMAKE_INSTALL_PREFIX}/include/fub/
    FILES_MATCHING PATTERN *.hpp
  )

################################################################################
# [ReactiveEulerSolver.SAMRAI.IdealGas]

  add_library(ReactiveEulerSolver.SAMRAI.IdealGas
    src/SAMRAI/ideal_gas/IdealGasEquation.cpp
    src/SAMRAI/ideal_gas/PerfectGasEquation.cpp
    src/SAMRAI/ideal_gas/FlameMasterKinetics.cpp
    src/SAMRAI/ideal_gas/TChemKinetics.cpp
    src/SAMRAI/ideal_gas/DimensionalSplitTimeIntegrator.cpp
    src/SAMRAI/ideal_gas/InitializeFromData.cpp
    src/SAMRAI/ideal_gas/HlleGodunovMethod.cpp
    src/SAMRAI/ideal_gas/KineticSourceTerm.cpp
    src/SAMRAI/ideal_gas/KineticDriver.cpp
    src/SAMRAI/ideal_gas/boundary_condition/ReflectiveBoundary.cpp
    src/SAMRAI/ideal_gas/boundary_condition/ConstantBoundary.cpp
    src/SAMRAI/ideal_gas/boundary_condition/IsentropicExpansionBoundary.cpp
    src/SAMRAI/ideal_gas/boundary_condition/PressureValveBoundary.cpp
    src/SAMRAI/ideal_gas/boundary_condition/MassflowBoundary.cpp)
  target_link_libraries(ReactiveEulerSolver.SAMRAI.IdealGas 
    PUBLIC
      ReactiveEulerSolver.SAMRAI
      ReactiveEulerSolver.IdealGas)
  target_include_directories(ReactiveEulerSolver.SAMRAI.IdealGas
    PRIVATE 
      $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>)

  install(
    TARGETS ReactiveEulerSolver.SAMRAI.IdealGas
    EXPORT ReactiveEulerSolver
    ARCHIVE DESTINATION ${CMAKE_INSTALL_PREFIX}/lib
    LIBRARY DESTINATION ${CMAKE_INSTALL_PREFIX}/lib
  )
  install(
    DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/include/fub/solver
    DESTINATION ${CMAKE_INSTALL_PREFIX}/include/fub/
    FILES_MATCHING PATTERN *.hpp
  )
  install(
    DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/include/fub/output
    DESTINATION ${CMAKE_INSTALL_PREFIX}/include/fub/
    FILES_MATCHING PATTERN *.hpp
  )
################################################################################
# Pyhton Bindings

  option(FUB_WITH_PYTHON_MODULE "Compile and install the python bindings of this solver." ON)
  if (FUB_WITH_PYTHON_MODULE)
    add_subdirectory(third_party/pybind11)
    pybind11_add_module(PyReactiveEulerSolver src/python/PyReactiveEulerSolver.cpp)
    target_link_libraries(PyReactiveEulerSolver PUBLIC ReactiveEulerSolver.SAMRAI.IdealGas)
  endif()
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

################################################################################
# Matlab Bindings

option(FUB_MEX_COMPONENTS "Compile and install the Matlab plugin for the ideal gas kinetics." ON)
if (FUB_MEX_COMPONENTS) 
  find_package(Matlab 
    REQUIRED COMPONENTS 
      MEX_COMPILER
      DATAARRAY_LIBRARY)

  if (TBB_FOUND) 
    matlab_add_mex(
      NAME FlamemasterKinetics_ 
      SRC src/matlab/FlamemasterKinetics_.cpp
      LINK_TO ${Matlab_DATAARRAY_LIBRARY} ReactiveEulerSolver.IdealGas TBB::tbb)
  else()
    matlab_add_mex(
      NAME FlamemasterKinetics_ 
      SRC src/matlab/FlamemasterKinetics_.cpp
      LINK_TO 
        ${Matlab_DATAARRAY_LIBRARY}
        ReactiveEulerSolver.IdealGas
        AddOpenMPFlag)
  endif()

  matlab_add_mex(
    NAME TChemKinetics_ 
    SRC src/matlab/TChemKinetics_.cpp
    LINK_TO ${Matlab_DATAARRAY_LIBRARY} ReactiveEulerSolver.IdealGas)

  install(
    TARGETS 
      FlamemasterKinetics_
      TChemKinetics_ 
    LIBRARY DESTINATION
      ${CMAKE_INSTALL_PREFIX}/lib)
  install(
    FILES
      src/matlab/FlamemasterKinetics2.m
      src/matlab/TChemKinetics.m
    DESTINATION
      ${CMAKE_INSTALL_PREFIX}/lib)
endif()

# Export everything into a ReactiveEulerSolver namespace for external cmake usage.

install(
  EXPORT ReactiveEulerSolver
  NAMESPACE ReactiveEulerSolver
  DESTINATION ${CMAKE_INSTALL_PREFIX}/cmake)