cmake_minimum_required(VERSION 3.13)
list(APPEND CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake")
project(HyperbolicSolver CXX Fortran)

# Detect if we use conan to build this project and its dependencies.
# If yes, we include the neccessary build information, i.e. include paths and
# linker options.
if (EXISTS ${CMAKE_BINARY_DIR}/conanbuildinfo.cmake)
  include(${CMAKE_BINARY_DIR}/conanbuildinfo.cmake)
  conan_basic_setup()
endif()

# Try to find SAMRAI (optional)
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
  EXPORT HyperbolicSolver
  ARCHIVE DESTINATION ${CMAKE_INSTALL_PREFIX}/lib
  LIBRARY DESTINATION ${CMAKE_INSTALL_PREFIX}/lib
)

# Write information about found libraries into fub/config.hpp
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/include/fub/config.hpp.in 
               ${CMAKE_BINARY_DIR}/include/fub/config.hpp)

################################################################################
# [HyperbolicSolver.Core]

add_library(HyperbolicSolver.Core INTERFACE)
target_include_directories(HyperbolicSolver.Core 
  INTERFACE
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
    $<BUILD_INTERFACE:${CMAKE_BINARY_DIR}/include>
    $<INSTALL_INTERFACE:${CMAKE_INSTALL_PREFIX}/include>)
target_link_libraries(HyperbolicSolver.Core INTERFACE Boost::boost)
target_compile_features(HyperbolicSolver.Core INTERFACE cxx_std_17)

install(
  TARGETS HyperbolicSolver.Core
  EXPORT HyperbolicSolver
  ARCHIVE DESTINATION ${CMAKE_INSTALL_PREFIX}/lib
  LIBRARY DESTINATION ${CMAKE_INSTALL_PREFIX}/lib
)

################################################################################
# [HyperbolicSolver.OdeSolver]

add_library(HyperbolicSolver.OdeSolver 
  src/ode_solver/RadauSolver.cpp
  src/ode_solver/OdeSolverFactory.cpp)

# If we found Sundials::cvode then we include its sources to the target.
if (Sundials_FOUND AND FUB_WITH_SUNDIALS)
  target_sources(HyperbolicSolver.OdeSolver 
    PRIVATE 
      src/ode_solver/CVodeSolver.cpp)
  target_link_libraries(HyperbolicSolver.OdeSolver PUBLIC Sundials::cvode)
endif()

target_link_libraries(HyperbolicSolver.OdeSolver
  PUBLIC 
    HyperbolicSolver.Core
    RadauSolverFortran)
target_include_directories(HyperbolicSolver.OdeSolver
  PRIVATE
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>)
install(
  TARGETS HyperbolicSolver.OdeSolver
  EXPORT HyperbolicSolver
  ARCHIVE DESTINATION ${CMAKE_INSTALL_PREFIX}/lib
  LIBRARY DESTINATION ${CMAKE_INSTALL_PREFIX}/lib
)

find_package(fmt REQUIRED)

add_library(HyperbolicSolver.Solver
  src/equations/Advection.cpp
  src/equations/Burgers.cpp
  src/equations/ShallowWater.cpp
  src/equations/PerfectGas.cpp
  src/CartesianCoordinates.cpp
  src/split_method/GodunovSplitting.cpp
  src/split_method/StrangSplitting.cpp
  src/RunSimulation.cpp)
target_link_libraries(HyperbolicSolver.Solver
  PUBLIC
    HyperbolicSolver.Core
    fmt::fmt-header-only
    Boost::boost
    Eigen3::Eigen)
install(
  TARGETS HyperbolicSolver.Solver
  EXPORT HyperbolicSolver
  ARCHIVE DESTINATION ${CMAKE_INSTALL_PREFIX}/lib
  LIBRARY DESTINATION ${CMAKE_INSTALL_PREFIX}/lib
)

################################################################################
# [HyperbolicSolver.AMReX]

if (AMReX_FOUND AND TARGET AMReX::amrex)
  add_library(HyperbolicSolver.AMReX
    src/AMReX/ScopeGuard.cpp
    src/AMReX/CartesianGridGeometry.cpp
    src/AMReX/PatchHierarchy.cpp
    src/AMReX/GriddingAlgorithm.cpp
    src/AMReX/HyperbolicSplitIntegratorContext.cpp)
  target_include_directories(HyperbolicSolver.AMReX 
    PUBLIC 
      $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
      $<INSTALL_INTERFACE:${CMAKE_INSTALL_PREFIX}/include>
    PRIVATE 
      ${CMAKE_CURRENT_SOURCE_DIR}
  )
  target_link_libraries(HyperbolicSolver.AMReX
    PUBLIC
      HyperbolicSolver.Solver
      AMReX::amrex
      Boost::boost
      Eigen3::Eigen)
endif()

################################################################################
# [HyperbolicSolver.SAMRAI]

if (SAMRAI_FOUND)
  add_library(HyperbolicSolver.SAMRAI
    src/SAMRAI/HyperbolicSplitIntegratorContext.cpp
    src/SAMRAI/ScopeGuard.cpp
    src/SAMRAI/utility.cpp
    src/SAMRAI/CartesianPatchHierarchy.cpp
    src/SAMRAI/GriddingAlgorithm.cpp
    src/SAMRAI/ViewPatch.cpp
    src/SAMRAI/RegisterVariables.cpp)
  target_include_directories(HyperbolicSolver.SAMRAI 
    PUBLIC 
      $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
      $<INSTALL_INTERFACE:${CMAKE_INSTALL_PREFIX}/include>
    PRIVATE 
      ${CMAKE_CURRENT_SOURCE_DIR}
  )
  target_link_libraries(HyperbolicSolver.SAMRAI
    PUBLIC
      HyperbolicSolver.Solver
      SAMRAI_appu SAMRAI_xfer SAMRAI_mesh SAMRAI_geom
      SAMRAI_pdat SAMRAI_hier SAMRAI_tbox)
  target_compile_features(HyperbolicSolver.SAMRAI PUBLIC cxx_std_17)
  install(
    TARGETS HyperbolicSolver.SAMRAI
    EXPORT HyperbolicSolver
    ARCHIVE DESTINATION ${CMAKE_INSTALL_PREFIX}/lib
    LIBRARY DESTINATION ${CMAKE_INSTALL_PREFIX}/lib
  )
  install(
    DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/include/fub/SAMRAI
    DESTINATION ${CMAKE_INSTALL_PREFIX}/include/fub/
    FILES_MATCHING PATTERN *.hpp
  )
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
  EXPORT HyperbolicSolver
  NAMESPACE HyperbolicSolver
  DESTINATION ${CMAKE_INSTALL_PREFIX}/cmake)
