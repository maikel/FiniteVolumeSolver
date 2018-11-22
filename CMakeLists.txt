cmake_minimum_required(VERSION 3.12)
list(APPEND CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake")
project(HyperbolicSolver CXX Fortran)

find_package(SAMRAI QUIET)
find_package(Eigen3 3.3 REQUIRED NO_MODULE)
find_package(Boost 1.61 REQUIRED COMPONENTS filesystem fiber)
find_package(Sundials QUIET)

add_library(Catch2 INTERFACE)
target_include_directories(Catch2 INTERFACE third_party/Catch2/single_include)

add_subdirectory(third_party/TChem2)

add_library(RadauSolverFortran
  src/solver/ode_solver/decsol.f
  src/solver/ode_solver/dc_decsol.f
  src/solver/ode_solver/radau.f)

add_library(HyperbolicSolver.OdeSolver
  src/solver/ode_solver/Radau.cpp
  src/ode_solver/RadauSolver.cpp)
if (Sundials_FOUND AND TARGET Sundials::cvode)
  target_sources(HyperbolicSolver.OdeSolver 
    PRIVATE src/ode_solver/CVodeSolver.cpp)
  target_link_libraries(HyperbolicSolver.OdeSolver PUBLIC Sundials::cvode)
endif()

target_link_libraries(HyperbolicSolver.OdeSolver
  PUBLIC 
    RadauSolverFortran)
target_include_directories(HyperbolicSolver.OdeSolver PUBLIC include/)
target_include_directories(HyperbolicSolver.OdeSolver PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})
target_compile_features(HyperbolicSolver.OdeSolver PUBLIC cxx_std_14)

add_library(HyperbolicSolver.IdealGas
  src/ideal_gas/FlameMasterReactor.cpp
  src/ideal_gas/TChemReactor.cpp
  src/ideal_gas/mechanism/Burke2012.cpp
  src/ideal_gas/mechanism/Gri30.cpp
  src/ideal_gas/mechanism/Zhao2008Dme.cpp)
target_link_libraries(HyperbolicSolver.IdealGas 
  PUBLIC
    HyperbolicSolver.OdeSolver
    Boost::boost
    TChem)
target_include_directories(HyperbolicSolver.IdealGas
  PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})
target_compile_features(HyperbolicSolver.IdealGas PUBLIC cxx_std_14)

if (SAMRAI_FOUND)
  add_library(HyperbolicSolver.SAMRAI
    src/SAMRAI/ScopeGuard.cpp
    src/SAMRAI/utility.cpp
    src/geometry/Halfspace.cpp
    src/output/GnuplotWriter.cpp
    src/tagging/GradientDetector.cpp
    src/initial_data/RiemannProblem.cpp
    src/solver/DimensionalSplitTimeIntegrator.cpp
    src/solver/DimensionalSplitSystemSolver.cpp
    src/solver/HyperbolicSystemSourceSolver.cpp
    src/solver/SourceTermIntegrator.cpp
    src/solver/GodunovSplitting.cpp)
  target_include_directories(HyperbolicSolver.SAMRAI PUBLIC include/)
  target_include_directories(HyperbolicSolver.SAMRAI PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})

  target_link_libraries(HyperbolicSolver.SAMRAI
    PUBLIC
      SAMRAI_appu SAMRAI_mesh SAMRAI_geom 
      Boost::boost
      Boost::filesystem
      Eigen3::Eigen)
  target_compile_features(HyperbolicSolver.SAMRAI PUBLIC cxx_std_14)

  add_library(HyperbolicSolver.SAMRAI.IdealGas
    src/ideal_gas/IdealGasEquation.cpp
    src/ideal_gas/PerfectGasEquation.cpp
    src/ideal_gas/FlameMasterKinetics.cpp
    src/ideal_gas/TChemKinetics.cpp
    src/ideal_gas/ForwardEulerTimeIntegrator.cpp
    src/ideal_gas/HlleRiemannSolver.cpp
    src/ideal_gas/KineticSourceTerm.cpp
    src/ideal_gas/boundary_condition/ReflectiveCondition.cpp)
  target_link_libraries(HyperbolicSolver.SAMRAI.IdealGas 
    PUBLIC
      HyperbolicSolver.SAMRAI
      HyperbolicSolver.IdealGas)
  target_include_directories(HyperbolicSolver.SAMRAI.IdealGas
    PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})
endif()

option(FUB_EXAMPLES "Build Examples." ON)
if (FUB_EXAMPLES)
  add_subdirectory(examples)
endif()

include(CTest)
if (BUILD_TESTING)
  add_subdirectory(tests)
endif()

option(FUB_MEX_COMPONENTS "Compile and install the Matlab plugin for the ideal gas kinetics." ON)
if (FUB_MEX_COMPONENTS) 
  find_package(Matlab 
    REQUIRED COMPONENTS 
      MEX_COMPILER
      ENGINE_LIBRARY
      DATAARRAY_LIBRARY)

  matlab_add_mex(
    NAME FlamemasterKinetics_ 
    SRC src/matlab/FlamemasterKinetics_.cpp
    LINK_TO ${Matlab_DATAARRAY_LIBRARY} HyperbolicSolver.IdealGas)

  matlab_add_mex(
    NAME TChemKinetics_ 
    SRC src/matlab/TChemKinetics_.cpp
    LINK_TO ${Matlab_DATAARRAY_LIBRARY} HyperbolicSolver.IdealGas)

  install(
    TARGETS FlamemasterKinetics_ TChemKinetics_ 
    LIBRARY DESTINATION ${CMAKE_INSTALL_PREFIX}/lib)
endif()