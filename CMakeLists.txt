cmake_minimum_required(VERSION 3.10)

project(HyperbolicSolver CXX Fortran)

find_package(SAMRAI REQUIRED)
find_package(Eigen3 3.3 REQUIRED NO_MODULE)
find_package(Boost 1.54 REQUIRED COMPONENTS filesystem)

add_library(Catch2 INTERFACE)
target_include_directories(Catch2 INTERFACE third_party/Catch2/single_include)

add_subdirectory(third_party/TChem2)

add_library(RadauSolverFortran
  src/solver/ode_solver/decsol.f
  src/solver/ode_solver/dc_decsol.f
  src/solver/ode_solver/radau.f)

add_library(HyperbolicSolver
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
  src/solver/GodunovSplitting.cpp
  src/solver/ode_solver/Radau.cpp)

target_include_directories(HyperbolicSolver PUBLIC include/)
target_include_directories(HyperbolicSolver PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})

target_link_libraries(HyperbolicSolver 
  PUBLIC 
    RadauSolverFortran 
    SAMRAI_appu SAMRAI_mesh SAMRAI_geom 
    Boost::boost
    Boost::filesystem
    Eigen3::Eigen)

target_compile_features(HyperbolicSolver PUBLIC cxx_std_14)


add_library(HyperbolicSolver.IdealGas
  src/ideal_gas/IdealGasEquation.cpp
  src/ideal_gas/PerfectGasEquation.cpp
  src/ideal_gas/Reactor.cpp
  src/ideal_gas/FlameMasterKinetics.cpp
  src/ideal_gas/mechanism/Burke2012.cpp
  src/ideal_gas/mechanism/Gri30.cpp
  src/ideal_gas/TChemKinetics.cpp
  src/ideal_gas/ForwardEulerTimeIntegrator.cpp
  src/ideal_gas/HlleRiemannSolver.cpp
  src/ideal_gas/KineticSourceTerm.cpp
  src/ideal_gas/boundary_condition/ReflectiveCondition.cpp)
target_link_libraries(HyperbolicSolver.IdealGas PUBLIC HyperbolicSolver TChem)
target_include_directories(HyperbolicSolver.IdealGas PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})

add_subdirectory(examples)

include(CTest)
if (BUILD_TESTING)
  add_subdirectory(tests)
endif()
