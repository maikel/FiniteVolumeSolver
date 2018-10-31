cmake_minimum_required(VERSION 3.10)

project(HyperbolicSolver CXX Fortran)

find_package(SAMRAI REQUIRED)
find_package(Eigen3 3.3 REQUIRED NO_MODULE)
find_package(Boost REQUIRED COMPONENTS container)

add_library(Catch2 INTERFACE)
target_include_directories(Catch2 INTERFACE third_party/Catch2/single_include)

add_library(RadauSolverFortran
  src/solver/ode_solver/decsol.f
  src/solver/ode_solver/dc_decsol.f
  src/solver/ode_solver/radau.f)

add_library(HyberbolicSolver
  src/samrai/ScopeGuard.cpp
  src/samrai/utility.cpp
  src/solver/DimensionalSplitTimeIntegrator.cpp
  src/solver/HyperbolicSystemSolver.cpp
  src/solver/GodunovSplitting.cpp
  src/solver/ode_solver/Radau.cpp
  src/solver/euler/Reactor.cpp
  src/solver/euler/mechanism/Burke2012.cpp
  src/solver/euler/IdealGas.cpp
  src/solver/euler/ForwardEulerTimeIntegrator.cpp
  src/solver/euler/HlleRiemannSolver.cpp)

target_include_directories(HyberbolicSolver PUBLIC include/)
target_include_directories(HyberbolicSolver PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})

target_link_libraries(HyberbolicSolver 
  PUBLIC 
    RadauSolverFortran 
    SAMRAI_appu SAMRAI_mesh SAMRAI_geom 
    Eigen3::Eigen 
    Boost::container)

    target_compile_features(HyberbolicSolver PUBLIC cxx_std_17)

add_subdirectory(examples)

include(CTest)
if (BUILD_TESTING)
  add_subdirectory(tests)
endif()
