cmake_minimum_required(VERSION 3.10)

project(HyperbolicSolver CXX)

find_package(SAMRAI REQUIRED)
find_package(Eigen3 3.3 REQUIRED NO_MODULE)
find_package(Boost REQUIRED COMPONENTS container)

add_library(Catch2 INTERFACE)
target_include_directories(Catch2 INTERFACE third_party/Catch2/single_include)

add_library(HyberbolicSolver
  src/samrai/ScopeGuard.cpp
  src/samrai/utility.cpp
  src/solver/DimensionalSplitIntegrator.cpp
  src/solver/euler/IdealGas.cpp
  src/solver/euler/IdealGasSplitIntegrator.cpp
  src/solver/euler/hlle_signal_velocities.cpp
  src/solver/euler/hll_flux_method.cpp)
target_include_directories(HyberbolicSolver PUBLIC include/)
target_include_directories(HyberbolicSolver PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})
target_link_libraries(HyberbolicSolver PUBLIC SAMRAI_appu SAMRAI_mesh SAMRAI_geom Eigen3::Eigen Boost::container)
target_compile_features(HyberbolicSolver PUBLIC cxx_std_17)

add_subdirectory(tests)
