// Copyright (c) 2020 Maikel Nadolski
// Copyright (c) 2019 Stefan Vater
//
// Permission is hereby granted, free of charge, to any person obtaining a copy
// of this software and associated documentation files (the "Software"), to deal
// in the Software without restriction, including without limitation the rights
// to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
// copies of the Software, and to permit persons to whom the Software is
// furnished to do so, subject to the following conditions:
//
// The above copyright notice and this permission notice shall be included in
// all copies or substantial portions of the Software.
//
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
// IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
// AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
// OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
// SOFTWARE.

#ifndef FUB_AMREX_BK19_SOLVER_HPP
#define FUB_AMREX_BK19_SOLVER_HPP

#include "fub/AMReX/CompressibleAdvectionIntegratorContext.hpp"
#include "fub/equations/CompressibleAdvection.hpp"

#include "fub/AMReX/MLMG/MLNodeHelmholtz.hpp"
#include "fub/ext/ProgramOptions.hpp"
#include "fub/solver/DimensionalSplitLevelIntegrator.hpp"
#include "fub/solver/NoSubcycleSolver.hpp"
#include "fub/split_method/SplittingMethod.hpp"

#include "AMReX_MLMG.H"

namespace fub::amrex {

struct BK19SolverOptions {
  BK19SolverOptions() = default;
  BK19SolverOptions(const ProgramOptions& map);

  void Print(SeverityLogger& log);

  bool do_initial_projection = true;
  double mlmg_tolerance_rel = 1.0e-4;
  double mlmg_tolerance_abs = -1.0;
  int mlmg_max_iter = 100;
  int mlmg_verbose = 0;
  double bottom_tolerance_rel = 1.0e-4;
  double bottom_tolerance_abs = -1.0;
  int bottom_max_iter = 20;
  int bottom_verbose = 0;
  int always_use_bnorm = 0;
  std::string prefix = "BK19LevelIntegrator";
  bool output_between_steps = false;
};

struct BK19PhysicalParameters {
  BK19PhysicalParameters() = default;
  BK19PhysicalParameters(const ProgramOptions& map);

  void Print(SeverityLogger& log);

  /// Specific gas constant
  double R_gas{287.4};

  /// Heat capacity ratio
  double gamma{1.4};

  /// Heat capacity at constant pressure
  double c_p{1006.0};

  /// Gravitational acceleration
  double g{10.0}; //  [m / s^2]

  /// Coriolis parameter in beta plane
  double f{0.0};
  std::array<double, 3> k_vect{0.0, 0.0, 1.0};

  double alpha_p{1.0};
  double Msq{0.0};
};

template <int Rank, int VelocityRank = Rank> class BK19Solver {
public:
  using Equation = CompressibleAdvection<Rank, VelocityRank>;
  using IndexMapping = ::fub::IndexMapping<Equation>;

  using SplittingMethod = ::fub::AnySplitMethod;
  using IntegratorContext = CompressibleAdvectionIntegratorContext;
  using AdvectionLevelIntegrator =
      DimensionalSplitLevelIntegrator<Rank, IntegratorContext, SplittingMethod>;
  using AdvectionSolver = NoSubcycleSolver<AdvectionLevelIntegrator>;
  using LinearOperator = ::amrex::MLNodeHelmholtz;

  ////////////////////////////////////////////////////////////////////////////
  // Constructor

  BK19Solver(const Equation& equation, AdvectionSolver advection,
             std::shared_ptr<::amrex::MLNodeHelmholtz> linop,
             const BK19PhysicalParameters& physical_parameters,
             const BK19SolverOptions& options = BK19SolverOptions());

  ////////////////////////////////////////////////////////////////////////////
  // Accessors

  AdvectionSolver& GetAdvectionSolver() noexcept { return advection_; }

  const AdvectionSolver& GetAdvectionSolver() const noexcept {
    return advection_;
  }

  const BK19PhysicalParameters& GetPhysicalParameters() const noexcept {
    return physical_parameters_;
  }

  const BK19SolverOptions& GetSolverOptions() const noexcept {
    return options_;
  }

  const std::shared_ptr<LinearOperator> GetLinearOperator() const noexcept {
    return lin_op_;
  }

  const std::shared_ptr<GriddingAlgorithm>&
  GetGriddingAlgorithm() const noexcept {
    return advection_.GetGriddingAlgorithm();
  }

  const Equation& GetEquation() const noexcept {
    return equation_;
  }

  ////////////////////////////////////////////////////////////////////////////
  // Solver interface

  auto&& GetCounterRegistry() const noexcept { return advection_.GetCounterRegistry(); }
  auto GetCycles() const noexcept { return advection_.GetCycles(); }
  auto GetTimePoint() const noexcept { return advection_.GetTimePoint(); }
  void PreAdvanceHierarchy() { advection_.PreAdvanceHierarchy(); }
  void PostAdvanceHierarchy(Duration dt) { advection_.PostAdvanceHierarchy(dt); }
  void ResetHierarchyConfiguration(std::shared_ptr<GriddingAlgorithm> grid) {
    advection_.ResetHierarchyConfiguration(std::move(grid));
  }

  ////////////////////////////////////////////////////////////////////////////
  // Modifiers

  Duration ComputeStableDt() { return advection_.ComputeStableDt(); }

  Result<void, TimeStepTooLarge> AdvanceHierarchy(Duration dt);

  void DoInitialProjection();

  void RecomputeAdvectiveFluxes();

  void FillAllGhostLayers();

private:
  Equation equation_;
  BK19PhysicalParameters physical_parameters_;
  BK19SolverOptions options_;
  AdvectionSolver advection_;
  std::shared_ptr<LinearOperator> lin_op_;
};

extern template class BK19Solver<2, 2>;

} // namespace fub::amrex

#endif