// Copyright (c) 2019 Maikel Nadolski
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

#ifndef FUB_BK19_LEVEL_INTEGRATOR_HPP
#define FUB_BK19_LEVEL_INTEGRATOR_HPP

#include "fub/AMReX/MLMG/MLNodeHelmholtz.hpp"
#include "fub/AMReX/CompressibleAdvectionIntegratorContext.hpp"
#include "fub/equations/CompressibleAdvection.hpp"
#include "fub/ext/Eigen.hpp"
#include "fub/ext/ProgramOptions.hpp"
#include "fub/solver/DimensionalSplitLevelIntegrator.hpp"

#include <AMReX_MLMG.H>

#include <AMReX_PlotFileUtil.H>
#include <AMReX_VisMF.H>

namespace fub::amrex {

void RecomputeAdvectiveFluxes(
    const IndexMapping<CompressibleAdvection<2>>& index,
    std::array<::amrex::MultiFab, 2>& Pv_faces, ::amrex::MultiFab& Pv_cells,
    const ::amrex::MultiFab& scratch, const ::amrex::Periodicity& periodicity);

struct BK19LevelIntegratorOptions {
  BK19LevelIntegratorOptions() = default;
  BK19LevelIntegratorOptions(const ProgramOptions& map);

  template <typename Log> void Print(Log& log);

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

  double alpha_p{1.0};
  double Msq{0.0};
};

class BK19LevelIntegrator
    : private DimensionalSplitLevelIntegrator<
          AMREX_SPACEDIM, CompressibleAdvectionIntegratorContext, AnySplitMethod> {
public:
  static constexpr int Rank = AMREX_SPACEDIM;

  using Coordinates = Eigen::Matrix<double, Rank, 1>;
  using Equation = CompressibleAdvection<Rank>;
  using Complete = ::fub::Complete<Equation>;
  using Conservative = ::fub::Conservative<Equation>;
  using SplittingMethod = ::fub::AnySplitMethod;

  using AdvectionSolver =
      DimensionalSplitLevelIntegrator<Rank, CompressibleAdvectionIntegratorContext,
                                      SplittingMethod>;

  using AdvectionSolver::ApplyFluxCorrection;
  using AdvectionSolver::CoarsenConservatively;
  using AdvectionSolver::CompleteFromCons;
  using AdvectionSolver::ComputeStableDt;
  using AdvectionSolver::CopyDataToScratch;
  using AdvectionSolver::CopyScratchToData;
  using AdvectionSolver::GetContext;
  using AdvectionSolver::GetCounterRegistry;
  using AdvectionSolver::GetCycles;
  using AdvectionSolver::GetGriddingAlgorithm;
  using AdvectionSolver::GetMpiCommunicator;
  using AdvectionSolver::GetRatioToCoarserLevel;
  using AdvectionSolver::GetTimePoint;
  using AdvectionSolver::LevelExists;
  using AdvectionSolver::PostAdvanceHierarchy;
  using AdvectionSolver::PostAdvanceLevel;
  using AdvectionSolver::PreAdvanceHierarchy;
  using AdvectionSolver::PreAdvanceLevel;
  using AdvectionSolver::ResetCoarseFineFluxes;
  using AdvectionSolver::ResetHierarchyConfiguration;

  AdvectionSolver& GetAdvection() { return *this; }
  const AdvectionSolver& GetAdvection() const { return *this; }

  BK19LevelIntegrator(
      const CompressibleAdvection<Rank>& equation, AdvectionSolver advection,
      std::shared_ptr<::amrex::MLNodeHelmholtz> linop,
      const BK19PhysicalParameters& physical_parameters,
      const BK19LevelIntegratorOptions& options = BK19LevelIntegratorOptions());

  void ResetPatchHierarchy(std::shared_ptr<GriddingAlgorithm> grid);

  Result<void, TimeStepTooLarge>
  AdvanceLevelNonRecursively(int level, Duration dt,
                             std::pair<int, int> subcycle);

  void InitialProjection(
                        int level,
                        Duration dt,
                        std::array<double,2> U0
                        );


private:
  BK19PhysicalParameters phys_param_;
  BK19LevelIntegratorOptions options_;
  CompressibleAdvection<Rank> equation_;
  fub::IndexMapping<fub::CompressibleAdvection<2>> index_;
  std::shared_ptr<::amrex::MLNodeHelmholtz> lin_op_;
};

void WriteRawField(const std::string& path, const std::string& name,
                   const ::amrex::MultiFab& data, int level);

struct WriteBK19Plotfile {
  std::string plotfilename{};

  void operator()(const CompressibleAdvectionIntegratorContext& context) const;
  void operator()(const GriddingAlgorithm& grid) const;
};

template <typename Log> void BK19LevelIntegratorOptions::Print(Log& log) {
  BOOST_LOG(log) << fmt::format(" - mlmg_tolerance_rel = {}",
                                mlmg_tolerance_rel);
  BOOST_LOG(log) << fmt::format(" - mlmg_tolerance_abs = {}",
                                mlmg_tolerance_abs);
  BOOST_LOG(log) << fmt::format(" - mlmg_max_iter = {}", mlmg_max_iter);
  BOOST_LOG(log) << fmt::format(" - mlmg_verbose = {}", mlmg_verbose);
  BOOST_LOG(log) << fmt::format(" - bottom_tolerance_rel = {}",
                                bottom_tolerance_rel);
  BOOST_LOG(log) << fmt::format(" - bottom_tolerance_abs = {}",
                                bottom_tolerance_abs);
  BOOST_LOG(log) << fmt::format(" - bottom_max_iter = {}", bottom_max_iter);
  BOOST_LOG(log) << fmt::format(" - bottom_verbose = {}", bottom_verbose);
  BOOST_LOG(log) << fmt::format(" - always_use_bnorm = {}", always_use_bnorm);
  BOOST_LOG(log) << fmt::format(" - prefix = {}", prefix);
  BOOST_LOG(log) << fmt::format(" - output_between_steps = {}",
                                output_between_steps);
}

} // namespace fub::amrex

#endif
