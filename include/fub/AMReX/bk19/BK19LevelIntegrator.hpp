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

#include "fub/AMReX/MLMG/MLNodeHelmDualCstVel.hpp"
#include "fub/AMReX/bk19/BK19IntegratorContext.hpp"
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

class BK19LevelIntegrator
    : private DimensionalSplitLevelIntegrator<AMREX_SPACEDIM,
                                              BK19IntegratorContext, AnySplitMethod> {
public:
  static constexpr int Rank = AMREX_SPACEDIM;

  using Coordinates = Eigen::Matrix<double, Rank, 1>;
  using Equation = CompressibleAdvection<Rank>;
  using Complete = ::fub::Complete<Equation>;
  using Conservative = ::fub::Conservative<Equation>;
  using SplittingMethod = ::fub::AnySplitMethod;

  using AdvectionSolver =
      DimensionalSplitLevelIntegrator<Rank, BK19IntegratorContext, SplittingMethod>;

  using AdvectionSolver::ApplyFluxCorrection;
  using AdvectionSolver::CoarsenConservatively;
  using AdvectionSolver::CompleteFromCons;
  using AdvectionSolver::ComputeStableDt;
  using AdvectionSolver::CopyDataToScratch;
  using AdvectionSolver::CopyScratchToData;
  using AdvectionSolver::GetContext;
  using AdvectionSolver::GetCycles;
  using AdvectionSolver::GetGriddingAlgorithm;
  using AdvectionSolver::GetCounterRegistry;
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
      std::shared_ptr<::amrex::MLNodeHelmDualCstVel> linop,
      const BK19LevelIntegratorOptions& options = BK19LevelIntegratorOptions());

  void ResetPatchHierarchy(std::shared_ptr<GriddingAlgorithm> grid);

  Result<void, TimeStepTooLarge>
  AdvanceLevelNonRecursively(int level, Duration dt,
                             std::pair<int, int> subcycle);

private:
  BK19LevelIntegratorOptions options_;
  CompressibleAdvection<Rank> equation_;
  fub::IndexMapping<fub::CompressibleAdvection<2>> index_;
  std::shared_ptr<::amrex::MLNodeHelmDualCstVel> lin_op_;
};

void WriteRawField(const std::string& path, const std::string& name,
                   const ::amrex::MultiFab& data, int level);

void WriteAdvectiveFluxes(const std::string& path,
                          const BK19AdvectiveFluxes& pv, int level);

struct WriteBK19Plotfile {
  std::string plotfilename{};

  void operator()(const BK19IntegratorContext& context) const;
  void operator()(const GriddingAlgorithm& grid) const;
};

template <typename Log> void BK19LevelIntegratorOptions::Print(Log& log) {
  BOOST_LOG(log) << fmt::format(" - mlmg_tolerance_rel = {}", mlmg_tolerance_rel);
  BOOST_LOG(log) << fmt::format(" - mlmg_tolerance_abs = {}", mlmg_tolerance_abs);
  BOOST_LOG(log) << fmt::format(" - mlmg_max_iter = {}", mlmg_max_iter);
  BOOST_LOG(log) << fmt::format(" - mlmg_verbose = {}", mlmg_verbose);
  BOOST_LOG(log) << fmt::format(" - bottom_tolerance_rel = {}", bottom_tolerance_rel);
  BOOST_LOG(log) << fmt::format(" - bottom_tolerance_abs = {}", bottom_tolerance_abs);
  BOOST_LOG(log) << fmt::format(" - bottom_max_iter = {}", bottom_max_iter);
  BOOST_LOG(log) << fmt::format(" - bottom_verbose = {}", bottom_verbose);
  BOOST_LOG(log) << fmt::format(" - always_use_bnorm = {}", always_use_bnorm);
  BOOST_LOG(log) << fmt::format(" - prefix = {}", prefix);
  BOOST_LOG(log) << fmt::format(" - output_between_steps = {}", output_between_steps);
}

} // namespace fub::amrex

#endif
