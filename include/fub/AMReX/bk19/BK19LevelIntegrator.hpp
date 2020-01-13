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

#include "fub/AMReX/MLMG/MLNodeHelmDualCstVel.hpp"
#include "fub/AMReX/bk19/BK19IntegratorContext.hpp"
#include "fub/equations/CompressibleAdvection.hpp"
#include "fub/ext/Eigen.hpp"
#include "fub/solver/DimensionalSplitLevelIntegrator.hpp"

#include <AMReX_MLMG.H>

#include <AMReX_PlotFileUtil.H>
#include <AMReX_VisMF.H>

namespace fub::amrex {

class BK19LevelIntegrator
    : private DimensionalSplitLevelIntegrator<AMREX_SPACEDIM,
                                              BK19IntegratorContext> {
public:
  static constexpr int Rank = AMREX_SPACEDIM;

  using Coordinates = Eigen::Matrix<double, Rank, 1>;
  using Equation = CompressibleAdvection<Rank>;
  using Complete = ::fub::Complete<Equation>;
  using Conservative = ::fub::Conservative<Equation>;

  using AdvectionSolver =
      DimensionalSplitLevelIntegrator<Rank, BK19IntegratorContext>;

  using AdvectionSolver::ApplyFluxCorrection;
  using AdvectionSolver::CoarsenConservatively;
  using AdvectionSolver::CompleteFromCons;
  using AdvectionSolver::ComputeStableDt;
  using AdvectionSolver::CopyDataToScratch;
  using AdvectionSolver::CopyScratchToData;
  using AdvectionSolver::GetContext;
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

  BK19LevelIntegrator(const CompressibleAdvection<Rank>& equation,
                      AdvectionSolver advection,
                      std::shared_ptr<::amrex::MLNodeHelmDualCstVel> linop);

  void ResetPatchHierarchy(std::shared_ptr<GriddingAlgorithm> grid);

  Result<void, TimeStepTooLarge>
  AdvanceLevelNonRecursively(int level, Duration dt,
                             std::pair<int, int> subcycle);

private:
  CompressibleAdvection<Rank> equation_;
  fub::IndexMapping<fub::CompressibleAdvection<2>> index_;
  std::shared_ptr<::amrex::MLNodeHelmDualCstVel> lin_op_;
  std::shared_ptr<::amrex::MLMG> nodal_solver_;
};

struct WriteBK19Plotfile {
  std::string plotfilename{};
void operator()(const fub::amrex::GriddingAlgorithm& grid) const
{
  using Equation = CompressibleAdvection<2>;
  fub::CompressibleAdvection<2> equation{};
  const fub::amrex::PatchHierarchy& hier = grid.GetPatchHierarchy();
  std::string name = fmt::format("{}/plt{:09}", plotfilename, grid.GetCycles());
  const int nlevels = hier.GetNumberOfLevels();
  const double time_point = hier.GetTimePoint().count();
  FUB_ASSERT(nlevels >= 0);
  std::size_t size = static_cast<std::size_t>(nlevels);
  ::amrex::Vector<const ::amrex::MultiFab*> mf(size);
  ::amrex::Vector<const ::amrex::MultiFab*> mfnodes(size);
  ::amrex::Vector<::amrex::Geometry> geoms(size);
  ::amrex::Vector<int> level_steps(size);
  ::amrex::Vector<::amrex::IntVect> ref_ratio(size);
  for (std::size_t i = 0; i < size; ++i) {
    mf[i] = &hier.GetPatchLevel(static_cast<int>(i)).data;
    mfnodes[i] = hier.GetPatchLevel(static_cast<int>(i)).nodes.get();
    geoms[i] = hier.GetGeometry(static_cast<int>(i));
    level_steps[i] = static_cast<int>(hier.GetCycles(static_cast<int>(i)));
    ref_ratio[i] = hier.GetRatioToCoarserLevel(static_cast<int>(i));
  }
  using Traits = StateTraits<Complete<Equation>>;
  constexpr auto names = Traits::names;
  const auto depths = Depths<Complete<Equation>>(equation);
  const std::size_t n_names =
      std::tuple_size<remove_cvref_t<decltype(names)>>::value;
  ::amrex::Vector<std::string> varnames;
  varnames.reserve(n_names);
  boost::mp11::tuple_for_each(Zip(names, StateToTuple(depths)), [&](auto xs) {
    const int ncomp = std::get<1>(xs);
    if (ncomp == 1) {
      varnames.push_back(std::get<0>(xs));
    } else {
      for (int i = 0; i < ncomp; ++i) {
        varnames.push_back(fmt::format("{}_{}", std::get<0>(xs), i));
      }
    }
  });

  ::amrex::Vector<std::string> rfs {"raw_fields"};
  ::amrex::WriteMultiLevelPlotfile(name, nlevels, mf, varnames, geoms,
                                   time_point, level_steps, ref_ratio,
                                   "HyperCLaw-V1.1", "Level_", "Cell", rfs);

  // write nodal raw fields
  ::amrex::VisMF::Header::Version plotfile_headerversion = ::amrex::VisMF::Header::Version_v1;
  ::amrex::VisMF::SetHeaderVersion(plotfile_headerversion);
  const std::string raw_pltname = name + "/" + rfs[0];
  for (int lev = 0; lev < nlevels; ++lev) {
    ::amrex::VisMF::Write(*mfnodes[lev],
      ::amrex::MultiFabFileFullPrefix(lev, raw_pltname, "Level_", "pi"));
  }

}
};

} // namespace fub::amrex
