// Copyright (c) 2019 Maikel Nadolski
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

#include "fub/AMReX/AxialSourceTerm.hpp"
#include "fub/AMReX/ForEachFab.hpp"
#include "fub/AMReX/ForEachIndex.hpp"
#include "fub/equations/IdealGasMix.hpp"

namespace fub::amrex {

AxialSourceTerm::AxialSourceTerm(const IdealGasMix<1>& eq,
                                 std::function<double(double)> diameter,
                                 const std::shared_ptr<GriddingAlgorithm>& grid)
    : diameter_(std::move(diameter)), equation_(eq) {
  ResetHierarchyConfiguration(grid);
}

AxialSourceTerm::AxialSourceTerm(const AxialSourceTerm& other)
    : diameter_(other.diameter_), equation_(other.equation_) {
  Ax_.reserve(other.Ax_.size());
  for (const ::amrex::MultiFab& other_Ax : other.Ax_) {
    ::amrex::BoxArray box_array = other_Ax.boxArray();
    ::amrex::DistributionMapping distribution_map = other_Ax.DistributionMap();
    const int ncomp = other_Ax.nComp();
    const int ngrow = other_Ax.nGrow();
    ::amrex::MultiFab& Ax =
        Ax_.emplace_back(box_array, distribution_map, ncomp, ngrow);
    Ax.copy(other_Ax);
  }
}

AxialSourceTerm& AxialSourceTerm::operator=(const AxialSourceTerm& other) {
  AxialSourceTerm tmp(other);
  return (*this = std::move(tmp));
}

namespace {
std::vector<::amrex::MultiFab>
ComputeDiameters(const amrex::GriddingAlgorithm& grid,
                 const std::function<double(double)>& diameter) {
  const PatchHierarchy& hier = grid.GetPatchHierarchy();
  const int nlevel = hier.GetNumberOfLevels();
  std::vector<::amrex::MultiFab> ds{};
  ds.reserve(static_cast<std::size_t>(nlevel));
  for (int ilvl = 0; ilvl < nlevel; ++ilvl) {
    const ::amrex::MultiFab& data = hier.GetPatchLevel(ilvl).data;
    const ::amrex::Geometry& geom = hier.GetGeometry(ilvl);
    const double dx = geom.CellSize(0);
    ::amrex::MultiFab& mf =
        ds.emplace_back(data.boxArray(), data.DistributionMap(), 2, 0);
    ForEachFab(execution::openmp, mf, [&](const ::amrex::MFIter& mfi) {
      ::amrex::FArrayBox& fab = mf[mfi];
      ForEachIndex(fab.box(), [&](auto... is) {
        ::amrex::IntVect i{int(is)...};
        const double xM = geom.CellCenter(i[0], 0);
        const double xR = geom.CellCenter(i[0] + 1, 0);
        fab(i, 0) = diameter(xM + 0.5 * dx);
        fab(i, 1) = (diameter(xR) - diameter(xM)) / dx;
      });
    });
  }
  return ds;
}
} // namespace

void AxialSourceTerm::ResetHierarchyConfiguration(
    const std::shared_ptr<amrex::GriddingAlgorithm>& grid) {
  Ax_ = ComputeDiameters(*grid, diameter_);
}

Duration AxialSourceTerm::ComputeStableDt(int /* level */) {
  return Duration(std::numeric_limits<double>::max());
}

Result<void, TimeStepTooLarge>
AxialSourceTerm::AdvanceLevel(IntegratorContext& simulation_data, int level,
                              Duration time_step_size) {
  ResetHierarchyConfiguration(simulation_data.GetGriddingAlgorithm());
  ::amrex::MultiFab& data = simulation_data.GetScratch(level);
  const double dt = time_step_size.count();
  Complete<IdealGasMix<1>> state(equation_);
  ForEachFab(data, [&](const ::amrex::MFIter& mfi) {
    ::amrex::Box box = mfi.tilebox();
    ::amrex::FArrayBox& fab = data[mfi];
    const ::amrex::FArrayBox& Ax = Ax_[static_cast<std::size_t>(level)][mfi];
    auto states = MakeView<Complete<IdealGasMix<1>>>(fab, equation_, box);
    FlameMasterReactor& reactor = equation_.GetReactor();
    ForEachIndex(Box<0>(states), [&](std::ptrdiff_t i) {
      const double dAdx = Ax({AMREX_D_DECL(int(i), 0, 0)}, 1);
      const double A = Ax({AMREX_D_DECL(int(i), 0, 0)}, 0);
      Load(state, states, {i});
      const double Rspec = state.c_p - state.c_p / state.gamma;
      const double u = state.momentum[0] / state.density;
      const double u_dAdx_A = u * dAdx / A;
      const double gamma = state.gamma;
      const double rho = state.density * std::exp(-u_dAdx_A * dt);
      const double p = state.pressure * std::exp(-gamma * u_dAdx_A * dt);
      const double T = p / (Rspec * rho);
      reactor.SetDensity(rho);
      reactor.SetMassFractions(state.species);
      reactor.SetTemperature(T);
      equation_.CompleteFromReactor(state, Eigen::Array<double, 1, 1>{u});
      Store(states, state, {i});
    });
  });
  return boost::outcome_v2::success();
}

} // namespace fub::amrex
