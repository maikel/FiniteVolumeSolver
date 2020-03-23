// Copyright (c) 2020 Stefan Vater
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

#include "fub/AMReX/ViscositySourceTerm.hpp"
#include "fub/AMReX/ForEachFab.hpp"
#include "fub/AMReX/ForEachIndex.hpp"
#include "fub/equations/PerfectGas.hpp"

namespace fub::amrex {

template <int Rank>
ViscositySourceTerm<Rank>::ViscositySourceTerm(const PerfectGas<Rank>& eq,
                                 const std::shared_ptr<GriddingAlgorithm>& grid)
    : equation_(eq) {
  ResetHierarchyConfiguration(grid);

  int finest_level = grid->finestLevel();

  ::amrex::LPInfo info_solve;
  info_solve.setMaxCoarseningLevel(m_mg_max_coarsening_level);
  ::amrex::LPInfo info_apply;
  info_apply.setMaxCoarseningLevel(0);

}

// ViscositySourceTerm::ViscositySourceTerm(const ViscositySourceTerm& other)
//     : equation_(other.equation_) {
//   Ax_.reserve(other.Ax_.size());
//   for (const ::amrex::MultiFab& other_Ax : other.Ax_) {
//     ::amrex::BoxArray box_array = other_Ax.boxArray();
//     ::amrex::DistributionMapping distribution_map = other_Ax.DistributionMap();
//     const int ncomp = other_Ax.nComp();
//     const int ngrow = other_Ax.nGrow();
//     ::amrex::MultiFab& Ax =
//         Ax_.emplace_back(box_array, distribution_map, ncomp, ngrow);
//     Ax.copy(other_Ax);
//   }
// }
//
// ViscositySourceTerm& ViscositySourceTerm::operator=(const ViscositySourceTerm& other) {
//   ViscositySourceTerm tmp(other);
//   return (*this = std::move(tmp));
// }

namespace {
} // namespace

template <int Rank>
void ViscositySourceTerm<Rank>::ResetHierarchyConfiguration(
    const std::shared_ptr<amrex::GriddingAlgorithm>& grid) {
}

template <int Rank>
Duration ViscositySourceTerm<Rank>::ComputeStableDt(int /* level */) {
  return Duration(std::numeric_limits<double>::max());
}

template <int Rank>
Result<void, TimeStepTooLarge>
ViscositySourceTerm<Rank>::AdvanceLevel(IntegratorContext& simulation_data, int level,
                              Duration time_step_size) {
  ::amrex::MultiFab& data = simulation_data.GetScratch(level);
  const double dt = time_step_size.count();

  return boost::outcome_v2::success();
}

template class ViscositySourceTerm<1>;
template class ViscositySourceTerm<2>;
template class ViscositySourceTerm<3>;

} // namespace fub::amrex
