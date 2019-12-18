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

#ifndef FUB_AMREX_BK19_GRIDDING_ALGORITHM_HPP
#define FUB_AMREX_BK19_GRIDDING_ALGORITHM_HPP

#include "fub/AMReX/GriddingAlgorithm.hpp"

namespace fub {
namespace amrex {

class BK19GriddingAlgorithm : public GriddingAlgorithm {
public:
  static constexpr int Rank = AMREX_SPACEDIM;

  // construct a new object

  BK19GriddingAlgorithm(PatchHierarchy hier, InitialData initial_data,
                        Tagging tagging, int pv_ghost_cell_width);

  BK19GriddingAlgorithm(PatchHierarchy hier, InitialData initial_data,
                        Tagging tagging, BoundaryCondition boundary,
                        int pv_ghost_cell_width);

  // copy and move constructors

  BK19GriddingAlgorithm(const BK19GriddingAlgorithm& other);
  BK19GriddingAlgorithm& operator=(const BK19GriddingAlgorithm& other);

  BK19GriddingAlgorithm(BK19GriddingAlgorithm&&) noexcept;
  BK19GriddingAlgorithm& operator=(BK19GriddingAlgorithm&&) noexcept;

  ::amrex::MultiFab& GetPvs(int level, Direction dir);
  const ::amrex::MultiFab& GetPvs(int level, Direction dir) const;

private:
  void AllocatePvs();
  void AllocatePvs(int level);
  void AllocatePvs(int level, Direction dir);

  void MakeNewLevelFromScratch(
      int level, double time_point, const ::amrex::BoxArray& box_array,
      const ::amrex::DistributionMapping& distribution_mapping) override;

  void MakeNewLevelFromCoarse(
      int level, double time_point, const ::amrex::BoxArray& box_array,
      const ::amrex::DistributionMapping& distribution_mapping) override;

  void RemakeLevel(
      int level, double time_point, const ::amrex::BoxArray& box_array,
      const ::amrex::DistributionMapping& distribution_mapping) override;

  void ClearLevel(int level) override;

  int pv_ghost_cell_width_;
  std::array<std::vector<::amrex::MultiFab>, Rank> Pvs_{};
};

} // namespace amrex
} // namespace fub

#endif
