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

#ifndef FUB_AMREX_COUPLED_GRIDDING_ALGORITHM_HPP
#define FUB_AMREX_COUPLED_GRIDDING_ALGORITHM_HPP

#include "fub/AMReX/multi_block/MultiBlockBoundary.hpp"

#include "fub/AMReX/GriddingAlgorithm.hpp"
#include "fub/AMReX/cutcell/GriddingAlgorithm.hpp"

#include "fub/core/span.hpp"

#include "fub/equations/ideal_gas_mix/FlameMasterReactor.hpp"

namespace fub {
namespace amrex {

class MultiBlockGriddingAlgorithm {
public:
  MultiBlockGriddingAlgorithm(
      FlameMasterReactor reactor,
      std::vector<std::shared_ptr<GriddingAlgorithm>> tubes,
      std::vector<std::shared_ptr<cutcell::GriddingAlgorithm>> plena,
      std::vector<BlockConnection> connectivity);

  MultiBlockGriddingAlgorithm(
      FlameMasterReactor reactor,
      std::vector<std::shared_ptr<GriddingAlgorithm>> tubes,
      std::vector<std::shared_ptr<cutcell::GriddingAlgorithm>> plena,
      std::vector<BlockConnection> connectivity,
      std::vector<std::shared_ptr<PressureValve>> valves);

  MultiBlockGriddingAlgorithm(const MultiBlockGriddingAlgorithm& other);
  MultiBlockGriddingAlgorithm&
  operator=(const MultiBlockGriddingAlgorithm& other);

  MultiBlockGriddingAlgorithm(MultiBlockGriddingAlgorithm&& other) noexcept =
      default;

  MultiBlockGriddingAlgorithm&
  operator=(MultiBlockGriddingAlgorithm&& other) noexcept = default;

  [[nodiscard]] span<const std::shared_ptr<GriddingAlgorithm>> GetTubes() const
      noexcept;
  [[nodiscard]] span<const std::shared_ptr<cutcell::GriddingAlgorithm>>
  GetPlena() const noexcept;

  [[nodiscard]] span<const BlockConnection> GetConnectivity() const noexcept;
  [[nodiscard]] span<MultiBlockBoundary> GetBoundaries(int level = 0) noexcept {
    return boundaries_[static_cast<std::size_t>(level)];
  }

  void RegridAllFinerLevels(int which_level);

private:
  FlameMasterReactor reactor_;
  std::vector<std::shared_ptr<GriddingAlgorithm>> tubes_;
  std::vector<std::shared_ptr<cutcell::GriddingAlgorithm>> plena_;
  std::vector<BlockConnection> connectivity_;
  std::vector<std::vector<MultiBlockBoundary>> boundaries_;
};

} // namespace amrex
} // namespace fub

#endif // FUB_AMREX_COUPLED_GRIDDING_ALGORITHM_HPP
