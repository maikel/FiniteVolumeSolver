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

#include "fub/AMReX/GriddingAlgorithm.hpp"
#include "fub/AMReX/cutcell/GriddingAlgorithm.hpp"

#include "fub/core/span.hpp"

#include "fub/equations/ideal_gas_mix/FlameMasterReactor.hpp"

namespace fub {
namespace amrex {

struct BlockEntry {
  std::size_t id;
  ::amrex::Box mirror_box;
};

struct BlockConnection {
  BlockEntry tube;
  BlockEntry plenum;
  Direction direction;
  int side;
};

class CoupledBoundaryCondition;

class CoupledGriddingAlgorithm {
public:
  CoupledGriddingAlgorithm(FlameMasterReactor reactor,
                           std::vector<GriddingAlgorithm> tubes,
                           std::vector<cutcell::GriddingAlgorithm> plena,
                           std::vector<BlockConnection> connectivity);

  CoupledGriddingAlgorithm(const CoupledGriddingAlgorithm& other);
  CoupledGriddingAlgorithm& operator=(const CoupledGriddingAlgorithm& other);

  CoupledGriddingAlgorithm(CoupledGriddingAlgorithm&& other) noexcept;
  CoupledGriddingAlgorithm& operator=(CoupledGriddingAlgorithm&& other) noexcept;

  span<GriddingAlgorithm> GetTubes() noexcept;
  span<const GriddingAlgorithm> GetTubes() const noexcept;

  span<cutcell::GriddingAlgorithm> GetPlena() noexcept;
  span<const cutcell::GriddingAlgorithm> GetPlena() const noexcept;

  span<const BlockConnection> GetConnectivity() const noexcept;

  CoupledBoundaryCondition& GetBoundaryCondition(int level) noexcept;

private:
  ideal_gas_mix::Mechanism mechanism_;
  std::vector<std::shared_ptr<GriddingAlgorithm>> tube_;
  std::vector<std::shared_ptr<cutcell::GriddingAlgorithm>> plenum_;
  std::vector<BlockConnection> connectivity_;
  std::vector<CoupledBoundaryCondition> boundary_condition_;
};

} // namespace amrex
} // namespace fub

#endif // FUB_AMREX_COUPLED_GRIDDING_ALGORITHM_HPP
