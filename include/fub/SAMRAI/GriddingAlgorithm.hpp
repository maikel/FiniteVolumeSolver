// Copyright (c) 2018 Maikel Nadolski
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

#ifndef FUB_SAMRAI_GRIDDING_ALGORITHM_HPP
#define FUB_SAMRAI_GRIDDING_ALGORITHM_HPP

#include "fub/SAMRAI/PatchHierarchy.hpp"
#include "fub/SAMRAI/RegisterVariables.hpp"
#include "fub/SAMRAI/ViewPatch.hpp"

#include "fub/AnyBoundaryCondition.hpp"
#include "fub/AnyInitialData.hpp"
#include "fub/AnyTaggingMethod.hpp"

#include "fub/Duration.hpp"

#include <SAMRAI/mesh/GriddingAlgorithm.h>

namespace fub {
namespace samrai {
class GriddingAlgorithm;
}

template <> struct GridTraits<samrai::GriddingAlgorithm> {
  using PatchLevel = std::shared_ptr<SAMRAI::hier::PatchLevel>;
  using TagDataHandle = int;
  using DataReference = SAMRAI::hier::Patch&;
};

namespace samrai {
using AnyBoundaryCondition = ::fub::AnyBoundaryCondition<GriddingAlgorithm>;
using AnyInitialData = ::fub::AnyInitialData<GriddingAlgorithm>;
using AnyTaggingMethod = ::fub::AnyTaggingMethod<GriddingAlgorithm>;

class GriddingAlgorithm {
public:
  GriddingAlgorithm(PatchHierarchy hier, AnyInitialData initial_data,
                    AnyTaggingMethod tagging, std::vector<int> tag_buffer);

  GriddingAlgorithm(const GriddingAlgorithm& ga);
  GriddingAlgorithm& operator=(const GriddingAlgorithm& ga) {
    GriddingAlgorithm tmp(ga);
    std::swap(*this, tmp);
    return *this;
  }

  GriddingAlgorithm(GriddingAlgorithm&& ph) = default;
  GriddingAlgorithm& operator=(GriddingAlgorithm&& ph) = default;

  void RegridAllFinerLevels(int level_num);

  void InitializeHierarchy(Duration initial_time = Duration(),
                           std::ptrdiff_t initial_cycle = 0);

  const PatchHierarchy& GetPatchHierarchy() const noexcept;
  PatchHierarchy& GetPatchHierarchy() noexcept;

  const std::vector<int>& GetTagBuffer() const noexcept;

  const AnyInitialData& GetInitialData() const noexcept;

  const AnyTaggingMethod& GetTaggingMethod() const noexcept;
  
  const AnyBoundaryCondition& GetBoundaryCondition() const noexcept;
  AnyBoundaryCondition& GetBoundaryCondition() noexcept;

private:
  PatchHierarchy hierarchy_;
  AnyInitialData initial_data_{};
  AnyTaggingMethod tagging_method_{};
  AnyBoundaryCondition boundary_condition_{};

  std::vector<int> tag_buffer_;
  std::shared_ptr<SAMRAI::xfer::RefineAlgorithm> refine_data_algorithm_{};
  std::shared_ptr<SAMRAI::mesh::GriddingAlgorithm> algorithm_{};
};

} // namespace samrai
} // namespace fub

#endif
