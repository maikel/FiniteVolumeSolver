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
#include "fub/SAMRAI/InitialData.hpp"
#include "fub/SAMRAI/RegisterVariables.hpp"
#include "fub/SAMRAI/Tagging.hpp"
#include "fub/SAMRAI/ViewPatch.hpp"

#include "fub/Duration.hpp"

#include <SAMRAI/mesh/GriddingAlgorithm.h>

namespace fub {
namespace samrai {

struct GriddingAlgorithm {
  GriddingAlgorithm(PatchHierarchy hier, InitialData initial_data,
                    Tagging tagging, std::vector<int> tag_buffer);

  GriddingAlgorithm(const GriddingAlgorithm& ga);
  GriddingAlgorithm& operator=(const GriddingAlgorithm& ga) {
    GriddingAlgorithm tmp(ga);
    return (*this = std::move(tmp));
  }

  GriddingAlgorithm(GriddingAlgorithm&& ph) = default;
  GriddingAlgorithm& operator=(GriddingAlgorithm&& ph) = default;

  void RegridAllFinerLevels(int level_num, int cycle, Duration time_point);

  void InitializeHierarchy(Duration initial_time = Duration(),
                           std::ptrdiff_t initial_cycle = 0);

  const PatchHierarchy& GetPatchHierarchy() const noexcept;
  PatchHierarchy& GetPatchHierarchy() noexcept;

  const InitialData& GetInitialData() const noexcept;
  InitialData& GetInitialData() noexcept;

  const Tagging& GetTagging() const noexcept;
  Tagging& GetTagging() noexcept;

  const DataDescription& GetDataDescription() const noexcept;
  DataDescription& GetDataDescription() noexcept;

  const std::vector<int>& GetTagBuffer() const noexcept;
  std::vector<int>& GetTagBuffer() noexcept;

  PatchHierarchy hierarchy_;
  InitialData initial_data_;
  Tagging tagging_;
  DataDescription id_set_;
  std::vector<int> tag_buffer_;
  std::shared_ptr<SAMRAI::xfer::RefineAlgorithm> refine_data_algorithm_{};
  std::shared_ptr<SAMRAI::mesh::GriddingAlgorithm> algorithm_;
};

} // namespace samrai
} // namespace fub

#endif
