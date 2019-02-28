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

#include "fub/grid/SAMRAI/CartesianPatchHierarchy.hpp"
#include "fub/grid/SAMRAI/InitialData.hpp"
#include "fub/grid/SAMRAI/RegisterVariables.hpp"
#include "fub/grid/SAMRAI/Tagging.hpp"
#include "fub/grid/SAMRAI/ViewPatch.hpp"

#include "fub/Duration.hpp"

#include <SAMRAI/mesh/GriddingAlgorithm.h>

namespace fub {
namespace samrai {

struct GriddingAlgorithm {
  GriddingAlgorithm(const std::shared_ptr<SAMRAI::hier::PatchHierarchy>& hier,
                    DataDescription description, InitialData initial_data,
                    Tagging tagging, std::vector<int> tag_buffer);

  void RegridAllFinerLevels(int level_num, int cycle, Duration time_point);

  void InitializeHierarchy(Duration initial_time = Duration(),
                           std::ptrdiff_t initial_cycle = 0);

  const std::shared_ptr<SAMRAI::hier::PatchHierarchy>& GetPatchHierarchy() const
      noexcept;

  std::shared_ptr<SAMRAI::hier::PatchHierarchy> hierarchy;
  InitialData initial_data;
  Tagging tagging;
  DataDescription id_set;
  std::vector<int> tag_buffer;
  std::shared_ptr<SAMRAI::xfer::RefineAlgorithm> refine_data_algorithm_{};
  std::shared_ptr<SAMRAI::mesh::GriddingAlgorithm> algorithm;
};

} // namespace samrai
} // namespace fub

#endif