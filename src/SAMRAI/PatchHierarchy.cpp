// Copyright (c) 2019 Maikel Nadolski
// Copyright (c) 2019 Patrick Denzler
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

#include "fub/SAMRAI/PatchHierarchy.hpp"

namespace fub::samrai {

SAMRAI::hier::ComponentSelector SelectComponents(span<const int> data_ids) {
  SAMRAI::hier::ComponentSelector selector;
  for (int id : data_ids) {
    selector.setFlag(id);
  }
  return selector;
}

const PatchHierarchyOptions& PatchHierarchy::GetOptions() const noexcept {
  return options_;
}

const std::shared_ptr<SAMRAI::hier::PatchHierarchy>&
PatchHierarchy::GetNative() const noexcept {
  return hierarchy_;
}

const DataDescription& PatchHierarchy::GetDataDescription() const noexcept {
  return data_desc_;
}

PatchHierarchy::PatchHierarchy(
    DataDescription dd,
    std::shared_ptr<SAMRAI::geom::CartesianGridGeometry> geom,
    PatchHierarchyOptions hier_opts)
    : data_desc_{std::move(dd)}, options_{std::move(hier_opts)} {
  hierarchy_ = std::make_shared<SAMRAI::hier::PatchHierarchy>(
      boost::uuids::to_string(boost::uuids::random_generator()()), geom);
  hierarchy_->setMaxNumberOfLevels(options_.max_number_of_levels);
  for (int i = 1; i < options_.max_number_of_levels; ++i) {
    hierarchy_->setRatioToCoarserLevel(options_.refine_ratio, i);
    hierarchy_->setSmallestPatchSize(SAMRAI::hier::IntVector(hierarchy_->getDim(), 8), i);
  }
}

PatchHierarchy::PatchHierarchy(const PatchHierarchy& ph)
    : data_desc_{ph.data_desc_}, options_{ph.options_} {

  hierarchy_ = std::make_shared<SAMRAI::hier::PatchHierarchy>(MakeUniqueName(), ph.GetNative()->getGridGeometry());

  hierarchy_->setMaxNumberOfLevels(ph.GetNative()->getMaxNumberOfLevels());
  for (int i = 1; i < options_.max_number_of_levels; ++i) {
    hierarchy_->setRatioToCoarserLevel(
        ph.GetNative()->getRatioToCoarserLevel(i), i);
  }
  for (int i = 0; i < ph.GetNative()->getMaxNumberOfLevels(); ++i) {
    if (ph.GetNative()->levelExists(i)) {
      hierarchy_->makeNewPatchLevel(
          i, *ph.GetNative()->getPatchLevel(i)->getBoxLevel());
      const SAMRAI::hier::PatchLevel& old_level = *ph.GetNative()->getPatchLevel(i);
      SAMRAI::hier::PatchLevel& new_level = *hierarchy_->getPatchLevel(i);
      new_level.allocatePatchData(SelectComponents(data_desc_.data_ids));
      FUB_ASSERT(new_level.getLocalNumberOfPatches() == old_level.getLocalNumberOfPatches());
      const std::size_t n_patches = new_level.getLocalNumberOfPatches();
      for (std::size_t p = 0; p < n_patches; ++p) {
        SAMRAI::hier::Patch& new_patch = *new_level.getPatch(p);
        const SAMRAI::hier::Patch& old_patch = *old_level.getPatch(p);
        for (int data_id : data_desc_.data_ids) {
          new_patch.getPatchData(data_id)->copy(*old_patch.getPatchData(data_id));
        }
      }
    }
  }
}

} // namespace fub::samrai