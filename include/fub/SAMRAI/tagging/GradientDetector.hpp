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

#ifndef FUB_TAGGING_GRADIENT_DETECTOR_HPP
#define FUB_TAGGING_GRADIENT_DETECTOR_HPP

#include "fub/SAMRAI/tagging/TaggingMethod.hpp"
#include "fub/SAMRAI/Direction.hpp"

#include "SAMRAI/tbox/Dimension.h"
#include "SAMRAI/hier/Patch.h"
#include "SAMRAI/hier/IntVector.h"

#include <vector>

namespace fub {

class GradientDetector : public TaggingMethod {
public:
  GradientDetector() = default;
  GradientDetector(std::vector<int> patch_data_ids,
                   std::vector<double> threshholds)
      : patch_data_ids_{std::move(patch_data_ids)}, threshholds_{std::move(
                                                        threshholds)} {}

  void fillTagsForRefinement(const SAMRAI::hier::Patch& patch, int tag_id,
                             double time_point) const override;

  SAMRAI::hier::IntVector
  GetStencilWidth(const SAMRAI::tbox::Dimension& dim) const override {
    return SAMRAI::hier::IntVector::getOne(dim);
  }

  void addGradientThreshhold(int patch_data_id, double threshhold) {
    patch_data_ids_.push_back(patch_data_id);
    try {
      threshholds_.push_back(threshhold);
    } catch (...) {
      patch_data_ids_.pop_back();
      throw;
    }
  }

private:
  std::vector<int> patch_data_ids_{};
  std::vector<double> threshholds_{};
};

} // namespace fub

#endif