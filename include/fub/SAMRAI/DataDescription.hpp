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

#ifndef FUB_SAMRAI_DATA_DESCRIPTION_HPP
#define FUB_SAMRAI_DATA_DESCRIPTION_HPP

#include "fub/core/span.hpp"

#include <SAMRAI/hier/Variable.h>

#include <memory>
#include <vector>

namespace fub {
namespace samrai {

class DataDescription {
public:
  DataDescription(std::vector<int> state, std::vector<int> conservative);

  span<const int> state() const noexcept { return state_ids_; }
  span<const std::shared_ptr<SAMRAI::hier::Variable>> state_variables() const
      noexcept {
    return state_variables_;
  };

  span<const int> conservative() const noexcept { return conservative_ids_; }
  span<const std::shared_ptr<SAMRAI::hier::Variable>>
  conservative_variables() const noexcept {
    return cons_variables_;
  };

private:
  std::vector<int> state_ids_;
  std::vector<std::shared_ptr<SAMRAI::hier::Variable>> state_variables_;

  std::vector<int> conservative_ids_;
  std::vector<std::shared_ptr<SAMRAI::hier::Variable>> cons_variables_;
};

} // namespace samrai
} // namespace fub

#endif