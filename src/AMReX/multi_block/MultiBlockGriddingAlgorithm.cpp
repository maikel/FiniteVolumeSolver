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

#include "fub/AMReX/multi_block/MultiBlockGriddingAlgorithm.hpp"

namespace fub::amrex {

MultiBlockGriddingAlgorithm::MultiBlockGriddingAlgorithm(
    FlameMasterReactor reactor,
    std::vector<std::shared_ptr<GriddingAlgorithm>> tubes,
    std::vector<std::shared_ptr<cutcell::GriddingAlgorithm>> plena,
    std::vector<BlockConnection> connectivity)
    : reactor_{std::move(reactor)}, tubes_{std::move(tubes)}, plena_{std::move(
                                                                  plena)},
      connectivity_(std::move(connectivity)) {
  boundaries_.reserve(connectivity_.size());
  for (const BlockConnection& conn : connectivity_) {
    boundaries_.emplace_back(*this, conn, 3, reactor_);
  }
}

MultiBlockGriddingAlgorithm::MultiBlockGriddingAlgorithm(
    const MultiBlockGriddingAlgorithm& other)
    : reactor_{other.reactor_}, tubes_(other.tubes_.size()),
      plena_(other.plena_.size()), connectivity_(other.connectivity_),
      boundaries_(other.boundaries_) {
  for (std::size_t i = 0; i < tubes_.size(); ++i) {
    tubes_[i] = std::make_shared<GriddingAlgorithm>(*other.tubes_[i]);
  }
  for (std::size_t i = 0; i < plena_.size(); ++i) {
    plena_[i] = std::make_shared<cutcell::GriddingAlgorithm>(*other.plena_[i]);
  }
}

MultiBlockGriddingAlgorithm& MultiBlockGriddingAlgorithm::
operator=(const MultiBlockGriddingAlgorithm& other) {
  MultiBlockGriddingAlgorithm tmp(other);
  *this = std::move(tmp);
  return *this;
}

span<const std::shared_ptr<GriddingAlgorithm>>
MultiBlockGriddingAlgorithm::GetTubes() const noexcept {
  return tubes_;
}

span<const std::shared_ptr<cutcell::GriddingAlgorithm>>
MultiBlockGriddingAlgorithm::GetPlena() const noexcept {
  return plena_;
}

span<const BlockConnection> MultiBlockGriddingAlgorithm::GetConnectivity() const
    noexcept {
  return connectivity_;
}

void MultiBlockGriddingAlgorithm::RegridAllFinerLevels(int which_level) {
  for (const std::shared_ptr<GriddingAlgorithm>& tube : tubes_) {
    tube->RegridAllFinerlevels(which_level);
  }
  for (const std::shared_ptr<cutcell::GriddingAlgorithm>& plenum : plena_) {
    plenum->RegridAllFinerlevels(which_level);
  }
}

} // namespace fub::amrex
