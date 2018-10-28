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

#ifndef FUB_SAMRAI_UTILITY_HPP
#define FUB_SAMRAI_UTILITY_HPP

#include "fub/core/function_ref.hpp"

#include <array>
#include <memory>

#include "SAMRAI/geom/CartesianPatchGeometry.h"
#include "SAMRAI/hier/Box.h"
#include "SAMRAI/hier/Index.h"
#include "SAMRAI/hier/Patch.h"
#include "SAMRAI/hier/PatchHierarchy.h"

#include "boost/container/static_vector.hpp"

namespace fub {
template <typename F>
F forEachPatchLevel(const SAMRAI::hier::PatchHierarchy& hierarchy, F function) {
  const int max_level_number = hierarchy.getNumberOfLevels();
  for (int level_number = 0; level_number < max_level_number; ++level_number) {
    const SAMRAI::hier::PatchLevel& level =
        *hierarchy.getPatchLevel(level_number);
    function(level);
  }
  return function;
}

template <typename F>
F forEachPatch(const SAMRAI::hier::PatchHierarchy& hierarchy, F function) {
  forEachPatchLevel(hierarchy, [&](const SAMRAI::hier::PatchLevel& level) {
    for (const std::shared_ptr<SAMRAI::hier::Patch>& patch : level) {
      function(*patch);
    }
  });
  return function;
}

struct IndexRange {
  SAMRAI::hier::Index lower;
  SAMRAI::hier::Index upper;
};

using Coordinates =
    boost::container::static_vector<double, SAMRAI_MAXIMUM_DIMENSION>;

struct CoordinateRange {
  Coordinates lower;
  Coordinates upper;
};

std::shared_ptr<SAMRAI::hier::PatchHierarchy>
makeCartesianPatchHierarchy(const IndexRange& ir, const CoordinateRange& cr);

Coordinates
computeCellCoordinates(const SAMRAI::geom::CartesianPatchGeometry& geometry,
                       const SAMRAI::hier::Box& box,
                       const SAMRAI::hier::Index& index);

SAMRAI::geom::CartesianPatchGeometry*
getCartesianPatchGeometry(const SAMRAI::hier::Patch& patch);

} // namespace fub

#endif