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

#ifndef FUB_SAMRAI_CARTESIAN_PATCH_HIERARCHY_HPP
#define FUB_SAMRAI_CARTESIAN_PATCH_HIERARCHY_HPP

#include "fub/CartesianCoordinates.hpp"

#include <SAMRAI/hier/PatchHierarchy.h>
#include <SAMRAI/geom/CartesianPatchGeometry.h>

#include <array>
#include <memory>

namespace SAMRAI {
namespace geom {

class CartesianPatchGeometry;

}
}

namespace fub {
namespace samrai {

struct CoordinatesRange {
  std::array<double, 3> lower;
  std::array<double, 3> upper;
};

std::shared_ptr<SAMRAI::hier::PatchHierarchy>
CartesianPatchHierarchy(const SAMRAI::hier::Box& box,
                        const CoordinatesRange& x_up);

SAMRAI::geom::CartesianPatchGeometry*
GetCartesianPatchGeometry(const SAMRAI::hier::Patch& patch);

CartesianCoordinates GetCartesianCoordinates(const SAMRAI::hier::Patch& patch);

struct CartesianPatchCoordinates {
  explicit CartesianPatchCoordinates(const SAMRAI::hier::Patch& patch);

  std::array<double, 3> operator()(const SAMRAI::hier::Index& index) const;

  SAMRAI::geom::CartesianPatchGeometry* geometry_;
  SAMRAI::hier::Box box_;
};

} // namespace samrai
} // namespace fub

#endif