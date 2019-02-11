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

#include "fub/SAMRAI/CartesianPatchHierarchy.hpp"
#include "fub/Eigen.hpp"

#include <SAMRAI/geom/CartesianGridGeometry.h>
#include <SAMRAI/geom/CartesianPatchGeometry.h>

#include <array>

namespace fub {
namespace samrai {

std::shared_ptr<SAMRAI::hier::PatchHierarchy>
CartesianPatchHierarchy(const SAMRAI::hier::Box& box,
                        const CoordinatesRange& cr) {
  SAMRAI::hier::BoxContainer domain{box};
  return std::make_shared<SAMRAI::hier::PatchHierarchy>(
      "Hierarchy", std::make_shared<SAMRAI::geom::CartesianGridGeometry>(
                       "Geometry", cr.lower.data(), cr.upper.data(), domain));
}

namespace {
std::array<double, 3>
ComputeCellCoordinates(const SAMRAI::geom::CartesianPatchGeometry& geometry,
                       const SAMRAI::hier::Box& box,
                       const SAMRAI::hier::Index& index) {
  FUB_ASSERT(geometry.getDim().getValue() == box.getDim().getValue());
  FUB_ASSERT(index.getDim().getValue() == box.getDim().getValue());
  Eigen::Map<const Eigen::Array<int, 3, 1>> i_lo(&box.lower()[0]);
  Eigen::Map<const Eigen::Array<int, 3, 1>> i_up(&box.upper()[0]);
  Eigen::Map<const Eigen::Array<int, 3, 1>> i(&index[0]);
  auto&& i_loc = i - i_lo;
  const int size = geometry.getDim().getValue();
  auto&& di = i_up - i_lo;
  auto&& lambda = i_loc.cast<double>() / di.cast<double>();
  Eigen::Map<const Eigen::Array<double, 3, 1>> x_lo(geometry.getXLower(), size);
  Eigen::Map<const Eigen::Array<double, 3, 1>> x_up(geometry.getXUpper(), size);
  Eigen::Map<const Eigen::Array<double, 3, 1>> dx(geometry.getDx());
  auto&& dx_2 = dx / 2;
  std::array<double, 3> buffer;
  Eigen::Map<Eigen::Array<double, 3, 1>> result(buffer.data(), buffer.size());
  result = (x_lo + dx_2) * (1.0 - lambda) + (x_up - dx_2) * lambda;
  return buffer;
}
} // namespace

SAMRAI::geom::CartesianPatchGeometry*
GetCartesianPatchGeometry(const SAMRAI::hier::Patch& patch) {
  return static_cast<SAMRAI::geom::CartesianPatchGeometry*>(
      patch.getPatchGeometry().get());
}

CartesianCoordinates GetCartesianCoordinates(const SAMRAI::hier::Patch& patch) {
  auto geom = GetCartesianPatchGeometry(patch);
  Eigen::Vector3d lower = Eigen::Vector3d::Map(geom->getXLower());
  Eigen::Vector3d upper = Eigen::Vector3d::Map(geom->getXUpper());
  std::array<std::ptrdiff_t, 3> extents;
  extents.fill(1);
  for (int i = 0; i < patch.getDim().getValue(); ++i) {
    extents[i] = patch.getBox().numberCells(i);
  }
  return CartesianCoordinates(lower, upper, DynamicExtents<3>(extents));
}

CartesianPatchCoordinates::CartesianPatchCoordinates(
    const SAMRAI::hier::Patch& patch)
    : geometry_{GetCartesianPatchGeometry(patch)}, box_{patch.getBox()} {}

std::array<double, 3> CartesianPatchCoordinates::
operator()(const SAMRAI::hier::Index& index) const {
  return ComputeCellCoordinates(*geometry_, box_, index);
}

} // namespace samrai
} // namespace fub