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

#include "fub/SAMRAI/utility.hpp"
#include "fub/core/assert.hpp"

#include "SAMRAI/geom/CartesianGridGeometry.h"
#include "SAMRAI/hier/BlockId.h"
#include "SAMRAI/hier/BoxContainer.h"

#include "SAMRAI/mesh/CascadePartitioner.h"
#include "SAMRAI/mesh/GriddingAlgorithm.h"
#include "SAMRAI/mesh/StandardTagAndInitialize.h"
#include "SAMRAI/mesh/TileClustering.h"

#include "Eigen/Eigen"

namespace fub {

std::shared_ptr<SAMRAI::hier::PatchHierarchy>
MakeCartesianPatchHierarchy(const IndexRange& ir, const CoordinateRange& cr) {
  SAMRAI::hier::BoxContainer domain{
      SAMRAI::hier::Box(ir.lower, ir.upper, SAMRAI::hier::BlockId(0))};
  return std::make_shared<SAMRAI::hier::PatchHierarchy>(
      "Hierarchy", std::make_shared<SAMRAI::geom::CartesianGridGeometry>(
                       "Geometry", cr.lower.data(), cr.upper.data(), domain));
}



namespace {
using ArrayXi = Eigen::Array<int, SAMRAI_MAXIMUM_DIMENSION, 1>;
using ArrayXd = Eigen::Array<double, SAMRAI_MAXIMUM_DIMENSION, 1>;

Eigen::Map<const ArrayXi> asEigen(const SAMRAI::hier::Index& index) {
  return Eigen::Map<const ArrayXi>{&index[0]};
}
} // namespace

Coordinates
ComputeCellCoordinates(const SAMRAI::geom::CartesianPatchGeometry& geometry,
                       const SAMRAI::hier::Box& box,
                       const SAMRAI::hier::Index& index) {
  FUB_ASSERT(geometry.getDim().getValue() == box.getDim().getValue());
  FUB_ASSERT(index.getDim().getValue() == box.getDim().getValue());
  Eigen::Map<const ArrayXi> i_lo = asEigen(box.lower());
  Eigen::Map<const ArrayXi> i_up = asEigen(box.upper());
  auto&& i_loc = asEigen(index) - i_lo;
  const int size = geometry.getDim().getValue();
  auto&& di = i_up - i_lo;
  auto&& lambda = i_loc.cast<double>() / di.cast<double>();
  Eigen::Map<const ArrayXd> x_lo(geometry.getXLower());
  Eigen::Map<const ArrayXd> x_up(geometry.getXUpper());
  Eigen::Map<const ArrayXd> dx(geometry.getDx());
  auto&& dx_2 = dx / 2;
  Coordinates buffer(size);
  Eigen::Map<ArrayXd> result(buffer.data());
  result = (x_lo + dx_2) * (1.0 - lambda) + (x_up - dx_2) * lambda;
  return buffer;
}

SAMRAI::geom::CartesianPatchGeometry*
GetCartesianPatchGeometry(const SAMRAI::hier::Patch& patch) {
  return static_cast<SAMRAI::geom::CartesianPatchGeometry*>(
      patch.getPatchGeometry().get());
}

} // namespace fub