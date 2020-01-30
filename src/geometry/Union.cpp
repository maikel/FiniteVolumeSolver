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

#include "fub/geometry/Union.hpp"

#include <algorithm>
#include <limits>
#include <numeric>

namespace fub {

template <std::size_t Rank>
PolymorphicUnion<Rank>::PolymorphicUnion(
    const std::vector<PolymorphicGeometry<Rank>>& geoms)
    : geoms_(geoms) {}

template <std::size_t Rank>
std::unique_ptr<Geometry<Rank>> PolymorphicUnion<Rank>::Clone() const {
  return std::make_unique<PolymorphicUnion<Rank>>(*this);
}

template <std::size_t Rank>
double PolymorphicUnion<Rank>::ComputeDistanceTo(
    const std::array<double, Rank>& x) const {
  return std::accumulate(
      geoms_.begin(), geoms_.end(), std::numeric_limits<double>::lowest(),
      [x](double max, const PolymorphicGeometry<Rank>& geom) {
        return std::max(max, geom.ComputeDistanceTo(x));
      });
}

template class PolymorphicUnion<2>;
template class PolymorphicUnion<3>;

} // namespace fub
