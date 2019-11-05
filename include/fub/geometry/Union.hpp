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

#ifndef FUB_GEOMETRY_UNION_HPP
#define FUB_GEOMETRY_UNION_HPP

#include "fub/geometry/PolymorphicGeometry.hpp"

#include <vector>
#include <memory>

namespace fub {

template <std::size_t Rank>
class PolymorphicUnion : public Geometry<Rank> {
public:
  PolymorphicUnion(const std::vector<PolymorphicGeometry<Rank>>& geoms);

  std::unique_ptr<Geometry<Rank>> Clone() const override;

  double ComputeDistanceTo(const std::array<double, Rank>& x) const override;

private:
  std::vector<PolymorphicGeometry<Rank>> geoms_; 
};

extern template class PolymorphicUnion<2>;
extern template class PolymorphicUnion<3>;

}

#endif