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

#ifndef FUB_GEOMETRY_POLYMORPHIC_GEOMETRY_HPP
#define FUB_GEOMETRY_POLYMORPHIC_GEOMETRY_HPP

#include "fub/core/type_traits.hpp"
#include "fub/geometry/Geometry.hpp"

#include <memory>

namespace fub {

class PolymorphicGeometry : public Geometry {
public:
  PolymorphicGeometry(const PolymorphicGeometry& other)
      : base_{other.Clone()} {}

  template <
      typename G,
      std::enable_if_t<!std::is_same<G, PolymorphicGeometry>::value>* = nullptr>
  PolymorphicGeometry(const G& geometry)
      : base_{std::make_unique<G>(geometry)} {
  }

  double ComputeDistanceTo(const Coordinates& x) const override {
    return base_->ComputeDistanceTo(x);
  }

  std::unique_ptr<Geometry> Clone() const override { return base_->Clone(); }

  Geometry& Base() noexcept { return *base_; }
  const Geometry& Base() const noexcept { return *base_; }

private:
  std::unique_ptr<Geometry> base_;
};

} // namespace fub

#endif