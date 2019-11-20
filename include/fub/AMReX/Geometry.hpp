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

#ifndef FUB_AMREX_GEOMETRY_HPP
#define FUB_AMREX_GEOMETRY_HPP

#include <AMReX.H>

#include <array>

namespace fub::amrex {

template <typename Geom> class Geometry : private Geom {
public:
  Geometry(const Geom& base) : Geom(base) {}
  Geometry(Geom&& base) : Geom(std::move(base)) {}

  [[nodiscard]] Geom& Base() noexcept { return *this; }
  [[nodiscard]] const Geom& Base() const noexcept { return *this; }

  [[nodiscard]] double operator()(AMREX_D_DECL(double x, double y, double z)) const {
    return Geom::ComputeDistanceTo({AMREX_D_DECL(x, y, z)});
  }

  [[nodiscard]] double operator()(const std::array<double, AMREX_SPACEDIM>& coords) const {
    return Geom::ComputeDistanceTo(coords);
  }
};

} // namespace fub::amrex

#endif
