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

#ifndef FUB_BOUNDARY_CONDITION_REFLECTIVE_BOUNDARY_HPP
#define FUB_BOUNDARY_CONDITION_REFLECTIVE_BOUNDARY_HPP

#include "fub/State.hpp"
#include "fub/ext/Eigen.hpp"

namespace fub {
std::array<std::ptrdiff_t, 1> ReflectIndex(std::array<std::ptrdiff_t, 1> i,
                                           const IndexBox<1>& domain,
                                           Direction dir, int side);
std::array<std::ptrdiff_t, 2> ReflectIndex(std::array<std::ptrdiff_t, 2> i,
                                           const IndexBox<2>& domain,
                                           Direction dir, int side);
std::array<std::ptrdiff_t, 3> ReflectIndex(std::array<std::ptrdiff_t, 3> i,
                                           const IndexBox<3>& domain,
                                           Direction dir, int side);

template <typename Equation> class ReflectiveBoundary {
public:
  using Complete = ::fub::Complete<Equation>;
  static constexpr int Rank = Equation::Rank();

  ReflectiveBoundary(const Equation& equation);

  void FillBoundary(const View<Complete>& states,
                    const IndexBox<Rank>& box_to_fill, Direction dir, int side);

  const Equation& GetEquation() const noexcept { return equation_; }

private:
  Equation equation_;
  Complete state_{equation_};
  Complete reflected_{equation_};
};

template <typename Equation>
ReflectiveBoundary<Equation>::ReflectiveBoundary(const Equation& equation)
    : equation_{equation} {}

template <typename Equation>
void ReflectiveBoundary<Equation>::FillBoundary(
    const View<Complete>& states, const IndexBox<Rank>& box_to_fill,
    Direction dir, int side) {
  const Eigen::Matrix<double, Rank, 1> unit = UnitVector<Rank>(dir);
  ForEachIndex(box_to_fill, [&](auto... is) {
    std::array<std::ptrdiff_t, Rank> dest{is...};
    std::array<std::ptrdiff_t, Rank> src =
        ReflectIndex(dest, box_to_fill, dir, side);
    Load(state_, states, src);
    Reflect(reflected_, state_, unit, equation_);
    Store(states, reflected_, dest);
  });
}

} // namespace fub

#endif
