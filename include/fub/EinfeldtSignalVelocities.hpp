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

#ifndef FUB_EINFELDT_SIGNAL_VELOCITIES_HPP
#define FUB_EINFELDT_SIGNAL_VELOCITIES_HPP

#include "fub/Direction.hpp"

namespace fub {

/// \defgroup EinfeldtSignalVelocities
/// This is a customization point for equations which can define two signal
/// velocities for usage with the \a HllMethod.
///
/// Semantically the velocities should be formed by the roe averaves of two
/// states.
template <typename Equation> struct EinfeldtSignalVelocitiesImpl {
  using Complete = typename Equation::Complete;

  static std::array<double, 2> apply(const Equation& equation,
                                     const Complete& left,
                                     const Complete& right, Direction dir);
};

template <typename Equation> struct EinfeldtSignalVelocitiesFn {
  using Complete = typename Equation::Complete;

  std::array<double, 2> operator()(const Equation& equation,
                                   const Complete& left, const Complete& right,
                                   Direction dir) const {
    return EinfeldtSignalVelocitiesImpl<Equation>::apply(equation, left, right,
                                                         dir);
  }
};

template <typename Equation>
inline constexpr EinfeldtSignalVelocitiesFn<Equation>
    EinfeldtSignalVelocities{};

} // namespace fub

#endif