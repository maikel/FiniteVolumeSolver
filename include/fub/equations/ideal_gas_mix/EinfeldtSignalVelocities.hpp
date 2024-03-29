// Copyright (c) 2020 Maikel Nadolski
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

#ifndef FUB_EQUATONS_IDEAL_GAS_MIX_EINFELDT_SIGNALS_HPP
#define FUB_EQUATONS_IDEAL_GAS_MIX_EINFELDT_SIGNALS_HPP

#include "fub/equations/IdealGasMix.hpp"
#include "fub/EinfeldtSignalVelocities.hpp"

namespace fub {

template <int Dim> struct EinfeldtSignalVelocities<IdealGasMix<Dim>> {
  using Complete = typename IdealGasMix<Dim>::Complete;
  using CompleteArray = typename IdealGasMix<Dim>::CompleteArray;

  std::array<double, 2> operator()(const IdealGasMix<Dim>& equation,
                                   const Complete& left, const Complete& right,
                                   Direction dir) const noexcept;

  std::array<Array1d, 2> operator()(const IdealGasMix<Dim>& equation,
                                    const CompleteArray& left,
                                    const CompleteArray& right,
                                    Direction dir) const noexcept;

  std::array<Array1d, 2> operator()(const IdealGasMix<Dim>& equation,
                                    const CompleteArray& left,
                                    const CompleteArray& right,
                                    MaskArray mask,
                                    Direction dir) const noexcept;
};

extern template struct EinfeldtSignalVelocities<IdealGasMix<1>>;
extern template struct EinfeldtSignalVelocities<IdealGasMix<2>>;
extern template struct EinfeldtSignalVelocities<IdealGasMix<3>>;

}

#endif