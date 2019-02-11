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

#include "fub/Eigen.hpp"
#include "fub/ideal_gas/IdealGasEquation.hpp"
#include "fub/ideal_gas/flux_method/HlleGodunovMethod.hpp"

#include <boost/hana/members.hpp>

#include <iostream>

template <typename... Vars>
using Complete = fub::ideal_gas::IdealGasEquation::Complete<Vars...>;

template <typename... Vars>
using FluxState = fub::ideal_gas::IdealGasEquation::Conservative<Vars...>;

int main() {
  using namespace fub;
  using namespace fub::ideal_gas;

  constexpr int n = 2;

  Eigen::Tensor<double, 3> density(n + 2, n, n);
  Eigen::Tensor<double, 4> momentum(n + 2, n, n, 3);
  Eigen::Tensor<double, 3> energy(n + 2, n, n);
  Eigen::Tensor<double, 3> pressure(n + 2, n, n);
  Eigen::Tensor<double, 3> temperature(n + 2, n, n);
  Eigen::Tensor<double, 3> speed_of_sound(n + 2, n, n);
  Eigen::Tensor<double, 4> species(n + 2, n, n, n);

  density.setConstant(1.0);
  momentum.setConstant(0.5);
  energy.setConstant(1.0 + 0.5 * 0.5);
  pressure.setConstant(1.0 / 1.4);
  species.setConstant(1.0);
  speed_of_sound.setConstant(1.4);

  auto state = boost::hana::make<fub::StateTag<Complete>>(
      ToTensorSpan(density), ToTensorSpan(momentum), ToTensorSpan(energy),
      ToTensorSpan(pressure), ToTensorSpan(temperature),
      ToTensorSpan(speed_of_sound), ToTensorSpan(species));

  auto eigen_state = Transform(
      [](auto variable) { return ToEigenTensorMap(variable); }, state);

  auto left = SubState(eigen_state, Direction::X, 0, -1);
  auto right = SubState(eigen_state, Direction::X, 1, -1);
  auto fluxes = ::fub::ideal_gas::HlleGodunovMethod::ComputeFluxes(left, right);
  std::cout << fluxes.density << '\n';
  std::cout << fluxes.momentum << '\n';
  std::cout << fluxes.energy << '\n';
}
