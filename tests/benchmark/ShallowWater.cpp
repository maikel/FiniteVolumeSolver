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

#include "fub/equations/ShallowWater.hpp"

#include "fub/flux_method/GodunovMethod.hpp"
#include "fub/flux_method/MusclHancockMethod.hpp"

#include <benchmark/benchmark.h>

static void BM_GodunovMethod_X(benchmark::State& state) {
  int width = state.range(0);
  std::vector<double> heigth_flux(width * width);
  std::vector<double> momentum_flux((width + 1) * width * 2);
  std::vector<double> heigth((width + 1) * width);
  std::vector<double> momentum((width + 1) * width * 2);
  std::iota(heigth.begin(), heigth.end(), 0.0);
  std::transform(heigth.begin(), heigth.end(), heigth.begin(), [&](double x) {
    return 1 + 0.1 * std::sin(x * 2 * M_PI / width);
  });

  using namespace fub;
  ShallowWater equation{};
  GodunovMethod godunov{equation};

  View<Conservative<ShallowWater>> fluxes;
  fluxes.heigth = mdspan<double, 2>(heigth_flux.data(), width, width);
  fluxes.momentum = mdspan<double, 3>(momentum_flux.data(), width, width, 2);

  View<const Complete<ShallowWater>> states;
  states.heigth = mdspan<const double, 2>(heigth.data(), width + 1, width);
  states.momentum =
      mdspan<const double, 3>(momentum.data(), width + 1, width, 2);

  static_assert(IsView<View<Complete<ShallowWater>>>::value);
  for (auto _ : state) {
    godunov.ComputeNumericFluxes(fluxes, states, Duration(0.7), 1.0,
                                 Direction::X);
  }
}
BENCHMARK(BM_GodunovMethod_X)->Range(8, 512);

static void BM_GodunovMethod_Y(benchmark::State& state) {
  int width = state.range(0);
  std::vector<double> heigth_flux(width * width);
  std::vector<double> momentum_flux((width + 1) * width * 2);
  std::vector<double> heigth((width + 1) * width);
  std::vector<double> momentum((width + 1) * width * 2);
  std::iota(heigth.begin(), heigth.end(), 0.0);
  std::transform(heigth.begin(), heigth.end(), heigth.begin(), [&](double x) {
    return 1 + 0.1 * std::sin(x * 2 * M_PI / width);
  });

  using namespace fub;
  ShallowWater equation{};
  GodunovMethod godunov{equation};

  View<Conservative<ShallowWater>> fluxes;
  fluxes.heigth = mdspan<double, 2>(heigth_flux.data(), width, width);
  fluxes.momentum = mdspan<double, 3>(momentum_flux.data(), width, width, 2);

  View<const Complete<ShallowWater>> states;
  states.heigth = mdspan<const double, 2>(heigth.data(), width, width + 1);
  states.momentum =
      mdspan<const double, 3>(momentum.data(), width, width + 1, 2);

  static_assert(IsView<View<Complete<ShallowWater>>>::value);
  for (auto _ : state) {
    godunov.ComputeNumericFluxes(fluxes, states, Duration(0.7), 1.0,
                                 Direction::Y);
  }
}
BENCHMARK(BM_GodunovMethod_Y)->Range(8, 512);

static void BM_MusclHancockMethod_X(benchmark::State& state) {
  int width = state.range(0);
  std::vector<double> heigth_flux(width * width);
  std::vector<double> momentum_flux(width * width * 2);
  std::vector<double> heigth((width + 3) * width);
  std::vector<double> momentum((width + 3) * width * 2);
  std::iota(heigth.begin(), heigth.end(), 0.0);
  std::transform(heigth.begin(), heigth.end(), heigth.begin(), [&](double x) {
    return 1 + 0.1 * std::sin(x * 2 * M_PI / width);
  });

  using namespace fub;
  ShallowWater equation{};
  MusclHancockMethod method{equation};

  View<Conservative<ShallowWater>> fluxes;
  fluxes.heigth = mdspan<double, 2>(heigth_flux.data(), width, width);
  fluxes.momentum = mdspan<double, 3>(momentum_flux.data(), width, width, 2);

  View<Complete<ShallowWater>> states;
  states.heigth = mdspan<double, 2>(heigth.data(), width + 3, width);
  states.momentum = mdspan<double, 3>(momentum.data(), width + 3, width, 2);

  static_assert(IsView<View<Complete<ShallowWater>>>::value);
  for (auto _ : state) {
    method.ComputeNumericFluxes(fluxes, states, Duration(0.7), 1.0,
                                Direction::X);
  }
}
BENCHMARK(BM_MusclHancockMethod_X)->Range(8, 512);

static void BM_MusclHancockMethod_Y(benchmark::State& state) {
  int width = state.range(0);
  std::vector<double> heigth_flux(width * width);
  std::vector<double> momentum_flux(width * width * 2);
  std::vector<double> heigth((width + 3) * width);
  std::vector<double> momentum((width + 3) * width * 2);
  std::iota(heigth.begin(), heigth.end(), 0.0);
  std::transform(heigth.begin(), heigth.end(), heigth.begin(), [&](double x) {
    return 1 + 0.1 * std::sin(x * 2 * M_PI / width);
  });

  using namespace fub;
  ShallowWater equation{};
  MusclHancockMethod method{equation};

  View<Conservative<ShallowWater>> fluxes;
  fluxes.heigth = mdspan<double, 2>(heigth_flux.data(), width, width);
  fluxes.momentum = mdspan<double, 3>(momentum_flux.data(), width, width, 2);

  View<const Complete<ShallowWater>> states;
  states.heigth = mdspan<const double, 2>(heigth.data(), width, width + 3);
  states.momentum =
      mdspan<const double, 3>(momentum.data(), width, width + 3, 2);

  static_assert(IsView<View<Complete<ShallowWater>>>::value);
  for (auto _ : state) {
    method.ComputeNumericFluxes(fluxes, states, Duration(0.7), 1.0,
                                Direction::Y);
  }
}
BENCHMARK(BM_MusclHancockMethod_Y)->Range(8, 512);

BENCHMARK_MAIN();