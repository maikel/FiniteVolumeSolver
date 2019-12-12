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

#include "fub/flux_method/HllMethod.hpp"
#include "fub/flux_method/MusclHancockMethod.hpp"

#include <benchmark/benchmark.h>

static void BM_HLL_CompueNumericFlux(benchmark::State& state) {
  using namespace fub;
  ShallowWater equation{};
  Hll<ShallowWater, ShallowWaterSignalVelocities> method{
      equation, ShallowWaterSignalVelocities{}};

  std::array<Complete<ShallowWater>, 2> states;
  states[0].height = 1.0;
  states[0].momentum << 0.0, 0.0;
  states[1].height = 1.0;
  states[1].momentum << 0.0, 0.0;

  Conservative<ShallowWater> numeric_flux;
  for (auto _ : state) {
    for (int i = 0; i < kDefaultChunkSize; ++i) {
      method.ComputeNumericFlux(numeric_flux, states, Duration(1.0), 1.0,
                                Direction::X);
      benchmark::DoNotOptimize(numeric_flux);
    }
  }
}
BENCHMARK(BM_HLL_CompueNumericFlux);

static void BM_SIMD_HLL_CompueNumericFlux(benchmark::State& state) {
  using namespace fub;
  ShallowWater equation{};
  Hll<ShallowWater, ShallowWaterSignalVelocities> method{
      equation, ShallowWaterSignalVelocities{}};

  std::array<CompleteArray<ShallowWater>, 2> states;
  states[0].height = Array1d::Constant(1.0);
  states[0].momentum = Array2d::Constant(0.0);
  states[1].height = Array1d::Constant(1.0);
  states[1].momentum = Array2d::Constant(0.0);

  ConservativeArray<ShallowWater> numeric_flux;
  for (auto _ : state) {
    method.ComputeNumericFlux(numeric_flux, states, Duration(1.0), 1.0,
                              Direction::X);
    benchmark::DoNotOptimize(numeric_flux);
  }
}
BENCHMARK(BM_SIMD_HLL_CompueNumericFlux);

static void BM_GodunovMethod_X(benchmark::State& state) {
  int width = state.range(0);
  std::vector<double> height_flux(width * width);
  std::vector<double> momentum_flux((width + 1) * width * 2);
  std::vector<double> height((width + 1) * width);
  std::vector<double> momentum((width + 1) * width * 2);
  std::iota(height.begin(), height.end(), 0.0);
  std::transform(height.begin(), height.end(), height.begin(), [&](double x) {
    return 1 + 0.1 * std::sin(x * 2 * M_PI / width);
  });

  using namespace fub;
  ShallowWater equation{};
  HllMethod godunov{equation, ShallowWaterSignalVelocities{}};

  BasicView<Conservative<ShallowWater>> basic_fluxes;
  basic_fluxes.height = PatchDataView<double, 2>(
      mdspan<double, 2>(height_flux.data(), width, width), {});
  basic_fluxes.momentum = PatchDataView<double, 3>(
      mdspan<double, 3>(momentum_flux.data(), width, width, 2), {});
  View<Conservative<ShallowWater>> fluxes =
      Subview(basic_fluxes, {{}, {width, width}});

  BasicView<const Complete<ShallowWater>> basic_states;
  basic_states.height = PatchDataView<const double, 2>(
      mdspan<const double, 2>(height.data(), width + 1, width), {-1});
  basic_states.momentum = PatchDataView<const double, 3>(
      mdspan<const double, 3>(momentum.data(), width + 1, width, 2), {-1});
  View<const Complete<ShallowWater>> states =
      Subview(basic_states, {{-1, 0}, {width, width}});

  static_assert(IsView<View<Complete<ShallowWater>>>::value);
  for (auto _ : state) {
    godunov.ComputeNumericFluxes(fluxes, states, Duration(0.7),
                                 1.0, Direction::X);
  }
}
BENCHMARK(BM_GodunovMethod_X)->Range(8, 512);

static void BM_SIMD_GodunovMethod_X(benchmark::State& state) {
  int width = state.range(0);
  std::vector<double> height_flux(width * width);
  std::vector<double> momentum_flux((width + 1) * width * 2);
  std::vector<double> height((width + 1) * width);
  std::vector<double> momentum((width + 1) * width * 2);
  std::iota(height.begin(), height.end(), 0.0);
  std::transform(height.begin(), height.end(), height.begin(), [&](double x) {
    return 1 + 0.1 * std::sin(x * 2 * M_PI / width);
  });

  using namespace fub;
  ShallowWater equation{};
  HllMethod method{equation, ShallowWaterSignalVelocities{}};

  BasicView<Conservative<ShallowWater>> basic_fluxes;
  basic_fluxes.height = PatchDataView<double, 2>(
      mdspan<double, 2>(height_flux.data(), width, width), {});
  basic_fluxes.momentum = PatchDataView<double, 3>(
      mdspan<double, 3>(momentum_flux.data(), width, width, 2), {});
  View<Conservative<ShallowWater>> fluxes =
      Subview(basic_fluxes, {{}, {width, width}});

  BasicView<const Complete<ShallowWater>> basic_states;
  basic_states.height = PatchDataView<const double, 2>(
      mdspan<const double, 2>(height.data(), width + 1, width), {-1});
  basic_states.momentum = PatchDataView<const double, 3>(
      mdspan<const double, 3>(momentum.data(), width + 1, width, 2), {-1});
  View<const Complete<ShallowWater>> states =
      Subview(basic_states, {{-1, 0}, {width, width}});

  static_assert(IsView<View<Complete<ShallowWater>>>::value);
  for (auto _ : state) {
    method.ComputeNumericFluxes(execution::simd, fluxes, states,
                                Duration(0.7), 1.0, Direction::X);
  }
}
BENCHMARK(BM_SIMD_GodunovMethod_X)->Range(8, 512);

static void BM_SIMD_GodunovMethod_Y(benchmark::State& state) {
  int width = state.range(0);
  std::vector<double> height_flux(width * width);
  std::vector<double> momentum_flux((width + 1) * width * 2);
  std::vector<double> height((width + 1) * width);
  std::vector<double> momentum((width + 1) * width * 2);
  std::iota(height.begin(), height.end(), 0.0);
  std::transform(height.begin(), height.end(), height.begin(), [&](double x) {
    return 1 + 0.1 * std::sin(x * 2 * M_PI / width);
  });

  using namespace fub;
  ShallowWater equation{};
  HllMethod method{equation, ShallowWaterSignalVelocities{}};

  BasicView<Conservative<ShallowWater>> basic_fluxes;
  basic_fluxes.height = PatchDataView<double, 2>(
      mdspan<double, 2>(height_flux.data(), width, width), {});
  basic_fluxes.momentum = PatchDataView<double, 3>(
      mdspan<double, 3>(momentum_flux.data(), width, width, 2), {});
  View<Conservative<ShallowWater>> fluxes =
      Subview(basic_fluxes, {{}, {width, width}});

  BasicView<const Complete<ShallowWater>> basic_states;
  basic_states.height = PatchDataView<const double, 2>(
      mdspan<const double, 2>(height.data(), width, width + 1), {-1});
  basic_states.momentum = PatchDataView<const double, 3>(
      mdspan<const double, 3>(momentum.data(), width, width + 1, 2), {-1});
  View<const Complete<ShallowWater>> states =
      Subview(basic_states, {{-1, 0}, {width, width}});

  static_assert(IsView<View<Complete<ShallowWater>>>::value);
  for (auto _ : state) {
    method.ComputeNumericFluxes(execution::simd, fluxes, states,
                                Duration(0.7), 1.0, Direction::Y);
  }
}
BENCHMARK(BM_SIMD_GodunovMethod_Y)->Range(8, 512);

BENCHMARK_MAIN();
