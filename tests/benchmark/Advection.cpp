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

#include "fub/equations/Advection.hpp"

#include "fub/flux_method/GodunovMethod.hpp"
#include "fub/flux_method/MusclHancockMethod.hpp"

#include <benchmark/benchmark.h>

static void BM_GodunovMethod_X(benchmark::State& state) {
  int width = state.range(0);
  std::vector<double> mass_flux(width * width);
  std::vector<double> mass((width + 1) * width);
  std::iota(mass_flux.begin(), mass_flux.end(), 0.0);
  std::iota(mass.begin(), mass.end(), 0.0);

  using namespace fub;
  Advection2d equation{{0.4, 0.8}};
  GodunovMethod godunov{equation};

  View<Conservative<Advection2d>> fluxes;
  fluxes.mass = mdspan<double, 2>(mass_flux.data(), width, width);

  View<Complete<Advection2d>> states;
  states.mass = mdspan<double, 2>(mass.data(), width + 1, width);

  static_assert(IsView<View<Complete<Advection2d>>>::value);
  for (auto _ : state) {
    godunov.ComputeNumericFluxes(fluxes, states, Duration(0.7), 1.0,
                                 Direction::X);
  }
}
BENCHMARK(BM_GodunovMethod_X)->Range(8, 512);

static void BM_GodunovMethod_Y(benchmark::State& state) {
  int width = state.range(0);
  std::vector<double> mass_flux(width * width);
  std::vector<double> mass((width + 1) * width);
  std::iota(mass_flux.begin(), mass_flux.end(), 0.0);
  std::iota(mass.begin(), mass.end(), 0.0);

  using namespace fub;
  Advection2d equation{{0.4, 0.8}};
  GodunovMethod godunov{equation};

  View<Conservative<Advection2d>> fluxes;
  fluxes.mass = mdspan<double, 2>(mass_flux.data(), width, width);

  View<const Complete<Advection2d>> states;
  states.mass = mdspan<const double, 2>(mass.data(), width, width + 1);

  static_assert(IsView<View<Complete<Advection2d>>>::value);
  for (auto _ : state) {
    godunov.ComputeNumericFluxes(fluxes, states, Duration(0.7), 1.0,
                                 Direction::Y);
  }
}
BENCHMARK(BM_GodunovMethod_Y)->Range(8, 512);

static void BM_SIMD_GodunovMethod_X(benchmark::State& state) {
  int width = state.range(0);
  std::vector<double> mass_flux(width * width);
  std::vector<double> mass((width + 1) * width);
  std::iota(mass_flux.begin(), mass_flux.end(), 0.0);
  std::iota(mass.begin(), mass.end(), 0.0);

  using namespace fub;
  Advection2d equation{{0.4, 0.8}};
  GodunovMethod godunov{equation};

  View<Conservative<Advection2d>> fluxes;
  fluxes.mass = mdspan<double, 2>(mass_flux.data(), width, width);

  View<Complete<Advection2d>> states;
  states.mass = mdspan<double, 2>(mass.data(), width + 1, width);

  static_assert(IsView<View<Complete<Advection2d>>>::value);
  for (auto _ : state) {
    godunov.ComputeNumericFluxes(execution::simd, fluxes, states, Duration(0.7),
                                 1.0, Direction::X);
  }
}
BENCHMARK(BM_SIMD_GodunovMethod_X)->Range(8, 512);

static void BM_SIMD_GodunovMethod_Y(benchmark::State& state) {
  int width = state.range(0);
  std::vector<double> mass_flux(width * width);
  std::vector<double> mass((width + 1) * width);
  std::iota(mass_flux.begin(), mass_flux.end(), 0.0);
  std::iota(mass.begin(), mass.end(), 0.0);

  using namespace fub;
  Advection2d equation{{0.4, 0.8}};
  GodunovMethod godunov{equation};

  View<Conservative<Advection2d>> fluxes;
  fluxes.mass = mdspan<double, 2>(mass_flux.data(), width, width);

  View<Complete<Advection2d>> states;
  states.mass = mdspan<double, 2>(mass.data(), width, width + 1);

  static_assert(IsView<View<Complete<Advection2d>>>::value);
  for (auto _ : state) {
    godunov.ComputeNumericFluxes(execution::simd, fluxes, states, Duration(0.7),
                                 1.0, Direction::Y);
  }
}
BENCHMARK(BM_SIMD_GodunovMethod_Y)->Range(8, 512);

static void BM_MusclHancockMethod_X(benchmark::State& state) {
  int width = state.range(0);
  std::vector<double> mass_flux(width * width);
  std::vector<double> mass((width + 3) * width);
  std::iota(mass_flux.begin(), mass_flux.end(), 0.0);
  std::iota(mass.begin(), mass.end(), 0.0);

  using namespace fub;
  Advection2d equation{{0.4, 0.8}};
  MusclHancockMethod method{equation};

  View<Conservative<Advection2d>> fluxes;
  fluxes.mass = mdspan<double, 2>(mass_flux.data(), width, width);

  View<Complete<Advection2d>> states;
  states.mass = mdspan<double, 2>(mass.data(), width + 3, width);

  static_assert(IsView<View<Complete<Advection2d>>>::value);
  for (auto _ : state) {
    method.ComputeNumericFluxes(fluxes, states, Duration(0.7), 1.0,
                                Direction::X);
  }
}
BENCHMARK(BM_MusclHancockMethod_X)->Range(8, 512);

static void BM_MusclHancockMethod_Y(benchmark::State& state) {
  int width = state.range(0);
  std::vector<double> mass_flux(width * width);
  std::vector<double> mass((width + 3) * width);
  std::iota(mass_flux.begin(), mass_flux.end(), 0.0);
  std::iota(mass.begin(), mass.end(), 0.0);

  using namespace fub;
  Advection2d equation{{0.4, 0.8}};
  MusclHancockMethod method{equation};

  View<Conservative<Advection2d>> fluxes;
  fluxes.mass = mdspan<double, 2>(mass_flux.data(), width, width);

  View<const Complete<Advection2d>> states;
  states.mass = mdspan<const double, 2>(mass.data(), width, width + 3);

  static_assert(IsView<View<Complete<Advection2d>>>::value);
  for (auto _ : state) {
    method.ComputeNumericFluxes(fluxes, states, Duration(0.7), 1.0,
                                Direction::Y);
  }
}
BENCHMARK(BM_MusclHancockMethod_Y)->Range(8, 512);

static void BM_SIMD_MusclHancockMethod_X(benchmark::State& state) {
  int width = state.range(0);
  std::vector<double> mass_flux(width * width);
  std::vector<double> mass((width + 3) * width);
  std::iota(mass_flux.begin(), mass_flux.end(), 0.0);
  std::iota(mass.begin(), mass.end(), 0.0);

  using namespace fub;
  Advection2d equation{{0.4, 0.8}};
  MusclHancockMethod godunov{equation};

  View<Conservative<Advection2d>> fluxes;
  fluxes.mass = mdspan<double, 2>(mass_flux.data(), width, width);

  View<Complete<Advection2d>> states;
  states.mass = mdspan<double, 2>(mass.data(), width + 3, width);

  static_assert(IsView<View<Complete<Advection2d>>>::value);
  for (auto _ : state) {
    godunov.ComputeNumericFluxes(execution::simd, fluxes, states, Duration(0.7),
                                 1.0, Direction::X);
  }
}
BENCHMARK(BM_SIMD_MusclHancockMethod_X)->Range(8, 512);

static void BM_SIMD_MusclHancockMethod_Y(benchmark::State& state) {
  int width = state.range(0);
  std::vector<double> mass_flux(width * width);
  std::vector<double> mass((width + 3) * width);
  std::iota(mass_flux.begin(), mass_flux.end(), 0.0);
  std::iota(mass.begin(), mass.end(), 0.0);

  using namespace fub;
  Advection2d equation{{0.4, 0.8}};
  MusclHancockMethod godunov{equation};

  View<Conservative<Advection2d>> fluxes;
  fluxes.mass = mdspan<double, 2>(mass_flux.data(), width, width);

  View<Complete<Advection2d>> states;
  states.mass = mdspan<double, 2>(mass.data(), width, width + 3);

  static_assert(IsView<View<Complete<Advection2d>>>::value);
  for (auto _ : state) {
    godunov.ComputeNumericFluxes(execution::simd, fluxes, states, Duration(0.7),
                                 1.0, Direction::Y);
  }
}
BENCHMARK(BM_SIMD_MusclHancockMethod_Y)->Range(8, 512);

BENCHMARK_MAIN();
