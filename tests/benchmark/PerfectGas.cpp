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

#include "fub/equations/PerfectGas.hpp"

#include "fub/flux_method/HllMethod.hpp"
#include "fub/flux_method/MusclHancockMethod.hpp"

#include <benchmark/benchmark.h>

static void BM_HLLE_Method_X(benchmark::State& state) {
  int width = state.range(0);

  std::vector<double> density_flux(width * width * width);
  std::vector<double> momentum_flux(width * width * width * 3);
  std::vector<double> energy_flux(width * width * width);

  std::vector<double> density((width + 1) * width * width);
  std::vector<double> energy((width + 1) * width * width);
  std::vector<double> pressure((width + 1) * width * width);
  std::vector<double> speed_of_sound((width + 1) * width * width);
  std::vector<double> momentum((width + 1) * width * width * 3);
  std::fill(density.begin(), density.end(), 1.0);
  std::iota(pressure.begin(), pressure.end(), 0.0);
  std::transform(
      pressure.begin(), pressure.end(), pressure.begin(),
      [&](double x) { return 1 + 0.1 * std::sin(x * 2 * M_PI / width); });
  std::transform(pressure.begin(), pressure.end(), density.begin(),
                 speed_of_sound.begin(),
                 [](double p, double rho) { return std::sqrt(1.4 * p / rho); });

  using namespace fub;
  PerfectGas<3> equation{};
  HllMethod hlle{equation, EinfeldtSignalVelocities<PerfectGas<3>>{}};

  BasicView<Conservative<PerfectGas<3>>> basic_fluxes;
  basic_fluxes.density = PatchDataView<double, 3>(
      mdspan<double, 3>(density_flux.data(), width, width, width), {});
  basic_fluxes.momentum = PatchDataView<double, 4>(
      mdspan<double, 4>(momentum_flux.data(), width, width, width, 3), {});
  basic_fluxes.energy = PatchDataView<double, 3>(
      mdspan<double, 3>(energy_flux.data(), width, width, width), {});

  View<Conservative<PerfectGas<3>>> fluxes =
      Subview(basic_fluxes, {{}, {width, width, width}});

  BasicView<const Complete<PerfectGas<3>>> basic_states;
  basic_states.density = PatchDataView<const double, 3>(
      mdspan<const double, 3>(density.data(), width + 1, width, width), {-1});
  basic_states.momentum = PatchDataView<const double, 4>(
      mdspan<const double, 4>(momentum.data(), width + 1, width, width, 3),
      {-1});
  basic_states.energy = PatchDataView<const double, 3>(
      mdspan<const double, 3>(energy.data(), width + 1, width, width), {-1});
  basic_states.pressure = PatchDataView<const double, 3>(
      mdspan<const double, 3>(pressure.data(), width + 1, width, width), {-1});
  basic_states.speed_of_sound = PatchDataView<const double, 3>(
      mdspan<const double, 3>(speed_of_sound.data(), width + 1, width, width),
      {-1});
  View<const Complete<PerfectGas<3>>> states =
      Subview(basic_states, {{-1, 0, 0}, {width, width, width}});

  static_assert(IsView<View<Complete<PerfectGas<3>>>>::value);
  for (auto _ : state) {
    hlle.ComputeNumericFluxes(fluxes, states, Duration(0.7), 1.0, Direction::X);
  }
}
BENCHMARK(BM_HLLE_Method_X)->Arg(8)->Arg(16)->Arg(32);

static void BM_HLLE_Method_Y(benchmark::State& state) {
  int width = state.range(0);

  std::vector<double> density_flux(width * width * width);
  std::vector<double> momentum_flux(width * width * width * 3);
  std::vector<double> energy_flux(width * width * width);

  std::vector<double> density((width + 1) * width * width);
  std::vector<double> energy((width + 1) * width * width);
  std::vector<double> pressure((width + 1) * width * width);
  std::vector<double> speed_of_sound((width + 1) * width * width);
  std::vector<double> momentum((width + 1) * width * width * 3);
  std::fill(density.begin(), density.end(), 1.0);
  std::iota(pressure.begin(), pressure.end(), 0.0);
  std::transform(
      pressure.begin(), pressure.end(), pressure.begin(),
      [&](double x) { return 1 + 0.1 * std::sin(x * 2 * M_PI / width); });
  std::transform(pressure.begin(), pressure.end(), density.begin(),
                 speed_of_sound.begin(),
                 [](double p, double rho) { return std::sqrt(1.4 * p / rho); });

  using namespace fub;
  PerfectGas<3> equation{};
  HllMethod hlle{equation, EinfeldtSignalVelocities<PerfectGas<3>>{}};

  BasicView<Conservative<PerfectGas<3>>> basic_fluxes;
  basic_fluxes.density = PatchDataView<double, 3>(
      mdspan<double, 3>(density_flux.data(), width, width, width), {});
  basic_fluxes.momentum = PatchDataView<double, 4>(
      mdspan<double, 4>(momentum_flux.data(), width, width, width, 3), {});
  basic_fluxes.energy = PatchDataView<double, 3>(
      mdspan<double, 3>(energy_flux.data(), width, width, width), {});

  View<Conservative<PerfectGas<3>>> fluxes =
      Subview(basic_fluxes, {{}, {width, width, width}});

  BasicView<const Complete<PerfectGas<3>>> basic_states;
  basic_states.density = PatchDataView<const double, 3>(
      mdspan<const double, 3>(density.data(), width, width + 1, width),
      {0, -1});
  basic_states.momentum = PatchDataView<const double, 4>(
      mdspan<const double, 4>(momentum.data(), width, width + 1, width, 3),
      {0, -1});
  basic_states.energy = PatchDataView<const double, 3>(
      mdspan<const double, 3>(energy.data(), width, width + 1, width), {0, -1});
  basic_states.pressure = PatchDataView<const double, 3>(
      mdspan<const double, 3>(pressure.data(), width, width + 1, width),
      {0, -1});
  basic_states.speed_of_sound = PatchDataView<const double, 3>(
      mdspan<const double, 3>(speed_of_sound.data(), width, width + 1, width),
      {0, -1});
  View<const Complete<PerfectGas<3>>> states =
      Subview(basic_states, {{0, -1, 0}, {width, width, width}});

  static_assert(IsView<View<Complete<PerfectGas<3>>>>::value);
  for (auto _ : state) {
    hlle.ComputeNumericFluxes(fluxes, states, Duration(0.7), 1.0, Direction::Y);
  }
}
BENCHMARK(BM_HLLE_Method_Y)->Arg(8)->Arg(16)->Arg(32);

static void BM_SIMD_HLLE_Method_X(benchmark::State& state) {
  int width = state.range(0);

  std::vector<double> density_flux(width * width * width);
  std::vector<double> momentum_flux(width * width * width * 3);
  std::vector<double> energy_flux(width * width * width);

  std::vector<double> density((width + 1) * width * width);
  std::vector<double> energy((width + 1) * width * width);
  std::vector<double> pressure((width + 1) * width * width);
  std::vector<double> speed_of_sound((width + 1) * width * width);
  std::vector<double> momentum((width + 1) * width * width * 3);
  std::fill(density.begin(), density.end(), 1.0);
  std::iota(pressure.begin(), pressure.end(), 0.0);
  std::transform(
      pressure.begin(), pressure.end(), pressure.begin(),
      [&](double x) { return 1 + 0.1 * std::sin(x * 2 * M_PI / width); });
  std::transform(pressure.begin(), pressure.end(), density.begin(),
                 speed_of_sound.begin(),
                 [](double p, double rho) { return std::sqrt(1.4 * p / rho); });

  using namespace fub;
  PerfectGas<3> equation{};
  HllMethod hlle{equation, EinfeldtSignalVelocities<PerfectGas<3>>{}};

  BasicView<Conservative<PerfectGas<3>>> basic_fluxes;
  basic_fluxes.density = PatchDataView<double, 3>(
      mdspan<double, 3>(density_flux.data(), width, width, width), {});
  basic_fluxes.momentum = PatchDataView<double, 4>(
      mdspan<double, 4>(momentum_flux.data(), width, width, width, 3), {});
  basic_fluxes.energy = PatchDataView<double, 3>(
      mdspan<double, 3>(energy_flux.data(), width, width, width), {});

  View<Conservative<PerfectGas<3>>> fluxes =
      Subview(basic_fluxes, {{}, {width, width, width}});

  BasicView<const Complete<PerfectGas<3>>> basic_states;
  basic_states.density = PatchDataView<const double, 3>(
      mdspan<const double, 3>(density.data(), width + 1, width, width), {-1});
  basic_states.momentum = PatchDataView<const double, 4>(
      mdspan<const double, 4>(momentum.data(), width + 1, width, width, 3),
      {-1});
  basic_states.energy = PatchDataView<const double, 3>(
      mdspan<const double, 3>(energy.data(), width + 1, width, width), {-1});
  basic_states.pressure = PatchDataView<const double, 3>(
      mdspan<const double, 3>(pressure.data(), width + 1, width, width), {-1});
  basic_states.speed_of_sound = PatchDataView<const double, 3>(
      mdspan<const double, 3>(speed_of_sound.data(), width + 1, width, width),
      {-1});
  View<const Complete<PerfectGas<3>>> states =
      Subview(basic_states, {{-1, 0, 0}, {width, width, width}});

  static_assert(IsView<View<Complete<PerfectGas<3>>>>::value);
  for (auto _ : state) {
    hlle.ComputeNumericFluxes(execution::simd, fluxes, states, Duration(0.7),
                              1.0, Direction::X);
  }
}
BENCHMARK(BM_SIMD_HLLE_Method_X)->Arg(8)->Arg(16)->Arg(32);

static void BM_SIMD_HLLE_Method_Y(benchmark::State& state) {
  int width = state.range(0);

  std::vector<double> density_flux(width * width * width);
  std::vector<double> momentum_flux(width * width * width * 3);
  std::vector<double> energy_flux(width * width * width);

  std::vector<double> density((width + 1) * width * width);
  std::vector<double> energy((width + 1) * width * width);
  std::vector<double> pressure((width + 1) * width * width);
  std::vector<double> speed_of_sound((width + 1) * width * width);
  std::vector<double> momentum((width + 1) * width * width * 3);
  std::fill(density.begin(), density.end(), 1.0);
  std::iota(pressure.begin(), pressure.end(), 0.0);
  std::transform(
      pressure.begin(), pressure.end(), pressure.begin(),
      [&](double x) { return 1 + 0.1 * std::sin(x * 2 * M_PI / width); });
  std::transform(pressure.begin(), pressure.end(), density.begin(),
                 speed_of_sound.begin(),
                 [](double p, double rho) { return std::sqrt(1.4 * p / rho); });

  using namespace fub;
  PerfectGas<3> equation{};
  HllMethod hlle{equation, EinfeldtSignalVelocities<PerfectGas<3>>{}};

  BasicView<Conservative<PerfectGas<3>>> basic_fluxes;
  basic_fluxes.density = PatchDataView<double, 3>(
      mdspan<double, 3>(density_flux.data(), width, width, width), {});
  basic_fluxes.momentum = PatchDataView<double, 4>(
      mdspan<double, 4>(momentum_flux.data(), width, width, width, 3), {});
  basic_fluxes.energy = PatchDataView<double, 3>(
      mdspan<double, 3>(energy_flux.data(), width, width, width), {});

  View<Conservative<PerfectGas<3>>> fluxes =
      Subview(basic_fluxes, {{}, {width, width, width}});

  BasicView<const Complete<PerfectGas<3>>> basic_states;
  basic_states.density = PatchDataView<const double, 3>(
      mdspan<const double, 3>(density.data(), width, width + 1, width),
      {0, -1});
  basic_states.momentum = PatchDataView<const double, 4>(
      mdspan<const double, 4>(momentum.data(), width, width + 1, width, 3),
      {0, -1});
  basic_states.energy = PatchDataView<const double, 3>(
      mdspan<const double, 3>(energy.data(), width, width + 1, width), {0, -1});
  basic_states.pressure = PatchDataView<const double, 3>(
      mdspan<const double, 3>(pressure.data(), width, width + 1, width),
      {0, -1});
  basic_states.speed_of_sound = PatchDataView<const double, 3>(
      mdspan<const double, 3>(speed_of_sound.data(), width, width + 1, width),
      {0, -1});
  View<const Complete<PerfectGas<3>>> states =
      Subview(basic_states, {{0, -1, 0}, {width, width, width}});

  static_assert(IsView<View<Complete<PerfectGas<3>>>>::value);
  for (auto _ : state) {
    hlle.ComputeNumericFluxes(execution::simd, fluxes, states, Duration(0.7),
                              1.0, Direction::Y);
  }
}
BENCHMARK(BM_SIMD_HLLE_Method_Y)->Arg(8)->Arg(16)->Arg(32);

BENCHMARK_MAIN();
