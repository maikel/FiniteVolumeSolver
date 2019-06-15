#include "fub/HyperbolicPatchIntegrator.hpp"

#include <benchmark/benchmark.h>

#include <algorithm>
#include <random>

template <typename T, typename... Is>
fub::PatchDataView<T, sizeof...(Is)> MakePdv(fub::span<T> data, Is... extents) {
  fub::mdspan<T, sizeof...(Is)> mdspan(data.data(), extents...);
  return fub::PatchDataView<T, sizeof...(Is)>{mdspan, {}};
}

static void Sequential(benchmark::State& state) {
  const int width = state.range(0);
  const int components = 3;

  std::vector<double> next(width * width * width * components);
  std::vector<double> prev(width * width * width * components);
  std::vector<double> flux((width + 1) * width * width * components);

  static std::default_random_engine generator;
  static std::uniform_real_distribution<> dis(0.0, 1.0);
  std::generate(prev.begin(), prev.end(), [] { return dis(generator); });
  std::generate(flux.begin(), flux.end(), [] { return dis(generator); });

  fub::IndexBox<4> cells{{}, {width, width, width, components}};
  fub::IndexBox<4> faces{{}, {width + 1, width, width, components}};

  auto n =
      MakePdv(fub::span{next}, width, width, width, components).Subview(cells);
  auto p =
      MakePdv(fub::span{prev}, width, width, width, components).Subview(cells);
  auto f = MakePdv(fub::span{flux}, width + 1, width, width, components)
               .Subview(faces);

  auto integrator = fub::HyperbolicPatchIntegrator(fub::execution::seq);

  for (auto _ : state) {
    integrator.UpdateConservatively(n, p, f, fub::Duration(0.4), 0.1,
                                    fub::Direction::X);
  }
}
BENCHMARK(Sequential)->RangeMultiplier(2)->Range(2, 64);

static void Simd(benchmark::State& state) {
  const int width = state.range(0);
  const int components = 3;

  std::vector<double> next(width * width * width * components);
  std::vector<double> prev(width * width * width * components);
  std::vector<double> flux((width + 1) * width * width * components);

  static std::default_random_engine generator;
  static std::uniform_real_distribution<> dis(0.0, 1.0);
  std::generate(prev.begin(), prev.end(), [] { return dis(generator); });
  std::generate(flux.begin(), flux.end(), [] { return dis(generator); });

  fub::IndexBox<4> cells{{}, {width, width, width, components}};
  fub::IndexBox<4> faces{{}, {width + 1, width, width, components}};

  auto n =
      MakePdv(fub::span{next}, width, width, width, components).Subview(cells);
  auto p =
      MakePdv(fub::span{prev}, width, width, width, components).Subview(cells);
  auto f = MakePdv(fub::span{flux}, width + 1, width, width, components)
               .Subview(faces);

  auto integrator = fub::HyperbolicPatchIntegrator(fub::execution::simd);

  for (auto _ : state) {
    integrator.UpdateConservatively(n, p, f, fub::Duration(0.4), 0.1,
                                    fub::Direction::X);
  }
}
BENCHMARK(Simd)->RangeMultiplier(2)->Range(2, 64);

BENCHMARK_MAIN();