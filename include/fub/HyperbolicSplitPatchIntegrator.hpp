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

#ifndef FUB_SOLVER_DIMENSIONAL_SPLIT_HYPERBOLIC_TIME_INTEGRATOR_HPP
#define FUB_SOLVER_DIMENSIONAL_SPLIT_HYPERBOLIC_TIME_INTEGRATOR_HPP

#include "fub/Direction.hpp"
#include "fub/Duration.hpp"
#include "fub/Equation.hpp"
#include "fub/Execution.hpp"
#include "fub/ForEach.hpp"

namespace fub {

template <typename Eq> class HyperbolicSplitPatchIntegrator {
public:
  using Equation = Eq;
  using Conservative = ::fub::Conservative<Equation>;
  using ConservativeArray = ::fub::ConservativeArray<Equation>;

  HyperbolicSplitPatchIntegrator(const Equation& eq) : equation{eq} {}

  const Equation& GetEquation() const noexcept { return equation; }

  void UpdateConservatively(const View<Conservative>& next,
                            const View<const Conservative>& fluxes,
                            const View<const Conservative>& prev, Direction dir,
                            Duration dt, double dx);

  void UpdateConservatively(execution::SimdTag, const View<Conservative>& next,
                            const View<const Conservative>& fluxes,
                            const View<const Conservative>& prev, Direction dir,
                            Duration dt, double dx);

private:
  Equation equation;

  /// @{
  /// States used as a buffer for intermediate loads and stores
  Conservative next_state{equation};
  Conservative prev_state{equation};
  Conservative flux_left{equation};
  Conservative flux_right{equation};
  /// @}
};

template <typename Equation>
void HyperbolicSplitPatchIntegrator<Equation>::UpdateConservatively(
    const View<Conservative>& next, const View<const Conservative>& fluxes,
    const View<const Conservative>& prev, Direction dir, Duration dt,
    double dx) {
  const double lambda = dt.count() / dx;
  constexpr int Rank = Equation::Rank();
  constexpr std::size_t sRank = static_cast<std::size_t>(Rank);
  ForEachIndex(Shrink(Box<0>(fluxes), dir, {0, 1}), [&](const auto... is) {
    const std::array<std::ptrdiff_t, sRank> face_left{is...};
    const std::array<std::ptrdiff_t, sRank> face_right =
        Shift(face_left, dir, 1);
    // Load fluxes left and right at index
    Load(flux_left, fluxes, face_left);
    Load(flux_right, fluxes, face_right);
    // Load state at index
    const std::array<std::ptrdiff_t, sRank> cell = face_left;
    Load(prev_state, prev, cell);
    // Do The computation
    ForEachVariable(
        [lambda](auto&& next, auto prev, auto flux_left, auto flux_right) {
          next = prev + lambda * (flux_left - flux_right);
        },
        next_state, prev_state, flux_left, flux_right);
    // Store the result at index
    Store(next, next_state, cell);
  });
}

template <typename Equation>
void HyperbolicSplitPatchIntegrator<Equation>::UpdateConservatively(
    execution::SimdTag, const View<Conservative>& next,
    const View<const Conservative>& fluxes,
    const View<const Conservative>& prev, Direction dir, Duration dt,
    double dx) {
  const double lambda = dt.count() / dx;
  constexpr int Rank = Equation::Rank();
  ForEachVariable(
      [lambda,
       dir](const PatchDataView<double, Rank, layout_stride>& next,
            const PatchDataView<const double, Rank, layout_stride>& prev,
            const PatchDataView<const double, Rank, layout_stride>& flux) {
        const std::ptrdiff_t next_extent = next.Extents().extent(0);
        const std::ptrdiff_t next_stride = next.Stride(1);
        const std::ptrdiff_t prev_stride = prev.Stride(1);
        const std::ptrdiff_t flux_stride = flux.Stride(1);
        const std::ptrdiff_t flux_dir = flux.Stride(int(dir));
        double* n = next.Span().begin();
        double* end = next.Span().end();
        const double* p = prev.Span().begin();
        const double* f = flux.Span().begin();
        while (end - n >= next_extent) {
          using Row = Eigen::Map<Eigen::Array<double, Eigen::Dynamic, 1>>;
          Eigen::Map<Row> nextv(n, next_extent);
          Eigen::Map<const Row> prevv(p, next_extent);
          Eigen::Map<const Row> fL(f, next_extent);
          Eigen::Map<const Row> fR(f + flux_dir, next_extent);
          nextv = prevv + lambda * (fL - fR);
          n += next_stride;
          p += prev_stride;
          f += flux_stride;
        }
      },
      next, prev, fluxes);
}

} // namespace fub

#endif
