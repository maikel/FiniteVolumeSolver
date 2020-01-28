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

#include "fub/HyperbolicPatchIntegrator.hpp"
#include "fub/ForEach.hpp"
#include "fub/StateRow.hpp"
#include "fub/ext/Vc.hpp"

#include <tuple>

namespace fub {
namespace {
template <typename... Pointers>
void AdvanceBy(std::ptrdiff_t n, Pointers&... ps) {
  ((ps += n), ...);
}

void UpdateConservatively_SIMD(span<double> nexts, span<const double> prevs,
                               span<const double> fluxes_left,
                               span<const double> fluxes_right, double dt_dx) {
  constexpr auto size = static_cast<std::ptrdiff_t>(Vc::Vector<double>::size());
  Vc::Vector<double> dt_over_dx(dt_dx);
  double* out = nexts.begin();
  double* end = nexts.end();
  const double* in = prevs.begin();
  const double* fluxL = fluxes_left.begin();
  const double* fluxR = fluxes_right.begin();
  std::ptrdiff_t n = end - out;
  while (n >= size) {
    const Vc::Vector<double> prev(in, Vc::Unaligned);
    const Vc::Vector<double> fL(fluxL, Vc::Unaligned);
    const Vc::Vector<double> fR(fluxR, Vc::Unaligned);
    const Vc::Vector<double> next = prev + dt_over_dx * (fL - fR);
    next.store(out, Vc::Unaligned);
    AdvanceBy(size, out, in, fluxL, fluxR);
    n = end - out;
  }
  const auto mask = Vc::Vector<double>([](int i) { return i; }) < int(n);
  const Vc::Vector<double> prev = mask_load(in, mask);
  const Vc::Vector<double> fL = mask_load(fluxL, mask);
  const Vc::Vector<double> fR = mask_load(fluxR, mask);
  const Vc::Vector<double> next = prev + dt_over_dx * (fL - fR);
  next.store(out, mask, Vc::Unaligned);
}

template <int Rank>
void UpdateConservatively_(execution::SimdTag,
                           const StridedDataView<double, Rank>& next,
                           const StridedDataView<const double, Rank>& prev,
                           const StridedDataView<const double, Rank>& fluxes,
                           Duration dt, double dx, Direction dir) {
  for (int comp = 0; comp < next.Extent(Rank - 1); ++comp) {
    const StridedDataView<double, Rank - 1> n = SliceLast(next, comp);
    const StridedDataView<const double, Rank - 1> p = SliceLast(prev, comp);
    const IndexBox<Rank> facesL = next.Box();
    const IndexBox<Rank> facesR = Grow(facesL, dir, {-1, 1});
    const StridedDataView<const double, Rank - 1> fL =
        SliceLast(fluxes.Subview(facesL), comp);
    const StridedDataView<const double, Rank - 1> fR =
        SliceLast(fluxes.Subview(facesR), comp);
    const double dt_dx = dt.count() / dx;
    ForEachRow(std::tuple{n, p, fL, fR},
               [dt_dx](span<double> n, span<const double> p,
                       span<const double> fL, span<const double> fR) {
                 UpdateConservatively_SIMD(n, p, fL, fR, dt_dx);
               });
  }
}

template <int Rank>
void UpdateConservatively_(execution::SequentialTag,
                           const StridedDataView<double, Rank>& next,
                           const StridedDataView<const double, Rank>& prev,
                           const StridedDataView<const double, Rank>& fluxes,
                           Duration dt, double dx, Direction dir) {
  double dt_dx = dt.count() / dx;
  ForEachIndex(next.Box(), [&](auto... is) {
    std::array<std::ptrdiff_t, sizeof...(is)> cell{is...};
    std::array<std::ptrdiff_t, sizeof...(is)> fL{is...};
    std::array<std::ptrdiff_t, sizeof...(is)> fR = Shift(fL, dir, 1);
    const double U_prev = prev(cell);
    const double f_left = fluxes(fL);
    const double f_right = fluxes(fR);
    const double U_next = U_prev + dt_dx * (f_left - f_right);
    next(cell) = U_next;
  });
}

template <int Rank>
void UpdateConservatively_(execution::OpenMpTag,
                           const StridedDataView<double, Rank>& next,
                           const StridedDataView<const double, Rank>& prev,
                           const StridedDataView<const double, Rank>& fluxes,
                           Duration dt, double dx, Direction dir) {
  return UpdateConservatively_(execution::seq, next, prev, fluxes, dt, dx, dir);
}

template <int Rank>
void UpdateConservatively_(execution::OpenMpSimdTag,
                           const StridedDataView<double, Rank>& next,
                           const StridedDataView<const double, Rank>& prev,
                           const StridedDataView<const double, Rank>& fluxes,
                           Duration dt, double dx, Direction dir) {
  return UpdateConservatively_(execution::simd, next, prev, fluxes, dt, dx,
                               dir);
}
} // namespace

template <typename Tag>
HyperbolicPatchIntegrator<Tag>::HyperbolicPatchIntegrator(Tag) {}

template <typename Tag>
void HyperbolicPatchIntegrator<Tag>::UpdateConservatively(
    const StridedDataView<double, 2>& next,
    const StridedDataView<const double, 2>& prev,
    const StridedDataView<const double, 2>& fluxes, Duration dt, double dx,
    Direction dir) {
  UpdateConservatively_(Tag(), next, prev, fluxes, dt, dx, dir);
}

template <typename Tag>
void HyperbolicPatchIntegrator<Tag>::UpdateConservatively(
    const StridedDataView<double, 3>& next,
    const StridedDataView<const double, 3>& prev,
    const StridedDataView<const double, 3>& fluxes, Duration dt, double dx,
    Direction dir) {
  UpdateConservatively_(Tag(), next, prev, fluxes, dt, dx, dir);
}

template <typename Tag>
void HyperbolicPatchIntegrator<Tag>::UpdateConservatively(
    const StridedDataView<double, 4>& next,
    const StridedDataView<const double, 4>& prev,
    const StridedDataView<const double, 4>& fluxes, Duration dt, double dx,
    Direction dir) {
  UpdateConservatively_(Tag(), next, prev, fluxes, dt, dx, dir);
}

template struct HyperbolicPatchIntegrator<execution::SimdTag>;
template struct HyperbolicPatchIntegrator<execution::SequentialTag>;
template struct HyperbolicPatchIntegrator<execution::OpenMpTag>;
template struct HyperbolicPatchIntegrator<execution::OpenMpSimdTag>;

} // namespace fub
