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
  const Vc::Vector<double> prev(in, mask);
  const Vc::Vector<double> fL(fluxL, mask);
  const Vc::Vector<double> fR(fluxR, mask);
  const Vc::Vector<double> next = prev + dt_over_dx * (fL - fR);
  next.store(out, mask, Vc::Unaligned);
}
} // namespace

template <typename Tag>
HyperbolicPatchIntegrator::HyperbolicPatchIntegrator(Tag) {}

template <>
void HyperbolicPatchIntegrator<exeuction::SimdTag>::UpdateConservatively(
    const PatchDataView<double, 2>& next,
    const PatchDataView<const double, 2>& prev,
    const PatchDataView<const double, 2>& fluxes, Duration dt, double dx,
    Direction dir) {
  for (int comp = 0; comp < next.Extent(1); ++comp) {
    const PatchDataView<double, 1> n = Slice(next, all, comp);
    const PatchDataView<double, 1> p = Slice(prev, all, comp);
    const IndexBox<1> facesL = n.Box();
    const IndexBox<1> facesR = Grow(cells, dir, {-1, 1});
    const PatchDataView<double, 1> fL =
        Slice(fluxes.Subview(facesL), all, comp);
    const PatchDataView<double, 1> fR =
        Slice(fluxes.Subview(facesR), all, comp);
    const double dt_dx = dt.count() / dx;
    UpdateConservatively_SIMD(n.Span(), p.Span(), fL.Span(), fR.Span(), dt_dx);
  }
}

} // namespace fub

} // namespace fub
