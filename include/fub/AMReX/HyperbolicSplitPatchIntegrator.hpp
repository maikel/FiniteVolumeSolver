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

#ifndef FUB_AMREX_HYPERBOLIC_SPLIT_PATCH_INTEGRATOR_HPP
#define FUB_AMREX_HYPERBOLIC_SPLIT_PATCH_INTEGRATOR_HPP

#include "fub/AMReX/PatchHandle.hpp"
#include "fub/Direction.hpp"
#include "fub/Duration.hpp"
#include "fub/State.hpp"

namespace fub {
namespace amrex {

template <typename Base> struct HyperbolicSplitPatchIntegrator : public Base {
  using Equation = typename Base::Equation;
  using Conservative = ::fub::Conservative<Equation>;
  using Complete = ::fub::Complete<Equation>;

  static constexpr int Rank = Equation::Rank();

  HyperbolicSplitPatchIntegrator(const Base& base) : Base(base) {}

  template <typename Context>
  void UpdateConservatively(Context& context, PatchHandle patch, Direction dir,
                            Duration dt) {
    const Equation& eq = Base::GetEquation();
    const int d = static_cast<int>(dir);

    IndexBox<Rank> cells =
        Grow(AsIndexBox<Rank>(patch.iterator->tilebox()), dir, {1, 1});
    IndexBox<Rank> faces =
        Grow(AsIndexBox<Rank>(patch.iterator->nodaltilebox(d)), dir, {1, 1});

    View<Conservative> scratch =
        Subview(AsCons(MakeView<BasicView<Complete>>(
                    context.GetScratch(patch, dir), eq)),
                cells);
    const double dx = context.GetDx(patch, dir);
    BasicView<const Conservative> basic_fluxes = AsConst(
        MakeView<BasicView<Conservative>>(context.GetFluxes(patch, dir), eq));
    View<const Conservative> fluxes = Subview(basic_fluxes, faces);
    Base::UpdateConservatively(scratch, fluxes, AsConst(scratch), dir, dt, dx);
  }
};

} // namespace amrex
} // namespace fub

#endif
