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

#ifndef FUB_AMREX_FLUX_METHOD_HPP
#define FUB_AMREX_FLUX_METHOD_HPP

#include "fub/Direction.hpp"
#include "fub/Duration.hpp"
#include "fub/State.hpp"
#include "fub/grid/AMReX/PatchHandle.hpp"
#include "fub/grid/AMReX/ViewFArrayBox.hpp"

namespace fub {
namespace amrex {

template <typename Base> struct FluxMethod : public Base {
  using Equation = typename Base::Equation;
  using Conservative = ::fub::Conservative<Equation>;
  using Complete = ::fub::Complete<Equation>;

  static const int Rank = Equation::Rank();

  FluxMethod(const Base& base) : Base(base) {}

  template <typename Context>
  double ComputeStableDt(Context& context, PatchHandle patch, Direction dir) {
    const Equation& equation = Base::GetEquation();
    const double dx = context.GetDx(patch, dir);
    BasicView<const Complete> scratch = AsConst(MakeView<BasicView<Complete>>(
        context.GetScratch(patch, dir), equation));
    static constexpr int Rank = Equation::Rank();
    const int gcw = context.GetGhostCellWidth(patch, dir);
    const IndexBox<Rank> tilebox =
        Grow(AsIndexBox<Rank>(patch.iterator->tilebox()), dir, {gcw, gcw});
    View<const Complete> subscratch = Subview(scratch, tilebox);
    return Base::ComputeStableDt(subscratch, dx, dir);
  }

  template <typename Context>
  void ComputeNumericFluxes(Context& context, PatchHandle patch, Direction dir,
                            Duration dt) {
    const int d = static_cast<int>(dir);
    const int gcw = context.GetGhostCellWidth(patch, dir);

    ::amrex::IntVect gcws{};
    gcws[d] = gcw;
    const IndexBox<Rank> tilebox_cells =
        AsIndexBox<Rank>(patch.iterator->growntilebox(gcws));

    gcws[d] = 1;
    const IndexBox<Rank> tilebox_faces =
        AsIndexBox<Rank>(patch.iterator->grownnodaltilebox(d, gcws));

    const double dx = context.GetDx(patch, dir);

    View<const Complete> scratch =
        Subview(AsConst(MakeView<BasicView<Complete>>(
                    context.GetScratch(patch, dir), Base::GetEquation())),
                tilebox_cells);

    View<Conservative> fluxes =
        Subview(MakeView<BasicView<Conservative>>(context.GetFluxes(patch, dir),
                                                  Base::GetEquation()),
                tilebox_faces);

    Base::ComputeNumericFluxes(fluxes, scratch, dir, dt, dx);
  }
};

} // namespace amrex
} // namespace fub

#endif
