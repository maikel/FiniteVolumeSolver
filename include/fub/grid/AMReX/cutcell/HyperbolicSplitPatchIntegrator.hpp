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

#include "fub/Direction.hpp"
#include "fub/Duration.hpp"
#include "fub/State.hpp"
#include "fub/grid/AMReX/PatchHandle.hpp"

namespace fub {
namespace amrex {
namespace cutcell {

template <typename Base> struct HyperbolicSplitPatchIntegrator : public Base {
  using Equation = typename Base::Equation;
  using Conservative = ::fub::Conservative<Equation>;
  using Complete = ::fub::Complete<Equation>;

  static constexpr int Rank = Equation::Rank();

  HyperbolicSplitPatchIntegrator(const Base& base) : Base(base) {}

  template <typename Context>
  void UpdateConservatively(Context& context, PatchHandle patch, Direction dir,
                            Duration dt) {
    const int gcw = context.GetGhostCellWidth(patch, dir);
    ::amrex::FabType type = context.GetCutCellPatchType(patch, gcw);
    if (type == ::amrex::FabType::covered) {
      return;
    }
    const Equation& equation = Base::GetEquation();
    const double dx = context.GetDx(patch, dir);
    const int d = static_cast<int>(dir);
    const ::amrex::IntVect gcws = ::amrex::IntVect::TheDimensionVector(d);
    const auto tilebox_cells = AsIndexBox(patch.iterator->growntilebox(gcws));
    const auto tilebox_faces =
        AsIndexBox(patch.iterator->grownnodaltilebox(d, gcws));

    StridedView<Conservative> scratch = Subview(
        MakeView<View<Complete>>(context.GetScratch(patch, dir), equation),
        tilebox_cells);
    StridedView<const Conservative> regular_fluxes = Subview(
        MakeView<View<Conservative>>(context.GetFluxes(patch, dir), equation),
        tilebox_faces);
    if (type == ::amrex::FabType::regular) {
      Base::UpdateConservatively(scratch, AsConst(scratch), regular_fluxes, dir,
                                 dt, dx);
    } else {
      FUB_ASSERT(type == ::amrex::FabType::singlevalued);
      CutCellData<Rank> cutcell_data = context.GetCutCellData(patch, dir);
      StridedView<const Conservative> stabilized_fluxes =
          Subview(MakeView<View<Conservative>>(
                      context.GetStabilizedFluxes(patch, dir), equation),
                  tilebox_faces);
      StridedView<const Conservative> boundary_fluxes =
          Subview(MakeView<View<Conservative>>(
                      context.GetBoundaryFluxes(patch, dir), equation),
                  tilebox_cells);
      Base::UpdateConservatively(scratch, AsConst(scratch), stabilized_fluxes,
                                 regular_fluxes, boundary_fluxes, cutcell_data,
                                 dir, dt, dx);
    }
  }
};

} // namespace cutcell
} // namespace amrex
} // namespace fub

#endif
