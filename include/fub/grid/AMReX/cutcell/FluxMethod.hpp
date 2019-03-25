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

#ifndef FUB_AMREX_CUTCELL_FLUX_METHOD_HPP
#define FUB_AMREX_CUTCELL_FLUX_METHOD_HPP

#include "fub/Direction.hpp"
#include "fub/Duration.hpp"
#include "fub/State.hpp"
#include "fub/grid/AMReX/PatchHandle.hpp"
#include "fub/grid/AMReX/ViewFArrayBox.hpp"

namespace fub {
namespace amrex {
namespace cutcell {

template <typename Base> struct FluxMethod : public Base {
  using Equation = typename Base::Equation;
  using Conservative = ::fub::Conservative<Equation>;
  using Complete = ::fub::Complete<Equation>;

  static constexpr int Rank = Equation::Rank();

  FluxMethod(const Base& base) : Base(base) {}

  const Equation& GetEquation() const { return Base::GetEquation(); }

  template <typename Context>
  double ComputeStableDt(Context& context, PatchHandle patch, Direction dir) {
    ::amrex::FabType type = context.GetCutCellPatchType(patch);
    if (type == ::amrex::FabType::covered) {
      return std::numeric_limits<double>::infinity();
    }
    const Equation& equation = Base::GetEquation();
    const double dx = context.GetDx(patch, dir);
    View<const Complete> scratch =
        MakeView<View<Complete>>(context.GetScratch(patch, dir), equation);
    if (type == ::amrex::FabType::regular) {
      return Base::ComputeStableDt(scratch, dx, dir);
    } else if (type == ::amrex::FabType::singlevalued) {
      CutCellData<Rank> cutcell_data = context.GetCutCellData(patch, dir);
      return Base::ComputeStableDt(scratch, cutcell_data, dir, dx);
    }
    return std::numeric_limits<double>::infinity();
  }

  template <typename Context>
  void ComputeNumericFluxes(Context& context, PatchHandle patch, Direction dir,
                            Duration dt) {
    ::amrex::FabType type = context.GetCutCellPatchType(
        patch, context.GetGhostCellWidth(patch, dir));
    if (type == ::amrex::FabType::covered) {
      return;
    }
    const Equation& equation = Base::GetEquation();
    const double dx = context.GetDx(patch, dir);
    View<const Complete> scratch =
        MakeView<View<Complete>>(context.GetScratch(patch, dir), equation);
    View<Conservative> regular_fluxes =
        MakeView<View<Conservative>>(context.GetFluxes(patch, dir), equation);

    if (type == ::amrex::FabType::regular) {
      Base::ComputeNumericFluxes(regular_fluxes, scratch, dir, dt, dx);
    } else if (type == ::amrex::FabType::singlevalued) {
      View<Conservative> boundary_fluxes = MakeView<View<Conservative>>(
          context.GetBoundaryFluxes(patch, dir), equation);
      CutCellData<Rank> cutcell_data = context.GetCutCellData(patch, dir);
      Base::ComputeBoundaryFluxes(boundary_fluxes, scratch, cutcell_data, dir,
                                  dt, dx);

      View<Conservative> stabilized_fluxes = MakeView<View<Conservative>>(
          context.GetStabilizedFluxes(patch, dir), equation);

      Base::ComputeCutCellFluxes(stabilized_fluxes, regular_fluxes,
                                 AsConst(boundary_fluxes), scratch,
                                 cutcell_data, dir, dt, dx);
    }
  }

  Complete ref_;
};

} // namespace cutcell
} // namespace amrex
} // namespace fub

#endif
