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
#include "fub/grid/AMReX/cutcell/PatchHierarchy.hpp"

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

  template <typename Context> void PreAdvanceHierarchy(Context& context) {
    const Equation& equation = Base::GetEquation();
    const PatchHierarchy& hierarchy = *context.GetPatchHierarchy();
    const int nlevels = hierarchy.GetNumberOfLevels();
    const int gcw = Base::GetStencilWidth() + 1;
    const int n_components = hierarchy.GetDataDescription().n_state_components;
    for (int level_num = 0; level_num < nlevels; ++level_num) {
      const PatchLevel& level = hierarchy.GetPatchLevel(level_num);
      ::amrex::MultiFab datas(level.data.boxArray(),
                              level.data.DistributionMap(), n_components, gcw,
                              ::amrex::MFInfo(), level.data.Factory());
      context.FillGhostLayer(datas, level_num);
      context.ForEachPatch(level_num, [&](PatchHandle patch) {
        ::amrex::FabType type = context.GetCutCellPatchType(patch, gcw);
        if (type == ::amrex::FabType::singlevalued) {
          CutCellData<Rank> cc_data =
              context.GetCutCellData(patch, Direction::X);
          const IndexBox<Rank> cells =
              AsIndexBox(patch.iterator->growntilebox(gcw));
          View<const Complete> data = AsConst(
              Subview(MakeView<BasicView<Complete>>(
                          MakePatchDataView(datas[*patch.iterator]), equation),
                      cells));
          View<Complete> reference_states =
              Subview(MakeView<BasicView<Complete>>(
                          context.GetReferenceStates(patch), equation),
                      cells);
          Base::PreAdvanceHierarchy(reference_states, data, cc_data);
        }
      });
    }
  }

  template <typename Context>
  double ComputeStableDt(Context& context, PatchHandle patch, Direction dir) {
    const int gcw = context.GetGhostCellWidth(patch, dir);
    ::amrex::FabType type = context.GetCutCellPatchType(patch, gcw);
    if (type == ::amrex::FabType::covered) {
      return std::numeric_limits<double>::infinity();
    }
    const Equation& equation = Base::GetEquation();
    const double dx = context.GetDx(patch, dir);
    const int d = static_cast<int>(dir);
    ::amrex::IntVect gcws{};
    gcws[d] = gcw;
    const auto tilebox_cells = AsIndexBox(patch.iterator->growntilebox(gcws));
    View<const Complete> scratch = AsConst(Subview(
        MakeView<BasicView<Complete>>(context.GetScratch(patch, dir), equation),
        tilebox_cells));
    //    BasicView<const Complete> scratch =
    //        AsConst(MakeView<BasicView<Complete>>(context.GetScratch(patch,
    //        dir), equation));
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
    const int gcw = context.GetGhostCellWidth(patch, dir);
    ::amrex::FabType type = context.GetCutCellPatchType(patch, gcw);
    if (type == ::amrex::FabType::covered) {
      return;
    }
    const Equation& equation = Base::GetEquation();
    const double dx = context.GetDx(patch, dir);

    const int d = static_cast<int>(dir);
    ::amrex::IntVect gcws{};
    gcws[d] = gcw;
    const auto tilebox_cells = AsIndexBox(patch.iterator->growntilebox(gcws));

    gcws[d] = 1;
    const auto tilebox_faces =
        AsIndexBox(patch.iterator->grownnodaltilebox(d, gcws));

    View<const Complete> scratch = AsConst(Subview(
        MakeView<BasicView<Complete>>(context.GetScratch(patch, dir), equation),
        tilebox_cells));

    View<Conservative> regular_fluxes =
        Subview(MakeView<BasicView<Conservative>>(context.GetFluxes(patch, dir),
                                                  equation),
                tilebox_faces);

    if (type == ::amrex::FabType::regular) {
      Base::ComputeNumericFluxes(regular_fluxes, scratch, dir, dt, dx);
    } else if (type == ::amrex::FabType::singlevalued) {
      View<Conservative> boundary_fluxes =
          Subview(MakeView<BasicView<Conservative>>(
                      context.GetBoundaryFluxes(patch, dir), equation),
                  tilebox_cells);

      View<const Complete> reference_states =
          AsConst(Subview(MakeView<BasicView<Complete>>(
                              context.GetReferenceStates(patch), equation),
                          tilebox_cells));

      CutCellData<Rank> cutcell_data = context.GetCutCellData(patch, dir);
      Base::ComputeBoundaryFluxes(boundary_fluxes, scratch, reference_states,
                                  cutcell_data, dir, dt, dx);

      View<Conservative> stabilized_fluxes =
          Subview(MakeView<BasicView<Conservative>>(
                      context.GetStabilizedFluxes(patch, dir), equation),
                  tilebox_faces);

      View<Conservative> shielded_left_fluxes =
          Subview(MakeView<BasicView<Conservative>>(
                      context.GetShieldedLeftFluxes(patch, dir), equation),
                  tilebox_faces);

      View<Conservative> shielded_right_fluxes =
          Subview(MakeView<BasicView<Conservative>>(
                      context.GetShieldedRightFluxes(patch, dir), equation),
                  tilebox_faces);

      View<Conservative> doubly_shielded_fluxes =
          Subview(MakeView<BasicView<Conservative>>(
                      context.GetDoublyShieldedFluxes(patch, dir), equation),
                  tilebox_faces);

      Base::ComputeCutCellFluxes(
          stabilized_fluxes, regular_fluxes, shielded_left_fluxes,
          shielded_right_fluxes, doubly_shielded_fluxes,
          AsConst(boundary_fluxes), scratch, cutcell_data, dir, dt, dx);
    }
  }

  Complete ref_;
};

} // namespace cutcell
} // namespace amrex
} // namespace fub

#endif
