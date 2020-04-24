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

#ifndef FUB_SAMRAI_FLUX_METHOD_HPP
#define FUB_SAMRAI_FLUX_METHOD_HPP

#include "fub/SAMRAI/IntegratorContext.hpp"
#include "fub/SAMRAI/ViewPatch.hpp"

#include <limits>
#include <numeric>

namespace fub::samrai {
/// This is a wrapper class which dispatches a given base method object and
/// dispatches SAMRAI typed patches.
///
/// The base method is expected to act on View objects of equation states.
template <typename Tag, typename BaseMethod>
class FluxMethodAdapter : private BaseMethod {
public:
  using Equation =
      std::decay_t<decltype(std::declval<const BaseMethod&>().GetEquation())>;
  using Conservative = ::fub::Conservative<Equation>;
  using Complete = ::fub::Complete<Equation>;

  FluxMethodAdapter(Tag, const BaseMethod& base) : BaseMethod(base) {}
  FluxMethodAdapter(Tag, BaseMethod&& base) : BaseMethod(std::move(base)) {}

  using BaseMethod::GetEquation;

  using BaseMethod::GetStencilWidth;

  /// Extracts the state variables patch data views including its ghost cells
  /// and compute a stable time step size in specified direction and refinement
  /// level.
  ///
  /// \pre This method expects the ghost cells to be filled with proper boundary
  /// data.
  Duration ComputeStableDt(IntegratorContext& context, int level,
                           Direction dir);

  void ComputeNumericFluxes(IntegratorContext& context, int level, Duration dt,
                            Direction dir);

  void ComputeNumericFluxes(span<SAMRAI::pdat::SideData<double>*> fluxes,
                            span<SAMRAI::pdat::CellData<double> const*> cells,
                            double dx, Duration dt, Direction dir);

  Duration ComputeStableDt(span<SAMRAI::pdat::CellData<double> const*> data,
                           double dx, Direction dir);
};

template <typename Tag, typename BaseMethod>
Duration
FluxMethodAdapter<Tag, BaseMethod>::ComputeStableDt(IntegratorContext& context,
                                             int level, Direction dir) {
  const SAMRAI::hier::PatchLevel& patch_level = context.GetPatchLevel(level);
  span<const int> data_ids = context.GetScratchIds();
  const SAMRAI::geom::CartesianGridGeometry& geom = context.GetGeometry(level);
  const int dir_value = static_cast<int>(dir);
  const double dx = geom.getDx()[dir_value];
  std::vector<const SAMRAI::pdat::CellData<double>*> patch_data(
      data_ids.size());
  return std::accumulate(patch_level.begin(), patch_level.end(),
                         Duration(std::numeric_limits<double>::max()),
                         [&](Duration min_dt, auto&& patch) {
                           GetPatchData(span{patch_data}, *patch, data_ids);
                           return std::min(
                               min_dt, ComputeStableDt(patch_data, dx, dir));
                         });
}

template <typename Tag, typename BaseMethod>
Duration FluxMethodAdapter<Tag, BaseMethod>::ComputeStableDt(
    span<SAMRAI::pdat::CellData<double> const*> data, double dx,
    Direction dir) {
  constexpr int Rank = Equation::Rank();
  auto states =
      MakeView<const Complete>(data, BaseMethod::GetEquation(),
                               AsIndexBox<Rank>(data[0]->getGhostBox()));
  return Duration(BaseMethod::ComputeStableDt(Tag(), states, dx, dir));
}

template <typename Tag, typename BaseMethod>
void FluxMethodAdapter<Tag, BaseMethod>::ComputeNumericFluxes(
    fub::samrai::IntegratorContext& context, int level, fub::Duration dt,
    fub::Direction dir) {
  const SAMRAI::geom::CartesianGridGeometry& geom = context.GetGeometry(level);
  const int dir_value = static_cast<int>(dir);
  const double dx = geom.getDx()[dir_value];
  const SAMRAI::hier::PatchLevel& patch_level = context.GetPatchLevel(level);
  span<const int> scratch_ids = context.GetScratchIds();
  std::vector<const SAMRAI::pdat::CellData<double>*> scratch_data(
      scratch_ids.size());
  span<const int> flux_ids = context.GetFluxIds();
  std::vector<SAMRAI::pdat::SideData<double>*> flux_data(flux_ids.size());
  for (std::shared_ptr<SAMRAI::hier::Patch> patch : patch_level) {
    GetPatchData(span{flux_data}, *patch, flux_ids);
    GetPatchData(span{scratch_data}, *patch, scratch_ids);
    ComputeNumericFluxes(flux_data, scratch_data, dx, dt, dir);
  }
}

template <typename Tag, typename BaseMethod>
void FluxMethodAdapter<Tag, BaseMethod>::ComputeNumericFluxes(
    span<SAMRAI::pdat::SideData<double>*> fluxes,
    span<SAMRAI::pdat::CellData<double> const*> cells, double dx, Duration dt,
    Direction dir) {
  constexpr int Rank = Equation::Rank();
  auto cell_box = AsIndexBox<Rank>(cells[0]->getGhostBox());
  const int gcw = GetStencilWidth();
  auto face_box = Shrink(cell_box, dir, {gcw, gcw});
  auto flux_view = MakeView<Conservative>(fluxes, GetEquation(), dir, face_box);
  auto scratch_view = MakeView<const Complete>(cells, GetEquation(), cell_box);
  BaseMethod::ComputeNumericFluxes(flux_view, scratch_view, dt, dx, dir);
}

} // namespace fub::samrai

#endif
