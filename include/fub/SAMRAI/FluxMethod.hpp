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

namespace fub::samrai {
/// This is a wrapper class which dispatches a given base method object and
/// dispatches SAMRAI typed patches.
///
/// The base method is expected to act on View objects of equation states.
template <typename Tag, typename BaseMethod>
class FluxMethod : private BaseMethod {
public:
  using Equation = std::decay_t<decltype(std::declval<const BaseMethod&>().GetEquation())>;

  FluxMethod(Tag, const BaseMethod& base) : BaseMethod(base) {}
  FluxMethod(Tag, BaseMethod&& base) : BaseMethod(std::move(base))) {}

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

  void ComputeNumericFluxes(span<SAMRAI::pdat::SideData*> fluxes,
                            span<SAMRAI::pdat::CellData const*> cells,
                            const SAMRAI::geom::CartesianGridGeometry& geom,
                            Duration dt, Direction dir);

  Duration ComputeStableDt(span<SAMRAI::pdat::CellData const*> data,
                           const SAMRAI::geom::CartesianGridGeometry& geom,
                           Direction dir);
};

template <typename Tag, typename BaseMethod>
Duration
FluxMethod<Tag, BaseMethod>::ComputeStableDt(IntegratorContext& context,
                                             int level, Direction dir) {
  const SAMRAI::hier::PatchLevel& scratch = context.GetScratch(level);
  span<const int> data_ids = context.GetScratchIds();
  const SAMRAI::geom::CartesianGridGeometry& geom = context.GetGeometry(level);
  const int dir_value = static_cast<int>(dir);
  const double dx = geom.getDx()[dir_value];
  std::vector<const SAMRAI::pdat::CellData*> patch_data(data_ids.size());
  return std::accumulate(scratch.begin(), scratch.end(), Duration(std::numeric_limits<double>::max()), [&](Duration min_dt, auto&& patch) {
    GetPatchData(patch_data, *patch, data_ids);
    return std::min(min_dt, ComputeStableDt(patch_data, dx, dir));
  });
}

template <typename Tag, typename BaseMethod>
Duration FluxMethod<Tag, BaseMethod>::ComputeStableDt(
    span<SAMRAI::pdat::CellData const*> data,
    double dx, Direction dir) {
  auto states = MakeView<const Complete<Equation>>(data, equation, data[0].getGhostBox());
  return BaseMethod::ComputeStableDt(Tag(), states, dx, dir);
}

} // namespace fub::samrai

#endif