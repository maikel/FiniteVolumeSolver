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

#include "fub/SAMRAI/HyperbolicTimeIntegrator.hpp"
#include "fub/HyperbolicPatchIntegrator.hpp"
#include <range/v3/view/zip.hpp>

namespace fub::samrai {

void HyperbolicTimeIntegrator::UpdateConservatively(IntegratorContext& context, int level,
                                          Duration dt, Direction dir) {
  SAMRAI::hier::PatchLevel& patch_level = context.GetPatchLevel(level);
  span<const int> flux_ids = context.GetFluxIds();
  span<const int> scratch_ids = context.GetScratchIds();
  const int n_conservative_variables =
      context.GetPatchHierarchy().GetDataDescription().n_cons_variables;
  const int dir_v = static_cast<int>(dir);
  const double dx = context.GetGeometry(level).getDx()[dir_v];
  for (std::shared_ptr<SAMRAI::hier::Patch> patch : patch_level) {
    constexpr int Rank = SAMRAI_MAXIMUM_DIMENSION;
    for (auto [flux, scratch] : ranges::views::zip(flux_ids, scratch_ids)) {
      SAMRAI::pdat::CellData<double>* state_data =
          dynamic_cast<SAMRAI::pdat::CellData<double>*>(
              patch->getPatchData(scratch).get());
      const SAMRAI::pdat::SideData<double>* flux_data =
          dynamic_cast<SAMRAI::pdat::SideData<double>*>(
              patch->getPatchData(flux).get());
      FUB_ASSERT(state_data && flux_data);
      IndexBox<Rank + 1> cells = Embed<Rank + 1>(
          Shrink(AsIndexBox<Rank>(state_data->getGhostBox()), dir, {1, 1}),
          {0, n_conservative_variables});
      IndexBox<Rank + 1> faces = Grow(cells, dir, {0, 1});
      auto next = MakePatchDataView<Rank + 1>(state_data->getArrayData())
                      .Subview(cells);
      auto prev = next;
      auto fluxes = MakePatchDataView<Rank + 1>(flux_data->getArrayData(dir_v))
                        .Subview(faces);
      using simd = execution::SimdTag;
      HyperbolicPatchIntegrator<simd>::UpdateConservatively(next, prev, fluxes,
                                                            dt, dx, dir);
    }
  }
}

} // namespace fub::samrai