// Copyright (c) 2019 Maikel Nadolski
// Copyright (c) 2019 Patrick Denzler
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

#ifndef FUB_SAMRAI_TAGGING_GRADIENT_DETECTOR_HPP
#define FUB_SAMRAI_TAGGING_GRADIENT_DETECTOR_HPP

#include "fub/tagging/GradientDetector.hpp"

#include "fub/SAMRAI/GriddingAlgorithm.hpp"
#include <SAMRAI/hier/Box.h>

namespace fub::samrai {

template <typename Equation, typename... Projections> class GradientDetector {
public:
  GradientDetector(const Equation& equation,
                   const std::pair<Projections, double>&... projs);

  void TagCellsForRefinement(GriddingAlgorithm& gridding, int level, int tag_id,
                             Duration time_point);

private:
  fub::GradientDetector<Equation, Projections...> detector_;
};

template <typename Eq, typename... Ps>
GradientDetector(const Eq& eq, const std::pair<Ps, double>&... ps)
    ->GradientDetector<Eq, Ps...>;

template <typename Equation, typename... Projections>
GradientDetector<Equation, Projections...>::GradientDetector(
    const Equation& equation, const std::pair<Projections, double>&... projs)
    : detector_(::fub::GradientDetector(equation, projs...)) {}

template <typename Equation, typename... Projections>
void GradientDetector<Equation, Projections...>::TagCellsForRefinement(
    GriddingAlgorithm& gridding, int level, int tag_id,
    Duration /* time_point */) {
  const fub::samrai::PatchHierarchy& patch_hierarchy = gridding.GetPatchHierarchy();
  const std::vector<int>& data_ids =
      patch_hierarchy.GetDataDescription().data_ids;
  std::vector<SAMRAI::pdat::CellData<double>*> datas(data_ids.size());

  for (const std::shared_ptr<SAMRAI::hier::Patch>& patch :
       *patch_hierarchy.GetNative()->getPatchLevel(level)) {
    // Get States
    std::transform(data_ids.begin(), data_ids.end(), datas.begin(),
                   [&](int id) -> SAMRAI::pdat::CellData<double>* {
                     return static_cast<SAMRAI::pdat::CellData<double>*>(
                         patch->getPatchData(id).get());
                   });
    fub::BasicView basic_view = fub::samrai::MakeView<fub::Complete<Equation>>(datas, detector_.GetEquation());
    fub::View<const fub::Complete<Equation>> statesview = AsConst(Subview(basic_view, fub::Box<0>(basic_view)));

    // Get Tags
    SAMRAI::pdat::CellData<int>& tags =
        *static_cast<SAMRAI::pdat::CellData<int>*>( // NOLINT
            patch->getPatchData(tag_id).get());

    auto tagsview = MakePatchDataView<Equation::Rank()>(tags.getArrayData()).Subview(Box<0>(statesview));

    detector_.TagCellsForRefinement(tagsview, statesview);
  }
}
} // namespace fub::samrai

#endif // FUB_GRADIENT_DETECTOR_HPP
