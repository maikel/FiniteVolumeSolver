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

#ifndef FUB_AMREX_TAGGING_GRADIENT_DETECTOR_HPP
#define FUB_AMREX_TAGGING_GRADIENT_DETECTOR_HPP

#include "fub/AMReX/Tagging.hpp"
#include "fub/tagging/GradientDetector.hpp"

namespace fub::amrex {

template <typename Equation, typename... Projections>
class GradientDetector
    : public ::fub::GradientDetector<Equation, Projections...> {
public:
  GradientDetector(const Equation& equation,
                   const std::pair<Projections, double>&... projs);

  void TagCellsForRefinement(::amrex::TagBoxArray& tags, Duration, int level,
                             GriddingAlgorithm& gridding);
};

template <typename Equation, typename... Projections>
GradientDetector<Equation, Projections...>::GradientDetector(
    const Equation& equation, const std::pair<Projections, double>&... projs)
    : ::fub::GradientDetector<Equation, Projections...>{equation, projs...} {}

template <typename Equation, typename... Projections>
void GradientDetector<Equation, Projections...>::TagCellsForRefinement(
    ::amrex::TagBoxArray& tags, Duration, int level,
    GriddingAlgorithm& gridding) {
  const PatchHierarchy& hierarchy = gridding.GetPatchHierarchy();
  ::amrex::BoxArray ba = hierarchy.GetPatchLevel(level).box_array;
  ::amrex::DistributionMapping dm =
      hierarchy.GetPatchLevel(level).distribution_mapping;
  const ::amrex::MultiFab& data = hierarchy.GetPatchLevel(level).data;
  ::amrex::MultiFab scratch(ba, dm, data.nComp(), 1, ::amrex::MFInfo(),
                            data.Factory());
  gridding.FillMultiFabFromLevel(scratch, level);
#if defined(_OPENMP) && defined(AMREX_USE_OMP)
#pragma omp parallel
#endif
  for (::amrex::MFIter mfi(scratch, true); mfi.isValid(); ++mfi) {
    ::amrex::Box grown = ::amrex::grow(mfi.tilebox(), 1);
    View<const Complete<Equation>> states = MakeView<const Complete<Equation>>(
        scratch[mfi],
        ::fub::GradientDetector<Equation, Projections...>::equation_, grown);
    PatchDataView<char, AMREX_SPACEDIM, layout_stride> tagsview =
        MakePatchDataView(tags[mfi], 0)
            .Subview(AsIndexBox<AMREX_SPACEDIM>(mfi.tilebox()));
    ::fub::GradientDetector<Equation, Projections...>::TagCellsForRefinement(
        tagsview, states);
  }
}

} // namespace fub::amrex

#endif
