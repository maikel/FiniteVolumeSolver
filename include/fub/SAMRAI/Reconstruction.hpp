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

#ifndef FUB_SAMRAI_RECONSTRUCTION_HPP
#define FUB_SAMRAI_RECONSTRUCTION_HPP

#include "fub/CompleteFromCons.hpp"
#include "fub/SAMRAI/IntegratorContext.hpp"
#include "fub/SAMRAI/ViewPatch.hpp"

namespace fub::samrai {

template <typename Tag, typename Equation_> class Reconstruction {
public:
  using Equation = Equation_;
  using ExecutionTag = Tag;

  static constexpr int Rank = Equation::Rank();

  Reconstruction(Tag, const Equation& eq) : rec_{eq} {}

  void CompleteFromCons(SAMRAI::hier::Patch& dest_patch,
                        span<const int> dest_ids,
                        const SAMRAI::hier::Patch& src_patch,
                        span<const int> src_ids);

  void CompleteFromCons(IntegratorContext& context, int level,
                        [[maybe_unused]] Duration time_step_size);

private:
  CompleteFromConsFn<Equation> rec_;
};

template <typename Tag, typename Equation>
void Reconstruction<Tag, Equation>::CompleteFromCons(
    SAMRAI::hier::Patch& dest_patch, span<const int> dest_ids,
    const SAMRAI::hier::Patch& src_patch, span<const int> src_ids) {
  Equation& equation = rec_.equation_;
  std::vector<SAMRAI::pdat::CellData<double>*> dest_data(dest_ids.size());
  GetPatchData(span{dest_data}, dest_patch, dest_ids);
  std::vector<SAMRAI::pdat::CellData<double>*> src_data(src_ids.size());
  GetPatchData(span{src_data}, src_patch, src_ids);
  IndexBox<Rank> box = AsIndexBox<Rank>(dest_data[0]->getGhostBox());
  View<Complete<Equation>> complete = MakeView<Complete<Equation>>(
      span{dest_data}, equation, box);
  View<const Conservative<Equation>> conservative =
      MakeView<const Conservative<Equation>>(src_data, equation, box);
  rec_.CompleteFromCons(Tag(), complete, conservative);
}

template <typename Tag, typename Equation>
void Reconstruction<Tag, Equation>::CompleteFromCons(
    IntegratorContext& context, int level,
    [[maybe_unused]] Duration time_step_size) {
  SAMRAI::hier::PatchLevel& scratch = context.GetScratch(level);
  const int n_cons = static_cast<int>(context.GetFluxIds().size());
  for (const std::shared_ptr<SAMRAI::hier::Patch>& patch : scratch) {
    CompleteFromCons(*patch, context.GetScratchIds(), *patch,
                     context.GetScratchIds().first(n_cons));
  }
}

} // namespace fub::samrai

#endif