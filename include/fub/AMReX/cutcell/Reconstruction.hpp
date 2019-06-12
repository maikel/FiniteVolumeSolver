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

#ifndef FUB_AMREX_CUTCELL_RECONSTRUCTION_HPP
#define FUB_AMREX_CUTCELL_RECONSTRUCTION_HPP

#include "fub/CompleteFromCons.hpp"
#include "fub/CutCellData.hpp"
#include "fub/ForEach.hpp"
#include "fub/State.hpp"

namespace fub::amrex::cutcell {

template <typename Equation> struct ScalarReconstructStates {
  using Complete = ::fub::Complete<Equation>;
  using Conservative = ::fub::Conservative<Equation>;
  static constexpr int Rank = Equation::Rank();

  explicit ScalarReconstructStates(const Equation& eq) : kernel_{eq} {}

  void CompleteFromCons(IntegratorContext& context, int level,
                        [[maybe_unused]] Duration dt, Direction dir) {
    ::amrex::MultiFab dest = context.GetData(level);
    const ::amrex::MultiFab src = context.GetScratch(level, dir);

    ForEachMFIter(execution::openmp, dest, [&](::amrex::MFIter& mfi) {
      const ::amrex::Box box = mfi.tilebox();
      auto state = MakeView<Complete>(dest[mfi], *equation_, box);
      auto scratch = MakeView<const Consevative>(src[mfi], *equation_, box);
      ::amrex::FabType type = context.GetCutCellPatchType(mfi);
      if (type == ::amrex::FabType::regular) {
        kernel_->CompleteFromCons(state, scratch);
      } else if (type == ::amrex::FabType::singlevalued) {
        ::amrex::TagBox flags = context.GetFlags(mfi);
        kernel_->CompleteFromCons(state, scratch, flags);
      }
    });
  }

  struct Kernel {
    void CompleteFromCons(const View<Complete>& dest,
                          const View<const Conservative>& src) {
      ForEachIndex(Box<0>(dest), [&](auto... is) {
        Load(cons_, src, {is...});
        CompleteFromCons(eq, complete_, cons_);
        Store(complete_, dest, {is...});
      });
    }

    void CompleteFromCons(const View<Complete>& dest,
                          const View<const Conservative>& src,
                          const ::amrex::TagBox& cutcell_flag) {
      ForEachIndex(Box<0>(dest), [&](auto... is) {
        if (cutcell_flag({int(is)...}) != 2) {
          Load(cons_, src, {is...});
          CompleteFromCons(eq, complete_, cons_);
          Store(complete_, dest, {is...});
        }
      });
    }

    Equation equation_;
    Complete complete_{equation_};
    Conservative cons_{equation_};
  };
  OmpLocal<Kernel> kernel_;
};

} // namespace fub::amrex::cutcell

#endif
