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

#ifndef FUB_AMREX_RECONSTRUCTION_HPP
#define FUB_AMREX_RECONSTRUCTION_HPP

#include "fub/CompleteFromCons.hpp"
#include "fub/Execution.hpp"
#include "fub/State.hpp"

#include "fub/AMReX/ForEachFab.hpp"
#include "fub/AMReX/IntegratorContext.hpp"

#include <AMReX_MultiFab.H>

namespace fub::amrex {

template <typename Tag, typename Equation> class Reconstruction {
public:
  Reconstruction(const Equation& eq) : Reconstruction(Tag(), eq) {}

  Reconstruction(Tag, const Equation& eq) : rec_{} {
    rec_ = Local<Tag, CompleteFromConsFn<Equation>>{
        CompleteFromConsFn<Equation>{eq}};
  }

  void CompleteFromCons(::amrex::MultiFab& dest, const ::amrex::MultiFab& src) {
    ForEachFab(Tag(), dest, [&](const ::amrex::MFIter& mfi) {
      Equation& eq = rec_->equation_;
      View<Complete<Equation>> complete =
          MakeView<Complete<Equation>>(dest[mfi], eq, mfi.growntilebox());
      View<const Conservative<Equation>> conservative =
          MakeView<const Conservative<Equation>>(src[mfi], eq,
                                                 mfi.growntilebox());
      rec_->CompleteFromCons(Tag(), complete, conservative);
    });
  }

  void CompleteFromCons(IntegratorContext& context, int level, Duration) {
    ::amrex::MultiFab& dest = context.GetScratch(level);
    const ::amrex::MultiFab& src = context.GetScratch(level);
    CompleteFromCons(dest, src);
  }

private:
  Local<Tag, CompleteFromConsFn<Equation>> rec_;
};

template <typename Equation>
Reconstruction(const Equation& eq)
    -> Reconstruction<execution::OpenMpSimdTag, Equation>;

template <typename Tag, typename Equation>
Reconstruction(Tag, const Equation& eq)->Reconstruction<Tag, Equation>;

struct NoReconstruction {
  void CompleteFromCons([[maybe_unused]] IntegratorContext& context,
                        [[maybe_unused]] int level,
                        [[maybe_unused]] Duration time_step_size) {}
};

} // namespace fub::amrex

#endif
