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
#include "fub/Execution.hpp"
#include "fub/ForEach.hpp"
#include "fub/State.hpp"

namespace fub::amrex::cutcell {
namespace detail {
template <typename Equation, bool IsSimd = true> struct ReconstructionKernel;
}

template <typename Tag, typename Equation> struct Reconstruction {
  static constexpr int Rank = Equation::Rank();
  static constexpr bool IsSimd = std::is_base_of<execution::SimdTag, Tag>();

  explicit Reconstruction(const Equation& eq);

  void CompleteFromCons(IntegratorContext& context, int level,
                        [[maybe_unused]] Duration dt, Direction dir);

  Local<Tag, detail::ReconstructionKernel<Equation, IsSimd>> kernel_;
};

////////////////////////////////////////////////////////////////////////////////
//                                                              Implementation

namespace detail {
// Scalar Case
template <typename Equation, bool IsSimd> struct ReconstructionKernel {
  Equation equation_;
  Complete<Equation> complete_(equation_);
  Conservative<Equation> cons_(equation_);

  void CompleteFromCons(const View<Complete>& dest,
                        const View<const Conservative>& src) {
    ForEachIndex(Box<0>(dest), [&](auto... is) {
      Load(cons_, src, {is...});
      ::fub::CompleteFromCons(eq, complete_, cons_);
      Store(complete_, dest, {is...});
    });
  }

  void
  CompleteFromCons(const View<Complete>& dest,
                   const View<const Conservative>& src,
                   const PatchDataView<double, Rank, layout_stride>& volume) {
    ForEachIndex(Box<0>(dest), [&](auto... is) {
      if (volume(is...)) {
        Load(cons_, src, {is...});
        ::fub::CompleteFromCons(eq, complete_, cons_);
        Store(complete_, dest, {is...});
      }
    });
  }
};

// SIMD Case
template <typename Equation> struct ReconstructionKernel<Equation, true> {
  Equation equation_;
  CompleteArray<Equation> complete_(equation_);
  ConservativeArray<Equation> cons_(equation_);

  void CompleteFromCons(const View<Complete>& dest,
                        const View<const Conservative>& src) {
    ForEachRow(std::tuple{dest, src}, [](const Row<Complete>& dest,
                                         const Row<const Conservative>& src) {
      ViewPointer<const Conservative> first = Begin(src);
      ViewPointer<const Conservative> last = End(src);
      ViewPointer<Complete> out = Begin(dest);
      constexpr std::ptrdiff_t size =
          std::ptrdiff_t(CompleteArray<Equation>::Size());
      std::ptrdiff_t n = get<0>(last) - get<0>(first);
      while (n >= size) {
        Load(cons_, first);
        ::fub::CompleteFromCons(eq, complete_, cons_);
        Store(out, complete_);
        Advance(first, size);
        Advance(out, size);
        n = get<0>(last) - get<0>(first);
      }
      LoadN(cons_, first, N);
      ::fub::CompleteFromCons(eq, complete_, cons_);
      StoreN(out, complete_, N);
    });
  }

  void
  CompleteFromCons(const View<Complete>& dest,
                   const View<const Conservative>& src,
                   const PatchDataView<double, Rank, layout_stride>& volume) {
    ForEachIndex(Box<0>(dest), [&](auto... is) {
      if (volume(is...)) {
        Load(cons_, src, {is...});
        ::fub::CompleteFromCons(eq, complete_, cons_);
        Store(complete_, dest, {is...});
      }
    });
  }
};
} // namespace detail

template <typename Tag, typename Equation>
Reconstruction<Tag, Equation>::Reconstruction(const Equation& eq)
    : kernel_(ReconstructionKernel<Equation, IsSimd>{eq}) {}

template <typename Tag, typename Equation>
void Reconstruction<Tag, Equation>::CompleteFromCons(IntegratorContext& context,
                                                     int level, Duration dt,
                                                     Direction dir) {
  ::amrex::MultiFab dest = context.GetData(level);
  ::amrex::MultiFab volumes = context.GetEmbeddedBoundary(level).getVolFrac();
  const ::amrex::MultiFab src = context.GetScratch(level, dir);

  ForEachMFIter(Tag(), dest, [&](::amrex::MFIter& mfi) {
    const ::amrex::Box box = mfi.tilebox();
    auto state = MakeView<Complete>(dest[mfi], *equation_, box);
    auto scratch = MakeView<const Consevative>(src[mfi], *equation_, box);
    auto volume = MakePatchDataView(volumes[mfi], 0);
    ::amrex::FabType type = context.GetCutCellPatchType(mfi);
    if (type == ::amrex::FabType::regular) {
      kernel_->CompleteFromCons(state, scratch);
    } else if (type == ::amrex::FabType::singlevalued) {
      ::amrex::TagBox flags = context.GetVolumes()[mfi];
      kernel_->CompleteFromCons(state, scratch, volume);
    }
  });
}

} // namespace fub::amrex::cutcell

#endif
