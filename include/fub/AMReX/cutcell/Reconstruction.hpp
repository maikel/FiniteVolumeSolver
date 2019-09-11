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

  explicit Reconstruction(const Tag& tag, const Equation& eq);

  void CompleteFromCons(IntegratorContext& context, int level,
                        [[maybe_unused]] Duration dt);

  Local<Tag, detail::ReconstructionKernel<Equation, IsSimd>> kernel_;
};

////////////////////////////////////////////////////////////////////////////////
//                                                              Implementation

namespace detail {
// Scalar Case
template <typename Equation, bool IsSimd> struct ReconstructionKernel {
  Equation equation_;
  Complete<Equation> complete_{equation_};
  Conservative<Equation> cons_{equation_};

  static constexpr int Rank = Equation::Rank();

  void CompleteFromCons(const View<Complete<Equation>>& dest,
                        const View<const Conservative<Equation>>& src) {
    ForEachIndex(Box<0>(dest), [&](auto... is) {
      Load(cons_, src, {is...});
      ::fub::CompleteFromCons(equation_, complete_, cons_);
      Store(dest, complete_, {is...});
    });
  }

  void CompleteFromCons(
      const View<Complete<Equation>>& dest,
      const View<const Conservative<Equation>>& src,
      const PatchDataView<const double, Rank, layout_stride>& volume) {
    ForEachIndex(Box<0>(dest), [&](auto... is) {
      const double alpha = volume(is...);
      if (alpha > 0.0) {
        Load(cons_, src, {is...});
        ::fub::CompleteFromCons(equation_, complete_, cons_);
        Store(dest, complete_, {is...});
      }
    });
  }
};

// SIMD Case
template <typename Equation> struct ReconstructionKernel<Equation, true> {
  Equation equation_;
  CompleteArray<Equation> complete_{equation_};
  ConservativeArray<Equation> cons_{equation_};

  static constexpr int Rank = Equation::Rank();

  void CompleteFromCons(const View<Complete<Equation>>& dest,
                        const View<const Conservative<Equation>>& src) {
    ForEachRow(std::tuple{dest, src},
               [&](const Row<Complete<Equation>>& dest,
                   const Row<const Conservative<Equation>>& src) {
                 ViewPointer<const Conservative<Equation>> first = Begin(src);
                 ViewPointer<const Conservative<Equation>> last = End(src);
                 ViewPointer<Complete<Equation>> out = Begin(dest);
                 constexpr std::ptrdiff_t size =
                     static_cast<std::ptrdiff_t>(kDefaultChunkSize);
                 std::ptrdiff_t n = get<0>(last) - get<0>(first);
                 while (n >= size) {
                   Load(cons_, first);
                   ::fub::CompleteFromCons(equation_, complete_, cons_);
                   Store(out, complete_);
                   Advance(first, size);
                   Advance(out, size);
                   n = get<0>(last) - get<0>(first);
                 }
                 LoadN(cons_, first, static_cast<int>(n));
                 ::fub::CompleteFromCons(equation_, complete_, cons_);
                 StoreN(out, complete_, static_cast<int>(n));
               });
  }

  void CompleteFromCons(
      const View<Complete<Equation>>& dest,
      const View<const Conservative<Equation>>& src,
      const PatchDataView<const double, Rank, layout_stride>& volume) {
    ForEachRow(std::tuple{dest, src, volume},
               [&](const Row<Complete<Equation>>& dest,
                   const Row<const Conservative<Equation>>& src,
                   span<const double> volume) {
                 ViewPointer<const Conservative<Equation>> first = Begin(src);
                 ViewPointer<const Conservative<Equation>> last = End(src);
                 ViewPointer<Complete<Equation>> out = Begin(dest);
                 Array1d alpha = Array1d::Zero();
                 constexpr std::ptrdiff_t size =
                     static_cast<std::ptrdiff_t>(kDefaultChunkSize);
                 std::ptrdiff_t n = get<0>(last) - get<0>(first);
                 while (n >= size) {
                   Load(cons_, first);
                   alpha = Array1d::Map(volume.data());
                   MaskArray mask = alpha > 0.0;
                   ::fub::CompleteFromCons(equation_, complete_, cons_, mask);
                   Store(out, complete_);
                   Advance(first, size);
                   Advance(out, size);
                   volume = volume.subspan(size);
                   n = get<0>(last) - get<0>(first);
                 }
                 LoadN(cons_, first, static_cast<int>(n));
                 std::copy_n(volume.data(), n, alpha.data());
                 std::fill_n(alpha.data() + n, size - n, 0.0);
                 MaskArray mask = alpha > 0.0;
                 ::fub::CompleteFromCons(equation_, complete_, cons_, mask);
                 StoreN(out, complete_, static_cast<int>(n));
               });
  }
};
} // namespace detail

template <typename Tag, typename Equation>
Reconstruction<Tag, Equation>::Reconstruction(const Tag&, const Equation& eq)
    : kernel_() {

  kernel_ = Local<Tag, detail::ReconstructionKernel<Equation, IsSimd>>{detail::ReconstructionKernel<Equation, IsSimd>{eq}};
    }

template <typename Tag, typename Equation>
void Reconstruction<Tag, Equation>::CompleteFromCons(IntegratorContext& context,
                                                     int level, Duration) {
  ::amrex::MultiFab& dest = context.GetData(level);
  const ::amrex::MultiFab& volumes =
      context.GetEmbeddedBoundary(level).getVolFrac();
  const ::amrex::MultiFab& src = context.GetScratch(level);

  ForEachFab(Tag(), dest, [&](const ::amrex::MFIter& mfi) {
    const ::amrex::Box box = mfi.tilebox();
    auto state =
        MakeView<Complete<Equation>>(dest[mfi], kernel_->equation_, box);
    auto scratch = MakeView<const Conservative<Equation>>(
        src[mfi], kernel_->equation_, box);
    ::amrex::FabType type =
        context.GetEmbeddedBoundary(level).getMultiEBCellFlagFab()[mfi].getType(
            box);
    if (type == ::amrex::FabType::regular) {
      kernel_->CompleteFromCons(state, scratch);
    } else if (type == ::amrex::FabType::singlevalued) {
      auto volume = MakePatchDataView(volumes[mfi], 0, box);
      kernel_->CompleteFromCons(state, scratch, volume);
    }
  });
}

} // namespace fub::amrex::cutcell

#endif
