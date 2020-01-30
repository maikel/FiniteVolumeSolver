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

#ifndef FUB_AMREX_FOR_EACH_FAB_HPP
#define FUB_AMREX_FOR_EACH_FAB_HPP

#include "fub/Execution.hpp"

#include <AMReX_MultiFab.H>

namespace fub::amrex {

template <typename Tag, typename F>
void ForEachFab(Tag, const ::amrex::FabArrayBase& fabarray, F function) {
  for (::amrex::MFIter mfi(fabarray); mfi.isValid(); ++mfi) {
    function(mfi);
  }
}

template <typename Tag, typename F>
void ForEachFab(Tag, const ::amrex::BoxArray& ba,
                const ::amrex::DistributionMapping& dm, F function) {
  for (::amrex::MFIter mfi(ba, dm); mfi.isValid(); ++mfi) {
    function(mfi);
  }
}

template <typename F>
void ForEachFab(execution::OpenMpTag, const ::amrex::FabArrayBase& fabarray,
                F function) {
#if defined(_OPENMP) && defined(AMREX_USE_OMP)
#pragma omp parallel
#endif
  for (::amrex::MFIter mfi(fabarray, true); mfi.isValid(); ++mfi) {
    function(mfi);
  }
}

template <typename F>
void ForEachFab(execution::OpenMpTag, const ::amrex::BoxArray& ba,
                const ::amrex::DistributionMapping& dm, F function) {
#if defined(_OPENMP) && defined(AMREX_USE_OMP)
#pragma omp parallel
#endif
  for (::amrex::MFIter mfi(ba, dm, true); mfi.isValid(); ++mfi) {
    function(mfi);
  }
}

template <typename F>
void ForEachFab(execution::OpenMpSimdTag, const ::amrex::FabArrayBase& fabarray,
                F function) {
  ForEachFab(execution::openmp, fabarray, std::move(function));
}

template <typename F>
void ForEachFab(const ::amrex::FabArrayBase& fabarray, F function) {
  ForEachFab(execution::seq, fabarray, std::move(function));
}

template <typename F>
void ForEachFab(execution::OpenMpSimdTag, const ::amrex::BoxArray& ba,
                const ::amrex::DistributionMapping& dm, F function) {
  ForEachFab(execution::openmp, ba, dm, std::move(function));
}

template <typename F>
void ForEachFab(const ::amrex::BoxArray& ba,
                const ::amrex::DistributionMapping& dm, F function) {
  ForEachFab(execution::seq, ba, dm, std::move(function));
}

} // namespace fub::amrex

#endif
