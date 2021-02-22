// Copyright (c) 2020 Maikel Nadolski
// Copyright (c) 2020 Christian Zenker
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

#ifndef FUB_AMREX_MULTIFAB_UTILITIES_HPP
#define FUB_AMREX_MULTIFAB_UTILITIES_HPP

#include "fub/Direction.hpp"
#include <AMReX_MultiFab.H>

namespace fub::amrex {

void Realloc(::amrex::MultiFab& mf, const ::amrex::BoxArray& ba,
             const ::amrex::DistributionMapping& dm, int ncomp,
             const ::amrex::IntVect& ngrow,
             const ::amrex::MFInfo& mf_info = ::amrex::MFInfo(),
             const ::amrex::FabFactory<::amrex::FArrayBox>& factory =
                 ::amrex::FArrayBoxFactory());

void ReallocLike(::amrex::MultiFab& mf, const ::amrex::MultiFab& other);

::amrex::Box GetInnerBox(const ::amrex::Box& box, int side, Direction dir,
                         int width);

double GetMeanValueInBox(const ::amrex::MultiFab& data, const ::amrex::Box& box,
                         int component);

template <typename GriddingAlgorithm>
::amrex::Box GetFineBoxAtLevel(const ::amrex::Box& coarse_box, const GriddingAlgorithm& grid, int level)
{
  ::amrex::Box box = coarse_box;
  for (int ilvl = 0; ilvl < level; ++ilvl) {
      box.refine(grid.GetPatchHierarchy().GetRatioToCoarserLevel(ilvl + 1));
  }
  return box;
}

::amrex::MultiFab CopyMultiFab(const ::amrex::MultiFab& other);
} // namespace fub::amrex

#endif