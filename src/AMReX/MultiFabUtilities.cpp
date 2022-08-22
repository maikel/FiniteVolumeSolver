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

#include "fub/AMReX/MultiFabUtilities.hpp"

#include "fub/AMReX/ForEachFab.hpp"
#include "fub/Execution.hpp"
#include "fub/ext/omp.hpp"

namespace fub::amrex {

void Realloc(::amrex::MultiFab& mf, const ::amrex::BoxArray& ba,
             const ::amrex::DistributionMapping& dm, int ncomp,
             const ::amrex::IntVect& ngrow, const ::amrex::MFInfo& mf_info,
             const ::amrex::FabFactory<::amrex::FArrayBox>& factory) {
  if (mf.boxArray() != ba || mf.DistributionMap() != dm ||
      mf.nComp() != ncomp || mf.nGrowVect() != ngrow) {
    mf.define(ba, dm, ncomp, ngrow, mf_info, factory);
  }
}

void ReallocLike(::amrex::MultiFab& mf, const ::amrex::MultiFab& other) {
  if (other.ok()) {
    Realloc(mf, other.boxArray(), other.DistributionMap(), other.nComp(),
            other.nGrowVect(), ::amrex::MFInfo().SetArena(other.arena()),
            other.Factory());
  }
}

::amrex::MultiFab CopyMultiFab(const ::amrex::MultiFab& other) {
  ::amrex::MultiFab mf(other.boxArray(), other.DistributionMap(), other.nComp(),
                       other.nGrowVect(),
                       ::amrex::MFInfo().SetArena(other.arena()),
                       other.Factory());
  ::amrex::MultiFab::Copy(mf, other, 0, 0, other.nComp(), other.nGrowVect());
  return mf;
}

/// \brief get the inner box
/// \param width number of cells
::amrex::Box GetInnerBox(const ::amrex::Box& box, int side, Direction dir,
                         int width) {
  const int dir_v = static_cast<int>(dir);
  const ::amrex::IntVect lower = box.smallEnd();
  const ::amrex::IntVect upper = box.bigEnd();
  if (side == 0) {
    const ::amrex::IntVect lower_new = lower;
    ::amrex::IntVect upper_new = upper;
    upper_new[dir_v] = lower_new[dir_v] + width - 1;
    return ::amrex::Box(lower_new, upper_new);
  }
  const ::amrex::IntVect upper_new = upper;
  ::amrex::IntVect lower_new = lower;
  lower_new[dir_v] = upper_new[dir_v] - width + 1;
  return ::amrex::Box(lower_new, upper_new);
}

double GetMeanValueInBox(const ::amrex::MultiFab& data, const ::amrex::Box& box,
                         int component) {
  const auto volume = box.volume();
  OmpLocal<double> omp_local_value{double{0.0}};
  ForEachFab(execution::openmp, data, [&](const ::amrex::MFIter& mfi) {
    const ::amrex::FArrayBox& fab = data[mfi];
    const ::amrex::Box subbox = box & fab.box();
    *omp_local_value += fab.sum(subbox, component) / volume;
  });
  const double local_value = omp_local_value.Accumulate();
  double average_value = 0.0;
  ::MPI_Allreduce(&local_value, &average_value, 1, MPI_DOUBLE, MPI_SUM,
                  ::amrex::ParallelContext::CommunicatorAll());
  return average_value;
}

} // namespace fub::amrex