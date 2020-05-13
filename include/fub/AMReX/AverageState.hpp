// Copyright (c) 2020 Maikel Nadolski
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

#ifndef FUB_AMREX_AVERAGE_STATE_HPP
#define FUB_AMREX_AVERAGE_STATE_HPP

#include <AMReX.H>

#include "fub/AMReX/ForEachFab.hpp"
#include "fub/AMReX/ViewFArrayBox.hpp"

#include "fub/ForEach.hpp"

#include <mpi.h>

namespace fub::amrex {

template <typename Equation>
void AverageState(Complete<Equation>& state, const ::amrex::MultiFab& mf,
                  const ::amrex::Geometry&, const ::amrex::Box& box) {
  const auto ncomp = static_cast<std::size_t>(mf.nComp());
  std::vector<double> state_buffer(ncomp);
  const std::ptrdiff_t num_cells = box.numPts();
  ForEachFab(execution::seq, mf, [&](const ::amrex::MFIter& mfi) {
    IndexBox<AMREX_SPACEDIM> section =
        AsIndexBox<AMREX_SPACEDIM>(box & mfi.tilebox());
    const ::amrex::FArrayBox& data = mf[mfi];
    for (std::size_t comp = 0; comp < ncomp; ++comp) {
      ForEachIndex(section, [&](auto... is) {
        ::amrex::IntVect index{static_cast<int>(is)...};
        state_buffer[comp] += data(index, static_cast<int>(comp)) /
                              static_cast<double>(num_cells);
      });
    }
  });
  std::vector<double> global_state_buffer(ncomp);
  MPI_Allreduce(state_buffer.data(), global_state_buffer.data(),
                static_cast<int>(ncomp), MPI_DOUBLE, MPI_SUM,
                ::amrex::ParallelDescriptor::Communicator());
  CopyFromBuffer(state, global_state_buffer);
}

template <typename Equation>
void AverageState(Conservative<Equation>& state, const ::amrex::MultiFab& mf,
                  const ::amrex::Geometry&, const ::amrex::Box& box) {
  const auto ncomp = static_cast<std::size_t>(mf.nComp());
  std::vector<double> state_buffer(ncomp);
  const std::ptrdiff_t num_cells = box.numPts();
  ForEachFab(mf, [&](const ::amrex::MFIter& mfi) {
    IndexBox<AMREX_SPACEDIM> section =
        AsIndexBox<AMREX_SPACEDIM>(box & mfi.tilebox());
    const ::amrex::FArrayBox& data = mf[mfi];
    for (std::size_t comp = 0; comp < ncomp; ++comp) {
      ForEachIndex(section, [&](auto... is) {
        ::amrex::IntVect index{static_cast<int>(is)...};
        state_buffer[comp] += data(index, static_cast<int>(comp)) /
                              static_cast<double>(num_cells);
      });
    }
  });
  std::vector<double> global_state_buffer(ncomp);
  MPI_Allreduce(state_buffer.data(), global_state_buffer.data(),
                static_cast<int>(ncomp), MPI_DOUBLE, MPI_SUM,
                ::amrex::ParallelDescriptor::Communicator());
  CopyFromBuffer(state, global_state_buffer);
}

} // namespace fub::amrex

#endif