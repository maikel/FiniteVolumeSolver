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

#include "fub/grid/AMReX/CutCellLevelIntegrator.hpp"

namespace fub {
namespace amrex {

void CutCellLevelIntegrator::CutCellUpdateConservatively(
    const ::amrex::MFIter& mfi, Direction dir, Duration dt) {
  const double dx = GetDx(dir);
  const double lambda = dt.count() / dx;
  const ::amrex::FArrayBox& flux_b = GetBoundaryFlux(mfi, dir);
  const ::amrex::FArrayBox& flux = GetFluxes(mfi, dir);
  const ::amrex::FArrayBox& states = GetScratch(mfi, dir);
  const ::amrex::FArrayBox& beta = GetFaceFractions(mfi, dir);
  for (BoxIterator bit(states.box()); bit.ok(); ++bit) {
    states(i) += lambda * (beta(i) * flux(i) - beta(i + 1) * flux(i + 1) +
                           (beta(i + 1) - beta(i)) * flux_bdry(i));
  }
}

void CutCellLevelIntegrator::UpdateConservatively(int level, Direction dir,
                                                  Duration dt) {
  const auto& flags = GetCutCellFlags(level);
  for (::amrex::MFIter mfi(states); mfi.isValid(); ++mfi) {
    const ::amrex::Box& box = mfi.tilebox();
    ::amrex::FabType type = flags[mfi].getType(box);
    if (type == ::amrex::FabType::regular) {
      RegularUpdateConservatively(mfi, dir, dt);
    } else if (type == ::amrex::FabType::singlevalued) {
      CutCellUpdateConservatively(mfi, dir, dt);
    }
  }
}

} // namespace amrex
} // namespace fub