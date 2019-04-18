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

#ifndef FUB_GRID_AMREX_BOUNDARY_CONDITION_HPP
#define FUB_GRID_AMREX_BOUNDARY_CONDITION_HPP

#include "fub/Direction.hpp"
#include "fub/Duration.hpp"
#include "fub/PatchDataView.hpp"
#include "fub/core/function_ref.hpp"
#include "fub/core/mdspan.hpp"
#include "fub/grid/AMReX/PatchHandle.hpp"

#include <AMReX_FillPatchUtil.H>
#include <AMReX_Interpolater.H>
#include <AMReX_MultiFabUtil.H>

namespace fub {
namespace amrex {

struct BoundaryCondition : public ::amrex::PhysBCFunctBase {
  using Function =
      function_ref<void(const PatchDataView<double, AMREX_SPACEDIM + 1>&,
                        PatchHandle, Location, int, Duration)>;

  BoundaryCondition(Function f, const ::amrex::Geometry& geom, int level);

  void FillBoundary(::amrex::MultiFab& mf, int, int, double time_point,
                    int) override;

  Function function_;
  ::amrex::Geometry geom_;
  int level_num_;
};

} // namespace amrex
} // namespace fub

#endif
