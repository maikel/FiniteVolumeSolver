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

#ifndef FUB_BOUNDARY_CONDITION_COUPLED
#define FUB_BOUNDARY_CONDITION_COUPLED

#include "fub/AMReX/coupled/CoupledGriddingAlgorithm.hpp"

#include "fub/AMReX/PatchHierarchy.hpp"
#include "fub/AMReX/cutcell/PatchHierarchy.hpp"
#include "fub/PatchDataView.hpp"
#include "fub/equations/IdealGasMix.hpp"

#include <AMReX.H>
#include <AMReX_MultiFab.H>

#include <vector>

namespace fub::amrex {

class CoupledBoundary {
public:
  static constexpr int Plenum_Rank = AMREX_SPACEDIM;
  static constexpr int Tube_Rank = 1;

  /// Constructs coupled boundary states by pre computing mirror and ghost
  /// states for each of the specified domains.
  ///
  /// This function might grow the specified mirror boxes to an extent which is
  /// required to fulfill the specified ghost cell width requirements.
  CoupledBoundary(CoupledGriddingAlgorithm& gridding,
                  const BlockConnection& connection, int gcw,
                  const FlameMasterReactor& reactor);

  /// Precompute Boundary states for each domain.
  ///
  /// Subsequent calls to FillXXXBoundary will use these new computed boundary
  /// states.
  ///
  /// \param[in] plenum  The higher dimensional patch hierarchy with geometry
  /// information. States here will be conservatively averaged and projected
  /// onto a one-dimesnional space.
  ///
  /// \param[in] tube The low dimensional tube data.
  void ComputeBoundaryData(const cutcell::PatchHierarchy& plenum,
                           const PatchHierarchy& tube);

  /// Assuming that mf represents a MultiFab living in the higher dimensional
  /// plenum simulation its ghost layer will be filled with data from the tube
  /// simulation.
  void FillBoundary(::amrex::MultiFab& mf, const ::amrex::Geometry& geom,
                    Duration time_point,
                    const cutcell::GriddingAlgorithm& gridding);

  /// Assuming that mf represents a MultiFab living in the one dimensional tube
  /// simulation its ghost layer will be filled with data from the plenum
  /// simulation.
  void FillBoundary(::amrex::MultiFab& mf, const ::amrex::Geometry& geom,
                    Duration time_point, const GriddingAlgorithm& gridding);

private:
  IdealGasMix<Plenum_Rank> plenum_equation_;
  IdealGasMix<Tube_Rank> tube_equation_;

  ::amrex::Box plenum_mirror_box_{};
  ::amrex::Box tube_mirror_box_{};

  ::amrex::FArrayBox plenum_mirror_data_{};
  ::amrex::FArrayBox tube_ghost_data_{};

  ::amrex::FArrayBox tube_mirror_data_{};
  ::amrex::FArrayBox plenum_ghost_data_{};

  Direction dir_{};
  int side_{};
};

//
// class CoupledBoundaryFunction {
// public:
//  CoupledBoundaryFunction(const cutcell::PatchHierarchy& plenum,
//                          const PatchHierarchy& tube, int gcw,
//                          const FlameMasterReactor& reactor)
//      : boundary_state_(
//            std::make_shared<CoupledBoundary>(plenum, tube, gcw, reactor)) {}
//
//  void operator()(const PatchDataView<double, AMREX_SPACEDIM + 1>& states,
//                  const cutcell::PatchHierarchy& plenum, PatchHandle patch,
//                  Location loc, int fill_width, Duration time_point) const {
//    boundary_state_->FillPlenumBoundary(states, plenum, patch, loc,
//    fill_width,
//                                        time_point);
//  }
//
//  void operator()(const PatchDataView<double, AMREX_SPACEDIM + 1>& states,
//                  const PatchHierarchy& plenum, PatchHandle patch, Location
//                  loc, int fill_width, Duration time_point) const {
//    boundary_state_->FillTubeBoundary(states, plenum, patch, loc, fill_width,
//                                      time_point);
//  }
//
//  const std::shared_ptr<CoupledBoundary>& GetSharedState() const noexcept {
//    return boundary_state_;
//  }
//
// private:
//  std::shared_ptr<CoupledBoundary> boundary_state_;
//};

} // namespace fub::amrex

#endif
