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

#ifndef FUB_AMREX_BOUNDARY_CONDITION_ISENTROPIC_PRESSURE_EXPANSION_HPP
#define FUB_AMREX_BOUNDARY_CONDITION_ISENTROPIC_PRESSURE_EXPANSION_HPP

#include "fub/AMReX/boundary_condition/IsentropicPressureBoundary.hpp"
#include "fub/equations/EulerEquation.hpp"

namespace fub::amrex {

inline std::array<std::ptrdiff_t, 1>
MapToSrc(const std::array<std::ptrdiff_t, 1>& dest,
         const ::amrex::Geometry& geom, int side, Direction dir) {
  const int boundary = (side == 0) * geom.Domain().smallEnd(int(dir)) +
                       (side == 1) * geom.Domain().bigEnd(int(dir));
  const auto distance = dest[std::size_t(dir)] - boundary;
  const int sign = int((distance > 0) - (distance < 0));
  std::array<std::ptrdiff_t, 1> src{dest};
  src[std::size_t(dir)] -= static_cast<int>(2 * distance - sign);
  //  src[std::size_t(dir)] = boundary;
  return src;
}

/// \ingroup BoundaryCondition
///
/// \brief This boundary models an isentropic pressure expansion for the
/// one-dimensional ideal gas equations for mixtures.
template <typename EulerEquation> class IsentropicPressureExpansion {
public:
  IsentropicPressureExpansion(const EulerEquation& eq,
                              const IsentropicPressureBoundaryOptions& options);

  IsentropicPressureExpansion(const EulerEquation& eq, double outer_pressure,
                              Direction dir, int side);

  void FillBoundary(::amrex::MultiFab& mf, const GriddingAlgorithm& gridding,
                    int level);

  void FillBoundary(::amrex::MultiFab& mf, const GriddingAlgorithm& gridding,
                    int level, Direction dir);

  void FillBoundary(::amrex::MultiFab& mf, const ::amrex::Geometry& geom);

private:
  EulerEquation equation_;
  double outer_pressure_;
  Direction dir_;
  int side_;
};

// this function seems to be incorrect... better not use it!
template <typename EulerEquation>
void ExpandState(EulerEquation& eq, Complete<EulerEquation>& dest,
                 const Complete<EulerEquation>& src, double pressure_dest,
                 double efficiency) {
  const auto old_velocity = euler::Velocity(eq, src);
  const double rhoE_kin = euler::KineticEnergy(src.density, src.momentum);
  static constexpr int N = EulerEquation::Rank();
  dest = src;
  dest.momentum = Array<double, N, 1>::Zero();
  dest.energy = src.energy - rhoE_kin;
  const double h_before = (dest.energy + dest.pressure) / dest.density;
  euler::SetIsentropicPressure(eq, dest, dest, pressure_dest);
  const double h_after = (dest.energy + dest.pressure) / dest.density;
  const double h_diff = h_after - h_before;
  const double e_kin_new =
      2.0 * (1.0 - efficiency) * h_diff + old_velocity.matrix().squaredNorm();
  FUB_ASSERT(e_kin_new >= 0);
  const Array<double, N, 1> u_border =
      std::sqrt(e_kin_new) * old_velocity.matrix().normalized();
  dest.momentum = dest.density * u_border;
  dest.energy = dest.energy + euler::KineticEnergy(dest.density, dest.momentum);
}

template <typename EulerEquation>
IsentropicPressureExpansion<EulerEquation>::IsentropicPressureExpansion(
    const EulerEquation& eq, const IsentropicPressureBoundaryOptions& options)
    : equation_{eq}, outer_pressure_{options.outer_pressure},
      dir_{options.direction}, side_{options.side} {}

template <typename EulerEquation>
IsentropicPressureExpansion<EulerEquation>::IsentropicPressureExpansion(
    const EulerEquation& eq, double outer_pressure, Direction dir, int side)
    : equation_{eq}, outer_pressure_{outer_pressure}, dir_{dir}, side_{side} {}

template <typename EulerEquation>
void IsentropicPressureExpansion<EulerEquation>::FillBoundary(
    ::amrex::MultiFab& mf, const GriddingAlgorithm& gridding, int level) {
  FillBoundary(mf, gridding.GetPatchHierarchy().GetGeometry(level));
}

template <typename EulerEquation>
void IsentropicPressureExpansion<EulerEquation>::FillBoundary(
    ::amrex::MultiFab& mf, const GriddingAlgorithm& gridding, int level,
    Direction dir) {
  if (dir == dir_) {
    FillBoundary(mf, gridding, level);
  }
}

template <typename EulerEquation>
void IsentropicPressureExpansion<EulerEquation>::FillBoundary(
    ::amrex::MultiFab& mf, const ::amrex::Geometry& geom) {
  const int ngrow = mf.nGrow(int(dir_));
  ::amrex::Box grown_box = geom.growNonPeriodicDomain(ngrow);
  ::amrex::BoxList boundaries =
      ::amrex::complementIn(grown_box, ::amrex::BoxList{geom.Domain()});
  Complete<EulerEquation> state{equation_};
  Complete<EulerEquation> expanded{equation_};
  if (boundaries.isEmpty()) {
    return;
  }
  ForEachFab(execution::seq, mf, [&](const ::amrex::MFIter& mfi) {
    ::amrex::FArrayBox& fab = mf[mfi];
    for (const ::amrex::Box& boundary : boundaries) {
      ::amrex::Box shifted =
          ::amrex::shift(boundary, int(dir_), GetSign(side_) * ngrow);
      if (!geom.Domain().intersects(shifted)) {
        continue;
      }
      ::amrex::Box box_to_fill = mfi.growntilebox() & boundary;
      if (!box_to_fill.isEmpty()) {
        auto states = MakeView<Complete<EulerEquation>>(fab, equation_,
                                                        mfi.growntilebox());
        ForEachIndex(AsIndexBox<EulerEquation::Rank()>(box_to_fill),
                     [&](auto... is) {
                       Index<EulerEquation::Rank()> dest{is...};
                       Index<EulerEquation::Rank()> src =
                           MapToSrc(dest, geom, side_, dir_);
                       Load(state, states, src);
                       euler::SetIsentropicPressure(equation_, expanded, state, outer_pressure_);
                      //  ExpandState(equation_, expanded, state, outer_pressure_, 1.0);
                       Store(states, expanded, dest);
                     });
      }
    }
  });
}

} // namespace fub::amrex

#endif