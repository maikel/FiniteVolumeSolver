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

#ifndef FUB_AMREX_CUTCELL_BOUNDARY_CONDITION_ISENTROPIC_PRESSURE_EXPANSION_HPP
#define FUB_AMREX_CUTCELL_BOUNDARY_CONDITION_ISENTROPIC_PRESSURE_EXPANSION_HPP

#include "fub/AMReX/boundary_condition/IsentropicPressureBoundary.hpp"
#include "fub/AMReX/cutcell/GriddingAlgorithm.hpp"
#include "fub/equations/EulerEquation.hpp"

namespace fub::amrex::cutcell {

inline std::array<std::ptrdiff_t, 2>
MapToSrc(const std::array<std::ptrdiff_t, 2>& dest,
         const ::amrex::Geometry& geom, int side, Direction dir) {
  const int boundary = (side == 0) * geom.Domain().smallEnd(int(dir)) +
                       (side == 1) * geom.Domain().bigEnd(int(dir));
  const auto distance = dest[std::size_t(dir)] - boundary;
  const int sign = int((distance > 0) - (distance < 0));
  std::array<std::ptrdiff_t, 2> src{dest};
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
  IsentropicPressureExpansion(
      const EulerEquation& eq,
      const fub::amrex::IsentropicPressureBoundaryOptions& options);

  void FillBoundary(::amrex::MultiFab& mf, const GriddingAlgorithm& gridding,
                    int level);

  void FillBoundary(::amrex::MultiFab& mf, const GriddingAlgorithm& gridding,
                    int level, Direction dir);

  void FillBoundary(::amrex::MultiFab& mf, const ::amrex::MultiFab& alphas,
                    const ::amrex::Geometry& geom);

private:
  EulerEquation equation_;
  fub::amrex::IsentropicPressureBoundaryOptions options_;
};

template <typename EulerEquation>
void ExpandState(EulerEquation& eq, Complete<EulerEquation>& dest,
                 const Complete<EulerEquation>& src, double pressure_dest,
                 double efficiency) {
  static constexpr int N = EulerEquation::Rank();
  const Array<double, N, 1> old_velocity = euler::Velocity(eq, src);
  dest = src;
  euler::SetVelocity(eq, dest, Array<double, N, 1>::Zero());
  const double h_before = euler::TotalEnthalpy(eq, dest);
  euler::SetIsentropicPressure(eq, dest, dest, pressure_dest);
  const double h_after = euler::TotalEnthalpy(eq, dest);
  const double h_diff = h_before - h_after;
  auto Sign = [](double x) { return (x >= 0) - (x < 0); };
  const double e_kin = efficiency * h_diff * 2.0 + old_velocity[0] * old_velocity[0];
  const double u_border0 = Sign(e_kin) * std::sqrt(std::abs(e_kin));
  Array<double, N, 1> u_border = old_velocity;
  u_border[0] = u_border0;
  euler::SetVelocity(eq, dest, u_border);
}

template <typename EulerEquation>
IsentropicPressureExpansion<EulerEquation>::IsentropicPressureExpansion(
    const EulerEquation& eq,
    const fub::amrex::IsentropicPressureBoundaryOptions& options)
    : equation_{eq}, options_{options} {}

template <typename EulerEquation>
void IsentropicPressureExpansion<EulerEquation>::FillBoundary(
    ::amrex::MultiFab& mf, const GriddingAlgorithm& grid, int level) {
  auto factory = grid.GetPatchHierarchy().GetEmbeddedBoundary(level);
  const ::amrex::MultiFab& alphas = factory->getVolFrac();
  FillBoundary(mf, alphas, grid.GetPatchHierarchy().GetGeometry(level));
}

template <typename EulerEquation>
void IsentropicPressureExpansion<EulerEquation>::FillBoundary(
    ::amrex::MultiFab& mf, const GriddingAlgorithm& gridding, int level,
    Direction dir) {
  if (dir == options_.direction) {
    FillBoundary(mf, gridding, level);
  }
}

template <typename EulerEquation>
void IsentropicPressureExpansion<EulerEquation>::FillBoundary(
    ::amrex::MultiFab& mf, const ::amrex::MultiFab& alphas,
    const ::amrex::Geometry& geom) {
  const int ngrow = mf.nGrow(int(options_.direction));
  ::amrex::Box grown_box = geom.growNonPeriodicDomain(ngrow);
  ::amrex::BoxList boundaries =
      ::amrex::complementIn(grown_box, ::amrex::BoxList{geom.Domain()});
  if (boundaries.isEmpty()) {
    return;
  }
  Complete<EulerEquation> state{equation_};
  Complete<EulerEquation> expanded{equation_};
  ForEachFab(execution::seq, mf, [&](const ::amrex::MFIter& mfi) {
    ::amrex::FArrayBox& fab = mf[mfi];
    const ::amrex::FArrayBox& alpha = alphas[mfi];
    for (const ::amrex::Box& boundary : boundaries) {
      ::amrex::Box shifted = ::amrex::shift(boundary, int(options_.direction),
                                            GetSign(options_.side) * ngrow);
      if (!geom.Domain().intersects(shifted)) {
        continue;
      }
      ::amrex::Box box_to_fill = mfi.growntilebox() & boundary;
      if (!box_to_fill.isEmpty()) {
        auto states = MakeView<Complete<EulerEquation>>(fab, equation_,
                                                        mfi.growntilebox());
        ForEachIndex(
            AsIndexBox<EulerEquation::Rank()>(box_to_fill), [&](auto... is) {
              Index<EulerEquation::Rank()> dest{is...};
              Index<EulerEquation::Rank()> src =
                  MapToSrc(dest, geom, options_.side, options_.direction);
              ::amrex::IntVect src_iv{
                  AMREX_D_DECL(int(src[0]), int(src[1]), int(src[2]))};
              ::amrex::IntVect iv{
                  AMREX_D_DECL(int(dest[0]), int(dest[1]), int(dest[2]))};
              if (alpha(iv) > 0.0 && alpha(src_iv) > 0.0) {
                Load(state, states, src);
                ExpandState(equation_, expanded, state, options_.outer_pressure,
                            options_.efficiency);
                Store(states, expanded, dest);
              }
            });
      }
    }
  });
}

} // namespace fub::amrex::cutcell

#endif