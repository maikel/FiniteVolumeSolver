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

#ifndef FUB_AMREX_CUTCELL_MASSFLOW_BOUNDARY2_HPP
#define FUB_AMREX_CUTCELL_MASSFLOW_BOUNDARY2_HPP

#include "fub/AMReX/cutcell/GriddingAlgorithm.hpp"
#include "fub/equations/EulerEquation.hpp"

#include "fub/Direction.hpp"

#include <AMReX.H>

namespace fub::amrex::cutcell {
/// \ingroup BoundaryCondition
///
struct MachnumberBoundaryOptions {
  MachnumberBoundaryOptions() = default;
  MachnumberBoundaryOptions(const ProgramOptions& options);

  std::string channel_name{"MachnumberBoundary"};
  ::amrex::Box boundary_section{};
  double required_mach_number = 1.0;
  Direction dir = Direction::X;
  int side = 0;
};

/// \ingroup BoundaryCondition
///
template <typename EulerEquation> class MachnumberBoundary {
public:
  using Complete = ::fub::Complete<EulerEquation>;

  MachnumberBoundary(const EulerEquation& eq,
                     const MachnumberBoundaryOptions& options);

  void FillBoundary(::amrex::MultiFab& mf, const GriddingAlgorithm& gridding,
                    int level);

  void FillBoundary(::amrex::MultiFab& mf, const ::amrex::MultiFab& alphas,
                    const ::amrex::Geometry& geom);

  void FillBoundary(::amrex::MultiFab& mf, const GriddingAlgorithm& gridding,
                    int level, Direction dir) {
    if (dir == options_.dir) {
      FillBoundary(mf, gridding, level);
    }
  }

  void ExpandState(Complete& expanded, const Complete& source);

private:
  EulerEquation equation_;
  MachnumberBoundaryOptions options_;
};

template <typename EulerEquation>
MachnumberBoundary<EulerEquation>::MachnumberBoundary(
    const EulerEquation& eq, const MachnumberBoundaryOptions& options)
    : equation_{eq}, options_{options} {}

template <typename EulerEquation>
void MachnumberBoundary<EulerEquation>::FillBoundary(
    ::amrex::MultiFab& mf, const GriddingAlgorithm& grid, int level) {
  auto factory = grid.GetPatchHierarchy().GetEmbeddedBoundary(level);
  const ::amrex::MultiFab& alphas = factory->getVolFrac();
  FillBoundary(mf, alphas, grid.GetPatchHierarchy().GetGeometry(level));
}

template <typename EulerEquation>
void MachnumberBoundary<EulerEquation>::FillBoundary(
    ::amrex::MultiFab& mf, const ::amrex::MultiFab& alphas,
    const ::amrex::Geometry& geom) {
  Complete state{equation_};
  Complete expanded{equation_};
  ForEachFab(execution::seq, mf, [&](const ::amrex::MFIter& mfi) {
    ::amrex::FArrayBox& fab = mf[mfi];
    const ::amrex::FArrayBox& alpha = alphas[mfi];
    ::amrex::Box box_to_fill =
          mfi.growntilebox() & options_.boundary_section;
      if (!box_to_fill.isEmpty()) {
        auto states = MakeView<Complete>(fab, equation_, mfi.growntilebox());
        ForEachIndex(
            AsIndexBox<EulerEquation::Rank()>(box_to_fill), [&](auto... is) {
              Index<EulerEquation::Rank()> dest{is...};
              Index<EulerEquation::Rank()> src =
                  MapToSrc(dest, geom, options_.side, options_.dir);
              ::amrex::IntVect src_iv{
                  AMREX_D_DECL(int(src[0]), int(src[1]), int(src[2]))};
              ::amrex::IntVect iv{
                  AMREX_D_DECL(int(dest[0]), int(dest[1]), int(dest[2]))};
              if (alpha(iv) > 0.0 && alpha(src_iv) > 0.0) {
                Load(state, states, src);
                ExpandState(expanded, state);
                Store(states, expanded, dest);
              }
            });
      }
  });
}

template <typename EulerEquation>
void MachnumberBoundary<EulerEquation>::ExpandState(Complete& expanded,
                                                    const Complete& source) {
  Primitive<EulerEquation> prim(equation_);
  PrimFromComplete(equation_, prim, source);
  const double rho = euler::Density(equation_, source);
  const double p = euler::Pressure(equation_, source);
  const double gamma = euler::Gamma(equation_, source);
  const double c = euler::SpeedOfSound(equation_, source);
  const double c2 = c*c;
  const int dir_v = static_cast<int>(options_.dir);
  const double u = prim.velocity[dir_v];
  const double u2 = u * u;
  const double Ma2 = u2 / c2;

  const double gammaMinus = gamma - 1.0;
  const double gammaPlus = gamma + 1.0;

  // compute total pressure and total density
  const double p_0 = p * std::pow(1.0 + 0.5 * Ma2 * gammaMinus, gamma / gammaMinus);
  const double rho_0 = rho * std::pow(1.0 + 0.5 * Ma2 * gammaMinus, 1.0 / gammaMinus);

  const double Ma_r = options_.required_mach_number;
  const double Ma2_r = Ma_r * Ma_r;
  const double laval2_n = Ma2_r * gammaPlus / (Ma2_r * gammaMinus + 2.0);

  // Use required laval number to determine new static values
  const double p_n = p_0 * std::pow(1.0 - gammaMinus / gammaPlus * laval2_n, gamma / gammaMinus);
  const double rho_n = rho_0 * std::pow(1.0 - gammaMinus / gammaPlus * laval2_n, 1.0 / gammaPlus);
  prim.density = rho_n;
  prim.pressure = p_n;

  const double c2_critical = (c2 + 0.5 * gammaMinus * u2) * 2.0 / gammaPlus;
  const double u2_n = laval2_n * c2_critical;
  prim.velocity[dir_v] = std::sqrt(u2_n);
  CompleteFromPrim(equation_, expanded, prim);
}

} // namespace fub::amrex::cutcell

#endif // !FUB_AMREX_CUTCELL_ISENTROPIC_PRESSURE_BOUNDARY_HPP
