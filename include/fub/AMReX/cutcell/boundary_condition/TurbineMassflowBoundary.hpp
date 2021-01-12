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

#ifndef FUB_AMREX_CUTCELL_TURBINE_MASS_FLOW_BOUNDARY_HPP
#define FUB_AMREX_CUTCELL_TURBINE_MASS_FLOW_BOUNDARY_HPP

#include "fub/AMReX/cutcell/GriddingAlgorithm.hpp"
#include "fub/equations/EulerEquation.hpp"

#include "fub/Direction.hpp"
#include "fub/NewtonIteration.hpp"
#include "fub/ext/Log.hpp"

#include <AMReX.H>

namespace fub::amrex::cutcell {
/// \ingroup BoundaryCondition
///
struct TurbineMassflowBoundaryOptions {
  TurbineMassflowBoundaryOptions() = default;
  TurbineMassflowBoundaryOptions(const ProgramOptions& options);

  void Print(SeverityLogger& log) const;

  std::string channel_name{"TurbineMassflowBoundary"};
  ::amrex::Box boundary_section{};
  double relative_surface_area = 1.0;
  double massflow_correlation = 0.06;
  Direction dir = Direction::X;
  int side = 0;
};

/// \ingroup BoundaryCondition
///
/// This is an outflow boundary condition that models the massflow condition of
/// a turbine machine.
///
/// The massflow is given by the relation
///
///         $$\dot{m} / A \cdot \frac{\sqrt{T_0}}{p_0} = \text{const}$$
///
/// Therefore, for given surface area $A$, total pressure $p_0$ and total
/// temperature $T_0$ one determines the required massflow $\dot m$ and
/// recomputes the static pressure and temperature values.
template <typename EulerEquation> class TurbineMassflowBoundary {
public:
  using Complete = ::fub::Complete<EulerEquation>;

  TurbineMassflowBoundary(const EulerEquation& eq,
                          const TurbineMassflowBoundaryOptions& options);

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

  void TransformState(Complete& expanded, const Complete& source);

private:
  EulerEquation equation_;
  TurbineMassflowBoundaryOptions options_;
};

template <typename EulerEquation>
TurbineMassflowBoundary<EulerEquation>::TurbineMassflowBoundary(
    const EulerEquation& eq, const TurbineMassflowBoundaryOptions& options)
    : equation_{eq}, options_{options} {}

template <typename EulerEquation>
void TurbineMassflowBoundary<EulerEquation>::FillBoundary(
    ::amrex::MultiFab& mf, const GriddingAlgorithm& grid, int level) {
  auto factory = grid.GetPatchHierarchy().GetEmbeddedBoundary(level);
  const ::amrex::MultiFab& alphas = factory->getVolFrac();
  FillBoundary(mf, alphas, grid.GetPatchHierarchy().GetGeometry(level));
}

template <typename EulerEquation>
void TurbineMassflowBoundary<EulerEquation>::FillBoundary(
    ::amrex::MultiFab& mf, const ::amrex::MultiFab& alphas,
    const ::amrex::Geometry& geom) {
  Complete state{equation_};
  Complete expanded{equation_};
  ForEachFab(execution::seq, mf, [&](const ::amrex::MFIter& mfi) {
    ::amrex::FArrayBox& fab = mf[mfi];
    const ::amrex::FArrayBox& alpha = alphas[mfi];
    ::amrex::Box box_to_fill = mfi.growntilebox() & options_.boundary_section;
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
              TransformState(expanded, state);
              Store(states, expanded, dest);
            }
          });
    }
  });
}

template <typename EulerEquation>
void TurbineMassflowBoundary<EulerEquation>::TransformState(
    Complete& dest, const Complete& source) {
  Primitive<EulerEquation> prim(equation_);
  PrimFromComplete(equation_, prim, source);
  const double rho = euler::Density(equation_, source);
  const double p = euler::Pressure(equation_, source);
  const double gamma = euler::Gamma(equation_, source);
  const double c = euler::SpeedOfSound(equation_, source);
  const double T = euler::Temperature(equation_, source);
  const double c2 = c * c;
  const int dir_v = static_cast<int>(options_.dir);
  const double u = prim.velocity[dir_v];
  const double u2 = u * u;
  const double Ma2 = u2 / c2;

  const double gammaMinus = gamma - 1.0;
  const double gammaMinusHalf = 0.5 * gammaMinus;
  const double gammaMinusInv = 1.0 / gammaMinus;
  const double gammaOverGammaMinus = gamma * gammaMinusInv;
  const double gammaPlus = gamma + 1.0;
  // const double gammaPlus = gamma + 1.0;

  // compute total quantities
  const double alpha_0 = 1.0 + Ma2 * gammaMinusHalf;
  const double p_0 = p * std::pow(alpha_0, gammaOverGammaMinus);
  const double rho_0 = rho * std::pow(alpha_0, gammaMinusInv);
  const double T_0 = T * alpha_0;
  const double c_02 = c2 + gammaMinusHalf * u2;
  // const double c_0 = std::sqrt(c_02);

  // const int enable_massflow = (p_0 > options_.pressure_threshold);
  const double required_massflow =
      options_.massflow_correlation * p_0 / std::sqrt(T_0);

  const double alpha = gammaMinusHalf / c_02;
  const double beta =
      required_massflow / (options_.relative_surface_area * rho_0);
  auto f = [alpha, beta, gammaMinusInv](double u_n) {
    return u_n * std::pow(1.0 - alpha * u_n * u_n, gammaMinusInv) - beta;
  };
  auto Df = [alpha, gammaMinusInv](double u_n) {
    const double x = 1.0 - alpha * u_n * u_n;
    return std::pow(x, gammaMinusInv) - 2.0 * alpha * u_n * u_n *
                                            gammaMinusInv *
                                            std::pow(x, gammaMinusInv - 1.0);
  };
  const double required_u = NewtonIteration(f, Df, u);
  // limit outflow velocity to a physically possible one
  const double critical_speed_of_sound = std::sqrt(2.0 * c_02 / (gamma + 1.0));
  const double required_laval = std::min(required_u / critical_speed_of_sound, 1.0);
  const double required_laval2 = required_laval * required_laval;
  
  const double alpha_n = 1.0 - gammaMinus / gammaPlus * required_laval2;
  prim.density = rho_0 * std::pow(alpha_n, gammaMinusInv);
  prim.pressure = p_0 * std::pow(alpha_n, gammaOverGammaMinus);
  prim.velocity[dir_v] = required_u;
  CompleteFromPrim(equation_, dest, prim);
}

} // namespace fub::amrex::cutcell

#endif // !FUB_AMREX_CUTCELL_ISENTROPIC_PRESSURE_BOUNDARY_HPP
