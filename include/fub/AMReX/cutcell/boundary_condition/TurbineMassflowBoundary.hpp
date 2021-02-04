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

/// \file
/// This file defines a Boundary Condition which emulates a Plenum-Turbine
/// Interaction.

#ifndef FUB_AMREX_CUTCELL_TURBINE_MASS_FLOW_BOUNDARY_HPP
#define FUB_AMREX_CUTCELL_TURBINE_MASS_FLOW_BOUNDARY_HPP

#include "fub/AMReX/AverageState.hpp"
#include "fub/AMReX/cutcell/GriddingAlgorithm.hpp"
#include "fub/AMReX/cutcell/boundary_condition/ConstantBoundary.hpp"
#include "fub/equations/EulerEquation.hpp"

#include "fub/Direction.hpp"
#include "fub/NewtonIteration.hpp"
#include "fub/ext/Log.hpp"

#include <AMReX.H>
#include <boost/log/utility/manipulators/add_value.hpp>
#include <boost/math/tools/roots.hpp>

namespace fub {
/// \brief Compute the new complete state from the old complete state to
/// enforce the required massflow. This is done by following the method from
/// Jirasek. TODO: cite references
struct RequireMassflow_Jirasek {
  template <typename EulerEquation>
  void operator()(EulerEquation& eq, Complete<EulerEquation>& expanded,
                  const Complete<EulerEquation>& source,
                  double required_massflow, double relative_surface_area,
                  Direction dir) const noexcept;
};

/// \brief Compute the new complete state from the old complete state to
/// enforce the required massflow. This is done by solving the exact Riemann
/// Problem as described in Toro. TODO: cite references
struct RequireMassflow_SolveExactRiemannProblem {
  template <typename EulerEquation>
  void operator()(EulerEquation& eq, Complete<EulerEquation>& expanded,
                  const Complete<EulerEquation>& source,
                  double required_massflow, double relative_surface_area,
                  Direction dir) const noexcept;
};
} // namespace fub

namespace fub::amrex::cutcell {
enum class TurbineMassflowMode {
  cellwise,
  average_inner_state,
  average_outer_state,
  average_massflow
};

/// \ingroup BoundaryCondition
///
struct TurbineMassflowBoundaryOptions {
  TurbineMassflowBoundaryOptions() = default;
  TurbineMassflowBoundaryOptions(const ProgramOptions& options);

  void Print(SeverityLogger& log) const;

  std::string channel_name{
      "TurbineMassflowBoundary"}; ///< the channel name for logging
  ::amrex::Box
      boundary_section{}; ///< the amrex box where the boundary is located
  TurbineMassflowMode mode =
      TurbineMassflowMode::cellwise; ///< the mode how the massflow is enforced
  ::amrex::Box coarse_average_mirror_box{};
  boost::uintmax_t max_iter{
      100}; ///< allowed maximum iterations for the bisection root fnding
  double relative_surface_area =
      1.0; ///< the relative surface area for the turbine
  /// the massflow correlation, this is given by the idealized compressor
  /// characteristics
  double massflow_correlation = 0.06;
  Direction dir =
      Direction::X; ///< the dimensional split direction which will be used
  int side = 0; ///< the side where the boundary condition is applied (0==left, 1==right)
};

/// \ingroup BoundaryCondition
///
/// This is an outflow boundary condition that models the massflow condition of
/// a turbine machine.
///
/// The massflow is given by the relation
///
///         \f$\dot{m} / A \cdot \frac{\sqrt{T_0}}{p_0} = \text{const}\f$
///
/// Therefore, for given surface area \f$A\f$, total pressure \f$p_0\f$ and
/// total temperature \f$T_0\f$ one determines the required massflow \f$\dot
/// m\f$ and recomputes the static pressure and temperature values.
template <typename EulerEquation,
          typename Transform = RequireMassflow_SolveExactRiemannProblem>
class TurbineMassflowBoundary {
public:
  using Complete = ::fub::Complete<EulerEquation>;
  using Conservative = ::fub::Conservative<EulerEquation>;

  TurbineMassflowBoundary(const EulerEquation& eq,
                          const TurbineMassflowBoundaryOptions& options);

  void FillBoundary(::amrex::MultiFab& mf, const GriddingAlgorithm& gridding,
                    int level);

  void FillBoundary(::amrex::MultiFab& mf, const ::amrex::MultiFab& alphas,
                    const ::amrex::Geometry& geom);

  void FillBoundary(::amrex::MultiFab& mf, const ::amrex::MultiFab& alphas,
                    const ::amrex::Geometry& geom, double required_massflow);

  void FillBoundary(::amrex::MultiFab& mf, const GriddingAlgorithm& gridding,
                    int level, Direction dir) {
    if (dir == options_.dir) {
      FillBoundary(mf, gridding, level);
    }
  }

  /// @{
  /// \brief Compute the new complete state from the old complete state to
  /// enforce the required massflow. This is done by solving the exact Riemann
  /// Problem as described in Toro or by following the method from Jirasek.
  /// TODO: cite references
  ///
  /// \param source the old complete state
  /// \param expanded the new complete state which enforces the required
  /// massflow condition
  void TransformState(Complete& expanded, const Complete& source);
  void TransformState(Complete& expanded, const Complete& source,
                      double required_massflow);
  /// @}

  /// \brief Returns the required mass flow for a given complete state.
  double ComputeRequiredMassflow(const Complete& state) const;

  /// \brief Returns the averaged required mass flow over all boundary cells.
  double AverageRequiredMassflow(::amrex::MultiFab& mf,
                                 const GriddingAlgorithm& grid, int level);

private:
  template <typename I>
  static I MapToSrc(I& dest, const ::amrex::Geometry& geom, int side,
                    Direction dir) noexcept {
    const int boundary = (side == 0) * geom.Domain().smallEnd(int(dir)) +
                         (side == 1) * geom.Domain().bigEnd(int(dir));
    I src{dest};
    src[int(dir)] = boundary;
    return src;
  }

  EulerEquation equation_;
  TurbineMassflowBoundaryOptions options_;
  Transform transform_{};
};

template <typename EulerEquation, typename Transform>
TurbineMassflowBoundary<EulerEquation, Transform>::TurbineMassflowBoundary(
    const EulerEquation& eq, const TurbineMassflowBoundaryOptions& options)
    : equation_{eq}, options_{options} {}

template <typename EulerEquation, typename Transform>
double
TurbineMassflowBoundary<EulerEquation, Transform>::AverageRequiredMassflow(
    ::amrex::MultiFab& mf, const GriddingAlgorithm& grid, int level) {
  Complete local_state{equation_};
  double local_average_massflow = 0.0;
  const double n =
      static_cast<double>(options_.coarse_average_mirror_box.numPts());
  AccumulateState(
      mf, options_.coarse_average_mirror_box, local_average_massflow,
      [&](double& average, span<const double> buffer) {
        CopyFromBuffer(local_state, buffer);
        const double local_massflow = ComputeRequiredMassflow(local_state);
        average += local_massflow / n;
      });
  double required_massflow = 0.0;
  ::MPI_Allreduce(&local_average_massflow, &required_massflow, 1, MPI_DOUBLE,
                  MPI_SUM, ::amrex::ParallelDescriptor::Communicator());
  return required_massflow;
}

template <typename EulerEquation, typename Transform>
void TurbineMassflowBoundary<EulerEquation, Transform>::FillBoundary(
    ::amrex::MultiFab& mf, const GriddingAlgorithm& grid, int level) {
  auto factory = grid.GetPatchHierarchy().GetEmbeddedBoundary(level);
  SeverityLogger debug = GetLogger(boost::log::trivial::debug);
  BOOST_LOG_SCOPED_LOGGER_TAG(debug, "Channel", options_.channel_name);
  BOOST_LOG_SCOPED_LOGGER_TAG(debug, "Time", grid.GetTimePoint().count());
  const double required_massflow = AverageRequiredMassflow(mf, grid, level);
  BOOST_LOG(debug) << fmt::format("Required massflow: {}", required_massflow)
                   << boost::log::add_value("required_massflow",
                                            required_massflow);
  if (options_.mode == TurbineMassflowMode::average_inner_state) {
    Conservative cons{equation_};
    Complete state{equation_};
    Complete expanded{equation_};
    const ::amrex::Geometry& geom = grid.GetPatchHierarchy().GetGeometry(level);
    AverageState(cons, mf, geom, options_.coarse_average_mirror_box);
    CompleteFromCons(equation_, state, cons);
    BOOST_LOG(debug) << fmt::format(
        "Average inner state: {} kg/m3, {} m/s, {} Pa", state.density,
        state.momentum.matrix().norm() / state.density, state.pressure);
    TransformState(expanded, state);
    ConstantBoundary<EulerEquation> constant_boundary{
        options_.dir, options_.side, equation_, expanded};
    const double mdot = expanded.momentum[static_cast<int>(options_.dir)] *
                        options_.relative_surface_area;
    BOOST_LOG(debug) << fmt::format(
        "Homogenous ghost state: {} kg/m3, {} m/s, {} Pa => massflow: {} kg/s, "
        "{:12g} kg/(m2 s)",
        expanded.density, expanded.momentum.matrix().norm() / expanded.density,
        expanded.pressure, mdot,
        expanded.momentum[static_cast<int>(options_.dir)]);
    constant_boundary.FillBoundary(mf, grid, level);
  } else if (options_.mode == TurbineMassflowMode::average_massflow) {
    const ::amrex::MultiFab& alphas = factory->getVolFrac();
    FillBoundary(mf, alphas, grid.GetPatchHierarchy().GetGeometry(level),
                 required_massflow);
  } else if (options_.mode == TurbineMassflowMode::average_outer_state) {
    Conservative avg{equation_};
    Conservative cons{equation_};
    Complete state{equation_};
    Complete expanded{equation_};

    const double n =
        static_cast<double>(options_.coarse_average_mirror_box.numPts());
    AccumulateState(mf, options_.coarse_average_mirror_box, avg,
                    [&](Conservative& average, span<const double> buffer) {
                      CopyFromBuffer(cons, buffer);
                      CompleteFromCons(equation_, state, cons);
                      TransformState(expanded, state);
                      ForEachVariable(
                          [n](auto&& avg, auto&& q) { avg += q / n; }, average,
                          AsCons(expanded));
                    });
    const int n_comps =
        grid.GetPatchHierarchy().GetDataDescription().n_cons_components;
    std::vector<double> local_buffer(n_comps);
    CopyToBuffer(local_buffer, avg);
    std::vector<double> global_buffer(n_comps);
    MPI_Allreduce(local_buffer.data(), global_buffer.data(), n_comps,
                  MPI_DOUBLE, MPI_SUM,
                  ::amrex::ParallelDescriptor::Communicator());
    CopyFromBuffer(avg, global_buffer);
    CompleteFromCons(equation_, state, avg);
    const double mdot = state.momentum[static_cast<int>(options_.dir)] *
                        options_.relative_surface_area;
    BOOST_LOG(debug) << fmt::format(
        "Average ghost state: {} kg/m3, {} m/s, {} Pa => massflow: {} kg/s",
        state.density, state.momentum.matrix().norm() / state.density,
        state.pressure, mdot);
    ConstantBoundary<EulerEquation> constant_boundary{
        options_.dir, options_.side, equation_, state};
    constant_boundary.FillBoundary(mf, grid, level);
  } else if (options_.mode == TurbineMassflowMode::cellwise) {
    const ::amrex::MultiFab& alphas = factory->getVolFrac();
    FillBoundary(mf, alphas, grid.GetPatchHierarchy().GetGeometry(level));
  }
}

template <typename EulerEquation, typename Transform>
void TurbineMassflowBoundary<EulerEquation, Transform>::FillBoundary(
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

template <typename EulerEquation, typename Transform>
void TurbineMassflowBoundary<EulerEquation, Transform>::FillBoundary(
    ::amrex::MultiFab& mf, const ::amrex::MultiFab& alphas,
    const ::amrex::Geometry& geom, double required_massflow) {
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
              TransformState(expanded, state, required_massflow);
              Store(states, expanded, dest);
            }
          });
    }
  });
}

template <typename EulerEquation, typename Transform>
void TurbineMassflowBoundary<EulerEquation, Transform>::TransformState(
    Complete& dest, const Complete& source, double required_massflow) {
  transform_(equation_, dest, source, required_massflow,
             options_.relative_surface_area, options_.dir);
}

template <typename EulerEquation, typename Transform>
double
TurbineMassflowBoundary<EulerEquation, Transform>::ComputeRequiredMassflow(
    const Complete& state) const {
  const double p = euler::Pressure(equation_, state);
  const double c = euler::SpeedOfSound(equation_, state);
  const double T = euler::Temperature(equation_, state);
  const double c2 = c * c;
  const int dir_v = static_cast<int>(options_.dir);
  const double u = euler::Velocity(equation_, state, dir_v);
  const double u2 = u * u;
  const double Ma2 = u2 / c2;

  const double gamma = euler::Gamma(equation_, state);
  const double gm1 = gamma - 1.0;
  const double oogm1 = 1.0 / gm1;

  // compute total quantities
  const double alpha_0 = 1.0 + 0.5 * Ma2 * gm1;
  const double p_0 = p * std::pow(alpha_0, gamma * oogm1);
  const double T_0 = T * alpha_0;

  // const int enable_massflow = (p_0 > options_.pressure_threshold);
  const double required_massflow =
      options_.massflow_correlation * p_0 / std::sqrt(T_0);

  return required_massflow;
}

template <typename EulerEquation, typename Transform>
void TurbineMassflowBoundary<EulerEquation, Transform>::TransformState(
    Complete& dest, const Complete& source) {
  const double required_massflow = ComputeRequiredMassflow(source);
  std::invoke(transform_, equation_, dest, source, required_massflow,
              options_.relative_surface_area, options_.dir);
}

} // namespace fub::amrex::cutcell

namespace fub {
template <typename EulerEquation>
void RequireMassflow_Jirasek::
operator()(EulerEquation& eq, Complete<EulerEquation>& expanded,
           const Complete<EulerEquation>& source, double required_massflow,
           double relative_surface_area, Direction dir) const noexcept {
  Primitive<EulerEquation> prim(eq);
  PrimFromComplete(eq, prim, source);
  const double rho = euler::Density(eq, source);
  const double p = euler::Pressure(eq, source);
  const double gamma = euler::Gamma(eq, source);
  const double c = euler::SpeedOfSound(eq, source);
  const double T = euler::Temperature(eq, source);
  const double c2 = c * c;
  const int dir_v = static_cast<int>(dir);
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
  // const double required_massflow =
  //     options_.massflow_correlation * p_0 / std::sqrt(T_0);

  const double alpha = gammaMinusHalf / c_02;
  const double beta = required_massflow / (relative_surface_area * rho_0);
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
  const double required_laval =
      std::min(required_u / critical_speed_of_sound, 1.0);
  const double required_laval2 = required_laval * required_laval;

  const double alpha_n = 1.0 - gammaMinus / gammaPlus * required_laval2;
  prim.density = rho_0 * std::pow(alpha_n, gammaMinusInv);
  prim.pressure = p_0 * std::pow(alpha_n, gammaOverGammaMinus);
  prim.velocity[dir_v] = required_u;
  CompleteFromPrim(eq, expanded, prim);
}

template <typename EulerEquation>
void RequireMassflow_SolveExactRiemannProblem::
operator()(EulerEquation& eq, Complete<EulerEquation>& expanded,
           const Complete<EulerEquation>& source, double required_massflow,
           double relative_surface_area, Direction dir) const noexcept {
  Primitive<EulerEquation> prim(eq);
  PrimFromComplete(eq, prim, source);
  const double rho = euler::Density(eq, source);
  const double p = euler::Pressure(eq, source);
  const double c = euler::SpeedOfSound(eq, source);
  const int dir_v = static_cast<int>(dir);
  const double u = prim.velocity[dir_v];

  const double gamma = euler::Gamma(eq, source);
  const double gp1 = gamma + 1.0;
  const double gm1 = gamma - 1.0;
  const double oogm1 = 1.0 / gm1;
  const double gm1_over_gp1 = gm1 / gp1;
  // const double gammaPlus = gamma + 1.0;

  struct UserData {
    double rhoL;
    double uL;
    double pL;
    double aL;
    double gamma;
    double gm1;
    double gp1;
    double gm1_over_gp1;
    double oogm1;
    double rhou_star;
  } data;

  data.rhoL = rho;
  data.uL = u;
  data.pL = p;
  data.aL = c;
  data.gamma = gamma;
  data.gm1 = gm1;
  data.gp1 = gp1;
  data.gamma = gamma;
  data.gm1_over_gp1 = gm1_over_gp1;
  data.oogm1 = oogm1;
  data.rhou_star = required_massflow / relative_surface_area;

  // Toro p.122 f_L for the case: left wave is shock wave
  static auto f_shock = [](double p, const UserData& d) {
    const double A = 2.0 / (d.gp1 * d.rhoL);
    const double B = d.pL * d.gm1_over_gp1;
    const double ooQ = std::sqrt(A / (p + B));
    return (p - d.pL) * ooQ;
  };

  // Toro p.122 eq(4.21) u_L for the case: left wave is shock wave
  static auto u_shock = [](double p, const UserData& d) {
    return d.uL - f_shock(p, d);
  };

  // Toro p.122 eq(4.19) rho_L for the case: left wave is shock wave
  static auto rho_shock = [](double p, const UserData& d) {
    FUB_ASSERT(p > 0.0);
    const double p_rel = p / d.pL;
    const double nom = d.gm1_over_gp1 + p_rel;
    const double denom = d.gm1_over_gp1 * p_rel + 1.0;
    FUB_ASSERT(denom > 0.0);
    const double coeff = nom / denom;
    return d.rhoL * coeff;
  };

  // Toro p.123 f_L for the case: left wave is rarefaction wave
  static auto f_rarefaction = [](double p, const UserData& d) {
    const double pre_coeff = 2.0 * d.aL * d.oogm1;
    const double p_rel = p / d.pL;
    const double exponent = d.gm1 / 2.0 / d.gamma;
    return pre_coeff * (std::pow(p_rel, exponent) - 1.0);
  };

  // Toro p.123 eq(4.26) u_L for the case: left wave is rarefaction wave
  static auto u_rarefaction = [](double p, const UserData& d) {
    return d.uL - f_rarefaction(p, d);
  };

  // Toro p.123 eq(4.23) rho_L for the case: left wave is rarefaction wave
  static auto rho_rarefaction = [](double p, const UserData& d) {
    FUB_ASSERT(p > 0.0);
    const double p_rel = p / d.pL;
    const double exponent = 1.0 / d.gamma;
    return d.rhoL * std::pow(p_rel, exponent);
  };

  static auto delta_rhou_shock = [](double p, const UserData& d) {
    const double u_star = u_shock(p, d);
    const double rho_star = rho_shock(p, d);
    return rho_star * u_star - d.rhou_star;
  };

  static auto delta_rhou_rarefaction = [](double p, const UserData& d) {
    const double u_star = u_rarefaction(p, d);
    const double rho_star = rho_rarefaction(p, d);
    return rho_star * u_star - d.rhou_star;
  };

  auto delta_rhou = [&data](double p) {
    if (p < data.pL) {
      return delta_rhou_rarefaction(p, data);
    } else {
      return delta_rhou_shock(p, data);
    }
  };

  auto FindLowerAndUpperBounds = [&] {
    using Result = std::optional<std::pair<double, double>>;
    double p0 = data.pL;
    double upper = p0;
    double lower = p0;
    double drhou = delta_rhou(p0);
    int counter = 12;
    if (drhou < 0.0) {
      do {
        lower *= 0.5;
        if (counter-- == 0) {
          return Result{};
        }
        drhou = delta_rhou(lower);
      } while (drhou < 0.0);
    } else {
      do {
        upper *= 2.0;
        if (counter-- == 0) {
          return Result{};
        }
        drhou = delta_rhou(upper);
      } while (drhou > 0.0);
    }
    return Result{std::pair{lower, upper}};
  };

  // Attempt to find lower and upper bounds for the bisection
  std::optional bounds = FindLowerAndUpperBounds();
  if (bounds) {
    auto [lower, upper] = *bounds;
    const int digits = std::numeric_limits<double>::digits;
    const int get_digits = digits - 5;
    boost::math::tools::eps_tolerance<double> tol(get_digits);
    boost::uintmax_t it = 100; // TODO: options_.max_iter;
    auto [lo, hi] =
        boost::math::tools::bisect(delta_rhou, lower, upper, tol, it);
    const double p_star = lo + 0.5 * (hi - lo);
    const double u_star =
        p_star < data.pL ? u_rarefaction(p_star, data) : u_shock(p_star, data);
    const double rho_star = p_star < data.pL ? rho_rarefaction(p_star, data)
                                             : rho_shock(p_star, data);
    prim.pressure = p_star;
    prim.velocity[dir_v] = u_star;
    prim.density = rho_star;
    CompleteFromPrim(eq, expanded, prim);
  } else {
    throw std::runtime_error("No alternative is currently implemented.");
  }
}

} // namespace fub

#endif // !FUB_AMREX_CUTCELL_TURBINE_MASS_FLOW_BOUNDARY_HPP
