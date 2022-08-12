#ifndef FUB_EQUATIONS_REQUIRE_MASSFLOW_HPP
#define FUB_EQUATIONS_REQUIRE_MASSFLOW_HPP

#include "fub/equations/EulerEquation.hpp"
#include "fub/State.hpp"
#include "fub/Direction.hpp"

#include <boost/log/utility/manipulators/add_value.hpp>
#include <boost/math/tools/roots.hpp>

#include <optional>

namespace fub {
/// \brief Compute the new complete state from the old complete state to
/// enforce the required massflow. This is done by solving the exact Riemann
/// Problem as described in Toro. TODO: cite references
struct RequireMassflow_SolveExactRiemannProblem {
  template <typename EulerEquation>
  void operator()(EulerEquation& eq, Complete<EulerEquation>& expanded,
                  const Complete<EulerEquation>& source,
                  double required_massflow, double relative_surface_area,
                  Direction dir) const;

  // Compute a left state
  template <typename EulerEquation>
  void operator()(EulerEquation& eq, Primitive<EulerEquation>& left,
                  double rhoL_star, double u_star, double p_star, double pL, Direction dir)
  {
    const double gamma = eq.gamma;
    const double gp1 = gamma + 1.0;
    const double gm1 = gamma - 1.0;
    const double oogm1 = 1.0 / gm1;
    const double gm1_over_gp1 = gm1 / gp1;
    const double p_star_over_pL = p_star / pL;
    const double rhoL_shock = rhoL_star * (gm1_over_gp1 * p_star_over_pL + 1.0) / (gm1_over_gp1 + p_star_over_pL);
    const double rhoL_rarefaction = rhoL_star / std::pow(p_star_over_pL, 1.0 / gamma);
    
    // Toro p.122 f_L for the case: left wave is shock wave
    auto f_shock = [&](double p) {
      const double A = 2.0 / (gp1 * rhoL_shock);
      const double B = pL * gm1_over_gp1;
      const double ooQ = std::sqrt(A / (p + B));
      return (p - pL) * ooQ;
    };

    // Toro p.122 eq(4.21) u_L for the case: left wave is shock wave
    auto u_shock = [&](double p) {
      return u_star + f_shock(p);
    };

    // Toro p.123 f_L for the case: left wave is rarefaction wave
    auto f_rarefaction = [&](double p) {
      const double aL = std::sqrt(gamma * pL / rhoL_rarefaction);
      const double pre_coeff = 2.0 * aL * oogm1;
      const double p_rel = p / pL;
      const double exponent = gm1 / 2.0 / gamma;
      return pre_coeff * (std::pow(p_rel, exponent) - 1.0);
    };

    // Toro p.123 eq(4.26) u_L for the case: left wave is rarefaction wave
    auto u_rarefaction = [&](double p) {
      return u_star + f_rarefaction(p);
    };
    
    const double rhoL = p_star > pL ? rhoL_shock : rhoL_rarefaction;
    const double uL = p_star > pL ? u_shock(p_star) : u_rarefaction(p_star);

    left.density = rhoL;
    left.velocity[int(dir)] = uL;
    left.pressure = pL; 
  }

  template <typename EulerEquation>
  void operator()(EulerEquation& eq, Conservative<EulerEquation>& expanded,
                  const Conservative<EulerEquation>& source,
                  double required_massflow, double relative_surface_area,
                  Direction dir) const
  {
    Complete<EulerEquation> src(eq);
    CompleteFromCons(eq, src, source);

    Complete<EulerEquation> dest(eq);
    this->operator()(eq, dest, src, required_massflow, relative_surface_area, dir);

    ForEachVariable([&](auto& exp, const auto& d) {
      exp = d;
    }, expanded, dest);
  }
};

template <typename EulerEquation>
void RequireMassflow_SolveExactRiemannProblem::
operator()(EulerEquation& eq, Complete<EulerEquation>& expanded,
           const Complete<EulerEquation>& source, double required_massflow,
           double relative_surface_area, Direction dir) const {
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

#endif