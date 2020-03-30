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

#ifndef FUB_INITIAL_DATA_SHOCK_MACHNUMBER_HPP
#define FUB_INITIAL_DATA_SHOCK_MACHNUMBER_HPP

#include "fub/AMReX/cutcell/initial_data/RiemannProblem.hpp"

#include <array>
#include <limits>

namespace fub {
namespace amrex {
namespace cutcell {

template <typename Eq, typename Geometry>
class ShockMachnumber : private RiemannProblem<Eq, Geometry> {
public:
  using Equation = Eq;
  using Complete = fub::Complete<Equation>;
  using Conservative = fub::Conservative<Equation>;
  static constexpr int Rank = Eq::Rank();

  static Complete ComputePreShockState(const Equation& equation,
                                       const Complete& post_shock, double M_S,
                                       const Array<double, Rank, 1>& normal) {
    const double g = equation.gamma;
    const double gp1 = g + 1.0;
    const double gm1 = g - 1.0;
    const double M_post = equation.Machnumber(post_shock);
    const double M2 = (M_post - M_S) * (M_post - M_S);

    const double rho_c = gp1 * M2 / (gm1 * M2 + 2.0);
    const double rho_post = post_shock.density;
    const double rho_pre = rho_post * rho_c;

    const double p_c = (2.0 * g * M2 - gm1) / gp1;
    const double p_post = post_shock.pressure;
    const double p_pre = p_post * p_c;

    const double t = rho_post / rho_pre;
    const double a_post = post_shock.speed_of_sound;
    const double u_S = M_S * a_post;
    const double u_post = equation.Velocity(post_shock).matrix().norm();
    const double u_pre = u_S * (1.0 - t) + u_post * t;

    return equation.CompleteFromPrim(rho_pre, u_pre * normal, p_pre);
  }

  ShockMachnumber(const Eq& eq, const Geometry& geom,
                  const Complete& post_shock, double mach_number,
                  const Array<double, Rank, 1>& normal)
      : RiemannProblem<Eq, Geometry>(
            eq, geom, ComputePreShockState(eq, post_shock, mach_number, normal),
            post_shock) {}

  using RiemannProblem<Eq, Geometry>::InitializeData;

  const RiemannProblem<Eq, Geometry>& GetRiemannProblem() const noexcept {
    return *this;
  }
};

} // namespace cutcell
} // namespace amrex
} // namespace fub

#endif
