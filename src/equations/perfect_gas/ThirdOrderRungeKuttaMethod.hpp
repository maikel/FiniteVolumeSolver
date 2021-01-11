// Copyright (c) 2020 Maikel Nadolski
// Copyright (c) 2020 Abbas Gholami Poshtehani
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

#include "fub/equations/perfect_gas/ThirdOrderRungeKuttaMethod.hpp"

namespace fub::perfect_gas {

namespace {

double CellToEdge_(double u_ll, double u_l, double u_r, double u_rr) {
  // This is alpha_1 = (sqrt(7) + 1) / 4
  constexpr double alpha_1 = 0.9114378277661477;
  // This is alpha_2 = (sqrt(7) - 1) / 4
  constexpr double alpha_2 = 0.4114378277661477;
  return alpha_1 * (u_l + u_r) - alpha_2 * (u_ll + u_rr);
}

Array1d CellToEdge_(Array1d u_ll, Array1d u_l, Array1d u_r, Array1d u_rr) {
  // This is alpha_1 = (sqrt(7) + 1) / 4
  constexpr double alpha_1 = 0.9114378277661477;
  // This is alpha_2 = (sqrt(7) - 1) / 4
  constexpr double alpha_2 = 0.4114378277661477;
  Array1d result = alpha_1 * (u_l + u_r) - alpha_2 * (u_ll + u_rr);
  return result;
}


template <int Rank>
Conservative<PerfectGas<Rank>>
CellToEdge(span<const Complete<PerfectGas<Rank>>, 4> stencil) {
  Conservative<PerfectGas<Rank>> result{};

  result.density = CellToEdge_(stencil[0].density, stencil[1].density,
                               stencil[2].density, stencil[3].density);
  for (int i = 0; i < Rank; ++i) {
    result.momentum[i] =
        CellToEdge_(stencil[0].momentum[i], stencil[1].momentum[i],
                    stencil[2].momentum[i], stencil[3].momentum[i]);
  }
  result.energy = CellToEdge_(stencil[0].energy, stencil[1].energy,
                              stencil[2].energy, stencil[3].energy);

  return result;
}

template <int Rank>
ConservativeArray<PerfectGas<Rank>>
CellToEdge(span<const CompleteArray<PerfectGas<Rank>>, 4> stencil) {
  ConservativeArray<PerfectGas<Rank>> result{};

  result.density = CellToEdge_(stencil[0].density, stencil[1].density,
                               stencil[2].density, stencil[3].density);
  for (int i = 0; i < Rank; ++i) {
    result.momentum.row(i) =
        CellToEdge_(stencil[0].momentum.row(i), stencil[1].momentum.row(i),
                    stencil[2].momentum.row(i), stencil[3].momentum.row(i));
  }
  result.energy = CellToEdge_(stencil[0].energy, stencil[1].energy,
                              stencil[2].energy, stencil[3].energy);

  return result;
}

span<const CompleteArray<PerfectGas<1>>, 4>
Stencil(span<const CompleteArray<PerfectGas<1>>, 12> xs, int offset) {
  return span<const CompleteArray<PerfectGas<1>>, 4>(xs.subspan(offset, 4));
}

span<const Complete<PerfectGas<1>>, 4>
Stencil(span<const Complete<PerfectGas<1>>, 12> xs, int offset) {
  return span<const Complete<PerfectGas<1>>, 4>(xs.subspan(offset, 4));
}

auto Stencil2(span<const Complete<PerfectGas<1>>, 12> xs, int offset) {
  return span<const Complete<PerfectGas<1>>, 2>(xs.subspan(offset, 2));
}

auto Stencil2(span<const CompleteArray<PerfectGas<1>>, 12> xs, int offset) {
  return span<const CompleteArray<PerfectGas<1>>, 2>(xs.subspan(offset, 2));
}

Conservative<PerfectGas<1>>
ComputeDiffusion(const PerfectGas<1>& eq,
                 span<const Complete<PerfectGas<1>>, 2> q,
                 const Complete<PerfectGas<1>>& q_edge, double dx, double eta_0,
                 double k_0, int dir_v) {
  const double u_edge = eq.Velocity(q_edge)[dir_v];
  const double T_edge = eq.Temperature(q_edge);
  const double sqrt_T_edge = std::sqrt(T_edge);
  const double eta_edge = eta_0 * sqrt_T_edge;
  const double k_edge = k_0 * sqrt_T_edge;

  const double uL = eq.Velocity(q[0])[dir_v];
  const double uR = eq.Velocity(q[1])[dir_v];
  const double du_dx = (uR - uL) / dx;

  const double TL = eq.Temperature(q[0]);
  const double TR = eq.Temperature(q[1]);
  const double dT_dx = (TR - TL) / dx;
  constexpr double c1 = 4.0 / 3.0;

  Conservative<PerfectGas<1>> diffusion{};
  diffusion.density = 0.0;
  diffusion.momentum = c1 * eta_edge * du_dx;
  diffusion.energy = c1 * eta_edge * u_edge * du_dx + k_edge * dT_dx;
  return diffusion;
}

ConservativeArray<PerfectGas<1>>
ComputeDiffusion(const PerfectGas<1>& eq,
                 span<const CompleteArray<PerfectGas<1>>, 2> q,
                 const CompleteArray<PerfectGas<1>>& q_edge, double dx, double eta_0,
                 double k_0, int dir_v) {
  const Array1d u_edge = eq.Velocity(q_edge).row(dir_v);
  const Array1d T_edge = eq.Temperature(q_edge);
  const Array1d sqrt_T_edge = T_edge.sqrt();
  const Array1d eta_edge = eta_0 * sqrt_T_edge;
  const Array1d k_edge = k_0 * sqrt_T_edge;

  const Array1d uL = eq.Velocity(q[0]).row(dir_v);
  const Array1d uR = eq.Velocity(q[1]).row(dir_v);
  const Array1d du_dx = (uR - uL) / dx;

  const Array1d TL = eq.Temperature(q[0]);
  const Array1d TR = eq.Temperature(q[1]);
  const Array1d dT_dx = (TR - TL) / dx;
  constexpr double c1 = 4.0 / 3.0;

  ConservativeArray<PerfectGas<1>> diffusion{};
  diffusion.density = Array1d::Zero();
  diffusion.momentum.row(0) = c1 * eta_edge * du_dx;
  diffusion.energy = c1 * eta_edge * u_edge * du_dx + k_edge * dT_dx;
  return diffusion;
}

} // namespace

template <int Rank>
void ThirdOrderRungeKutta<Rank>::ComputeNumericFlux(
    Conservative& flux, span<const Complete, 12> u_n, Duration dt, double dx,
    Direction dir) {
  // F^n
  //              v
  //  | | | | | | | | | | | | |       u_n
  //
  //      | | | | | | | | |           u_n13
  //
  //          | | | | |               u_n23
  const double lambda = dt.count() / dx;
  std::array<Complete, 12> u_n13;
  std::array<Conservative, 13> f_n;
  const int dir_v = static_cast<int>(dir);
  for (int i = 2; i < 11; ++i) {
    Complete u_edge;
    equation_.CompleteFromCons(u_edge, CellToEdge(Stencil(u_n, i - 2)));
    equation_.Flux(f_n[i], u_edge, dir);
    const Conservative diffusion = ComputeDiffusion(
        equation_, Stencil2(u_n, i - 1), u_edge, dx, eta_0, k_0, dir_v);
    f_n[i].density -= diffusion.density;
    f_n[i].momentum -= diffusion.momentum;
    f_n[i].energy -= diffusion.energy;
  }
  for (int i = 2; i < 10; ++i) {
    u_n13[i].density =
        u_n[i].density + lambda * (f_n[i].density - f_n[i + 1].density);
    u_n13[i].momentum =
        u_n[i].momentum + lambda * (f_n[i].momentum - f_n[i + 1].momentum);
    u_n13[i].energy =
        u_n[i].energy + lambda * (f_n[i].energy - f_n[i + 1].energy);
    equation_.CompleteFromCons(u_n13[i], u_n13[i]);
  }

  // F^n+1/3
  std::array<Complete, 12> u_n23;
  std::array<Conservative, 13> f_n13;
  for (int i = 4; i < 9; ++i) {
    Complete u_edge;
    equation_.CompleteFromCons(u_edge, CellToEdge(Stencil(u_n13, i - 2)));
    equation_.Flux(f_n13[i], u_edge, dir);
    const Conservative diffusion = ComputeDiffusion(
        equation_, Stencil2(u_n13, i - 1), u_edge, dx, eta_0, k_0, dir_v);
    f_n13[i].density -= diffusion.density;
    f_n13[i].momentum -= diffusion.momentum;
    f_n13[i].energy -= diffusion.energy;
  }
  for (int i = 4; i < 8; ++i) {
    u_n23[i].density =
        0.75 * u_n[i].density + 0.25 * u_n13[i].density +
        0.25 * lambda * (f_n13[i].density - f_n13[i + 1].density);
    u_n23[i].momentum =
        0.75 * u_n[i].momentum + 0.25 * u_n13[i].momentum +
        0.25 * lambda * (f_n13[i].momentum - f_n13[i + 1].momentum);
    u_n23[i].energy = 0.75 * u_n[i].energy + 0.25 * u_n13[i].energy +
                      0.25 * lambda * (f_n13[i].energy - f_n13[i + 1].energy);
    equation_.CompleteFromCons(u_n23[i], u_n23[i]);
  }

  // F^n+2/3
  std::array<Conservative, 13> f_n23;
  for (int i = 6; i < 7; ++i) {
    Complete u_edge;
    equation_.CompleteFromCons(u_edge, CellToEdge(Stencil(u_n23, i - 2)));
    equation_.Flux(f_n23[i], u_edge, dir);
    const Conservative diffusion = ComputeDiffusion(
        equation_, Stencil2(u_n23, i - 1), u_edge, dx, eta_0, k_0, dir_v);
    f_n23[i].density -= diffusion.density;
    f_n23[i].momentum -= diffusion.momentum;
    f_n23[i].energy -= diffusion.energy;
  }

  constexpr double c1 = 1.0 / 6.0;
  constexpr double c2 = 2.0 / 3.0;
  flux.density =
      c1 * f_n[6].density + c1 * f_n13[6].density + c2 * f_n23[6].density;
  flux.momentum =
      c1 * f_n[6].momentum + c1 * f_n13[6].momentum + c2 * f_n23[6].momentum;
  flux.energy =
      c1 * f_n[6].energy + c1 * f_n13[6].energy + c2 * f_n23[6].energy;
}

template <int Rank>
void ThirdOrderRungeKutta<Rank>::ComputeNumericFlux(
    ConservativeArray& , Array1d,
    span<const CompleteArray, 12>,
    span<const Array1d, 12>, Duration, double,
    Direction) {
  throw std::runtime_error("Not implemented yet.");
}

template <int Rank>
void ThirdOrderRungeKutta<Rank>::ComputeNumericFlux(
    ConservativeArray& flux, span<const CompleteArray, 12> q_n, Duration dt,
    double dx, Direction dir) {
  const double lambda = dt.count() / dx;
  std::array<CompleteArray, 12> q_n13;
  std::array<ConservativeArray, 13> f_n;
  const int dir_v = static_cast<int>(dir);
  for (int i = 2; i < 11; ++i) {
    CompleteArray u_edge;
    equation_.CompleteFromCons(u_edge, CellToEdge(Stencil(q_n, i - 2)));
    equation_.Flux(f_n[i], u_edge, dir);
    const ConservativeArray diffusion = ComputeDiffusion(
        equation_, Stencil2(q_n, i - 1), u_edge, dx, eta_0, k_0, dir_v);
    f_n[i].density -= diffusion.density;
    f_n[i].momentum -= diffusion.momentum;
    f_n[i].energy -= diffusion.energy;
  }
  for (int i = 2; i < 10; ++i) {
    q_n13[i].density =
        q_n[i].density + lambda * (f_n[i].density - f_n[i + 1].density);
    q_n13[i].momentum =
        q_n[i].momentum + lambda * (f_n[i].momentum - f_n[i + 1].momentum);
    q_n13[i].energy =
        q_n[i].energy + lambda * (f_n[i].energy - f_n[i + 1].energy);
    equation_.CompleteFromCons(q_n13[i], q_n13[i]);
  }

  // F^n+1/3
  std::array<CompleteArray, 12> q_n23;
  std::array<ConservativeArray, 13> f_n13;
  for (int i = 4; i < 9; ++i) {
    CompleteArray u_edge;
    equation_.CompleteFromCons(u_edge, CellToEdge(Stencil(q_n13, i - 2)));
    equation_.Flux(f_n13[i], u_edge, dir);
    const ConservativeArray diffusion = ComputeDiffusion(
        equation_, Stencil2(q_n13, i - 1), u_edge, dx, eta_0, k_0, dir_v);
    f_n13[i].density -= diffusion.density;
    f_n13[i].momentum -= diffusion.momentum;
    f_n13[i].energy -= diffusion.energy;
  }
  for (int i = 4; i < 8; ++i) {
    q_n23[i].density =
        0.75 * q_n[i].density + 0.25 * q_n13[i].density +
        0.25 * lambda * (f_n13[i].density - f_n13[i + 1].density);
    q_n23[i].momentum =
        0.75 * q_n[i].momentum + 0.25 * q_n13[i].momentum +
        0.25 * lambda * (f_n13[i].momentum - f_n13[i + 1].momentum);
    q_n23[i].energy = 0.75 * q_n[i].energy + 0.25 * q_n13[i].energy +
                      0.25 * lambda * (f_n13[i].energy - f_n13[i + 1].energy);
    equation_.CompleteFromCons(q_n23[i], q_n23[i]);
  }

  // F^n+2/3
  std::array<ConservativeArray, 13> f_n23;
  for (int i = 6; i < 7; ++i) {
    CompleteArray u_edge;
    equation_.CompleteFromCons(u_edge, CellToEdge(Stencil(q_n23, i - 2)));
    equation_.Flux(f_n23[i], u_edge, dir);
    const ConservativeArray diffusion = ComputeDiffusion(
        equation_, Stencil2(q_n23, i - 1), u_edge, dx, eta_0, k_0, dir_v);
    f_n23[i].density -= diffusion.density;
    f_n23[i].momentum -= diffusion.momentum;
    f_n23[i].energy -= diffusion.energy;
  }

  constexpr double c1 = 1.0 / 6.0;
  constexpr double c2 = 2.0 / 3.0;
  flux.density =
      c1 * f_n[6].density + c1 * f_n13[6].density + c2 * f_n23[6].density;
  flux.momentum =
      c1 * f_n[6].momentum + c1 * f_n13[6].momentum + c2 * f_n23[6].momentum;
  flux.energy =
      c1 * f_n[6].energy + c1 * f_n13[6].energy + c2 * f_n23[6].energy;
}

template <int Rank>
Array1d ThirdOrderRungeKutta<Rank>::ComputeStableDt(
    span<const CompleteArray, 12> states, double dx, Direction dir) {
  const int dir_v = static_cast<int>(dir);
  const Array1d cL = states[5].speed_of_sound;
  const Array1d cR = states[6].speed_of_sound;
  const Array1d uL = states[5].momentum.row(dir_v) / states[5].density;
  const Array1d uR = states[6].momentum.row(dir_v) / states[6].density;
  return dx / (uL.abs() + cL).max(uR.abs() + cR);
}

template <int Rank>
Array1d ThirdOrderRungeKutta<Rank>::ComputeStableDt(
    span<const CompleteArray, 12> , Array1d,
    span<const Array1d, 2>, double, Direction) {
  throw std::runtime_error("Not implemented yet.");
}

template <int Rank>
double
ThirdOrderRungeKutta<Rank>::ComputeStableDt(span<const Complete, 12> states,
                                            double dx, Direction dir) {
  const int dir_v = static_cast<int>(dir);
  const double cL = states[5].speed_of_sound;
  const double cR = states[6].speed_of_sound;
  const double uL = states[5].momentum[dir_v] / states[5].density;
  const double uR = states[6].momentum[dir_v] / states[6].density;
  return dx / (std::max(std::abs(uL) + cL, std::abs(uR) + cR));
}

} // namespace fub::perfect_gas