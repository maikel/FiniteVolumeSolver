// Copyright (c) 2020 Maikel Nadolski
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

#include "fub/core/function_ref.hpp"
#include "fub/equations/ideal_gas_mix/MusclHancockCharactersticMethod.hpp"

namespace fub::ideal_gas {

using Limiter = std::function<double(double, double)>;

std::array<int, 3> MakeIndices(Direction dir) {
  const int ix = static_cast<int>(dir);
  const int iy = (ix + 1) % 3;
  const int iz = (iy + 1) % 3;
  return {ix, iy, iz};
}

BasicCharacteristics ComputeAmplitudes(const BasicPrimitives& left,
                                       const BasicPrimitives& right,
                                       double rhoc, double ooc2, int ix) {
  const double dp = right.pressure - left.pressure;
  const double du = right.velocity[ix] - left.velocity[ix];
  const double drho = right.density - left.density;
  BasicCharacteristics amplitudes;
  amplitudes.minus = dp - rhoc * du;
  amplitudes.zero = drho - ooc2 * dp;
  amplitudes.plus = dp + rhoc * du;
  return amplitudes;
}

template <int Rank>
Characteristics ComputeSlopes(Characteristics& limited,
                              span<const Complete<IdealGasMix<Rank>>, 3> states,
                              double rhoc, double ooc2, Limiter& limiter,
                              Direction dir) {
  const auto [ix, iy, iz] = MakeIndices(dir);
  BasicPrimitives pL(states[0]);
  BasicPrimitives pM(states[1]);
  BasicPrimitives pR(states[2]);
  BasicCharacteristics amplitudes_left =
      ComputeAmplitudes(pL, pM, rhoc, ooc2, ix);
  BasicCharacteristics amplitudes_right =
      ComputeAmplitudes(pM, pR, rhoc, ooc2, ix);
  limited.minus = limiter(amplitudes_left.minus, amplitudes_right.minus);
  limited.zero = limiter(amplitudes_left.zero, amplitudes_right.zero);
  limited.plus = limiter(amplitudes_left.plus, amplitudes_right.plus);

  const double dv_left = pM.velocity[iy] - pL.velocity[iy];
  const double dv_right = pR.velocity[iy] - pM.velocity[iy];
  limited.v = limiter(dv_left, dv_right);

  const double dw_left = pM.velocity[iz] - pL.velocity[iz];
  const double dw_right = pR.velocity[iz] - pM.velocity[iz];
  limited.w = limiter(dw_left, dw_right);

  const int nspec = limited.species.size();
  for (int s = 0; s < nspec; ++s) {
    const double dY_left = states[1].species[s] / states[1].density -
                           states[0].species[s] / states[0].density;
    const double dY_right = states[2].species[s] / states[2].density -
                            states[1].species[s] / states[1].density;
    limited.species[s] = limiter(dY_left, dY_right);
  }
  return limited;
}

// (rhou)_t + (rho u^2 + p)_x = 0

void Recombine(Primitives& diffs, const Characteristics& q, double rhoc,
               double ooc2, Direction dir) {
  diffs.pressure = 0.5 * (q.minus + q.plus);
  diffs.density = ooc2 * diffs.pressure + q.zero;

  const auto [ix, iy, iz] = MakeIndices(dir);
  diffs.velocity[ix] = 0.5 * (q.plus - q.minus) / rhoc;
  diffs.velocity[iy] = q.v;
  diffs.velocity[iz] = q.w;
  for (int s = 0; s < q.species.size(); ++s) {
    diffs.species[s] = q.species[s];
  }
}

void DoTimeStep(Characteristics& ampls, const Characteristics& slope,
                int is_right_side, double lambda, double u, double c) {
  // clang-format off
  int sign = 1;
  if (is_right_side == 1) {
    sign = -1;
  }
  ampls.minus = sign * 0.5 * slope.minus * (1.0 - sign * lambda * (u - c));
  ampls.zero  = sign * 0.5 * slope.zero  * (1.0 - sign * lambda *    u   );
  ampls.plus  = sign * 0.5 * slope.plus  * (1.0 - sign * lambda * (u + c));
  ampls.v     = sign * 0.5 * slope.v     * (1.0 - sign * lambda *    u   );
  ampls.w     = sign * 0.5 * slope.w     * (1.0 - sign * lambda *    u   );
  for (int s = 0; s < ampls.species.size(); ++s) {
    ampls.species[s] = sign * 0.5 * slope.species[s] * (1.0 - sign * lambda * u);
  }
  // clang-format on
}

template <int Rank>
void CompleteFromPrim(IdealGasMix<Rank>& eq, Complete<IdealGasMix<Rank>>& q,
                      const Primitives& w) {
  FlameMasterReactor& reactor = eq.GetReactor();
  reactor.SetMassFractions(w.species);
  reactor.SetDensity(w.density);
  const double M = reactor.GetMeanMolarMass();
  const double R = reactor.GetUniversalGasConstant();
  const double T = w.pressure * M / w.density / R;
  reactor.SetTemperature(T);

  Eigen::Array<double, Rank, 1> velocity;
  for (int d = 0; d < Rank; ++d) {
    velocity[d] = w.velocity[d];
  }
  eq.CompleteFromReactor(q, velocity);
}

template <int Rank>
void Reconstruct(Complete<IdealGasMix<Rank>>& rec, Primitives& diffs,
                 Characteristics& amplitudes, Characteristics& slopes,
                 int is_right_side, IdealGasMix<Rank>& eq,
                 span<const Complete<IdealGasMix<Rank>>, 3> states,
                 double lambda, Limiter& limiter, Direction dir) {
  const auto& q = states[1];
  const double rho = q.density;
  const double u = q.momentum[int(dir)] / rho;
  const double c = q.speed_of_sound;
  const double rhoc = rho * c;
  const double ooc2 = 1.0 / (c * c);
  ComputeSlopes(slopes, states, rhoc, ooc2, limiter, dir);
  DoTimeStep(amplitudes, slopes, is_right_side, lambda, u, c);
  Recombine(diffs, amplitudes, rhoc, ooc2, dir);
  diffs.density += rho;
  for (int d = 0; d < Rank; ++d) {
    diffs.velocity[d] += (q.momentum[d] / q.density);
  }
  diffs.pressure += q.pressure;
  for (int s = 0; s < q.species.size(); ++s) {
    diffs.species[s] += (q.species[s] / q.density);
  }
  CompleteFromPrim(eq, rec, diffs);
}

template <int Rank>
void MusclHancockCharacteristic<Rank>::ComputeNumericFlux(
    Conservative& flux, span<const Complete, 4> stencil, Duration dt, double dx,
    Direction dir) {
  const double dt_over_dx_half = 0.5 * dt.count() / dx;
  int left_side = 0;
  int right_side = 1;
  Reconstruct(reconstruction_[0], diffs_, amplitudes_, slopes_, left_side,
              equation_, stencil.template first<3>(), dt_over_dx_half, limiter_,
              dir);
  Reconstruct(reconstruction_[1], diffs_, amplitudes_, slopes_, right_side,
              equation_, stencil.template last<3>(), dt_over_dx_half, limiter_,
              dir);
  hlle_.ComputeNumericFlux(flux, reconstruction_, dt, dx, dir);
}

template <int Rank>
void MusclHancockCharacteristic<Rank>::ComputeNumericFlux(
    ConservativeArray& flux, span<const CompleteArray, 4> stencil, Duration dt,
    double dx, Direction dir) {}

template <int Rank>
void MusclHancockCharacteristic<Rank>::ComputeNumericFlux(
    ConservativeArray& flux, Array1d face_fractions,
    span<const CompleteArray, 4> stencil,
    span<const Array1d, 4> volume_fractions, Duration dt, double dx,
    Direction dir) {}

template <int Rank>
double MusclHancockCharacteristic<Rank>::ComputeStableDt(
    span<const Complete, 4> states, double dx, Direction dir) noexcept {
  return hlle_.ComputeStableDt(states.template subspan<1, 2>(), dx, dir);
}

template <int Rank>
Array1d MusclHancockCharacteristic<Rank>::ComputeStableDt(
    span<const CompleteArray, 4> states, Array1d face_fraction,
    span<const Array1d, 4> volume_fraction, double dx, Direction dir) noexcept {
  return {};
}

template <int Rank>
Array1d MusclHancockCharacteristic<Rank>::ComputeStableDt(
    span<const CompleteArray, 4> states, double dx, Direction dir) noexcept {
  return {};
}

} // namespace fub::ideal_gas