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
using ArrayLimiter = std::function<Array1d(Array1d, Array1d)>;

std::array<int, 3> MakeIndices(Direction dir) {
  const int ix = static_cast<int>(dir);
  const int iy = (ix + 1) % 3;
  const int iz = (iy + 1) % 3;
  return {ix, iy, iz};
}

BasicCharacteristics ComputeAmplitudes(const BasicPrimitives& left,
                                       const BasicPrimitives& right,
                                       double rhoc, double ooc2, int ix, int iy,
                                       int iz) {
  const double dp = right.pressure - left.pressure;
  const double du = right.velocity[ix] - left.velocity[ix];
  const double drho = right.density - left.density;
  BasicCharacteristics amplitudes;
  amplitudes.minus = dp - rhoc * du;
  amplitudes.zero = drho - ooc2 * dp;
  amplitudes.plus = dp + rhoc * du;
  amplitudes.v = right.velocity[iy] - left.velocity[iy];
  amplitudes.w = right.velocity[iz] - left.velocity[iz];
  return amplitudes;
}

BasicCharacteristicsArray ComputeAmplitudes(const BasicPrimitivesArray& left,
                                            const BasicPrimitivesArray& right,
                                            Array1d rhoc, Array1d ooc2, int ix,
                                            int iy, int iz) {
  const Array1d dp = right.pressure - left.pressure;
  const Array1d du = right.velocity.row(ix) - left.velocity.row(ix);
  const Array1d drho = right.density - left.density;
  BasicCharacteristicsArray amplitudes;
  amplitudes.minus = dp - rhoc * du;
  amplitudes.zero = drho - ooc2 * dp;
  amplitudes.plus = dp + rhoc * du;
  amplitudes.v = right.velocity.row(iy) - left.velocity.row(iy);
  amplitudes.w = right.velocity.row(iz) - left.velocity.row(iz);
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
      ComputeAmplitudes(pL, pM, rhoc, ooc2, ix, iy, iz);
  BasicCharacteristics amplitudes_right =
      ComputeAmplitudes(pM, pR, rhoc, ooc2, ix, iy, iz);
  limited.minus = limiter(amplitudes_left.minus, amplitudes_right.minus);
  limited.zero = limiter(amplitudes_left.zero, amplitudes_right.zero);
  limited.plus = limiter(amplitudes_left.plus, amplitudes_right.plus);
  limited.v = limiter(amplitudes_left.v, amplitudes_right.v);
  limited.w = limiter(amplitudes_left.w, amplitudes_right.w);

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

template <int Rank>
CharacteristicsArray
ComputeSlopes(CharacteristicsArray& limited,
              span<const CompleteArray<IdealGasMix<Rank>>, 3> states, Array1d rhoc,
              Array1d ooc2, ArrayLimiter& limiter, Direction dir) {
  const auto [ix, iy, iz] = MakeIndices(dir);
  BasicPrimitivesArray pL(states[0]);
  BasicPrimitivesArray pM(states[1]);
  BasicPrimitivesArray pR(states[2]);
  BasicCharacteristicsArray amplitudes_left =
      ComputeAmplitudes(pL, pM, rhoc, ooc2, ix, iy, iz);
  BasicCharacteristicsArray amplitudes_right =
      ComputeAmplitudes(pM, pR, rhoc, ooc2, ix, iy, iz);
  limited.minus = limiter(amplitudes_left.minus, amplitudes_right.minus);
  limited.zero = limiter(amplitudes_left.zero, amplitudes_right.zero);
  limited.plus = limiter(amplitudes_left.plus, amplitudes_right.plus);
  limited.v = limiter(amplitudes_left.v, amplitudes_right.v);
  limited.w = limiter(amplitudes_left.w, amplitudes_right.w);

  const int nspec = limited.species.rows();
  for (int s = 0; s < nspec; ++s) {
    const Array1d dY_left = states[1].species.row(s) / states[1].density -
                            states[0].species.row(s) / states[0].density;
    const Array1d dY_right = states[2].species.row(s) / states[2].density -
                             states[1].species.row(s) / states[1].density;
    limited.species.row(s) = limiter(dY_left, dY_right);
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

void Recombine(PrimitivesArray& diffs, const CharacteristicsArray& q,
               Array1d rhoc, Array1d ooc2, Direction dir) {
  diffs.pressure = 0.5 * (q.minus + q.plus);
  diffs.density = ooc2 * diffs.pressure + q.zero;

  const auto [ix, iy, iz] = MakeIndices(dir);
  diffs.velocity.row(ix) = 0.5 * (q.plus - q.minus) / rhoc;
  diffs.velocity.row(iy) = q.v;
  diffs.velocity.row(iz) = q.w;
  for (int s = 0; s < q.species.rows(); ++s) {
    diffs.species.row(s) = q.species.row(s);
  }
}

void DoTimeStep(Characteristics& ampls, const Characteristics& slope,
                int is_right_side, double lambda, double u, double c) {
  // clang-format off
  int sign = 1 - 2 * is_right_side;
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

void DoTimeStep(CharacteristicsArray& ampls, const CharacteristicsArray& slope,
                int is_right_side, Array1d lambda, Array1d u, Array1d c) {
  // clang-format off
  int sign = 1 - 2 * is_right_side;
  ampls.minus = sign * 0.5 * slope.minus * (1.0 - sign * lambda * (u - c));
  ampls.zero  = sign * 0.5 * slope.zero  * (1.0 - sign * lambda *    u   );
  ampls.plus  = sign * 0.5 * slope.plus  * (1.0 - sign * lambda * (u + c));
  ampls.v     = sign * 0.5 * slope.v     * (1.0 - sign * lambda *    u   );
  ampls.w     = sign * 0.5 * slope.w     * (1.0 - sign * lambda *    u   );
  for (int s = 0; s < ampls.species.rows(); ++s) {
    ampls.species.row(s) = sign * 0.5 * slope.species.row(s) * (1.0 - sign * lambda * u);
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
void CompleteFromPrim(IdealGasMix<Rank>& eq,
                      CompleteArray<IdealGasMix<Rank>>& q,
                      const PrimitivesArray& w) {
  FlameMasterReactor& reactor = eq.GetReactor();
  reactor.SetMassFractionsArray(w.species);
  reactor.SetDensityArray(w.density);
  const Array1d M = reactor.GetMeanMolarMassArray();
  const double R = reactor.GetUniversalGasConstant();
  const Array1d Rspec = R / M;
  const Array1d T = w.pressure / (Rspec * w.density);
  reactor.SetTemperatureArray(T);
  Array<double, Rank> velocity;
  for (int d = 0; d < Rank; ++d) {
    velocity.row(d) = w.velocity.row(d);
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
void Reconstruct(CompleteArray<IdealGasMix<Rank>>& rec, PrimitivesArray& diffs,
                 CharacteristicsArray& amplitudes, CharacteristicsArray& slopes,
                 int is_right_side, IdealGasMix<Rank>& eq,
                 span<const CompleteArray<IdealGasMix<Rank>>, 3> states,
                 Array1d lambda, ArrayLimiter& limiter, Direction dir) {
  const auto& q = states[1];
  const Array1d rho = q.density;
  const Array1d u = q.momentum.row(int(dir)) / rho;
  const Array1d c = q.speed_of_sound;
  const Array1d rhoc = rho * c;
  const Array1d ooc2 = 1.0 / (c * c);
  ComputeSlopes(slopes, states, rhoc, ooc2, limiter, dir);
  DoTimeStep(amplitudes, slopes, is_right_side, lambda, u, c);
  Recombine(diffs, amplitudes, rhoc, ooc2, dir);
  diffs.density += rho;
  for (int d = 0; d < Rank; ++d) {
    diffs.velocity.row(d) += (q.momentum.row(d) / q.density);
  }
  diffs.pressure += q.pressure;
  for (int s = 0; s < q.species.size(); ++s) {
    diffs.species.row(s) += (q.species.row(s) / q.density);
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
    double dx, Direction dir) {
  const Array1d dt_over_dx_half = Array1d::Constant(0.5 * dt.count() / dx);
  int left_side = 0;
  int right_side = 1;
  Reconstruct(reconstruction_array_[0], diffs_array_, amplitudes_array_,
              slopes_array_, left_side, equation_, stencil.template first<3>(),
              dt_over_dx_half, array_limiter_, dir);
  Reconstruct(reconstruction_array_[1], diffs_array_, amplitudes_array_,
              slopes_array_, right_side, equation_, stencil.template last<3>(),
              dt_over_dx_half, array_limiter_, dir);
  hlle_.ComputeNumericFlux(flux, reconstruction_array_, dt, dx, dir);
}

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
  return hlle_.ComputeStableDt(states.template subspan<1, 2>(), face_fraction,
                               volume_fraction.template subspan<1, 2>(), dx,
                               dir);
}

template <int Rank>
Array1d MusclHancockCharacteristic<Rank>::ComputeStableDt(
    span<const CompleteArray, 4> states, double dx, Direction dir) noexcept {
  return hlle_.ComputeStableDt(states.template subspan<1, 2>(), dx, dir);
}

} // namespace fub::ideal_gas