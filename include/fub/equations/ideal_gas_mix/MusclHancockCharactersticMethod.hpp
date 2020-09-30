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

#ifndef FUB_EQUATIONS_IDEAL_GAS_MIX_MUSCL_HANCOCK_PRIM_HPP
#define FUB_EQUATIONS_IDEAL_GAS_MIX_MUSCL_HANCOCK_PRIM_HPP

#include "fub/equations/IdealGasMix.hpp"
#include "fub/equations/ideal_gas_mix/HlleMethod.hpp"

namespace fub {
namespace ideal_gas {

struct BasicCharacteristics {
  double minus{};
  double zero{};
  double plus{};
};

struct Characteristics : BasicCharacteristics {
  double v{};
  double w{};
  Array<double, -1, 1> species;

  template <int Rank>
  Characteristics(const IdealGasMix<Rank>& eq)
      : species(eq.GetReactor().GetNSpecies()) {}
};

struct BasicPrimitives {
  double density{};
  Array<double, 3, 1> velocity{};
  double pressure{};

  BasicPrimitives(const Complete<IdealGasMix<3>>& q)
      : density{q.density},
        velocity(q.momentum / q.density), pressure{q.pressure} {}

  BasicPrimitives(const Complete<IdealGasMix<2>>& q)
      : density{q.density}, velocity{q.momentum[0] / q.density,
                                     q.momentum[1] / q.density, 0},
        pressure{q.pressure} {}

  BasicPrimitives(const Complete<IdealGasMix<1>>& q)
      : density{q.density}, velocity{q.momentum[0] / q.density, 0, 0},
        pressure{q.pressure} {}
};

struct Primitives : BasicPrimitives {
  template <int Rank>
  Primitives(const Complete<IdealGasMix<Rank>>& q)
      : BasicPrimitives(q), species(q.species.size()) {}

  Array<double, -1, 1> species;
};

struct VanLeerLimiter {
  double operator()(double sL, double sR) {
    double r = 0.0;
    if (sL * sR > 0.0) {
      r = sL / sR;
    }
    if (r < 0.0) {
      return 0.0;
    } else {
      return 0.5 * std::min(2 * r / (1 + r), 2 / (1 + r)) * (sL + sR);
    }
  }
};

/// \ingroup FluxMethod
///
/// This is a variation of the Muscl Hancock Method where the reconstruction at
/// the half time level is based on the primitive variables (p, u, T, Y) instead
/// of on conservative variables.
template <int Rank>
class MusclHancockCharacteristic {
public:
  using Equation = IdealGasMix<Rank>;
  using Complete = ::fub::Complete<Equation>;
  using Conservative = ::fub::Conservative<Equation>;
  using CompleteArray = ::fub::CompleteArray<Equation>;
  using ConservativeArray = ::fub::ConservativeArray<Equation>;

  explicit MusclHancockCharacteristic(const IdealGasMix<Rank>& equation)
      : equation_(equation), limiter_(VanLeerLimiter()) {}

  [[nodiscard]] static constexpr int GetStencilWidth() noexcept { return 2; }

  /// Returns a stable time step estimate based on HLL signal velocities.
  [[nodiscard]] double ComputeStableDt(span<const Complete, 4> states,
                                       double dx, Direction dir) noexcept;

  /// Returns an array of stable time step estimates based on HLL signal
  /// velocities.
  [[nodiscard]] Array1d ComputeStableDt(span<const CompleteArray, 4> states,
                                        Array1d face_fraction,
                                        span<const Array1d, 4> volume_fraction,
                                        double dx, Direction dir) noexcept;

  [[nodiscard]] Array1d ComputeStableDt(span<const CompleteArray, 4> states,
                                        double dx, Direction dir) noexcept;

  void ComputeNumericFlux(Conservative& flux, span<const Complete, 4> stencil,
                          Duration dt, double dx, Direction dir);

  void ComputeNumericFlux(ConservativeArray& flux,
                          span<const CompleteArray, 4> stencil, Duration dt,
                          double dx, Direction dir);

  void ComputeNumericFlux(ConservativeArray& flux, Array1d face_fractions,
                          span<const CompleteArray, 4> stencil,
                          span<const Array1d, 4> volume_fractions, Duration dt,
                          double dx, Direction dir);

  [[nodiscard]] const Equation& GetEquation() const noexcept {
    return equation_;
  }
  [[nodiscard]] Equation& GetEquation() noexcept { return equation_; }

private:
  IdealGasMix<Rank> equation_;
  std::array<Complete, 2> reconstruction_{Complete(equation_),
                                          Complete(equation_)};
  Primitives diffs_{reconstruction_[0]};
  Characteristics amplitudes_{equation_};
  Characteristics slopes_{equation_};
  std::function<double(double, double)> limiter_;
  Hll<IdealGasMix<Rank>, EinfeldtSignalVelocities<IdealGasMix<Rank>>>
      hlle_{equation_, EinfeldtSignalVelocities<IdealGasMix<Rank>>()};
};

/// \ingroup FluxMethod
template <int Rank>
using MusclHancockCharMethod =
    ::fub::FluxMethod<MusclHancockCharacteristic<Rank>>;

extern template class MusclHancockCharacteristic<1>;
extern template class MusclHancockCharacteristic<2>;
extern template class MusclHancockCharacteristic<3>;
} // namespace ideal_gas

extern template class FluxMethod<ideal_gas::MusclHancockCharacteristic<1>>;
extern template class FluxMethod<ideal_gas::MusclHancockCharacteristic<2>>;
extern template class FluxMethod<ideal_gas::MusclHancockCharacteristic<3>>;

} // namespace fub

#endif