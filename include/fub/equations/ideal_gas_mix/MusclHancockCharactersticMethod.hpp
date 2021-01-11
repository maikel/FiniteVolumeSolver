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

#ifndef FUB_EQUATIONS_IDEAL_GAS_MIX_MUSCL_HANCOCK_CHAR_HPP
#define FUB_EQUATIONS_IDEAL_GAS_MIX_MUSCL_HANCOCK_CHAR_HPP

#include "fub/equations/IdealGasMix.hpp"
#include "fub/equations/ideal_gas_mix/HlleMethod.hpp"
#include "fub/flux_method/MusclHancockMethod.hpp"

namespace fub {
namespace ideal_gas {

struct BasicCharacteristics {
  double minus{};
  double zero{};
  double plus{};
  double v{};
  double w{};
};

struct Characteristics : BasicCharacteristics {
  Array<double, -1, 1> species;

  template <int Rank>
  Characteristics(const IdealGasMix<Rank>& eq)
      : species(eq.GetReactor().GetNSpecies()) {}
};

struct BasicCharacteristicsArray {
  Array1d minus;
  Array1d zero;
  Array1d plus;
  Array1d v;
  Array1d w;
};

struct CharacteristicsArray : BasicCharacteristicsArray {
  ArrayXd species;

  template <int Rank>
  CharacteristicsArray(const IdealGasMix<Rank>& eq)
      : species(eq.GetReactor().GetNSpecies(), kDefaultChunkSize) {}
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

struct BasicPrimitivesArray {
  Array1d density{Array1d::Zero()};
  Array3d velocity{Array3d::Zero()};
  Array1d pressure{Array1d::Zero()};

  BasicPrimitivesArray(const CompleteArray<IdealGasMix<3>>& q)
      : density{q.density}, pressure{q.pressure} {
    for (int d = 0; d < 3; ++d) {
      velocity.row(d) = q.momentum.row(d) / q.density;
    }
  }

  BasicPrimitivesArray(const CompleteArray<IdealGasMix<2>>& q)
      : density{q.density}, pressure{q.pressure} {
    for (int d = 0; d < 2; ++d) {
      velocity.row(d) = q.momentum.row(d) / q.density;
    }
  }

  BasicPrimitivesArray(const CompleteArray<IdealGasMix<1>>& q)
      : density{q.density}, pressure{q.pressure} {
    for (int d = 0; d < 1; ++d) {
      velocity.row(d) = q.momentum.row(d) / q.density;
    }
  }
};

struct PrimitivesArray : BasicPrimitivesArray {
  template <int Rank>
  PrimitivesArray(const CompleteArray<IdealGasMix<Rank>>& q)
      : BasicPrimitivesArray(q), species(q.species.rows(), kDefaultChunkSize) {}

  ArrayXd species;
};

/// \ingroup FluxMethod
///
/// This is a variation of the Muscl Hancock Method where the reconstruction at
/// the half time level is based on the primitive variables (p, u, T, Y) instead
/// of on conservative variables.
template <int Rank> class MusclHancockCharacteristic {
public:
  using Equation = IdealGasMix<Rank>;
  using Complete = ::fub::Complete<Equation>;
  using Conservative = ::fub::Conservative<Equation>;
  using CompleteArray = ::fub::CompleteArray<Equation>;
  using ConservativeArray = ::fub::ConservativeArray<Equation>;

  explicit MusclHancockCharacteristic(const IdealGasMix<Rank>& equation)
      : equation_(equation), limiter_(VanLeer()), array_limiter_(VanLeer()) {}

  template <typename Limiter>
  MusclHancockCharacteristic(const IdealGasMix<Rank>& equation,
                             Limiter&& limiter)
      : equation_(equation), limiter_(std::forward<Limiter>(limiter)),
        array_limiter_(std::forward<Limiter>(limiter)) {}

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

  std::array<CompleteArray, 2> reconstruction_array_{CompleteArray(equation_),
                                                     CompleteArray(equation_)};
  PrimitivesArray diffs_array_{reconstruction_array_[0]};
  CharacteristicsArray amplitudes_array_{equation_};
  CharacteristicsArray slopes_array_{equation_};
  std::function<Array1d(Array1d, Array1d)> array_limiter_;

  using HlleMethod =
      Hll<IdealGasMix<Rank>, EinfeldtSignalVelocities<IdealGasMix<Rank>>>;
  HlleMethod hlle_{equation_};
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