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

#ifndef FUB_FLUX_METHOD_MUSCL_HANCOCK_METHOD2
#define FUB_FLUX_METHOD_MUSCL_HANCOCK_METHOD2

#include "fub/CompleteFromCons.hpp"
#include "fub/Equation.hpp"
#include "fub/core/span.hpp"

#include "fub/flux_method/FluxMethod.hpp"
#include "fub/flux_method/GodunovMethod.hpp"

#include "fub/flux_method/Gradient.hpp"
#include "fub/flux_method/Reconstruct.hpp"

namespace fub {
template <
    typename Equation, typename GradientMethodT = ConservativeGradient<Equation>,
    typename ReconstructionMethodT = ConservativeReconstruction<Equation>,
    typename BaseMethod = GodunovMethod<Equation, ExactRiemannSolver<Equation>>>
struct MusclHancock2 {
  using GradientMethod = GradientMethodT;
  using ReconstructionMethod = ReconstructionMethodT;

  using Complete = typename Equation::Complete;
  using Conservative = typename Equation::Conservative;
  using Gradient = typename GradientMethod::Gradient;

  using CompleteArray = ::fub::CompleteArray<Equation>;
  using ConservativeArray = ::fub::ConservativeArray<Equation>;
  using GradientArray = typename GradientMethod::GradientArray;

  explicit MusclHancock2(const Equation& eq) : equation_{eq} {}

  MusclHancock2(const Equation& eq, const BaseMethod& method)
      : equation_{eq}, flux_method_{method} {}

  static constexpr int GetStencilWidth() noexcept { return 2; }

  double ComputeStableDt(span<const Complete, 4> states, double dx,
                         Direction dir) noexcept {
    return flux_method_.ComputeStableDt(states.template subspan<1, 2>(), dx,
                                        dir);
  }

  Array1d ComputeStableDt(span<const CompleteArray, 4> states, double dx,
                          Direction dir) noexcept {
    return flux_method_.ComputeStableDt(states.template subspan<1, 2>(), dx,
                                        dir);
  }

  Array1d ComputeStableDt(span<const CompleteArray, 4> states,
                          Array1d face_fraction,
                          span<const Array1d, 4> volume_fraction, double dx,
                          Direction dir) {
    return flux_method_.ComputeStableDt(
        states.template subspan<1, 2>(), face_fraction,
        volume_fraction.template subspan<1, 2>(), dx, dir);
  }

  void ComputeNumericFlux(Conservative& flux, span<const Complete, 4> stencil,
                          Duration dt, double dx, Direction dir);

  void ComputeNumericFlux(Conservative& flux, span<const Complete, 2> stencil,
                          span<const Gradient, 2> gradients, Duration dt,
                          double dx, Direction dir);

  decltype(auto) ComputeNumericFlux(ConservativeArray& flux,
                     span<const CompleteArray, 4> stencil, Duration dt,
                     double dx, Direction dir);

  decltype(auto) ComputeNumericFlux(ConservativeArray& flux,
                     span<const CompleteArray, 2> stencil,
                     span<const GradientArray, 2> gradients, Duration dt,
                     double dx, Direction dir);

  void ComputeNumericFlux(ConservativeArray& flux, Array1d face_fractions,
                          span<const CompleteArray, 4> stencil,
                          span<Array1d, 4> volume_fractions, Duration dt,
                          double dx, Direction dir);

  void ComputeNumericFlux(ConservativeArray& flux, Array1d face_fractions,
                          span<const CompleteArray, 2> stencil,
                          span<const GradientArray, 2> gradient,
                          span<Array1d, 2> volume_fractions, Duration dt,
                          double dx, Direction dir);

  const Equation& GetEquation() const noexcept { return equation_; }
  Equation& GetEquation() noexcept { return equation_; }

  const BaseMethod& GetBaseMethod() const noexcept { return flux_method_; }

  const GradientMethod& GetGradientMethod() const noexcept { return gradient_method_; }
  GradientMethod& GetGradientMethod() noexcept { return gradient_method_; }

  const ReconstructionMethod& GetReconstruction() const noexcept { return reconstruction_method_; }
  ReconstructionMethod& GetReconstruction() noexcept { return reconstruction_method_; }

private:
  // These member variables control the behaviour of this method
  Equation equation_;
  GradientMethod gradient_method_{equation_};
  ReconstructionMethod reconstruction_method_{equation_};
  BaseMethod flux_method_{equation_};

  std::array<Gradient, 2> gradient_{Gradient{equation_}, Gradient{equation_}};
  std::array<Complete, 2> reconstruction_{Complete{equation_},
                                          Complete{equation_}};

  std::array<GradientArray, 2> gradient_array_{GradientArray{equation_},
                                               GradientArray{equation_}};
  std::array<CompleteArray, 2> reconstruction_array_{CompleteArray{equation_},
                                                     CompleteArray{equation_}};
};

template <typename Equation, typename GradientMethod,
          typename ReconstructionMethod, typename BaseMethod>
void MusclHancock2<Equation, GradientMethod, ReconstructionMethod, BaseMethod>::
    ComputeNumericFlux(Conservative& flux, span<const Complete, 4> stencil,
                       Duration dt, double dx, Direction dir) {
  gradient_method_.ComputeGradient(gradient_[0],
                                   stencil.template subspan<0, 3>(), dx, dir);
  gradient_method_.ComputeGradient(gradient_[1],
                                   stencil.template subspan<1, 3>(), dx, dir);
  ComputeNumericFlux(flux, stencil.template subspan<1, 2>(), gradient_, dt, dx,
                     dir);
}

template <typename Equation, typename GradientMethod,
          typename ReconstructionMethod, typename BaseMethod>
void MusclHancock2<Equation, GradientMethod, ReconstructionMethod, BaseMethod>::
    ComputeNumericFlux(Conservative& flux, span<const Complete, 2> stencil,
                       span<const Gradient, 2> gradients, Duration dt,
                       double dx, Direction dir) {
  reconstruction_method_.Reconstruct(reconstruction_[0], stencil[0],
                                     gradients[0], dt, dx, dir, Side::Upper);
  reconstruction_method_.Reconstruct(reconstruction_[1], stencil[1],
                                     gradients[1], dt, dx, dir, Side::Lower);
  flux_method_.ComputeNumericFlux(flux, reconstruction_, dt, dx, dir);
}

template <typename Equation, typename GradientMethod,
          typename ReconstructionMethod, typename BaseMethod>
decltype(auto)
MusclHancock2<Equation, GradientMethod, ReconstructionMethod, BaseMethod>::
    ComputeNumericFlux(ConservativeArray& flux,
                       span<const CompleteArray, 4> stencil, Duration dt,
                       double dx, Direction dir) {
  gradient_method_.ComputeGradient(gradient_array_[0],
                                   stencil.template subspan<0, 3>(), dx, dir);
  gradient_method_.ComputeGradient(gradient_array_[1],
                                   stencil.template subspan<1, 3>(), dx, dir);
  return ComputeNumericFlux(flux, stencil.template subspan<1, 2>(), gradient_array_,
                     dt, dx, dir);
}

template <typename Equation, typename GradientMethod,
          typename ReconstructionMethod, typename BaseMethod>
decltype(auto)
MusclHancock2<Equation, GradientMethod, ReconstructionMethod, BaseMethod>::
    ComputeNumericFlux(ConservativeArray& flux,
                       span<const CompleteArray, 2> stencil,
                       span<const GradientArray, 2> gradients, Duration dt,
                       double dx, Direction dir) {
  reconstruction_method_.Reconstruct(reconstruction_array_[0], stencil[0],
                                     gradients[0], dt, dx, dir, Side::Upper);
  reconstruction_method_.Reconstruct(reconstruction_array_[1], stencil[1],
                                     gradients[1], dt, dx, dir, Side::Lower);
  return flux_method_.ComputeNumericFlux(flux, reconstruction_array_, dt, dx, dir);
}

template <typename Equation, typename GradientMethod,
          typename ReconstructionMethod, typename BaseMethod>
void MusclHancock2<Equation, GradientMethod, ReconstructionMethod, BaseMethod>::
    ComputeNumericFlux(ConservativeArray& flux, Array1d face_fractions,
                       span<const CompleteArray, 4> stencil,
                       span<Array1d, 4> volume_fractions, Duration dt,
                       double dx, Direction dir) {
  gradient_method_.ComputeGradient(gradient_array_[0],
                                   stencil.template subspan<0, 3>(), dx, dir);
  gradient_method_.ComputeGradient(gradient_array_[1],
                                   stencil.template subspan<1, 3>(), dx, dir);
  MaskArray left_mask =
      (volume_fractions[0] > 0.0 && volume_fractions[1] > 0.0 &&
       volume_fractions[2] > 0.0);
  MaskArray right_mask =
      (volume_fractions[1] > 0.0 && volume_fractions[2] > 0.0 &&
       volume_fractions[3] > 0.0);
  ForEachComponent(
      [&](auto&& x, auto&& y) {
        x = left_mask.select(x, 0.0);
        y = right_mask.select(y, 0.0);
      },
      gradient_array_[0], gradient_array_[1]);
  ComputeNumericFlux(flux, face_fractions, stencil.template subspan<1, 2>(),
                     gradient_array_, volume_fractions.template subspan<1, 2>(),
                     dt, dx, dir);
}

template <typename Equation, typename GradientMethod,
          typename ReconstructionMethod, typename BaseMethod>
void MusclHancock2<Equation, GradientMethod, ReconstructionMethod, BaseMethod>::
    ComputeNumericFlux(ConservativeArray& flux, Array1d face_fractions,
                       span<const CompleteArray, 2> stencil,
                       span<const GradientArray, 2> gradients,
                       span<Array1d, 2> volume_fractions, Duration dt,
                       double dx, Direction dir) {
  reconstruction_method_.Reconstruct(reconstruction_array_[0], stencil[0],
                                     gradients[0], dt, dx, dir, Side::Upper);
  reconstruction_method_.Reconstruct(reconstruction_array_[1], stencil[1],
                                     gradients[1], dt, dx, dir, Side::Lower);
  flux_method_.ComputeNumericFlux(flux, face_fractions, reconstruction_array_,
                                  volume_fractions, dt, dx, dir);
  
  MaskArray mask = (face_fractions > 0);
  ForEachComponent(
      [&](auto&& f) { 
        f = mask.select(f, 0);
        FUB_ASSERT(!f.isNaN().any()); }, flux);
}

} // namespace fub

#endif
