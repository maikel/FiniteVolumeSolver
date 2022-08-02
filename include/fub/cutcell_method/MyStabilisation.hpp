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

#ifndef FUB_CUTCELL_METHOD_MY_STABILISATION_HPP
#define FUB_CUTCELL_METHOD_MY_STABILISATION_HPP

#include "fub/CutCellData.hpp"
#include "fub/Duration.hpp"
#include "fub/Equation.hpp"
#include "fub/ExactRiemannSolver.hpp"
#include "fub/ForEach.hpp"
#include "fub/State.hpp"
#include "fub/StateRow.hpp"
#include "fub/core/span.hpp"
#include "fub/core/tuple.hpp"
#include "fub/equations/RequireMassflow.hpp"

#include "fub/cutcell_method/AnyMdLimiter.hpp"
#include "fub/cutcell_method/HGridReconstruction.hpp"
#include "fub/cutcell_method/MdGradients.hpp"

#include <algorithm>

namespace fub {

void MyStab_ComputeStableFluxComponents(
    const PatchDataView<double, 2, layout_stride>& stabilised_fluxes,
    const PatchDataView<double, 2, layout_stride>& shielded_left_fluxes,
    const PatchDataView<double, 2, layout_stride>& shielded_right_fluxes,
    const PatchDataView<const double, 2, layout_stride>& regular_fluxes,
    const PatchDataView<const double, 2, layout_stride>& boundary_fluxes,
    const CutCellData<2>& geom, Duration dt, double dx, Direction dir);

void MyStab_ComputeStableFluxComponents(
    const PatchDataView<double, 3, layout_stride>& stabilised_fluxes,
    const PatchDataView<double, 3, layout_stride>& shielded_left_fluxes,
    const PatchDataView<double, 3, layout_stride>& shielded_right_fluxes,
    const PatchDataView<const double, 3, layout_stride>& regular_fluxes,
    const PatchDataView<const double, 3, layout_stride>& boundary_fluxes,
    const CutCellData<3>& geom, Duration dt, double dx, Direction dir);

template <typename EquationT, typename FluxMethod, typename HGridRec,
          typename RiemannSolver =
              ExactRiemannSolver<typename FluxMethod::Equation>>
class MyCutCellMethod : public FluxMethod {
public:
  using Equation = EquationT;

  // Typedefs
  using GradientMethod = typename FluxMethod::GradientMethod;
  using Gradient = typename FluxMethod::Gradient;
  using GradientArray = typename FluxMethod::GradientArray;
  using Conservative = ::fub::Conservative<Equation>;
  using ConservativeArray = ::fub::ConservativeArray<Equation>;
  using Complete = ::fub::Complete<Equation>;
  using CompleteArray = ::fub::CompleteArray<Equation>;

  // Static Variables

  static constexpr int Rank = Equation::Rank();
  static constexpr std::size_t sRank = static_cast<std::size_t>(Rank);
  static constexpr int StencilWidth = FluxMethod::GetStencilWidth();
  static_assert(StencilWidth > 0);
  static constexpr std::size_t StencilSize =
      static_cast<std::size_t>(2 * StencilWidth);

  template <typename T> using DataView = PatchDataView<T, Rank, layout_stride>;

  // Constructors

  /// Constructs a CutCell method from a given base flux method.
  ///
  /// This constructor uses a default constructed riemann problem solver.
  MyCutCellMethod(const Equation& equation, HGridRec rec,
                  AnyLimiter<Rank> limiter);
  MyCutCellMethod(const Equation& equation, HGridRec rec,
                  const FluxMethod& flux_method, AnyLimiter<Rank> limiter);

  using FluxMethod::ComputeStableDt;

  void ComputeBoundaryStates(const View<Complete>& references,
                             const View<const Complete>& states,
                             const span<View<const Gradient>, 3>& gradients,
                             const CutCellData<Rank>& cutcell_data,
                             const Coordinates<Rank>& dx);

  /// \todo compute stable dt inside of cutcells, i.e. in the reflection with
  /// their boundary state.
  double ComputeStableDt(const View<const Complete>& states,
                         const CutCellData<Rank>& cutcell_data, double dx,
                         Direction dir);

  void ComputeRegularFluxes(const View<Conservative>& fluxes,
                            const View<const Complete>& states,
                            const View<const Gradient>& gradient_x,
                            const View<const Gradient>& gradient_y,
                            const View<const Gradient>& gradient_z,
                            const CutCellData<Rank>& cutcell_data, Duration dt,
                            double dx, Direction dir);

  void ComputeCutCellFluxes(
      const View<Conservative>& stabilised_fluxes,
      const View<Conservative>& shielded_left_fluxes,
      const View<Conservative>& shielded_right_fluxes,
      const View<Conservative>& doubly_shielded_fluxes,
      const View<Conservative>& regular_fluxes,
      const View<Conservative>& boundary_fluxes,
      const View<const Complete>& boundary_reference_states,
      const View<const Gradient>& gradient_x,
      const View<const Gradient>& gradient_y,
      const View<const Gradient>& gradient_z,
      const View<const Complete>& states, const CutCellData<Rank>& geom,
      Duration dt, const Eigen::Matrix<double, Rank, 1>& dx, Direction dir);

  void ComputeGradients(const View<Gradient>& gradient_x,
                        const View<Gradient>& gradient_y,
                        const View<Gradient>& gradient_z,
                        const View<const Complete>& states,
                        const CutCellData<Rank>& geom,
                        const Coordinates<Rank>& dx) {
    md_gradient_method_.ComputeGradients(gradient_x, gradient_y, gradient_z, states, geom, dx);
  }

  static constexpr int GetStencilWidth() noexcept {
    return FluxMethod::GetStencilWidth() + 1;
  }

private:
  Equation equation_;
  MdGradients<GradientMethod> md_gradient_method_;
  HGridRec h_grid_reconstruction_;
  RiemannSolver riemann_solver_{equation_};

  std::array<Complete, 2> h_grid_eb_{};
  std::array<Gradient, 2> h_grid_eb_gradients_{};
  std::array<Complete, 2> h_grid_singly_shielded_{};
  std::array<Gradient, 2> h_grid_singly_shielded_gradients_{};
  std::array<Complete, 2> h_grid_regular_{};
  std::array<Gradient, 2> h_grid_regular_gradients_{};

  Conservative boundary_flux_{equation_};
  Conservative singly_shielded_flux_{equation_};
  Conservative regular_flux_{equation_};
  Conservative regular_flux_left_{equation_};
  Conservative regular_flux_right_{equation_};
  Conservative stable_flux_{equation_};

  Complete state_{equation_};
  Complete reflected_{equation_};
  Complete riemann_solution_{equation_};

  std::array<Gradient, 3> gradients_{Gradient(equation_), Gradient(equation_),
                                     Gradient(equation_)};
  Gradient gradient_{equation_};
  Gradient scratch_{equation_};
  std::array<CompleteArray, 4> stencil_array_{};
  std::array<GradientArray, 2> gradient_array_{};
  ConservativeArray numeric_flux_array_{equation_};
};

// IMPLEMENTATION

template <typename Equation, typename FluxMethod, typename HGridRec,
          typename RiemannSolver>
MyCutCellMethod<Equation, FluxMethod, HGridRec, RiemannSolver>::MyCutCellMethod(
    const Equation& eq, HGridRec hgrid_rec, AnyLimiter<Rank> limiter)
    : MyCutCellMethod(eq, std::move(hgrid_rec), FluxMethod(eq),
                      std::move(limiter)) {}

template <typename Equation, typename FluxMethod, typename HGridRec,
          typename RiemannSolver>
MyCutCellMethod<Equation, FluxMethod, HGridRec, RiemannSolver>::MyCutCellMethod(
    const Equation& eq, HGridRec hgrid_rec, const FluxMethod& flux_method,
    AnyLimiter<Rank> limiter)
    : FluxMethod(flux_method), equation_(eq),
      md_gradient_method_(FluxMethod::GetGradientMethod(), std::move(limiter)),
      h_grid_reconstruction_(std::move(hgrid_rec)) {
  h_grid_eb_.fill(Complete(equation_));
  h_grid_eb_gradients_.fill(Gradient(equation_));
  h_grid_singly_shielded_.fill(Complete(equation_));
  h_grid_singly_shielded_gradients_.fill(Gradient(equation_));
  h_grid_regular_.fill(Complete(equation_));
  h_grid_regular_gradients_.fill(Gradient(equation_));
  stencil_array_.fill(CompleteArray(equation_));
}

template <typename Equation, typename FluxMethod, typename HGridRec,
          typename RiemannSolver>
void MyCutCellMethod<Equation, FluxMethod, HGridRec, RiemannSolver>::
    ComputeBoundaryStates(const View<Complete>& references,
                          const View<const Complete>& states,
                          const span<View<const Gradient>, 3>& gradients,
                          const CutCellData<Rank>& geom,
                          const Coordinates<Rank>& h) {
  const Eigen::Matrix<double, Rank, 1> unit = UnitVector<Rank>(Direction::X);
  ForEachIndex(Box<0>(references), [&](auto... is) {
    Index<Rank> index{is...};
    if (IsCutCell(geom, index)) {
      Coordinates<Rank> xM = GetAbsoluteVolumeCentroid(geom, index, h);
      Coordinates<Rank> xB = GetAbsoluteBoundaryCentroid(geom, index, h);
      Coordinates<Rank> dx = xB - xM;
      Load(state_, states, index);
      Load(gradients_[0], gradients[0], index);
      Load(gradients_[1], gradients[1], index);
      Load(gradients_[2], gradients[2], index);
      span<const Gradient, Rank> grads{gradients_.data(), Rank};
      ApplyGradient(gradient_, grads, dx);
      StateFromComplete(equation_, scratch_, state_);
      scratch_ += gradient_;
      CompleteFromState(equation_, state_, scratch_);
      const Eigen::Matrix<double, Rank, 1> normal =
          GetBoundaryNormal(geom, index);
      Rotate(state_, state_, MakeRotation(normal, unit), equation_);
      Reflect(reflected_, state_, unit, equation_);
      riemann_solver_.SolveRiemannProblem(riemann_solution_, reflected_, state_,
                                          Direction::X);
      Rotate(riemann_solution_, riemann_solution_, MakeRotation(unit, normal),
             equation_);
      Reflect(reflected_, state_, normal, equation_);
      Store(references, riemann_solution_, index);
    }
  });
}

/// \todo compute stable dt inside of cutcells, i.e. in the reflection with
/// their boundary state.
template <typename Equation, typename FluxMethod, typename HGridRec,
          typename RiemannSolver>
double
MyCutCellMethod<Equation, FluxMethod, HGridRec, RiemannSolver>::ComputeStableDt(
    const View<const Complete>& states, const CutCellData<Rank>& geom,
    double dx, Direction dir) {
  double min_dt = std::numeric_limits<double>::infinity();
  static constexpr int kWidth = StencilWidth;
  IndexBox<Rank> cellbox = Box<0>(states);
  IndexBox<Rank> fluxbox = Shrink(cellbox, dir, {kWidth, kWidth - 1});
  View<const Complete> base = Subview(states, cellbox);
  using ArrayView = PatchDataView<const double, Rank, layout_stride>;
  ArrayView volumes = geom.volume_fractions.Subview(cellbox);
  const int d = static_cast<int>(dir);
  ArrayView faces = geom.face_fractions[d].Subview(fluxbox);
  std::array<View<const Complete>, StencilSize> stencil_views;
  std::array<ArrayView, StencilSize> stencil_volumes;
  for (std::size_t i = 0; i < StencilSize; ++i) {
    stencil_views[i] =
        Shrink(base, dir,
               {static_cast<std::ptrdiff_t>(i),
                static_cast<std::ptrdiff_t>(StencilSize - i) - 1});
    stencil_volumes[i] = volumes.Subview(
        Shrink(cellbox, dir,
               {static_cast<std::ptrdiff_t>(i),
                static_cast<std::ptrdiff_t>(StencilSize - i) - 1}));
  }
  std::tuple views = std::tuple_cat(std::tuple(faces), AsTuple(stencil_volumes),
                                    AsTuple(stencil_views));
  ForEachRow(views, [this, &min_dt, dx, dir](span<const double> faces,
                                             auto... rest) {
    std::tuple args{rest...};
    std::array<span<const double>, StencilSize> volumes =
        AsArray(Take<StencilSize>(args));
    std::array states = std::apply(
        [](const auto&... xs)
            -> std::array<ViewPointer<const Complete>, StencilSize> {
          return {Begin(xs)...};
        },
        Drop<StencilSize>(args));
    std::array<Array1d, StencilSize> alphas;
    alphas.fill(Array1d::Zero());
    Array1d betas = Array1d::Zero();
    int n = static_cast<int>(faces.size());
    while (n >= kDefaultChunkSize) {
      betas = Array1d::Map(faces.data());
      for (std::size_t i = 0; i < StencilSize; ++i) {
        Load(stencil_array_[i], states[i]);
        alphas[i] = Array1d::Map(volumes[i].data());
      }
      Array1d dts =
          FluxMethod::ComputeStableDt(stencil_array_, betas, alphas, dx, dir);
      min_dt = std::min(min_dt, dts.minCoeff());
      for (std::size_t i = 0; i < StencilSize; ++i) {
        Advance(states[i], kDefaultChunkSize);
        volumes[i] = volumes[i].subspan(kDefaultChunkSize);
      }
      faces = faces.subspan(kDefaultChunkSize);
      n = static_cast<int>(faces.size());
    }
    std::copy_n(faces.data(), n, betas.data());
    std::fill_n(betas.data() + n, kDefaultChunkSize - n, 0.0);
    for (std::size_t i = 0; i < StencilSize; ++i) {
      LoadN(stencil_array_[i], states[i], n);
      std::copy_n(volumes[i].data(), n, alphas[i].data());
      std::fill_n(alphas[i].data() + n, kDefaultChunkSize - n, 0.0);
    }
    Array1d dts =
        FluxMethod::ComputeStableDt(stencil_array_, betas, alphas, dx, dir);
    min_dt = std::min(min_dt, dts.minCoeff());
  });
  return min_dt;
}

template <typename Equation, typename FluxMethod, typename HGridRec,
          typename RiemannSolver>
void MyCutCellMethod<Equation, FluxMethod, HGridRec, RiemannSolver>::
    ComputeRegularFluxes(const View<Conservative>& fluxes,
                         const View<const Complete>& states,
                         const View<const Gradient>& gradient_x,
                         const View<const Gradient>& gradient_y,
                         const View<const Gradient>& gradient_z,
                         const CutCellData<Rank>& cutcell_data, Duration dt,
                         double dx, Direction dir) {
  IndexBox<Rank> fluxbox = Box<0>(fluxes);
  IndexBox<Rank> cellbox = Grow(fluxbox, dir, {1, 0});
  View<const Complete> base = Subview(states, cellbox);
  using ArrayView = PatchDataView<const double, Rank, layout_stride>;
  ArrayView volumes = cutcell_data.volume_fractions.Subview(cellbox);
  const int d = static_cast<int>(dir);
  std::array<View<const Gradient>, 3> grads{gradient_x, gradient_y, gradient_z};
  View<const Gradient> base_grads = Subview(grads[d], cellbox);
  ArrayView faces = cutcell_data.face_fractions[d].Subview(fluxbox);
  std::array<View<const Complete>, 2> stencil_views{};
  std::array<View<const Gradient>, 2> gradient_views{};
  std::array<ArrayView, 2> stencil_volumes{};
  for (std::size_t i = 0; i < 2; ++i) {
    stencil_views[i] = Shrink(base, dir,
                              {static_cast<std::ptrdiff_t>(i),
                               static_cast<std::ptrdiff_t>(2 - i) - 1});
    gradient_views[i] = Shrink(base_grads, dir,
                               {{static_cast<std::ptrdiff_t>(i),
                                 static_cast<std::ptrdiff_t>(2 - i) - 1}});
    stencil_volumes[i] =
        volumes.Subview(Shrink(cellbox, dir,
                               {static_cast<std::ptrdiff_t>(i),
                                static_cast<std::ptrdiff_t>(2 - i) - 1}));
  }
  std::tuple views =
      std::tuple_cat(std::tuple(fluxes, faces), AsTuple(stencil_volumes),
                     AsTuple(stencil_views), AsTuple(gradient_views));
  ForEachRow(views, [this, dt, dx, dir](const Row<Conservative>& fluxes,
                                        span<const double> faces,
                                        span<const double> volumeL,
                                        span<const double> volumeR,
                                        const Row<const Complete>& statesL,
                                        const Row<const Complete>& statesR,
                                        const Row<const Gradient>& gradsL,
                                        const Row<const Gradient>& gradsR) {
    ViewPointer fit = Begin(fluxes);
    ViewPointer fend = End(fluxes);
    std::array<span<const double>, 2> volumes{volumeL, volumeR};
    std::array states{Begin(statesL), Begin(statesR)};
    std::array grads{Begin(gradsL), Begin(gradsR)};
    std::array<Array1d, 2> alphas{Array1d::Zero(), Array1d::Zero()};
    Array1d betas = Array1d::Zero();
    int n = static_cast<int>(get<0>(fend) - get<0>(fit));
    while (n >= kDefaultChunkSize) {
      betas = Array1d::Map(faces.data());
      for (std::size_t i = 0; i < 2; ++i) {
        Load(stencil_array_[i], states[i]);
        Load(gradient_array_[i], grads[i]);
        alphas[i] = Array1d::Map(volumes[i].data());
      }
      FluxMethod::ComputeNumericFlux(
          numeric_flux_array_, betas,
          span{stencil_array_}.template subspan<0, 2>(), gradient_array_,
          alphas, dt, dx, dir);
      for (int i = 0; i < betas.size(); ++i) {
        ForEachComponent(
            [&](auto&& flux [[maybe_unused]]) {
              FUB_ASSERT(betas[i] == 0.0 ||
                         (betas[i] > 0.0 && !std::isnan(flux[i])));
            },
            numeric_flux_array_);
      }
      Store(fit, numeric_flux_array_);
      Advance(fit, kDefaultChunkSize);
      for (std::size_t i = 0; i < StencilSize; ++i) {
        Advance(states[i], kDefaultChunkSize);
        Advance(grads[i], kDefaultChunkSize);
        volumes[i] = volumes[i].subspan(kDefaultChunkSize);
      }
      faces = faces.subspan(kDefaultChunkSize);
      n = static_cast<int>(get<0>(fend) - get<0>(fit));
    }
    std::copy_n(faces.data(), n, betas.data());
    std::fill_n(betas.data() + n, kDefaultChunkSize - n, 0.0);
    for (std::size_t i = 0; i < 2; ++i) {
      LoadN(stencil_array_[i], states[i], n);
      LoadN(gradient_array_[i], grads[i], n);
      std::copy_n(volumes[i].data(), n, alphas[i].data());
      std::fill_n(alphas[i].data() + n, kDefaultChunkSize - n, 0.0);
    }
    FluxMethod::ComputeNumericFlux(
        numeric_flux_array_, betas,
        span{stencil_array_}.template subspan<0, 2>(), gradient_array_, alphas,
        dt, dx, dir);
    StoreN(fit, numeric_flux_array_, n);
  });
}

template <typename Equation, typename FluxMethod, typename HGridRec,
          typename RiemannSolver>
void MyCutCellMethod<Equation, FluxMethod, HGridRec, RiemannSolver>::
    ComputeCutCellFluxes(
        const View<Conservative>& stabilised_fluxes,
        const View<Conservative>& shielded_left_fluxes,
        const View<Conservative>& shielded_right_fluxes,
        const View<Conservative>& /* doubly_shielded_fluxes */,
        const View<Conservative>& regular_fluxes,
        const View<Conservative>& boundary_fluxes,
        const View<const Complete>& boundary_reference_states,
        const View<const Gradient>& gradient_x,
        const View<const Gradient>& gradient_y,
        const View<const Gradient>& gradient_z,
        const View<const Complete>& states, const CutCellData<Rank>& geom,
        Duration dt, const Eigen::Matrix<double, Rank, 1>& dx, Direction dir) {

  // Iterate through each EB cell and then
  //
  //   (i) Compute h grid centered at EB and keep the slopes for the reflected
  //   states
  //  (ii) Compute boundary flux using the underlying flux method and its h grid
  // (iii) Compute h grid centered at singly shielded fluxes where boundary
  // states and slopes are interpolated from step (i)
  //  (iv) Compute singly shielded fluxes using the underlying flux method and
  //  its h grid
  const auto d = static_cast<std::size_t>(dir);
  const PatchDataView<const double, Rank>& betas = geom.face_fractions[d];
  const PatchDataView<const double, Rank>& betaUs =
      geom.unshielded_fractions[d];
  auto box = Box<0>(boundary_fluxes);
  ForEachIndex(box, [&](auto... is) {
    Index<Rank> cell{is...};
    Index<Rank> faceL = cell;
    Index<Rank> faceR = Shift(faceL, dir, 1);
    const double betaL = betas(faceL);
    const double betaR = betas(faceR);
    const double betaUsL = betaUs(faceL);
    const double betaUsR = betaUs(faceR);

    if (betaUsR > 0.0 && Contains(Box<0>(regular_fluxes), faceR)) {
      h_grid_reconstruction_.ReconstructRegularStencil(
          h_grid_regular_, h_grid_regular_gradients_, states, gradient_x,
          gradient_y, gradient_z, geom, faceR, dt, dx, dir);
      FluxMethod::ComputeNumericFlux(regular_flux_, h_grid_regular_,
                                     h_grid_regular_gradients_, dt, dx[d], dir);
      Store(regular_fluxes, regular_flux_, faceR);
    }

    if (betaUsL > 0.0 && Contains(Box<0>(regular_fluxes), faceL)) {
      h_grid_reconstruction_.ReconstructRegularStencil(
          h_grid_regular_, h_grid_regular_gradients_, states, gradient_x,
          gradient_y, gradient_z, geom, faceL, dt, dx, dir);
      FluxMethod::ComputeNumericFlux(regular_flux_, h_grid_regular_,
                                     h_grid_regular_gradients_, dt, dx[d], dir);
      Store(regular_fluxes, regular_flux_, faceL);
    }

    if (betaL == betaR) {
      return;
    }

    h_grid_reconstruction_.ReconstructEmbeddedBoundaryStencil(
        h_grid_eb_, h_grid_eb_gradients_, states, gradient_x, gradient_y,
        gradient_z, boundary_reference_states, geom, cell, dt, dx, dir);
    FluxMethod::ComputeNumericFlux(boundary_flux_, h_grid_eb_,
                                   h_grid_eb_gradients_, dt, dx[d], dir);
    Store(boundary_fluxes, boundary_flux_, cell);

    h_grid_reconstruction_.ReconstructSinglyShieldedStencil(
        h_grid_singly_shielded_, h_grid_singly_shielded_gradients_, h_grid_eb_,
        h_grid_eb_gradients_, states, gradient_x, gradient_y, gradient_z, geom,
        cell, dt, dx, dir);
    FluxMethod::ComputeNumericFlux(
        singly_shielded_flux_, h_grid_singly_shielded_,
        h_grid_singly_shielded_gradients_, dt, dx[d], dir);

    if (betaL < betaR) {
      if (Contains(Box<0>(shielded_left_fluxes), faceR)) {
        Store(shielded_left_fluxes, singly_shielded_flux_, faceR);
      }
    } else if (betaR < betaL) {
      if (Contains(Box<0>(shielded_right_fluxes), faceL)) {
        Store(shielded_right_fluxes, singly_shielded_flux_, faceL);
      }
    }
  });

  // Do the convex combination of regular and singly shielded fluxes
  ForEachComponent(
      [&](DataView<double> stabilised, DataView<double> shielded_left,
          DataView<double> shielded_right, DataView<const double> regular,
          DataView<const double> boundary) {
        MyStab_ComputeStableFluxComponents(stabilised, shielded_left,
                                           shielded_right, regular, boundary,
                                           geom, dt, dx[d], dir);
      },
      stabilised_fluxes, shielded_left_fluxes, shielded_right_fluxes,
      regular_fluxes, boundary_fluxes);
}

} // namespace fub

#endif
