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

#ifndef FUB_CUTCELL_METHOD_KBN_STABILISATION_HPP
#define FUB_CUTCELL_METHOD_KBN_STABILISATION_HPP

#include "fub/CutCellData.hpp"
#include "fub/Duration.hpp"
#include "fub/Equation.hpp"
#include "fub/ExactRiemannSolver.hpp"
#include "fub/ForEach.hpp"
#include "fub/State.hpp"
#include "fub/StateRow.hpp"
#include "fub/core/span.hpp"
#include "fub/core/tuple.hpp"

#include <algorithm>

namespace fub {

void ComputeStableFluxComponents(
    const PatchDataView<double, 2, layout_stride>& stabilised_fluxes,
    const PatchDataView<double, 2, layout_stride>& shielded_left_fluxes,
    const PatchDataView<double, 2, layout_stride>& shielded_right_fluxes,
    const PatchDataView<const double, 2, layout_stride>& regular_fluxes,
    const PatchDataView<const double, 2, layout_stride>& boundary_fluxes,
    const CutCellData<2>& geom, Duration dt, double dx, Direction dir);

void ComputeStableFluxComponents(
    const PatchDataView<double, 3, layout_stride>& stabilised_fluxes,
    const PatchDataView<double, 3, layout_stride>& shielded_left_fluxes,
    const PatchDataView<double, 3, layout_stride>& shielded_right_fluxes,
    const PatchDataView<const double, 3, layout_stride>& regular_fluxes,
    const PatchDataView<const double, 3, layout_stride>& boundary_fluxes,
    const CutCellData<3>& geom, Duration dt, double dx, Direction dir);

template <typename FM,
          typename RiemannSolver = ExactRiemannSolver<typename FM::Equation>>
class KbnCutCellMethod : public FM {
public:
  // Typedefs

  using Equation = typename FM::Equation;
  using Conservative = ::fub::Conservative<Equation>;
  using ConservativeArray = ::fub::ConservativeArray<Equation>;
  using Complete = ::fub::Complete<Equation>;
  using CompleteArray = ::fub::CompleteArray<Equation>;

  // Static Variables

  static constexpr int Rank = Equation::Rank();
  static constexpr std::size_t sRank = static_cast<std::size_t>(Rank);
  static constexpr int StencilWidth = FM::GetStencilWidth();
  static_assert(StencilWidth > 0);
  static constexpr std::size_t StencilSize =
      static_cast<std::size_t>(2 * StencilWidth);

  template <typename T> using DataView = PatchDataView<T, Rank, layout_stride>;

  // Constructors

  /// Constructs a CutCell method from a given base flux method.
  ///
  /// This constructor uses a default constructed riemann problem solver.
  explicit KbnCutCellMethod(const FM& fm);

  /// Constructs a CutCell Method from a specified base flux and riemann solver.
  KbnCutCellMethod(const FM& fm, const RiemannSolver& rs);

  /// This function computes a reference state for each cut-cell.
  ///
  /// These states will be used to compute the embedded boundary fluxes and
  /// need to be fixed for a whole split cycle.
  ///
  /// \param[out] references the destination view of references states which
  ///                        will be filled by this method.
  ///
  /// \param[in] states the source view of states which will used as a
  ///            reference state.
  void PreAdvanceHierarchy(const View<Complete>& references,
                           const View<const Complete>& states,
                           const CutCellData<Rank>& cutcell_data);

  /// This function can be used to compute a boundary flux for a given cut-cell.
  ///
  /// \param[out] flux the storage for the result.
  /// \param[in,out] state the state in the cut cell.
  /// \param[in] boundary_normal the boundary normal of the cutcell
  /// \param[in] dir the split direction for the flux
  /// \param[in] dt the current time step size
  /// \param[in] dx the cell width length
  void ComputeBoundaryFlux(
      Conservative& flux, Complete& state, const Complete& reference_state,
      const Complete& reference_mirror_state,
      const std::optional<Eigen::Vector2d>& inflow_boundary_normal,
      const Eigen::Matrix<double, Rank, 1>& boundary_normal, int mask,
      Direction dir, Duration dt, double dx);

  /// This function can be used to compute a boundary flux for all cut cells.
  void ComputeBoundaryFluxes(
      const View<Conservative>& boundary_fluxes,
      const View<const Complete>& states,
      const View<const Complete>& reference_states,
      const View<const Complete>& reference_mirror_states,
      const StridedDataView<int, Rank>& reference_masks,
      const std::optional<Eigen::Vector2d>& inflow_boundary_normal,
      const CutCellData<Rank>& cutcell_data, Duration dt, double dx,
      Direction dir);

  using FM::ComputeStableDt;

  /// \todo compute stable dt inside of cutcells, i.e. in the reflection with
  /// their boundary state.
  double ComputeStableDt(const View<const Complete>& states,
                         const CutCellData<Rank>& cutcell_data, double dx,
                         Direction dir);

  void ComputeRegularFluxes(const View<Conservative>& regular_fluxes,
                            const View<const Complete>& states,
                            const CutCellData<Rank>& cutcell_data, Duration dt,
                            double dx, Direction dir);

  void ComputeCutCellFluxes(const View<Conservative>& stabilised_fluxes,
                            const View<Conservative>& shielded_left_fluxes,
                            const View<Conservative>& shielded_right_fluxes,
                            const View<Conservative>& doubly_shielded_fluxes,
                            const View<Conservative>& boundary_fluxes,
                            const View<const Conservative>& regular_fluxes,
                            const View<const Complete>& states,
                            const CutCellData<Rank>& cutcell_data, Duration dt,
                            double dx, Direction dir);

private:
  std::array<Complete, StencilSize> stencil_{};
  Complete state_{FM::GetEquation()};
  Complete solution_{FM::GetEquation()};
  Complete reflected_{FM::GetEquation()};
  Complete reference_state_{FM::GetEquation()};
  Complete reference_mirror_state_{FM::GetEquation()};
  Conservative cutcell_flux_{FM::GetEquation()};
  Conservative regular_flux_{FM::GetEquation()};
  Conservative boundary_flux_left_{FM::GetEquation()};
  Conservative boundary_flux_right_{FM::GetEquation()};
  Conservative shielded_right_flux_{FM::GetEquation()};
  Conservative shielded_left_flux_{FM::GetEquation()};
  Conservative doubly_shielded_flux_{FM::GetEquation()};

  std::array<CompleteArray, StencilSize> stencil_array_{};
  ConservativeArray numeric_flux_array_{FM::GetEquation()};

  RiemannSolver riemann_solver_{FM::GetEquation()};
};

// IMPLEMENTATION

template <typename FM, typename RiemannSolver>
KbnCutCellMethod<FM, RiemannSolver>::KbnCutCellMethod(const FM& fm) : FM(fm) {
  stencil_.fill(Complete(FM::GetEquation()));
  stencil_array_.fill(CompleteArray(FM::GetEquation()));
}

template <typename FM, typename RiemannSolver>
KbnCutCellMethod<FM, RiemannSolver>::KbnCutCellMethod(const FM& fm,
                                                      const RiemannSolver& rs)
    : FM(fm), riemann_solver_(rs) {
  stencil_.fill(Complete(FM::GetEquation()));
  stencil_array_.fill(CompleteArray(FM::GetEquation()));
}

template <typename FM, typename RiemannSolver>
void KbnCutCellMethod<FM, RiemannSolver>::PreAdvanceHierarchy(
    const View<Complete>& references, const View<const Complete>& states,
    const CutCellData<Rank>& geom) {
  const Equation& equation = FM::GetEquation();
  const Eigen::Matrix<double, Rank, 1> unit = UnitVector<Rank>(Direction::X);
  ForEachIndex(Box<0>(states), [&](auto... is) {
    std::array<std::ptrdiff_t, sRank> cell{is...};
    if (IsCutCell(geom, cell)) {
      Load(state_, states, cell);
      if (state_.density > 0.0) {
        const Eigen::Matrix<double, Rank, 1> normal =
            GetBoundaryNormal(geom, cell);
        Rotate(state_, state_, MakeRotation(normal, unit), equation);
        Reflect(reflected_, state_, unit, equation);
        riemann_solver_.SolveRiemannProblem(solution_, reflected_, state_,
                                            Direction::X);
        Rotate(solution_, solution_, MakeRotation(unit, normal), equation);
        Store(references, solution_, cell);
      }
    }
  });
}

template <typename FM, typename RiemannSolver>
void KbnCutCellMethod<FM, RiemannSolver>::ComputeBoundaryFlux(
    Conservative& flux, Complete& state, const Complete& reference_state,
    const Complete& reference_mirror_state,
    const std::optional<Eigen::Vector2d>& inflow_boundary_normal,
    const Eigen::Matrix<double, Rank, 1>& boundary_normal, int mask,
    Direction dir, Duration /* dt */, double /* dx */) {
  const Equation& equation = FM::GetEquation();
  equation.Flux(flux, reference_state, dir);
  const Eigen::Matrix<double, Rank, 1> unit = UnitVector<Rank>(Direction::X);

  if (mask) {
    if (inflow_boundary_normal) {
      Rotate(state, state, MakeRotation(*inflow_boundary_normal, unit),
             equation);
      riemann_solver_.SolveRiemannProblem(solution_, reference_mirror_state,
                                          state, Direction::X);
      Rotate(solution_, solution_, MakeRotation(unit, *inflow_boundary_normal),
             equation);
    } else {
      Rotate(state, state, MakeRotation(boundary_normal, unit), equation);
      riemann_solver_.SolveRiemannProblem(solution_, reference_mirror_state,
                                          state, Direction::X);
      Rotate(solution_, solution_, MakeRotation(unit, boundary_normal),
             equation);
    }
  } else {
    // Rotate states such that the boundary is left and the state is right
    // Reflect state in the split direction
    Rotate(state, state, MakeRotation(boundary_normal, unit), equation);
    Reflect(reflected_, state, unit, equation);
    riemann_solver_.SolveRiemannProblem(solution_, reflected_, state,
                                        Direction::X);
    Rotate(solution_, solution_, MakeRotation(unit, boundary_normal), equation);
  }
  const int d = static_cast<int>(dir);
  const double u_advective =
      reference_state.momentum[d] / reference_state.density;
  const double u_solution = solution_.momentum[d] / solution_.density;
  flux.momentum[d] += solution_.pressure - reference_state.pressure;
  flux.energy +=
      u_solution * solution_.pressure - u_advective * reference_state.pressure;
}

template <typename FM, typename RiemannSolver>
void KbnCutCellMethod<FM, RiemannSolver>::ComputeBoundaryFluxes(
    const View<Conservative>& boundary_fluxes,
    const View<const Complete>& states,
    const View<const Complete>& reference_states,
    const View<const Complete>& reference_mirror_states,
    const StridedDataView<int, Rank>& masks,
    const std::optional<Eigen::Vector2d>& inflow_boundary_normal,
    const CutCellData<Rank>& geom, Duration dt, double dx, Direction dir) {
  FUB_ASSERT(Extents<0>(boundary_fluxes) == Extents<0>(states));
  ForEachIndex(Box<0>(reference_states), [&](auto... is) {
    if (IsCutCell(geom, {is...})) {
      // Get the state and the boundary normal in this cell.
      std::array<std::ptrdiff_t, sRank> cell{is...};
      Load(state_, states, cell);
      Load(reference_state_, reference_states, cell);
      Load(reference_mirror_state_, reference_mirror_states, cell);
      const Eigen::Matrix<double, Rank, 1> normal =
          GetBoundaryNormal(geom, cell);

      this->ComputeBoundaryFlux(boundary_flux_left_, state_, reference_state_,
                                reference_mirror_state_, inflow_boundary_normal,
                                normal, masks(cell), dir, dt, dx);

      // Store the result in our array
      Store(boundary_fluxes, boundary_flux_left_, cell);
    }
  });
}

/// \todo compute stable dt inside of cutcells, i.e. in the reflection with
/// their boundary state.
template <typename FM, typename RiemannSolver>
double KbnCutCellMethod<FM, RiemannSolver>::ComputeStableDt(
    const View<const Complete>& states, const CutCellData<Rank>& geom,
    double dx, Direction dir) {
  double min_dt = std::numeric_limits<double>::infinity();
  static constexpr int kWidth = FM::GetStencilWidth();
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
      Array1d dts = FM::ComputeStableDt(stencil_array_, betas, alphas, dx, dir);
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
    Array1d dts = FM::ComputeStableDt(stencil_array_, betas, alphas, dx, dir);
    min_dt = std::min(min_dt, dts.minCoeff());
  });
  return min_dt;
}

template <typename FM, typename RiemannSolver>
void KbnCutCellMethod<FM, RiemannSolver>::ComputeRegularFluxes(
    const View<Conservative>& fluxes, const View<const Complete>& states,
    const CutCellData<Rank>& cutcell_data, Duration dt, double dx,
    Direction dir) {
  IndexBox<Rank> fluxbox = Box<0>(fluxes);
  static constexpr int kWidth = FM::GetStencilWidth();
  IndexBox<Rank> cellbox = Grow(fluxbox, dir, {kWidth, kWidth - 1});
  View<const Complete> base = Subview(states, cellbox);
  using ArrayView = PatchDataView<const double, Rank, layout_stride>;
  ArrayView volumes = cutcell_data.volume_fractions.Subview(cellbox);
  const int d = static_cast<int>(dir);
  ArrayView faces = cutcell_data.face_fractions[d].Subview(fluxbox);
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
  std::tuple views =
      std::tuple_cat(std::tuple(fluxes, faces), AsTuple(stencil_volumes),
                     AsTuple(stencil_views));
  ForEachRow(views,
             [this, dt, dx, dir](const Row<Conservative>& fluxes,
                                 span<const double> faces, auto... rest) {
               ViewPointer fit = Begin(fluxes);
               ViewPointer fend = End(fluxes);
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
               int n = static_cast<int>(get<0>(fend) - get<0>(fit));
               while (n >= kDefaultChunkSize) {
                 betas = Array1d::Map(faces.data());
                 for (std::size_t i = 0; i < StencilSize; ++i) {
                   Load(stencil_array_[i], states[i]);
                   alphas[i] = Array1d::Map(volumes[i].data());
                 }
                 FM::ComputeNumericFlux(numeric_flux_array_, betas,
                                        stencil_array_, alphas, dt, dx, dir);
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
                   volumes[i] = volumes[i].subspan(kDefaultChunkSize);
                 }
                 faces = faces.subspan(kDefaultChunkSize);
                 n = static_cast<int>(get<0>(fend) - get<0>(fit));
               }
               std::copy_n(faces.data(), n, betas.data());
               std::fill_n(betas.data() + n, kDefaultChunkSize - n, 0.0);
               for (std::size_t i = 0; i < StencilSize; ++i) {
                 LoadN(stencil_array_[i], states[i], n);
                 std::copy_n(volumes[i].data(), n, alphas[i].data());
                 std::fill_n(alphas[i].data() + n, kDefaultChunkSize - n, 0.0);
               }
               FM::ComputeNumericFlux(numeric_flux_array_, betas,
                                      stencil_array_, alphas, dt, dx, dir);
               StoreN(fit, numeric_flux_array_, n);
             });
}

template <typename FM, typename RiemannSolver>
void KbnCutCellMethod<FM, RiemannSolver>::ComputeCutCellFluxes(
    const View<Conservative>& stabilised_fluxes,
    const View<Conservative>& shielded_left_fluxes,
    const View<Conservative>& shielded_right_fluxes,
    const View<Conservative>& /* doubly_shielded_fluxes */,
    const View<Conservative>& boundary_fluxes,
    const View<const Conservative>& regular_fluxes,
    const View<const Complete>& /* states */, const CutCellData<Rank>& geom,
    Duration dt, double dx, Direction dir) {
  ForEachComponent(
      [&](DataView<double> stabilised, DataView<double> shielded_left,
          DataView<double> shielded_right, DataView<const double> regular,
          DataView<const double> boundary) {
        ComputeStableFluxComponents(stabilised, shielded_left, shielded_right,
                                    regular, boundary, geom, dt, dx, dir);
      },
      stabilised_fluxes, shielded_left_fluxes, shielded_right_fluxes,
      regular_fluxes, AsConst(boundary_fluxes));
}

} // namespace fub

#endif
