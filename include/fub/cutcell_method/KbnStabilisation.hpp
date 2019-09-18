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
  void
  ComputeBoundaryFlux(Conservative& flux, Complete& state,
                      const Complete& reference_state,
                      const Eigen::Matrix<double, Rank, 1>& boundary_normal,
                      Direction dir, Duration dt, double dx);

  /// This function can be used to compute a boundary flux for all cut cells.
  void ComputeBoundaryFluxes(const View<Conservative>& boundary_fluxes,
                             const View<const Complete>& states,
                             const View<const Complete>& reference_states,
                             const CutCellData<Rank>& cutcell_data, Duration dt,
                             double dx, Direction dir);

  void
  ReflectCoveredStates(const std::array<int, StencilSize>& is_covered,
                       const std::array<std::ptrdiff_t, sRank> leftmost_cell,
                       const CutCellData<Rank>& cutcell_data, Direction dir);

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
                            const View<const Conservative>& regular_fluxes,

                            const View<const Conservative>& boundary_fluxes,
                            const View<const Complete>& states,
                            const CutCellData<Rank>& cutcell_data, Duration dt,
                            double dx, Direction dir);

private:
  std::array<Complete, StencilSize> stencil_{};
  Complete state_{FM::GetEquation()};
  Complete solution_{FM::GetEquation()};
  Complete reflected_{FM::GetEquation()};
  Complete reference_state_{FM::GetEquation()};
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
    const Eigen::Matrix<double, Rank, 1>& boundary_normal, Direction dir,
    Duration /* dt */, double /* dx */) {
  const Equation& equation = FM::GetEquation();
  const Eigen::Matrix<double, Rank, 1> unit = UnitVector<Rank>(Direction::X);

  // Rotate states such that the boundary is left and the state is right
  // Reflect state in the split direction
  Rotate(state, state, MakeRotation(boundary_normal, unit), equation);
  Reflect(reflected_, state, unit, equation);
  riemann_solver_.SolveRiemannProblem(solution_, reflected_, state,
                                      Direction::X);
  Rotate(solution_, solution_, MakeRotation(unit, boundary_normal), equation);

  const int d = static_cast<int>(dir);
  const double u_advective =
      reference_state.momentum[d] / reference_state.density;
  const double u_solution = solution_.momentum[d] / solution_.density;

  equation.Flux(flux, reference_state, dir);
  flux.momentum[d] += solution_.pressure - reference_state.pressure;
  flux.energy +=
      u_solution * solution_.pressure - u_advective * reference_state.pressure;
}

template <typename FM, typename RiemannSolver>
void KbnCutCellMethod<FM, RiemannSolver>::ComputeBoundaryFluxes(
    const View<Conservative>& boundary_fluxes,
    const View<const Complete>& states,
    const View<const Complete>& reference_states, const CutCellData<Rank>& geom,
    Duration dt, double dx, Direction dir) {
  FUB_ASSERT(Extents<0>(boundary_fluxes) == Extents<0>(states));
  ForEachIndex(Box<0>(reference_states), [&](auto... is) {
    if (IsCutCell(geom, {is...})) {
      // Get the state and the boundary normal in this cell.
      std::array<std::ptrdiff_t, sRank> cell{is...};
      Load(state_, states, cell);
      Load(reference_state_, reference_states, cell);
      const Eigen::Matrix<double, Rank, 1> normal =
          GetBoundaryNormal(geom, cell);

      this->ComputeBoundaryFlux(boundary_flux_left_, state_, reference_state_,
                                normal, dir, dt, dx);

      // Store the result in our array
      Store(boundary_fluxes, boundary_flux_left_, cell);
    }
  });
}

template <typename FM, typename RiemannSolver>
void KbnCutCellMethod<FM, RiemannSolver>::ReflectCoveredStates(
    const std::array<int, StencilSize>& is_covered,
    const std::array<std::ptrdiff_t, sRank> /* leftmost_cell */,
    const CutCellData<Rank>& /* cutcell_data */, Direction dir) {
  const Equation& equation = FM::GetEquation();
  if constexpr (StencilWidth == 1) {
    if (is_covered[0] == 2) {
      const Eigen::Matrix<double, Rank, 1> normal = UnitVector<Rank>(dir);
      Reflect(stencil_[0], stencil_[1], normal, equation);
    } else if (is_covered[1] == 2) {
      const Eigen::Matrix<double, Rank, 1> normal = -UnitVector<Rank>(dir);
      Reflect(stencil_[0], stencil_[1], normal, equation);
    }
  } else if constexpr (StencilWidth == 2) {
    // Do cases for | X | R | X | X |
    if (is_covered[1] == 0) {
      // Check for
      // | B | R | X | X | or
      // | R | R | X | X |
      if (is_covered[0] == 2) {
        const Eigen::Matrix<double, Rank, 1> normal = UnitVector<Rank>(dir);
        Reflect(stencil_[0], stencil_[1], normal, equation);
      }
      // Check for
      // | X | R | B | (R|SV) | or
      // | X | R | B | B |      or
      // | X | R | SV | B |
      // | X | R| R | B |
      if (is_covered[2] == 2) {
        const Eigen::Matrix<double, Rank, 1> normal = -UnitVector<Rank>(dir);
        Reflect(stencil_[2], stencil_[1], normal, equation);
        if (is_covered[3] == 2) {
          Reflect(stencil_[3], stencil_[0], normal, equation);
        }
      } else if (is_covered[2] == 1 && is_covered[3] == 2) {
        //          const Eigen::Matrix<double, Rank, 1> normal =
        //          GetBoundaryNormal(cutcell_data, Shift(leftmost_cell, dir,
        //          2)); Reflect(stencil_[3], stencil_[2], normal, equation);
        stencil_[3] = stencil_[2];
      } else if (is_covered[2] == 0 && is_covered[3] == 2) {
        const Eigen::Matrix<double, Rank, 1> normal = -UnitVector<Rank>(dir);
        Reflect(stencil_[3], stencil_[2], normal, equation);
      }
    } else {
      // Do cases for | X | X | R | X |
      FUB_ASSERT(is_covered[2] == 0);
      // Check for
      // | X | X | R | B | or
      // | X | X | R | (R|SV) |
      if (is_covered[3] == 2) {
        const Eigen::Matrix<double, Rank, 1> normal = -UnitVector<Rank>(dir);
        Reflect(stencil_[3], stencil_[2], normal, equation);
      }

      if (is_covered[1] == 2) {
        const Eigen::Matrix<double, Rank, 1> normal = UnitVector<Rank>(dir);
        Reflect(stencil_[1], stencil_[2], normal, equation);
        if (is_covered[0] == 2) {
          Reflect(stencil_[0], stencil_[3], normal, equation);
        }
      } else if (is_covered[0] == 2 && is_covered[1] == 1) {
        //          const Eigen::Matrix<double, Rank, 1> normal =
        //          GetBoundaryNormal(cutcell_data, Shift(leftmost_cell, dir,
        //          1)); Reflect(stencil_[0], stencil_[1], normal, equation);
        stencil_[0] = stencil_[1];
      }
      //        else if (is_covered[1] == 0 && is_covered[0] == 1) {
      //          const Eigen::Matrix<double, Rank, 1> normal =
      //          UnitVector<Rank>(dir); Reflect(stencil_[0], stencil_[1],
      //          normal, equation);
      //        }
    }
  } else {
    FUB_ASSERT(false);
    //      static_assert(false, "ReflectCoveredStates is not implemented for
    //      this stencil width");
  }
}

/// \todo compute stable dt inside of cutcells, i.e. in the reflection with
/// their boundary state.
template <typename FM, typename RiemannSolver>
double KbnCutCellMethod<FM, RiemannSolver>::ComputeStableDt(
    const View<const Complete>& states, const CutCellData<Rank>& geom,
    double dx, Direction dir) {
  double min_dt = std::numeric_limits<double>::infinity();
  auto&& vols = geom.volume_fractions;
  const Eigen::Matrix<double, Rank, 1> unit = UnitVector<Rank>(Direction::X);
  ForEachIndex(
      Shrink(Box<0>(states), dir, {0, StencilSize}), [&](const auto... is) {
        using Index = std::array<std::ptrdiff_t, sizeof...(is)>;
        const Index cell0{is...};
        std::array<int, StencilSize> is_covered{};
        if (IsCutCell(geom, cell0)) {
          Load(state_, states, cell0);
          const Eigen::Matrix<double, Rank, 1> normal =
              GetBoundaryNormal(geom, cell0);
          Rotate(state_, state_, MakeRotation(normal, unit), FM::GetEquation());
          Reflect(reflected_, state_, unit, FM::GetEquation());
          auto half = std::fill_n(stencil_.begin(), StencilWidth, reflected_);
          std::fill_n(half, StencilWidth, state_);
          const double dt = FM::ComputeStableDt(stencil_, dx, Direction::X);
          min_dt = std::min(dt, min_dt);
        }
        for (std::size_t i = 0; i < stencil_.size(); ++i) {
          const Index cell = Shift(cell0, dir, static_cast<int>(i));
          Load(stencil_[i], states, cell);
          if (IsCutCell(geom, cell)) {
            is_covered[i] = 1;
          } else if (vols(cell) == 1.0) {
            is_covered[i] = 0;
          } else {
            is_covered[i] = 2;
          }
        }
        const std::size_t iL = StencilWidth - 1;
        const std::size_t iR = StencilWidth;
        if (is_covered[iL] == 0 || is_covered[iR] == 0) {
          ReflectCoveredStates(is_covered, cell0, geom, dir);
          double dt = FM::ComputeStableDt(stencil_, dx, dir);
          min_dt = std::min(dt, min_dt);
        }
        if (is_covered[iL] == 1 && is_covered[iR] == 1) {
          for (std::size_t i = 0; i < iL; ++i) {
            stencil_[i] = stencil_[iL];
          }
          for (std::size_t i = iR + 1; i < StencilSize; ++i) {
            stencil_[i] = stencil_[iR];
          }
          double dt = FM::ComputeStableDt(stencil_, dx, dir);
          min_dt = std::min(dt, min_dt);
        }
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
                   ForEachComponent([&](auto&& flux) {
                     FUB_ASSERT(betas[i] == 0.0 || (betas[i] > 0.0 && !std::isnan(flux[i])));
                   }, numeric_flux_array_);
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
    const View<const Conservative>& regular_fluxes,
    const View<const Conservative>& boundary_fluxes,
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
      regular_fluxes, boundary_fluxes);
}

} // namespace fub

#endif
