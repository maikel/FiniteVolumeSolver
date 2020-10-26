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

template <int Rank> using Coordinates = Eigen::Matrix<double, Rank, 1>;

template <int Rank>
Coordinates<Rank>
ComputeReflectedCoordinates(const Coordinates<Rank>& offset,
                            const Coordinates<Rank>& boundary_normal) {
  return offset - offset.dot(boundary_normal) * boundary_normal;
}

template <int Rank>
Index<Rank> RelativeCellIndex(const Coordinates<Rank>& x,
                              const Coordinates<Rank>& dx) {
  Index<Rank> i{};
  Coordinates<Rank> rel_x = x / dx;
  for (int d = 0; d < Rank; ++d) {
    i[d] = std::floor(rel_x[d]);
  }
  return i;
}

template <int Rank>
Index<Rank + 1> EmbedIndex(const Index<Rank>& index, Direction dir) {
  return std::apply(
      [dir](auto... is) {
        return Index<Rank>{is..., static_cast<int>(dir)};
      },
      index);
}

template <int Rank> struct AuxiliaryReconstructionData {
  std::array<Index<Rank>, 3> sources;
  std::array<double, 3> start;
  std::array<double, 3> end;
  Coordinates<Rank> direction;
  Coordinates<Rank> origin;
};

template <typename State, int Rank>
void ApplyGradient(State& u, span<const State, Rank> grad,
                   const Eigen::Matrix<double, Rank, 1>& x) noexcept {
  std::apply(
      [&x, &u](auto&&... grad_x) {
        ForEachComponent(
            [&x](double& u, auto&&... grad_x) {
              Coordinates<Rank> grad{grad_x...};
              u = grad.dot(x);
            },
            u, grad_x...);
      },
      grad);
}

template <typename Equation> struct ConservativeHGridReconstruction {
  using Conservative = ::fub::Conservative<Equation>;
  using ConservativeArray = ::fub::ConservativeArray<Equation>;
  using Complete = ::fub::Complete<Equation>;
  using CompleteArray = ::fub::CompleteArray<Equation>;

  static constexpr int Rank = Equation::Rank();

  ConservativeHGridReconstruction(const Equation& equation)
      : equation_(equation) {}

  // int_0^l (u0 + gradient * x) = l ( u_0 + l/2 * gradient )
  void IntegrateCellState(Conservative& integral,
                          const View<const Conservative>& states,
                          const View<const Conservative>& gradient_x,
                          const View<const Conservative>& gradient_y,
                          const View<const Conservative>& gradient_z,
                          const CutCellData<Rank>& geom,
                          const AuxiliaryReconstructionData<Rank>& aux_data) {
    // The first index starts in the cut cell itself
    const Index<Rank> cc_index = aux_data.sources[0];

    const Coordinates<Rank> normal = GetBoundaryNormal(cc_index);
    const Coordinates<Rank> x0 = GetBoundaryCentroid(geom, cc_index);
    const int dir_v = static_cast<int>(dir);
    const Coordinates<Rank> unit_vector =
        Eigen::Matrix<double, Rank, 1>::Unit(dir_v);
    const Coordinates<Rank> direction_vector =
        ReflectOnNormal(unit_vector, normal);

    Load(gradient[0], gradient_x, cc_index);
    Load(gradient[1], gradient_y, cc_index);
    Load(gradient[2], gradient_z, cc_index);
    Coordinates<Rank> xM = GetVolumeCentroid(geom, cc_index);
    ApplyGradient(gradient_0, gradient, x0 - xM);

    Load(state_0, states, cc_index);
    state_0 += gradient_0;

    ApplyGradient(gradient_0, gradient, direction_vector);

    // The second and third indices are neighboring cells

    Index<Rank> index = aux_data.sources[1];
    Coordinates<Rank> xStart = x0 + aux_data.start[1] * direction_vector;

    Load(gradient[0], gradient_x, index);
    Load(gradient[1], gradient_y, index);
    Load(gradient[2], gradient_z, index);
    Coordinates<Rank> xM = GetVolumeCentroid(geom, index);
    ApplyGradient(gradient_1, gradient, xStart - xM);

    Load(state_1, states, index);
    state_1 += gradient_1;

    ApplyGradient(gradient_1, gradient, direction_vector);

    index = aux_data.sources[2];
    xStart = x0 + aux_data.start[2] * direction_vector;

    Load(gradient[0], gradient_x, index);
    Load(gradient[1], gradient_y, index);
    Load(gradient[2], gradient_z, index);
    Coordinates<Rank> xM = GetVolumeCentroid(geom, index);
    ApplyGradient(gradient_2, gradient, xStart - xM);

    Load(state_2, states, index);
    state_2 += gradient_2;

    ApplyGradient(gradient_1, gradient, direction_vector);

    double length_0 = aux_data.end[0] - aux_data.start[0];
    double length_1 = aux_data.end[1] - aux_data.start[1];
    double length_2 = aux_data.end[2] - aux_data.start[2];

    ForEachVariable(
        [length_0, length_1, length_2](
            auto&& integral, auto&& state_0, auto&& state_1, auto&& state_2,
            auto&& gradient_0, auto&& gradient_1, auto&& gradient_2) {
          integral = length_0 * (state_0 + 0.5 * length_0 * gradient_0) +
                     length_1 * (state_1 + 0.5 * length_1 * gradient_1) +
                     length_2 * (state_2 + 0.5 * length_2 * gradient_2);
        },
        integral, state_0, state_1, state_2, gradient_0, gradient_1,
        gradient_2);
  }

  void ReconstructSinglyShieldedStencil(
      span<Complete, 4> h_grid_singly_shielded,
      span<const Complete, 4> h_grid_embedded_boundary,
      span<const Conservative, 2> limited_slopes,
      const View<const Complete>& states, const CutCellData<Rank>& cutcell_data,
      const Index<Rank>& cell, Duration /*dt*/, double /*dx*/, Direction dir) {
    const int d = static_cast<std::size_t>(dir);
    const Index<Rank + 1> cell_d = std::apply(
        [=](auto... is) {
          return Index<Rank + 1>{is..., d};
        },
        cell);
    const Index<Rank> face_L = cell;
    const Index<Rank> face_R = Shift(face_L, dir, 1);
    const Index<Rank> face_RR = Shift(face_L, dir, 2);
    const Index<Rank> face_LL = Shift(face_L, dir, -1);
    const double betaL = cutcell_data.face_fractions[d](face_L);
    const double betaR = cutcell_data.face_fractions[d](face_R);
    const double betaRR = cutcell_data.face_fractions[d](face_RR);
    const double betaLL = cutcell_data.face_fractions[d](face_LL);
    View<const Conservative> cons = AsCons(states);
    ////////////////////////////////////////////////////////////////////////
    // Case 1: get data   FROM LEFT
    //                    ^^^^^^^^^^
    if (betaL < betaR) {
      Load(h_grid_singly_shielded[2], states, Shift(cell, dir, 1));
      if (betaL < betaRR && cell[d] + 2 < Box<0>(states).upper[d]) {
        Load(h_grid_singly_shielded[3], states, Shift(cell, dir, 2));
      } else {
        h_grid_singly_shielded[3] = h_grid_singly_shielded[2];
      }
      Load(stencil[2], cons, Shift(cell, dir, 0));
      stencil[1] = h_grid_embedded_boundary[1];
      stencil[0] = h_grid_embedded_boundary[0];
      // Compute slopes, centroid is in local cell coordinates [-0.5, 0.5]
      const double d = 0.5 - cutcell_data.boundary_centeroids(cell_d);
      std::array<double, 3> xs{0.0, 1.0, 1.0 + 0.5 * (1.0 + d)};
      ComputeLimitedSlope(limited_slope_, span(stencil).template first<3>(),
                          xs);

      // Do the h-grid average for left states
      const double alpha = d;
      const double alpha_half = 0.5 * alpha;
      const double omalpha = 1.0 - alpha;
      ForEachComponent(
          [=](double& hL, double& hLL, double uCC, double uL, double uLL,
              double grad_uL, double grad_uLL) {
            hL = alpha * uCC + omalpha * (uL + alpha_half * grad_uL);
            hLL = alpha * uL + omalpha * (uLL + alpha_half * grad_uLL);
          },
          AsCons(h_grid_singly_shielded[1]), AsCons(h_grid_singly_shielded[0]),
          stencil[2], stencil[1], stencil[0], limited_slope_,
          limited_slopes[1]);
      CompleteFromCons(equation_, h_grid_singly_shielded[0],
                       h_grid_singly_shielded[0]);
      CompleteFromCons(equation_, h_grid_singly_shielded[1],
                       h_grid_singly_shielded[1]);
      ////////////////////////////////////////////////////////////////////////
      // Case 2: get data   FROM RIGHT
      //                    ^^^^^^^^^^
    } else if (betaR < betaL) {
      Load(h_grid_singly_shielded[1], states, Shift(cell, dir, -1));
      if (betaR < betaLL && Box<0>(states).lower[d] <= cell[d] - 2) {
        Load(h_grid_singly_shielded[0], states, Shift(cell, dir, -2));
      } else {
        h_grid_singly_shielded[0] = h_grid_singly_shielded[1];
      }
      Load(stencil[0], cons, Shift(cell, dir, 0));
      stencil[1] = h_grid_embedded_boundary[2];
      stencil[2] = h_grid_embedded_boundary[3];
      // Compute slopes
      const double d = 0.5 + cutcell_data.boundary_centeroids(cell_d);
      std::array<double, 3> xs{0.0, 0.5 * (1.0 + d), 1.0 + 0.5 * (1.0 + d)};
      ComputeLimitedSlope(limited_slope_, span(stencil).template first<3>(),
                          xs);

      // Do the h-grid average for right states
      const double alpha = d;
      const double alpha_half = 0.5 * alpha;
      const double omalpha = 1.0 - alpha;
      ForEachComponent(
          [=](double& hR, double& hRR, double uCC, double uR, double uRR,
              double grad_uR, double grad_uRR) {
            hR = alpha * uCC + omalpha * (uR - alpha_half * grad_uR);
            hRR = alpha * uR + omalpha * (uRR - alpha_half * grad_uRR);
          },
          AsCons(h_grid_singly_shielded[2]), AsCons(h_grid_singly_shielded[3]),
          stencil[0], stencil[1], stencil[2], limited_slope_,
          limited_slopes[0]);
      CompleteFromCons(equation_, h_grid_singly_shielded[2],
                       h_grid_singly_shielded[2]);
      CompleteFromCons(equation_, h_grid_singly_shielded[3],
                       h_grid_singly_shielded[3]);
    }
  }

  void ReconstructEmbeddedBoundaryStencil(
      span<Complete, 4> h_grid_embedded_boundary,
      span<Conservative, 2> limited_slopes, const View<const Complete>& states,
      const CutCellData<Rank>& cutcell_data, const Index<Rank>& cell,
      Duration /*dt*/, double /*dx*/, Direction dir) {
    const int dir_v = static_cast<int>(dir);
    const Index<Rank + 1> cell_d = std::apply(
        [=](auto... is) {
          return Index<Rank + 1>{is..., dir_v};
        },
        cell);
    const int d = static_cast<std::size_t>(dir);
    const Index<Rank> face_L = cell;
    const Index<Rank> face_R = Shift(face_L, dir, 1);
    const Index<Rank> face_RR = Shift(face_L, dir, 2);
    const Index<Rank> face_RRR = Shift(face_L, dir, 3);
    const Index<Rank> face_LL = Shift(face_L, dir, -1);
    const Index<Rank> face_LLL = Shift(face_L, dir, -2);
    const double betaL = cutcell_data.face_fractions[d](face_L);
    const double betaR = cutcell_data.face_fractions[d](face_R);
    const double betaRR = cutcell_data.face_fractions[d](face_RR);
    const double betaRRR = cutcell_data.face_fractions[d](face_RRR);
    const double betaLL = cutcell_data.face_fractions[d](face_LL);
    const double betaLLL = cutcell_data.face_fractions[d](face_LLL);
    View<const Conservative> cons = AsCons(states);
    ////////////////////////////////////////////////////////////////////////
    // Case 1: get data   FROM RIGHT
    //                    ^^^^^^^^^^
    if (betaL < betaR) {
      Load(stencil[0], cons, Shift(cell, dir, 0));
      Load(stencil[1], cons, Shift(cell, dir, 1));
      if (betaL < betaRR && cell[d] + 2 < Box<0>(states).upper[d]) {
        Load(stencil[2], cons, Shift(cell, dir, 2));
      } else {
        stencil[2] = stencil[1];
      }
      if (betaL < betaRR && betaL < betaRRR &&
          cell[d] + 3 < Box<0>(states).upper[d]) {
        Load(stencil[3], cons, Shift(cell, dir, 3));
      } else {
        stencil[3] = stencil[2];
      }

      // Compute slopes
      // const double dL = 0.5 - geom.centerL[face];
      // const double dR = 0.5 + geom.centerR[face];
      const double alphaL = 0.5 - cutcell_data.boundary_centeroids(cell_d);
      const double alphaR =
          betaRRR < betaL || betaRRR == 1.0
              ? 1.0
              : 0.5 +
                    std::apply(
                        [&](auto... is) {
                          return cutcell_data.boundary_centeroids(is..., dir_v);
                        },
                        Shift(cell, dir, +3));
      std::array<double, 4> xs{0.0, 0.5 * (1.0 + alphaL), 0.5 * (3.0 + alphaL),
                               0.5 * (4.0 + alphaL + alphaR)};
      ComputeLimitedSlope(limited_slopes[0], span(stencil).template first<3>(),
                          span(xs).template first<3>());
      ComputeLimitedSlope(limited_slopes[1], span(stencil).template last<3>(),
                          span(xs).template last<3>());

      // Do the h-grid average for right states
      const double alpha = alphaL;
      const double alpha_half = 0.5 * alpha;
      const double omalpha = 1.0 - alpha;
      ForEachComponent(
          [=](double& hR, double& hRR, double uCC, double uR, double uRR,
              double grad_uR, double grad_uRR) {
            hR = alpha * uCC + omalpha * (uR - alpha_half * grad_uR);
            hRR = alpha * uR + omalpha * (uRR - alpha_half * grad_uRR);
          },
          AsCons(h_grid_embedded_boundary[2]),
          AsCons(h_grid_embedded_boundary[3]), stencil[0], stencil[1],
          stencil[2], limited_slopes[0], limited_slopes[1]);
      CompleteFromCons(equation_, h_grid_embedded_boundary[2],
                       h_grid_embedded_boundary[2]);
      CompleteFromCons(equation_, h_grid_embedded_boundary[3],
                       h_grid_embedded_boundary[3]);

      // Impose embedded boundary conditions for left states
      Eigen::Array<double, Rank, 1> normal =
          GetBoundaryNormal(cutcell_data, cell);
      Reflect(h_grid_embedded_boundary[1], h_grid_embedded_boundary[2], -normal,
              equation_);
      Reflect(h_grid_embedded_boundary[0], h_grid_embedded_boundary[3], -normal,
              equation_);

      ////////////////////////////////////////////////////////////////////////
      // Case 2: get data   FROM LEFT
      //                    ^^^^^^^^^
    } else if (betaR < betaL) {
      Load(stencil[3], cons, Shift(cell, dir, 0));
      Load(stencil[2], cons, Shift(cell, dir, -1));
      if (betaR < betaLL && Box<0>(states).lower[d] <= cell[d] - 2) {
        Load(stencil[1], cons, Shift(cell, dir, -2));
      } else {
        stencil[1] = stencil[2];
      }
      if (betaR < betaLL && betaR < betaLLL &&
          Box<0>(states).lower[d] <= cell[d] - 3) {
        Load(stencil[0], cons, Shift(cell, dir, -3));
      } else {
        stencil[0] = stencil[1];
      }

      // Compute slopes
      // const double dL = 0.5 - geom.centerL[face];
      // const double dR = 0.5 + geom.centerR[face];
      const double alphaR = 0.5 + cutcell_data.boundary_centeroids(cell_d);
      const double alphaL =
          betaLLL < betaR || betaLLL == 1.0
              ? 1.0
              : 0.5 -
                    std::apply(
                        [&](auto... is) {
                          return cutcell_data.boundary_centeroids(is..., dir_v);
                        },
                        Shift(cell, dir, -3));
      std::array<double, 4> xs{0.0, 0.5 * (1.0 + alphaL), 0.5 * (3.0 + alphaL),
                               0.5 * (4.0 + alphaL + alphaR)};
      ComputeLimitedSlope(limited_slopes[0], span(stencil).template first<3>(),
                          span(xs).template first<3>());
      ComputeLimitedSlope(limited_slopes[1], span(stencil).template last<3>(),
                          span(xs).template last<3>());

      // Do the h-grid average for left states
      const double alpha = alphaR;
      const double alpha_half = 0.5 * alpha;
      const double omalpha = 1.0 - alpha;
      ForEachComponent(
          [=](double& hL, double& hLL, double uCC, double uL, double uLL,
              double grad_uL, double grad_uLL) {
            hL = alpha * uCC + omalpha * (uL + alpha_half * grad_uL);
            hLL = alpha * uL + omalpha * (uLL + alpha_half * grad_uLL);
          },
          AsCons(h_grid_embedded_boundary[1]),
          AsCons(h_grid_embedded_boundary[0]), stencil[3], stencil[2],
          stencil[1], limited_slopes[1], limited_slopes[0]);
      CompleteFromCons(equation_, h_grid_embedded_boundary[0],
                       h_grid_embedded_boundary[0]);
      CompleteFromCons(equation_, h_grid_embedded_boundary[1],
                       h_grid_embedded_boundary[1]);

      // Impose embedded boundary conditions for left states
      Eigen::Array<double, Rank, 1> normal =
          GetBoundaryNormal(cutcell_data, cell);
      Reflect(h_grid_embedded_boundary[2], h_grid_embedded_boundary[1], -normal,
              equation_);
      Reflect(h_grid_embedded_boundary[3], h_grid_embedded_boundary[0], -normal,
              equation_);

      ////////////////////////////////////////////////////////////////////////
      // Case 3:     NO RECONSTRUCTION    needed
      //             ^^^^^^^^^^^^^^^^^
    } else {
    }
  }

  void ComputeLimitedSlope(Conservative& slope, span<const Conservative, 3> u,
                           span<double, 3> x) const noexcept {
    const double hminus = x[1] - x[0];
    const double hplus = x[2] - x[1];
    ForEachComponent(
        [=](double& slope, double uL, double uM, double uR) {
          const double sL = (uM - uL) / hminus;
          const double sR = (uR - uM) / hplus;
          if (sL * sR <= 0.0) {
            slope = 0.0;
          } else {
            // slope = 0.0;
            slope = std::min(sL, sR);
          }
        },
        slope, u[0], u[1], u[2]);
  }

  Equation equation_;
  std::array<Conservative, 4> stencil{
      Conservative(equation_), Conservative(equation_), Conservative(equation_),
      Conservative(equation_)};
  Conservative limited_slope_{equation_};
};

template <typename EquationT, typename FluxMethod,
          typename HGridReconstruction =
              ConservativeHGridReconstruction<EquationT>>
class MyCutCellMethod : public FluxMethod {
public:
  using Equation = EquationT;

  // Typedefs
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
  explicit MyCutCellMethod(const Equation& equation);

  using FluxMethod::ComputeStableDt;

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
                            const CutCellData<Rank>& geom, Duration dt,
                            double dx, Direction dir);

private:
  Equation equation_;
  HGridReconstruction h_grid_reconstruction_;

  std::array<Complete, 4> h_grid_eb_{};
  std::array<Complete, 4> h_grid_singly_shielded_{};
  std::array<Conservative, 2> limited_slopes_{};
  Conservative boundary_flux_{equation_};
  Conservative shielded_right_flux_{equation_};
  Conservative shielded_left_flux_{equation_};

  std::array<CompleteArray, StencilSize> stencil_array_{};
  ConservativeArray numeric_flux_array_{equation_};
};

// IMPLEMENTATION

template <typename Equation, typename FluxMethod, typename HGridReconstruction>
MyCutCellMethod<Equation, FluxMethod, HGridReconstruction>::MyCutCellMethod(
    const Equation& eq)
    : FluxMethod(eq), equation_(eq), h_grid_reconstruction_(eq) {
  h_grid_eb_.fill(Complete(equation_));
  h_grid_singly_shielded_.fill(Complete(equation_));
  limited_slopes_.fill(Complete(equation_));
  stencil_array_.fill(CompleteArray(equation_));
}

/// \todo compute stable dt inside of cutcells, i.e. in the reflection with
/// their boundary state.
template <typename Equation, typename FluxMethod, typename HGridReconstruction>
double
MyCutCellMethod<Equation, FluxMethod, HGridReconstruction>::ComputeStableDt(
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

template <typename Equation, typename FluxMethod, typename HGridReconstruction>
void MyCutCellMethod<Equation, FluxMethod, HGridReconstruction>::
    ComputeRegularFluxes(const View<Conservative>& fluxes,
                         const View<const Complete>& states,
                         const CutCellData<Rank>& cutcell_data, Duration dt,
                         double dx, Direction dir) {
  IndexBox<Rank> fluxbox = Box<0>(fluxes);
  static constexpr int kWidth = FluxMethod::GetStencilWidth();
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
  ForEachRow(views, [this, dt, dx, dir](const Row<Conservative>& fluxes,
                                        span<const double> faces,
                                        auto... rest) {
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
      FluxMethod::ComputeNumericFlux(numeric_flux_array_, betas, stencil_array_,
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
    FluxMethod::ComputeNumericFlux(numeric_flux_array_, betas, stencil_array_,
                                   alphas, dt, dx, dir);
    StoreN(fit, numeric_flux_array_, n);
  });
}

template <typename Equation, typename FluxMethod, typename HGridReconstruction>
void MyCutCellMethod<Equation, FluxMethod, HGridReconstruction>::
    ComputeCutCellFluxes(const View<Conservative>& stabilised_fluxes,
                         const View<Conservative>& shielded_left_fluxes,
                         const View<Conservative>& shielded_right_fluxes,
                         const View<Conservative>& /* doubly_shielded_fluxes */,
                         const View<Conservative>& boundary_fluxes,
                         const View<const Conservative>& regular_fluxes,
                         const View<const Complete>& states,
                         const CutCellData<Rank>& geom, Duration dt, double dx,
                         Direction dir) {

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
  ForEachIndex(Shrink(Box<0>(states), dir, {1, 1}), [&](auto... is) {
    Index<Rank> cell{is...};
    Index<Rank> faceL = cell;
    Index<Rank> faceR = Shift(faceL, dir, 1);
    const double betaL = betas(faceL);
    const double betaR = betas(faceR);
    if (betaL < betaR) {
      h_grid_reconstruction_.ReconstructEmbeddedBoundaryStencil(
          h_grid_eb_, limited_slopes_, states, geom, cell, dt, dx, dir);
      FluxMethod::ComputeNumericFlux(boundary_flux_, h_grid_eb_, dt, dx, dir);
      h_grid_reconstruction_.ReconstructSinglyShieldedStencil(
          h_grid_singly_shielded_, h_grid_eb_, limited_slopes_, states, geom,
          cell, dt, dx, dir);
      FluxMethod::ComputeNumericFlux(shielded_left_flux_,
                                     h_grid_singly_shielded_, dt, dx, dir);
      Store(boundary_fluxes, boundary_flux_, cell);
      if (Contains(Box<0>(shielded_left_fluxes), faceR)) {
        Store(shielded_left_fluxes, shielded_left_flux_, faceR);
      }
    } else if (betaR < betaL) {
      h_grid_reconstruction_.ReconstructEmbeddedBoundaryStencil(
          h_grid_eb_, limited_slopes_, states, geom, cell, dt, dx, dir);
      FluxMethod::ComputeNumericFlux(boundary_flux_, h_grid_eb_, dt, dx, dir);
      h_grid_reconstruction_.ReconstructSinglyShieldedStencil(
          h_grid_singly_shielded_, h_grid_eb_, limited_slopes_, states, geom,
          cell, dt, dx, dir);
      FluxMethod::ComputeNumericFlux(shielded_right_flux_,
                                     h_grid_singly_shielded_, dt, dx, dir);
      Store(boundary_fluxes, boundary_flux_, cell);
      if (Contains(Box<0>(shielded_right_fluxes), faceL)) {
        Store(shielded_right_fluxes, shielded_right_flux_, faceL);
      }
    }
  });

  // Do the convex combiination of regular and singly shielded fluxes
  ForEachComponent(
      [&](DataView<double> stabilised, DataView<double> shielded_left,
          DataView<double> shielded_right, DataView<const double> regular,
          DataView<const double> boundary) {
        MyStab_ComputeStableFluxComponents(stabilised, shielded_left,
                                           shielded_right, regular, boundary,
                                           geom, dt, dx, dir);
      },
      stabilised_fluxes, shielded_left_fluxes, shielded_right_fluxes,
      regular_fluxes, boundary_fluxes);
}

} // namespace fub

#endif
