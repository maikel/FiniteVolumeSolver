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

template <int Rank>
Eigen::Matrix<double, Rank, 1>
GetAbsoluteVolumeCentroid(const CutCellData<Rank>& geom,
                          const Index<Rank>& index,
                          const Eigen::Matrix<double, Rank, 1>& dx) {
  Eigen::Matrix<double, Rank, 1> relative_xM = GetVolumeCentroid(geom, index);
  Eigen::Matrix<double, Rank, 1> offset;
  for (int i = 0; i < Rank; ++i) {
    offset[i] = static_cast<double>(index[i]);
  }
  Eigen::Matrix<double, Rank, 1> xM =
      (offset + relative_xM).array() * dx.array();
  return xM;
}

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

template <typename State, std::ptrdiff_t Rank>
void ApplyGradient(
    State& u, span<const State, Rank> grad,
    nodeduce_t<const Eigen::Matrix<double, Rank, 1>&> x) noexcept {
  if constexpr (Rank == 2) {
    ForEachComponent(
        [&x](double& u, auto&&... grad_x) {
          Coordinates<Rank> grad{grad_x...};
          u = grad.dot(x);
        },
        u, grad[0], grad[1]);
  } else if constexpr (Rank == 3) {
    ForEachComponent(
        [&x](double& u, auto&&... grad_x) {
          Coordinates<Rank> grad{grad_x...};
          u = grad.dot(x);
        },
        u, grad[0], grad[1], grad[2]);
  }
}

template <int Rank> struct BasicHGridReconstruction {
  void ComputeGradients(span<double, 2> gradient, span<const double, 4> states,
                        span<const Coordinates<Rank>, 4> x);

  void ComputeGradients(span<double, 2> gradient, span<const double, 5> states,
                        span<const Coordinates<Rank>, 5> x);

  void
  LimitGradients(const std::array<StridedDataView<double, Rank>, Rank>& grad_u,
                 StridedDataView<const double, Rank> u,
                 StridedDataView<const char, Rank> needs_limiter,
                 const CutCellData<Rank>& geom,
                 const Coordinates<Rank>& dx) const;

private:
  void LimitGradientsAtIndex(
      const std::array<StridedDataView<double, Rank>, Rank>& grad_u,
      StridedDataView<const double, Rank> u, const CutCellData<Rank>& geom,
      const Index<Rank>& index, const Coordinates<Rank>& dx) const;
};

extern template struct BasicHGridReconstruction<2>;

template <typename Equation>
struct ConservativeHGridReconstruction
    : BasicHGridReconstruction<Equation::Rank()> {
  using Base = BasicHGridReconstruction<Equation::Rank()>;

  using Conservative = ::fub::Conservative<Equation>;
  using ConservativeArray = ::fub::ConservativeArray<Equation>;
  using Complete = ::fub::Complete<Equation>;
  using CompleteArray = ::fub::CompleteArray<Equation>;

  static constexpr int Rank = Equation::Rank();

  struct AuxiliaryReconstructionData {
    static constexpr int kMaxSources = 6;

    std::array<Index<Rank>, kMaxSources> sources{};
    std::array<double, kMaxSources> start{};
    std::array<double, kMaxSources> end{};

    Coordinates<Rank> slope{};
    Coordinates<Rank> xB{};
  };

  ConservativeHGridReconstruction(const Equation& equation)
      : equation_(equation) {}

  void ComputeGradients(span<Conservative, 2> gradient,
                        span<const Conservative, 4> states,
                        span<const Coordinates<Rank>, 4> x) {
    ForEachComponent(
        [&](double& grad_x, double& grad_y, double uM, double u1, double u2,
            double u3) {
          std::array<double, 2> grad{0.0, 0.0};
          const std::array<double, 4> quantities{uM, u1, u2, u3};
          Base::ComputeGradients(grad, quantities, x);
          grad_x = grad[0];
          grad_y = grad[1];
        },
        gradient[0], gradient[1], states[0], states[1], states[2], states[3]);
  }

  void ComputeGradients(span<Conservative, 2> gradient,
                        span<const Conservative, 5> states,
                        span<const Coordinates<Rank>, 5> x) {
    ForEachComponent(
        [&](double& grad_x, double& grad_y, double uM, double u1, double u2,
            double u3, double u4) {
          std::array<double, 2> grad{0.0, 0.0};
          const std::array<double, 5> quantities{uM, u1, u2, u3, u4};
          Base::ComputeGradients(grad, quantities, x);
          grad_x = grad[0];
          grad_y = grad[1];
        },
        gradient[0], gradient[1], states[0], states[1], states[2], states[3],
        states[4]);
  }

  void ComputeGradients(const View<Conservative>& gradient_x,
                        const View<Conservative>& gradient_y,
                        const View<Conservative>& gradient_z,
                        const View<const Conservative>& states,
                        const StridedDataView<const char, Rank>& flags,
                        const CutCellData<Rank>& geom,
                        const Coordinates<Rank>& dx) {
    Conservative zero{equation_};
    if constexpr (Rank == 2) {
      std::array<Conservative, 2> gradient{equation_};
      std::array<Conservative, 5> u;
      u.fill(Conservative{equation_});
      const IndexBox<Rank> box =
          Shrink(Shrink(Box<0>(gradient_x), Direction::X, {1, 1}), Direction::Y,
                 {1, 1});
      ForEachIndex(box, [&](int i, int j) {
        //////////////////////////////////////////////////
        // All regular case
        if (geom.volume_fractions(i, j) == 1.0 &&
            geom.volume_fractions(i + 1, j) == 1.0 &&
            geom.volume_fractions(i - 1, j) == 1.0 &&
            geom.volume_fractions(i, j + 1) == 1.0 &&
            geom.volume_fractions(i, j - 1) == 1.0) {
          Load(u[0], AsCons(states), {i - 1, j});
          Load(u[1], AsCons(states), {i + 1, j});
          ForEachComponent(
              [&](double& gradient, double qL, double qR) {
                return gradient = 0.5 * (qR - qL) / dx[0];
              },
              gradient[0], u[0], u[1]);
          Load(u[2], AsCons(states), {i, j - 1});
          Load(u[3], AsCons(states), {i, j + 1});
          ForEachComponent(
              [&](double& gradient, double qL, double qR) {
                return gradient = 0.5 * (qR - qL) / dx[1];
              },
              gradient[1], u[2], u[3]);
          //////////////////////////////////////////////////
          // Cut-Cell case
        } else if (geom.volume_fractions(i, j) > 0.0) {
          std::array<Index<Rank + 1>, 4> edges{
              Index<Rank + 1>{i, j, 0}, Index<Rank + 1>{i + 1, j, 0},
              Index<Rank + 1>{i, j, 1}, Index<Rank + 1>{i, j + 1, 1}};
          std::array<double, 4> betas{};
          std::transform(edges.begin(), edges.end(), betas.begin(),
                         [&geom](const Index<Rank + 1>& ijd) {
                           return geom.face_fractions[ijd[2]](ijd[0], ijd[1]);
                         });
          std::array<Index<Rank>, 4> neighbors{
              Index<Rank>{i - 1, j}, Index<Rank>{i + 1, j},
              Index<Rank>{i, j - 1}, Index<Rank>{i, j + 1}};
          std::array<double, 4> alphas{};
          std::transform(neighbors.begin(), neighbors.end(), alphas.begin(),
                         [&geom](const Index<Rank>& ij) {
                           return geom.volume_fractions(ij);
                         });
          std::array<int, 4> is{0, 1, 2, 3};
          std::sort(is.begin(), is.end(), [&](int i, int j) {
            return betas[i] >= betas[j] && alphas[i] >= alphas[j];
          });
          if (betas[is[2]] <= 0.0) {
            // Set neighbors[is[2]] to corner between neighbors[is[0]] and
            // neighbors[is[1]]
            FUB_ASSERT(edges[is[0]][2] != edges[is[1]][2]);
            neighbors[is[2]] = Index<Rank>{neighbors[is[0]][edges[is[0]][2]],
                                           neighbors[is[1]][edges[is[1]][2]]};
          } else if (betas[is[3]] > 0.0) {
            FUB_ASSERT(betas[is[0]] > 0.0 && betas[is[2]] > 0.0 &&
                       betas[is[1]] > 0.0);
            Load(u[0], AsCons(states), {i, j});
            Load(u[1], AsCons(states), neighbors[is[0]]);
            Load(u[2], AsCons(states), neighbors[is[1]]);
            Load(u[3], AsCons(states), neighbors[is[2]]);
            Load(u[4], AsCons(states), neighbors[is[3]]);
            std::array<Coordinates<Rank>, 5> xM;
            xM[0] = GetAbsoluteVolumeCentroid(geom, {i, j}, dx);
            xM[1] = GetAbsoluteVolumeCentroid(geom, neighbors[is[0]], dx);
            xM[2] = GetAbsoluteVolumeCentroid(geom, neighbors[is[1]], dx);
            xM[3] = GetAbsoluteVolumeCentroid(geom, neighbors[is[2]], dx);
            xM[4] = GetAbsoluteVolumeCentroid(geom, neighbors[is[3]], dx);
            ComputeGradients(gradient, u, xM);
          } else if (betas[is[2]] > 0.0) {
            FUB_ASSERT(betas[is[0]] > 0.0 && betas[is[1]] > 0.0);
            Load(u[0], AsCons(states), {i, j});
            Load(u[1], AsCons(states), neighbors[is[0]]);
            Load(u[2], AsCons(states), neighbors[is[1]]);
            Load(u[3], AsCons(states), neighbors[is[2]]);
            std::array<Coordinates<Rank>, 4> xM;
            xM[0] = GetAbsoluteVolumeCentroid(geom, {i, j}, dx);
            xM[1] = GetAbsoluteVolumeCentroid(geom, neighbors[is[0]], dx);
            xM[2] = GetAbsoluteVolumeCentroid(geom, neighbors[is[1]], dx);
            xM[3] = GetAbsoluteVolumeCentroid(geom, neighbors[is[2]], dx);
            ComputeGradients(gradient, span(u).template subspan<0, 4>(), xM);
          }
        }
        Store(gradient_x, gradient[0], {i, j});
        Store(gradient_y, gradient[1], {i, j});
        Store(gradient_z, zero, {i, j});
      });
      ForEachComponent(
          [&](const StridedDataView<double, 2>& u_x,
              const StridedDataView<double, 2>& u_y,
              const StridedDataView<const double, 2>& u) {
            std::array<StridedDataView<double, 2>, 2> grad_u{u_x.Subview(box),
                                                             u_y.Subview(box)};
            Base::LimitGradients(grad_u, u.Subview(box), flags, geom, dx);
          },
          gradient_x, gradient_y, states);
    }
  }

  std::array<double, 2> FindAdmissibleInterval(const Coordinates<Rank>& xM,
                                               const Coordinates<Rank>& xB,
                                               const Coordinates<Rank>& slope,
                                               const Coordinates<Rank>& dx) {
    double lower = std::numeric_limits<double>::lowest();
    double upper = std::numeric_limits<double>::max();
    for (int d = 0; d < Rank; ++d) {
      if (slope[d] != 0) {
        // xM - (xB + slope * lambda) <= 0.5 * dx
        const double lambda_lower = (xM[d] - xB[d] - 0.5 * dx[d]) / slope[d];
        // (xB + slope * lambda) - xM <= 0.5 * dx
        const double lambda_upper = (0.5 * dx[d] + xM[d] - xB[d]) / slope[d];
        if (slope[d] > 0) {
          upper = std::min(upper, lambda_upper);
          lower = std::max(lower, lambda_lower);
        } else {
          upper = std::min(upper, lambda_lower);
          lower = std::max(lower, lambda_upper);
        }
      } else {
        lower = std::numeric_limits<double>::max();
        upper = std::numeric_limits<double>::lowest();
      }
    }
    return {lower, upper};
  }

  std::array<double, 2> Intersect(const std::array<double, 2>& i1,
                                  const std::array<double, 2>& i2) {
    return {std::max(i1[0], i2[0]), std::min(i1[1], i2[1])};
  }

  double Length(const std::array<double, 2>& interval) {
    return interval[1] - interval[0];
  }

  IndexBox<Rank> MakeIndexBox(const Index<Rank>& signs) noexcept {
    IndexBox<Rank> box{};
    for (int d = 0; d < Rank; ++d) {
      if (signs[d] >= 0) {
        box.upper[d] = 2 * signs[d] + 1;
      } else {
        box.upper[d] = 1;
        box.lower[d] = 2 * signs[d];
      }
    }
    return box;
  }

  AuxiliaryReconstructionData Sort(const AuxiliaryReconstructionData& aux_data,
                                   int count) {
    std::array<int, AuxiliaryReconstructionData::kMaxSources> indices;
    std::iota(indices.begin(), indices.begin() + count, 0);
    std::sort(indices.begin(), indices.begin() + count, [&](int i, int j) {
      return aux_data.start[i] < aux_data.start[j];
    });
    AuxiliaryReconstructionData sorted_aux{};
    sorted_aux.slope = aux_data.slope;
    sorted_aux.xB = aux_data.xB;
    for (int i = 0; i < count; ++i) {
      sorted_aux.sources[i] = aux_data.sources[indices[i]];
      sorted_aux.start[i] = aux_data.start[indices[i]];
      sorted_aux.end[i] = aux_data.end[indices[i]];
    }
    for (int i = count; i < AuxiliaryReconstructionData::kMaxSources; ++i) {
      sorted_aux.sources[i] = sorted_aux.sources[count - 1];
      sorted_aux.start[i] = sorted_aux.end[count - 1];
      sorted_aux.end[i] = sorted_aux.end[count - 1];
    }
    return sorted_aux;
  }

  std::array<AuxiliaryReconstructionData, 2>
  SplitAt(const AuxiliaryReconstructionData& aux, double dx) {
    std::array<AuxiliaryReconstructionData, 2> datas{};
    datas[0].slope = aux.slope;
    datas[0].xB = aux.xB;
    datas[1].slope = aux.slope;
    datas[1].xB = aux.xB;
    int i = 0;
    while (i < AuxiliaryReconstructionData::kMaxSources) {
      datas[0].sources[i] = aux.sources[i];
      datas[0].start[i] = aux.start[i];
      if (dx < aux.end[i]) {
        datas[0].end[i] = dx;
        break;
      } else {
        datas[0].end[i] = aux.end[i];
      }
      ++i;
    }
    datas[1].sources[0] = aux.sources[i];
    datas[1].start[0] = dx;
    datas[1].end[0] = aux.end[i];
    for (int j = i + 1; j < AuxiliaryReconstructionData::kMaxSources; ++j) {
      datas[1].sources[j - i] = aux.sources[j];
      datas[1].start[j - i] = aux.start[j];
      datas[1].end[j - i] = aux.end[j];
    }
    return datas;
  }

  AuxiliaryReconstructionData GetAuxiliaryData(const Index<Rank>& index,
                                               const CutCellData<Rank>& geom,
                                               const Coordinates<Rank>& dx,
                                               Direction dir) {
    const Coordinates<Rank> normal = GetBoundaryNormal(geom, index);
    const int dir_v = static_cast<int>(dir);
    const Coordinates<Rank> unit_vector =
        ((normal[dir_v] < 0) - (normal[dir_v] >= 0)) *
        Eigen::Matrix<double, Rank, 1>::Unit(dir_v);
    const Coordinates<Rank> slope =
        unit_vector - 2.0 * unit_vector.dot(normal) * normal;
    Index<Rank> signs{};
    for (int d = 0; d < Rank; ++d) {
      signs[d] = (slope[d] > 0.0) - (slope[d] < 0.0);
    }
    // Fill aux data for the cut-cell itself
    const Coordinates<Rank> xB =
        GetBoundaryCentroid(geom, index).array() * dx.array();
    AuxiliaryReconstructionData aux_data{};
    aux_data.slope = slope;
    aux_data.xB = xB;
    int count = 0;
    ForEachIndex(MakeIndexBox(signs), [&](auto... is) {
      Coordinates<Rank> x_i{double(is)...};
      const Coordinates<Rank> xM = x_i.array() * dx.array();
      const auto interval = FindAdmissibleInterval(xM, xB, slope, dx);
      const auto [a, b] = Intersect(interval, {0.0, 2.0 * dx[int(dir)]});
      if (b - a > 0.0) {
        FUB_ASSERT(count < AuxiliaryReconstructionData::kMaxSources);
        Index<Rank> i{is...};
        Index<Rank> global{};
        std::transform(i.begin(), i.end(), index.begin(), global.begin(),
                       std::plus<>());
        aux_data.sources[count] = global;
        aux_data.start[count] = a;
        aux_data.end[count] = b;
        ++count;
      }
    });
    FUB_ASSERT(count > 0 && count <= AuxiliaryReconstructionData::kMaxSources);
    AuxiliaryReconstructionData sorted_aux = Sort(aux_data, count);
    return sorted_aux;
  }

  double TotalLength(const AuxiliaryReconstructionData& aux) noexcept {
    double total = 0.0;
    for (int i = 0; i < AuxiliaryReconstructionData::kMaxSources; ++i) {
      total += aux.end[i] - aux.start[i];
    }
    return total;
  }

  // int_0^l (u0 + gradient * x) = l ( u_0 + l/2 * gradient )
  bool IntegrateCellState(Conservative& integral,
                          const View<const Conservative>& states,
                          const View<const Conservative>& gradient_x,
                          const View<const Conservative>& gradient_y,
                          const View<const Conservative>& gradient_z,
                          const CutCellData<Rank>& geom,
                          const AuxiliaryReconstructionData& aux_data,
                          const Coordinates<Rank>& dx) {
    ForEachComponent([](double& u) { u = 0.0; }, integral);
    const double total_length = TotalLength(aux_data);
    for (int i = 0; i < AuxiliaryReconstructionData::kMaxSources; ++i) {
      const Index<Rank> index = aux_data.sources[i];
      if (Contains(Box<0>(states), index) &&
          geom.volume_fractions(index) > 0.0) {
        Load(gradient_[0], gradient_x, index);
        Load(gradient_[1], gradient_y, index);
        Load(gradient_[2], gradient_z, index);
        span<const Conservative, Rank> grads{gradient_.data(), Rank};
        const Coordinates<Rank> xM = GetAbsoluteVolumeCentroid(geom, index, dx);
        const Coordinates<Rank> x0 =
            aux_data.xB + aux_data.start[i] * aux_data.slope;
        ApplyGradient(gradient_dir_, grads, x0 - xM);
        Load(state_, states, index);
        state_ += gradient_dir_;
        ApplyGradient(gradient_dir_, grads, aux_data.slope);
        const double length = aux_data.end[i] - aux_data.start[i];
        ForEachVariable(
            [length, total_length](auto&& integral, auto&& u_0, auto&& grad_u) {
              integral += length * (u_0 + 0.5 * length * grad_u) / total_length;
            },
            integral, state_, gradient_dir_);
      } else {
        return false;
      }
    }
    return true;
  }

  void ReconstructSinglyShieldedStencil(
      span<Complete, 4> h_grid_singly_shielded,
      span<const Complete, 4> /*h_grid_embedded_boundary*/,
      span<const Conservative, 2> /*limited_slopes*/,
      const View<const Complete>& states,
      const View<const Conservative>& gradient_x,
      const View<const Conservative>& gradient_y,
      const View<const Conservative>& gradient_z,
      const CutCellData<Rank>& cutcell_data, const Index<Rank>& cell,
      Duration /*dt*/, const Coordinates<Rank>& dx, Direction dir) {
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
      const double d = 0.5 - cutcell_data.boundary_centeroids(cell_d);

      AuxiliaryReconstructionData aux_data =
          GetAuxiliaryData(cell, cutcell_data, dx, dir);
      auto [aux_near, aux2] = SplitAt(aux_data, (1.0 - d) * dx[d]);
      auto [aux_far, aux_rest] = SplitAt(aux2, dx[d]);
      if (!IntegrateCellState(near_, states, gradient_x, gradient_y, gradient_z,
                              cutcell_data, aux_near, dx)) {
        Load(near_, AsCons(states), cell);
      }
      if (!IntegrateCellState(far_, states, gradient_x, gradient_y, gradient_z,
                              cutcell_data, aux_far, dx)) {
        far_ = near_;
      }

      // Compute slopes, centroid is in local cell coordinates [-0.5, 0.5]

      // Do the h-grid average for left states
      const double alpha = d;
      const double omalpha = 1.0 - alpha;
      ForEachComponent([=](double& hL, double uCC,
                           double uL) { hL = alpha * uCC + omalpha * uL; },
                       AsCons(h_grid_singly_shielded[1]), stencil[2], near_);
      CompleteFromCons(equation_, h_grid_singly_shielded[0], far_);
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

      const double d = 0.5 + cutcell_data.boundary_centeroids(cell_d);

      AuxiliaryReconstructionData aux_data =
          GetAuxiliaryData(cell, cutcell_data, dx, dir);
      auto [aux_near, aux2] = SplitAt(aux_data, (1.0 - d) * dx[d]);
      auto [aux_far, aux_rest] = SplitAt(aux2, dx[d]);
      if (!IntegrateCellState(near_, states, gradient_x, gradient_y, gradient_z,
                              cutcell_data, aux_near, dx)) {
        Load(near_, AsCons(states), cell);
      }
      if (!IntegrateCellState(far_, states, gradient_x, gradient_y, gradient_z,
                              cutcell_data, aux_far, dx)) {
        far_ = near_;
      }

      // Do the h-grid average for right states
      const double alpha = d;
      const double omalpha = 1.0 - alpha;
      ForEachComponent([=](double& hR, double uCC,
                           double uR) { hR = alpha * uCC + omalpha * uR; },
                       AsCons(h_grid_singly_shielded[2]), stencil[0], near_);
      CompleteFromCons(equation_, h_grid_singly_shielded[2],
                       h_grid_singly_shielded[2]);
      CompleteFromCons(equation_, h_grid_singly_shielded[3], far_);
    }
  }

  void ReconstructEmbeddedBoundaryStencil(
      span<Complete, 4> h_grid_embedded_boundary,
      span<Conservative, 2> limited_slopes, const View<const Complete>& states,
      const View<const Conservative>& gradient_x,
      const View<const Conservative>& gradient_y,
      const View<const Conservative>& gradient_z,
      const CutCellData<Rank>& cutcell_data, const Index<Rank>& cell,
      Duration /*dt*/, Eigen::Matrix<double, Rank, 1> dx, Direction dir) {
    AuxiliaryReconstructionData aux_data =
        GetAuxiliaryData(cell, cutcell_data, dx, dir);
    const int dir_v = static_cast<int>(dir);
    auto [aux_near, aux_far] = SplitAt(aux_data, dx[dir_v]);
    if (!IntegrateCellState(near_, states, gradient_x, gradient_y, gradient_z,
                            cutcell_data, aux_near, dx)) {
      Load(near_, AsCons(states), cell);
    }
    if (!IntegrateCellState(far_, states, gradient_x, gradient_y, gradient_z,
                            cutcell_data, aux_far, dx)) {
      far_ = near_;
    }
    const int d = static_cast<std::size_t>(dir);
    const Index<Rank + 1> cell_d = std::apply(
        [=](auto... is) {
          return Index<Rank + 1>{is..., dir_v};
        },
        cell);
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
      AsCons(h_grid_embedded_boundary[1]) = near_;
      AsCons(h_grid_embedded_boundary[0]) = far_;
      CompleteFromCons(equation_, h_grid_embedded_boundary[0],
                       h_grid_embedded_boundary[0]);
      CompleteFromCons(equation_, h_grid_embedded_boundary[1],
                       h_grid_embedded_boundary[1]);
      Reflect(h_grid_embedded_boundary[1], h_grid_embedded_boundary[1], normal,
              equation_);
      Reflect(h_grid_embedded_boundary[0], h_grid_embedded_boundary[0], normal,
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
      AsCons(h_grid_embedded_boundary[2]) = near_;
      AsCons(h_grid_embedded_boundary[3]) = far_;
      CompleteFromCons(equation_, h_grid_embedded_boundary[2],
                       h_grid_embedded_boundary[2]);
      CompleteFromCons(equation_, h_grid_embedded_boundary[3],
                       h_grid_embedded_boundary[3]);
      Reflect(h_grid_embedded_boundary[2], h_grid_embedded_boundary[2], normal,
              equation_);
      Reflect(h_grid_embedded_boundary[3], h_grid_embedded_boundary[3], normal,
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
            slope = 0.0;
            // slope = std::min(sL, sR);
          }
        },
        slope, u[0], u[1], u[2]);
  }

  Equation equation_;
  std::array<Conservative, 3> gradient_{Conservative(equation_),
                                        Conservative(equation_),
                                        Conservative(equation_)};
  std::array<Conservative, 4> stencil{
      Conservative(equation_), Conservative(equation_), Conservative(equation_),
      Conservative(equation_)};
  Conservative limited_slope_{equation_};
  Conservative state_{equation_};
  Conservative gradient_dir_{equation_};
  Conservative near_{equation_};
  Conservative far_{equation_};
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
                            const View<const Conservative>& gradient_x,
                            const View<const Conservative>& gradient_y,
                            const View<const Conservative>& gradient_z,
                            const View<const Conservative>& regular_fluxes,
                            const View<const Complete>& states,
                            const CutCellData<Rank>& geom, Duration dt,
                            const Eigen::Matrix<double, Rank, 1>& dx,
                            Direction dir);

  void ComputeGradients(const View<Conservative>& gradient_x,
                        const View<Conservative>& gradient_y,
                        const View<Conservative>& gradient_z,
                        const View<const Conservative>& states,
                        const StridedDataView<const char, Rank>& flags,
                        const CutCellData<Rank>& geom,
                        const Coordinates<Rank>& dx) {
    h_grid_reconstruction_.ComputeGradients(gradient_x, gradient_y, gradient_z,
                                            states, flags, geom, dx);
  }

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
                         const View<const Conservative>& gradient_x,
                         const View<const Conservative>& gradient_y,
                         const View<const Conservative>& gradient_z,
                         const View<const Conservative>& regular_fluxes,
                         const View<const Complete>& states,
                         const CutCellData<Rank>& geom, Duration dt,
                         const Eigen::Matrix<double, Rank, 1>& dx,
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
          h_grid_eb_, limited_slopes_, states, gradient_x, gradient_y,
          gradient_z, geom, cell, dt, dx, dir);
      FluxMethod::ComputeNumericFlux(boundary_flux_, h_grid_eb_, dt, dx[d],
                                     dir);
      h_grid_reconstruction_.ReconstructSinglyShieldedStencil(
          h_grid_singly_shielded_, h_grid_eb_, limited_slopes_, states,
          gradient_x, gradient_y, gradient_z, geom, cell, dt, dx, dir);
      FluxMethod::ComputeNumericFlux(shielded_left_flux_,
                                     h_grid_singly_shielded_, dt, dx[d], dir);
      Store(boundary_fluxes, boundary_flux_, cell);
      if (Contains(Box<0>(shielded_left_fluxes), faceR)) {
        Store(shielded_left_fluxes, shielded_left_flux_, faceR);
      }
    } else if (betaR < betaL) {
      h_grid_reconstruction_.ReconstructEmbeddedBoundaryStencil(
          h_grid_eb_, limited_slopes_, states, gradient_x, gradient_y,
          gradient_z, geom, cell, dt, dx, dir);
      FluxMethod::ComputeNumericFlux(boundary_flux_, h_grid_eb_, dt, dx[d],
                                     dir);
      h_grid_reconstruction_.ReconstructSinglyShieldedStencil(
          h_grid_singly_shielded_, h_grid_eb_, limited_slopes_, states,
          gradient_x, gradient_y, gradient_z, geom, cell, dt, dx, dir);
      FluxMethod::ComputeNumericFlux(shielded_right_flux_,
                                     h_grid_singly_shielded_, dt, dx[d], dir);
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
                                           geom, dt, dx[d], dir);
      },
      stabilised_fluxes, shielded_left_fluxes, shielded_right_fluxes,
      regular_fluxes, boundary_fluxes);
}

} // namespace fub

#endif
