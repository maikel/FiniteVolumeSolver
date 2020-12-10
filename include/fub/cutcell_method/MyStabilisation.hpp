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

inline Coordinates<2>
ComputeReflectedCoordinates(const Coordinates<2>& offset,
                            const Coordinates<2>& boundary_normal) {
  Coordinates<2> reflected =
      offset - 2.0 * offset.dot(boundary_normal) * boundary_normal;
  return reflected;
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
        return Index<Rank + 1>{is..., static_cast<int>(dir)};
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

    int n_sources;

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
        ////////////////////////////////////////////////
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
        } else if (geom.volume_fractions(i, j) == 0.0) {
          SetZero(gradient[0]);
          SetZero(gradient[1]);
        } else {
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
            FUB_ASSERT(betas[is[0]] > 0.0 && betas[is[1]] > 0.0);
            FUB_ASSERT(edges[is[0]][2] != edges[is[1]][2]);
            Index<Rank> corner{};
            corner[edges[is[0]][2]] = neighbors[is[0]][edges[is[0]][2]];
            corner[edges[is[1]][2]] = neighbors[is[1]][edges[is[1]][2]];
            neighbors[is[2]] = corner;
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
            // Load(u[0], AsCons(states), {i, j});
            // Load(u[1], AsCons(states), neighbors[is[0]]);
            // Load(u[2], AsCons(states), neighbors[is[1]]);
            // u[3] = u[1];
            // u[4] = u[2];
            // const Coordinates<Rank> xB =
            //     GetAbsoluteBoundaryCentroid(geom, {i, j}, dx);
            // const Coordinates<Rank> n = GetBoundaryNormal(geom, {i, j});
            // std::array<Coordinates<Rank>, 5> xM;
            // xM[0] = GetAbsoluteVolumeCentroid(geom, {i, j}, dx);
            // xM[1] = GetAbsoluteVolumeCentroid(geom, neighbors[is[0]], dx);
            // xM[2] = GetAbsoluteVolumeCentroid(geom, neighbors[is[1]], dx);
            // xM[3] = xB + ComputeReflectedCoordinates(xM[1] - xB, n);
            // xM[4] = xB + ComputeReflectedCoordinates(xM[2] - xB, n);
            // ComputeGradients(gradient, span(u), xM);
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
      } else if (std::abs(xM[d] - xB[d]) > 0.5 * dx[d]) {
        lower = std::numeric_limits<double>::max();
        upper = std::numeric_limits<double>::lowest();
      }
    }
    return {lower, upper};
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

  AuxiliaryReconstructionData Sort(const AuxiliaryReconstructionData& aux_data) {
    std::array<int, AuxiliaryReconstructionData::kMaxSources> indices;
    const int count = aux_data.n_sources;
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
    sorted_aux.n_sources = count;
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
      datas[0].n_sources = i + 1;
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
      datas[1].n_sources = aux.n_sources - datas[0].n_sources;
    }
    return datas;
  }

  void SetZero(Conservative& cons) {
    ForEachComponent([](auto&& u) { u = 0.0; }, cons);
  }

  AuxiliaryReconstructionData GetAuxiliaryData(const Index<Rank>& index,
                                               const CutCellData<Rank>& geom,
                                               const Coordinates<Rank>& dx,
                                               const Coordinates<Rank>& slope,
                                               double required_length) {
    Index<Rank> signs{};
    for (int d = 0; d < Rank; ++d) {
      signs[d] = (slope[d] > 0.0) - (slope[d] < 0.0);
    }
    // Fill aux data for the cut-cell itself
    const Coordinates<Rank> xB =
        GetBoundaryCentroid(geom, index).array() * dx.array();
    AuxiliaryReconstructionData aux_data{};
    aux_data.slope = slope;
    aux_data.xB = GetAbsoluteBoundaryCentroid(geom, index, dx);
    int count = 0;
    ForEachIndex(MakeIndexBox(signs), [&](auto... is) {
      Coordinates<Rank> x_i{double(is)...};
      const Coordinates<Rank> xM = x_i.array() * dx.array();
      const auto interval = FindAdmissibleInterval(xM, xB, slope, dx);
      const auto [a, b] = Intersect(interval, {0.0, required_length});
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
    aux_data.n_sources = count;
    AuxiliaryReconstructionData sorted_aux = Sort(aux_data);
    return sorted_aux;
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
    return GetAuxiliaryData(index, geom, dx, slope, dx[int(dir)]);
  }

  double TotalLength(const AuxiliaryReconstructionData& aux) noexcept {
    double total = 0.0;
    for (int i = 0; i < AuxiliaryReconstructionData::kMaxSources; ++i) {
      total += aux.end[i] - aux.start[i];
    }
    return total;
  }

  bool IntegrateCellState(Conservative& integral,
                          Conservative& integral_gradient,
                          const View<const Conservative>& states,
                          const View<const Conservative>& gradient_x,
                          const View<const Conservative>& gradient_y,
                          const View<const Conservative>& gradient_z,
                          const CutCellData<Rank>& geom,
                          const AuxiliaryReconstructionData& aux_data,
                          const Coordinates<Rank>& dx, double total_length) {
    for (int i = 0; i < aux_data.n_sources; ++i) {
      const Index<Rank> index = aux_data.sources[i];
      if (Contains(Box<0>(states), index) &&
          geom.volume_fractions(index) > 0.0) {
        Load(gradient_[0], gradient_x, index);
        Load(gradient_[1], gradient_y, index);
        Load(gradient_[2], gradient_z, index);
        span<const Conservative, Rank> grads{gradient_.data(), Rank};
        const Coordinates<Rank> xM = GetAbsoluteVolumeCentroid(geom, index, dx);
        const Coordinates<Rank> x0 =
            aux_data.xB +
            0.5 * (aux_data.start[i] + aux_data.end[i]) * aux_data.slope;
        const Coordinates<Rank> dx = x0 - xM;
        ApplyGradient(gradient_dir_, grads, dx);
        Load(state_, states, index);
        state_ += gradient_dir_;
        ApplyGradient(gradient_dir_, grads, aux_data.slope);
        const double length = aux_data.end[i] - aux_data.start[i];
        ForEachVariable(
            [length, total_length](auto&& u, auto&& grad_u, auto&& u_0,
                                   auto&& grad_u_0) {
              u += length / total_length * u_0;
              grad_u += length / total_length * grad_u_0;
            },
            integral, integral_gradient, state_, gradient_dir_);
      } else {
        return false;
      }
    }
    return true;
  }

  // int_0^l (u0 + gradient * x) = l ( u_0 + l/2 * gradient )
  bool IntegrateCellState(Conservative& integral,
                          Conservative& integral_gradient,
                          const View<const Conservative>& states,
                          const View<const Conservative>& gradient_x,
                          const View<const Conservative>& gradient_y,
                          const View<const Conservative>& gradient_z,
                          const CutCellData<Rank>& geom,
                          const AuxiliaryReconstructionData& aux_data,
                          const Coordinates<Rank>& dx) {
    const double total_length = TotalLength(aux_data);
    return IntegrateCellState(integral, integral_gradient, states, gradient_x,
                              gradient_y, gradient_z, geom, aux_data, dx,
                              total_length);
  }

  void ReconstructSinglyShieldedStencil(
      span<Complete, 2> h_grid_singly_shielded,
      span<Conservative, 2> h_grid_singly_shielded_gradients,
      const View<const Complete>& states,
      const View<const Conservative>& gradient_x,
      const View<const Conservative>& gradient_y,
      const View<const Conservative>& gradient_z,
      const CutCellData<Rank>& cutcell_data, const Index<Rank>& cell,
      Duration /*dt*/, const Coordinates<Rank>& dx, Direction dir) {
    AuxiliaryReconstructionData boundary_aux_data =
        GetAuxiliaryData(cell, cutcell_data, dx, dir);
    const Coordinates<Rank> normal = GetBoundaryNormal(cutcell_data, cell);
    const int dir_v = static_cast<int>(dir);
    const Coordinates<Rank> unit_vector =
        ((normal[dir_v] >= 0) - (normal[dir_v] < 0)) *
        Eigen::Matrix<double, Rank, 1>::Unit(dir_v);
    AuxiliaryReconstructionData interior_aux_data =
        GetAuxiliaryData(cell, cutcell_data, dx, unit_vector, 2.0 * dx[dir_v]);
    SetZero(boundary_state_);
    SetZero(boundary_gradient_);
    SetZero(interior_state_);
    SetZero(interior_gradient_);
    const int d = static_cast<std::size_t>(dir);
    const Index<Rank> face_L = cell;
    const Index<Rank> face_R = Shift(face_L, dir, 1);
    const double betaL = cutcell_data.face_fractions[d](face_L);
    const double betaR = cutcell_data.face_fractions[d](face_R);

    const Index<Rank + 1> cell_d = EmbedIndex<Rank>(cell, dir);

    if (betaL < betaR) {
      const double d = 0.5 - cutcell_data.boundary_centeroids(cell_d);

      auto [laux1, raux1] = SplitAt(boundary_aux_data, (1.0 - d) * dx[d]);
      auto [laux2, raux2] = SplitAt(interior_aux_data, d * dx[d]);
      auto [interior_aux, raux3] = SplitAt(raux2, (1.0 + d) * dx[d]);

      if (!IntegrateCellState(boundary_state_, boundary_gradient_, states,
                              gradient_x, gradient_y, gradient_z, cutcell_data,
                              laux1, dx, dx[d])) {
        Load(boundary_state_, AsCons(states), cell);
        SetZero(boundary_gradient_);
      } else {
        ForEachComponent([](auto&& u) { u *= -1.0; }, boundary_gradient_);
        if (!IntegrateCellState(boundary_state_, boundary_gradient_, states,
                                gradient_x, gradient_y, gradient_z,
                                cutcell_data, laux2, dx, dx[d])) {
          Load(boundary_state_, AsCons(states), cell);
          SetZero(boundary_gradient_);
        }
      }
      Reflect(boundary_state_, boundary_state_, normal, equation_);
      Reflect(boundary_gradient_, boundary_gradient_, normal, equation_);

      CompleteFromCons(equation_, h_grid_singly_shielded[0], boundary_state_);
      h_grid_singly_shielded_gradients[0] = boundary_gradient_;

      if (!IntegrateCellState(interior_state_, interior_gradient_, states,
                              gradient_x, gradient_y, gradient_z, cutcell_data,
                              interior_aux, dx)) {
        Load(interior_state_, AsCons(states), cell);
        SetZero(interior_gradient_);
      }

      CompleteFromCons(equation_, h_grid_singly_shielded[1], interior_state_);
      h_grid_singly_shielded_gradients[1] = interior_gradient_;

    } else if (betaR < betaL) {
      const double d = 0.5 + cutcell_data.boundary_centeroids(cell_d);

      auto [laux1, raux1] = SplitAt(boundary_aux_data, (1.0 - d) * dx[d]);
      auto [laux2, raux2] = SplitAt(interior_aux_data, d * dx[d]);
      auto [interior_aux, raux3] = SplitAt(raux2, (1.0 + d) * dx[d]);

      if (!IntegrateCellState(interior_state_, interior_gradient_, states,
                              gradient_x, gradient_y, gradient_z, cutcell_data,
                              interior_aux, dx)) {
        Load(interior_state_, AsCons(states), cell);
        SetZero(interior_gradient_);
      } else {
        ForEachComponent([](auto&& u) { u *= -1.0; }, interior_gradient_);
      }

      CompleteFromCons(equation_, h_grid_singly_shielded[0], interior_state_);
      h_grid_singly_shielded_gradients[0] = interior_gradient_;

      if (!IntegrateCellState(boundary_state_, boundary_gradient_, states,
                              gradient_x, gradient_y, gradient_z, cutcell_data,
                              laux2, dx, dx[d])) {
        Load(boundary_state_, AsCons(states), cell);
        SetZero(boundary_gradient_);
      } else {
        ForEachComponent([](auto&& u) { u *= -1.0; }, boundary_gradient_);
        if (!IntegrateCellState(boundary_state_, boundary_gradient_, states,
                                gradient_x, gradient_y, gradient_z,
                                cutcell_data, laux1, dx, dx[d])) {
          Load(boundary_state_, AsCons(states), cell);
          SetZero(boundary_gradient_);
        }
      }
      Reflect(boundary_state_, boundary_state_, normal, equation_);
      Reflect(boundary_gradient_, boundary_gradient_, normal, equation_);

      CompleteFromCons(equation_, h_grid_singly_shielded[1], boundary_state_);
      h_grid_singly_shielded_gradients[1] = boundary_gradient_;
    }
  }

  void ReconstructEmbeddedBoundaryStencil(
      span<Complete, 2> h_grid_embedded_boundary,
      span<Conservative, 2> h_grid_embedded_boundary_slopes,
      const View<const Complete>& states,
      const View<const Conservative>& gradient_x,
      const View<const Conservative>& gradient_y,
      const View<const Conservative>& gradient_z,
      const CutCellData<Rank>& cutcell_data, const Index<Rank>& cell,
      Duration /*dt*/, Eigen::Matrix<double, Rank, 1> dx, Direction dir) {
    AuxiliaryReconstructionData boundary_aux_data =
        GetAuxiliaryData(cell, cutcell_data, dx, dir);
    const Coordinates<Rank> normal = GetBoundaryNormal(cutcell_data, cell);
    const int dir_v = static_cast<int>(dir);
    const Coordinates<Rank> unit_vector =
        ((normal[dir_v] >= 0) - (normal[dir_v] < 0)) *
        Eigen::Matrix<double, Rank, 1>::Unit(dir_v);
    AuxiliaryReconstructionData interior_aux_data =
        GetAuxiliaryData(cell, cutcell_data, dx, unit_vector, dx[dir_v]);
    SetZero(boundary_state_);
    SetZero(boundary_gradient_);
    if (!IntegrateCellState(boundary_state_, boundary_gradient_, states,
                            gradient_x, gradient_y, gradient_z, cutcell_data,
                            boundary_aux_data, dx)) {
      Load(boundary_state_, AsCons(states), cell);
      SetZero(boundary_gradient_);
    }
    SetZero(interior_state_);
    SetZero(interior_gradient_);
    if (!IntegrateCellState(interior_state_, interior_gradient_, states,
                            gradient_x, gradient_y, gradient_z, cutcell_data,
                            interior_aux_data, dx)) {
      interior_state_ = boundary_state_;
      interior_gradient_ = boundary_gradient_;
    }
    const int d = static_cast<std::size_t>(dir);
    const Index<Rank> face_L = cell;
    const Index<Rank> face_R = Shift(face_L, dir, 1);
    const double betaL = cutcell_data.face_fractions[d](face_L);
    const double betaR = cutcell_data.face_fractions[d](face_R);

    Reflect(boundary_state_, boundary_state_, normal, equation_);
    Reflect(boundary_gradient_, boundary_gradient_, normal, equation_);

    if (betaL < betaR) {
      ForEachComponent([](auto&& u) { u *= -1.0; }, boundary_gradient_);

      CompleteFromCons(equation_, h_grid_embedded_boundary[0], boundary_state_);
      h_grid_embedded_boundary_slopes[0] = boundary_gradient_;
      CompleteFromCons(equation_, h_grid_embedded_boundary[1], interior_state_);
      h_grid_embedded_boundary_slopes[1] = interior_gradient_;

    } else if (betaR < betaL) {
      ForEachComponent([](auto&& u) { u *= -1.0; }, interior_gradient_);

      CompleteFromCons(equation_, h_grid_embedded_boundary[1], boundary_state_);
      h_grid_embedded_boundary_slopes[1] = boundary_gradient_;
      CompleteFromCons(equation_, h_grid_embedded_boundary[0], interior_state_);
      h_grid_embedded_boundary_slopes[0] = interior_gradient_;
    }
  }

  void ReconstructRegularStencil(span<Complete, 2> h_grid_regular,
                                 span<Conservative, 2> h_grid_regular_gradients,
                                 const View<const Complete>& states,
                                 const View<const Conservative>& gradient_x,
                                 const View<const Conservative>& gradient_y,
                                 const View<const Conservative>& gradient_z,
                                 const CutCellData<Rank>& geom,
                                 const Index<Rank>& face, Duration /*dt*/,
                                 Eigen::Matrix<double, Rank, 1> dx,
                                 Direction dir) {
    Index<Rank> iL = LeftTo(face, dir);
    Index<Rank> iR = RightTo(face, dir);

    const Coordinates<Rank> face_xM =
        GetUnshieldedCentroid(geom, face, dx, dir);
    const Coordinates<Rank> xL_us = Shift(face_xM, dir, -0.5 * dx[int(dir)]);
    const Coordinates<Rank> xR_us = Shift(face_xM, dir, +0.5 * dx[int(dir)]);
    const Coordinates<Rank> xL = GetAbsoluteVolumeCentroid(geom, iL, dx);
    const Coordinates<Rank> xR = GetAbsoluteVolumeCentroid(geom, iR, dx);
    Load(state_, AsCons(states), iL);
    Load(gradient_[0], gradient_x, iL);
    Load(gradient_[1], gradient_y, iL);
    Load(gradient_[2], gradient_z, iL);

    span<const Conservative, Rank> grads{gradient_.data(), Rank};

    const Coordinates<Rank> delta_xL = xL_us - xL;
    ApplyGradient(gradient_dir_, grads, delta_xL);
    state_ += gradient_dir_;
    CompleteFromCons(equation_, h_grid_regular[0], state_);
    h_grid_regular_gradients[0] = gradient_[int(dir)];

    Load(state_, AsCons(states), iR);
    Load(gradient_[0], gradient_x, iR);
    Load(gradient_[1], gradient_y, iR);
    Load(gradient_[2], gradient_z, iR);

    const Coordinates<Rank> delta_xR = xR_us - xR;
    ApplyGradient(gradient_dir_, grads, delta_xR);
    state_ += gradient_dir_;
    CompleteFromCons(equation_, h_grid_regular[1], state_);
    h_grid_regular_gradients[1] = gradient_[int(dir)];
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
  Conservative boundary_state_{equation_};
  Conservative boundary_gradient_{equation_};
  Conservative interior_state_{equation_};
  Conservative interior_gradient_{equation_};
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
  MyCutCellMethod(const Equation& equation, const FluxMethod& flux_method);

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
                            const View<Conservative>& regular_fluxes,
                            const View<Conservative>& boundary_fluxes,
                            const View<const Conservative>& gradient_x,
                            const View<const Conservative>& gradient_y,
                            const View<const Conservative>& gradient_z,
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

  std::array<Complete, 2> h_grid_eb_{};
  std::array<Conservative, 2> h_grid_eb_gradients_{};
  std::array<Complete, 2> h_grid_singly_shielded_{};
  std::array<Conservative, 2> h_grid_singly_shielded_gradients_{};
  std::array<Complete, 2> h_grid_regular_{};
  std::array<Conservative, 2> h_grid_regular_gradients_{};

  Conservative boundary_flux_{equation_};
  Conservative singly_shielded_flux_{equation_};
  Conservative regular_flux_{equation_};

  std::array<CompleteArray, StencilSize> stencil_array_{};
  ConservativeArray numeric_flux_array_{equation_};
};

template <typename Equation, typename FluxMethod>
MyCutCellMethod(const Equation&, const FluxMethod&) -> MyCutCellMethod<Equation, FluxMethod>;

// IMPLEMENTATION

template <typename Equation, typename FluxMethod, typename HGridReconstruction>
MyCutCellMethod<Equation, FluxMethod, HGridReconstruction>::MyCutCellMethod(
    const Equation& eq)
    : MyCutCellMethod(eq, FluxMethod(eq)) {}

template <typename Equation, typename FluxMethod, typename HGridReconstruction>
MyCutCellMethod<Equation, FluxMethod, HGridReconstruction>::MyCutCellMethod(
    const Equation& eq, const FluxMethod& flux_method)
    : FluxMethod(flux_method), equation_(eq), h_grid_reconstruction_(eq) {
  h_grid_eb_.fill(Complete(equation_));
  h_grid_eb_gradients_.fill(Conservative(equation_));
  h_grid_singly_shielded_.fill(Complete(equation_));
  h_grid_singly_shielded_gradients_.fill(Conservative(equation_));
  h_grid_regular_.fill(Complete(equation_));
  h_grid_regular_gradients_.fill(Conservative(equation_));
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
                         const View<Conservative>& regular_fluxes,
                         const View<Conservative>& boundary_fluxes,
                         const View<const Conservative>& gradient_x,
                         const View<const Conservative>& gradient_y,
                         const View<const Conservative>& gradient_z,
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
  const PatchDataView<const double, Rank>& betaUs =
      geom.unshielded_fractions[d];
  ForEachIndex(Shrink(Box<0>(states), dir, {1, 1}), [&](auto... is) {
    Index<Rank> cell{is...};
    Index<Rank> faceL = cell;
    Index<Rank> faceR = Shift(faceL, dir, 1);
    const double betaL = betas(faceL);
    const double betaR = betas(faceR);
    const double betaUsL = betaUs(faceL);
    const double betaUsR = betaUs(faceR);
    if (betaL < betaR) {
      h_grid_reconstruction_.ReconstructEmbeddedBoundaryStencil(
          h_grid_eb_, h_grid_eb_gradients_, states, gradient_x, gradient_y,
          gradient_z, geom, cell, dt, dx, dir);
      FluxMethod::ComputeNumericFlux(boundary_flux_, h_grid_eb_,
                                     h_grid_eb_gradients_, dt, dx[d], dir);
      Store(boundary_fluxes, boundary_flux_, cell);

      h_grid_reconstruction_.ReconstructSinglyShieldedStencil(
          h_grid_singly_shielded_, h_grid_singly_shielded_gradients_, states,
          gradient_x, gradient_y, gradient_z, geom, cell, dt, dx, dir);
      FluxMethod::ComputeNumericFlux(
          singly_shielded_flux_, h_grid_singly_shielded_,
          h_grid_singly_shielded_gradients_, dt, dx[d], dir);
      if (Contains(Box<0>(shielded_left_fluxes), faceR)) {
        Store(shielded_left_fluxes, singly_shielded_flux_, faceR);
        if (betaUsR > 0.0) {
          h_grid_reconstruction_.ReconstructRegularStencil(
              h_grid_regular_, h_grid_regular_gradients_, states, gradient_x,
              gradient_y, gradient_z, geom, faceR, dt, dx, dir);
          FluxMethod::ComputeNumericFlux(regular_flux_, h_grid_regular_,
                                         h_grid_regular_gradients_, dt, dx[d],
                                         dir);
          Store(regular_fluxes, regular_flux_, faceR);
        }
      }
    } else if (betaR < betaL) {
      h_grid_reconstruction_.ReconstructEmbeddedBoundaryStencil(
          h_grid_eb_, h_grid_eb_gradients_, states, gradient_x, gradient_y,
          gradient_z, geom, cell, dt, dx, dir);
      FluxMethod::ComputeNumericFlux(boundary_flux_, h_grid_eb_,
                                     h_grid_eb_gradients_, dt, dx[d], dir);
      Store(boundary_fluxes, boundary_flux_, cell);

      h_grid_reconstruction_.ReconstructSinglyShieldedStencil(
          h_grid_singly_shielded_, h_grid_singly_shielded_gradients_, states,
          gradient_x, gradient_y, gradient_z, geom, cell, dt, dx, dir);
      FluxMethod::ComputeNumericFlux(
          singly_shielded_flux_, h_grid_singly_shielded_,
          h_grid_singly_shielded_gradients_, dt, dx[d], dir);
      if (Contains(Box<0>(shielded_right_fluxes), faceL)) {
        Store(shielded_right_fluxes, singly_shielded_flux_, faceL);
        if (betaUsL > 0.0) {
          h_grid_reconstruction_.ReconstructRegularStencil(
              h_grid_regular_, h_grid_regular_gradients_, states, gradient_x,
              gradient_y, gradient_z, geom, faceL, dt, dx, dir);
          FluxMethod::ComputeNumericFlux(regular_flux_, h_grid_regular_,
                                         h_grid_regular_gradients_, dt, dx[d],
                                         dir);
          Store(regular_fluxes, regular_flux_, faceL);
        }
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
