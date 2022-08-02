#ifndef FUB_CUTCELL_MD_GRADIENTS_HPP
#define FUB_CUTCELL_MD_GRADIENTS_HPP

#include "fub/cutcell_method/AnyMdLimiter.hpp"

#include "fub/CutCellData.hpp"
#include "fub/Equation.hpp"
#include "fub/ForEach.hpp"
#include "fub/PatchDataView.hpp"
#include "fub/State.hpp"
#include "fub/core/span.hpp"
#include "fub/ext/Eigen.hpp"

namespace fub {
template <typename GradientMethod> class MdGradients {
public:
  using Equation = typename GradientMethod::Equation;
  using Gradient = typename GradientMethod::Gradient;
  using Complete = typename GradientMethod::Complete;

  static constexpr int Rank = Equation::Rank();

  MdGradients(GradientMethod gradient_method, AnyLimiter<Rank> limiter)
      : equation_(gradient_method.GetEquation()),
        gradient_method_(std::move(gradient_method)),
        limiter_(std::move(limiter)) {}

  void ComputeGradients(const View<Gradient>& gradient_x,
                        const View<Gradient>& gradient_y,
                        const View<Gradient>& gradient_z,
                        const View<const Complete>& states,
                        const CutCellData<Rank>& geom,
                        const Coordinates<Rank>& dx);

private:
  void ComputeGradients(span<Gradient, 2> gradient,
                        span<const Complete, 4> states,
                        span<const Coordinates<Rank>, 4> x);

  void ComputeGradients(span<Gradient, 2> gradient,
                        span<const Complete, 5> states,
                        span<const Coordinates<Rank>, 5> x);

  void ComputeGradients(span<double, 2> gradient, span<const double, 4> states,
                        span<const Coordinates<Rank>, 4> x);

  void ComputeGradients(span<double, 2> gradient, span<const double, 5> states,
                        span<const Coordinates<Rank>, 5> x);

  void
  LimitGradients(const std::array<StridedDataView<double, Rank>, Rank>& grad_u,
                 StridedDataView<const double, Rank> u,
                 const CutCellData<Rank>& geom,
                 const Coordinates<Rank>& dx) const;

private:
  Equation equation_;
  GradientMethod gradient_method_;
  AnyLimiter<Rank> limiter_;
};

template <typename GradientMethod>
void MdGradients<GradientMethod>::ComputeGradients(
    const View<Gradient>& gradient_x, const View<Gradient>& gradient_y,
    const View<Gradient>& gradient_z, const View<const Complete>& states,
    const CutCellData<Rank>& geom, const Coordinates<Rank>& dx) {
  Gradient zero{equation_};
  if constexpr (Rank == 2) {
    std::array<Gradient, 2> gradient{Gradient{equation_}, Gradient{equation_}};
    std::array<Complete, 5> u;
    u.fill(Complete{equation_});
    const IndexBox<Rank> box = Box<0>(gradient_x);
    FUB_ASSERT(box == Box<0>(gradient_y));
    FUB_ASSERT(box == Box<0>(gradient_z));
    ForEachIndex(box, [&](int i, int j) {
      ////////////////////////////////////////////////
      // All regular case
      if (geom.volume_fractions(i, j) == 1.0 &&
          geom.volume_fractions(i + 1, j) == 1.0 &&
          geom.volume_fractions(i - 1, j) == 1.0 &&
          geom.volume_fractions(i, j + 1) == 1.0 &&
          geom.volume_fractions(i, j - 1) == 1.0) {
        Load(u[0], states, {i - 1, j});
        Load(u[1], states, {i, j});
        Load(u[2], states, {i + 1, j});
        gradient_method_.ComputeGradient(
            gradient[0], span<const Complete>{u}.template subspan<0, 3>(),
            dx[0], Direction::X);

        Load(u[0], states, {i, j - 1});
        Load(u[1], states, {i, j});
        Load(u[2], states, {i, j + 1});
        gradient_method_.ComputeGradient(
            gradient[1], span<const Complete>{u}.template subspan<0, 3>(),
            dx[1], Direction::Y);
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
            Index<Rank>{i - 1, j}, Index<Rank>{i + 1, j}, Index<Rank>{i, j - 1},
            Index<Rank>{i, j + 1}};
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
          Load(u[0], states, {i, j});
          Load(u[1], states, neighbors[is[0]]);
          Load(u[2], states, neighbors[is[1]]);
          Load(u[3], states, neighbors[is[2]]);
          std::array<Coordinates<Rank>, 4> xM;
          xM[0] = GetAbsoluteVolumeCentroid(geom, {i, j}, dx);
          xM[1] = GetAbsoluteVolumeCentroid(geom, neighbors[is[0]], dx);
          xM[2] = GetAbsoluteVolumeCentroid(geom, neighbors[is[1]], dx);
          xM[3] = GetAbsoluteVolumeCentroid(geom, neighbors[is[2]], dx);
          ComputeGradients(
              gradient, fub::span<const Complete>(u).template subspan<0, 4>(),
              xM);
        } else if (betas[is[3]] > 0.0) {
          FUB_ASSERT(betas[is[0]] > 0.0 && betas[is[2]] > 0.0 &&
                     betas[is[1]] > 0.0);
          Load(u[0], states, {i, j});
          Load(u[1], states, neighbors[is[0]]);
          Load(u[2], states, neighbors[is[1]]);
          Load(u[3], states, neighbors[is[2]]);
          Load(u[4], states, neighbors[is[3]]);
          std::array<Coordinates<Rank>, 5> xM;
          xM[0] = GetAbsoluteVolumeCentroid(geom, {i, j}, dx);
          xM[1] = GetAbsoluteVolumeCentroid(geom, neighbors[is[0]], dx);
          xM[2] = GetAbsoluteVolumeCentroid(geom, neighbors[is[1]], dx);
          xM[3] = GetAbsoluteVolumeCentroid(geom, neighbors[is[2]], dx);
          xM[4] = GetAbsoluteVolumeCentroid(geom, neighbors[is[3]], dx);
          ComputeGradients(gradient, u, xM);
        } else if (betas[is[2]] > 0.0) {
          FUB_ASSERT(betas[is[0]] > 0.0 && betas[is[1]] > 0.0);
          Load(u[0], states, {i, j});
          Load(u[1], states, neighbors[is[0]]);
          Load(u[2], states, neighbors[is[1]]);
          Load(u[3], states, neighbors[is[2]]);
          std::array<Coordinates<Rank>, 4> xM;
          xM[0] = GetAbsoluteVolumeCentroid(geom, {i, j}, dx);
          xM[1] = GetAbsoluteVolumeCentroid(geom, neighbors[is[0]], dx);
          xM[2] = GetAbsoluteVolumeCentroid(geom, neighbors[is[1]], dx);
          xM[3] = GetAbsoluteVolumeCentroid(geom, neighbors[is[2]], dx);
          ComputeGradients(
              gradient, fub::span<const Complete>(u).template subspan<0, 4>(),
              xM);
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
          LimitGradients(grad_u, u.Subview(box), geom, dx);
        },
        gradient_x, gradient_y, states);
  }
}

template <typename GradientMethod>
void MdGradients<GradientMethod>::ComputeGradients(
    span<Gradient, 2> gradient, span<const Complete, 4> states,
    span<const Coordinates<Rank>, 4> x) {
  std::array<Gradient, 4> states_as_gradient_type;
  states_as_gradient_type.fill(Gradient(equation_));
  for (int i = 0; i < states.size(); ++i) {
    StateFromComplete(equation_, states_as_gradient_type[i], states[size_t(i)]);
  }
  ForEachComponent(
      [&](double& grad_x, double& grad_y, double uM, double u1, double u2,
          double u3) {
        std::array<double, 2> grad{0.0, 0.0};
        const std::array<double, 4> quantities{uM, u1, u2, u3};
        ComputeGradients(grad, quantities, x);
        grad_x = grad[0];
        grad_y = grad[1];
      },
      gradient[0], gradient[1], states_as_gradient_type[0],
      states_as_gradient_type[1], states_as_gradient_type[2],
      states_as_gradient_type[3]);
}

template <typename GradientMethod>
void MdGradients<GradientMethod>::ComputeGradients(
    span<Gradient, 2> gradient, span<const Complete, 5> states,
    span<const Coordinates<Rank>, 5> x) {
  std::array<Gradient, 4> states_as_gradient_type;
  states_as_gradient_type.fill(Gradient(equation_));
  for (int i = 0; i < states.size(); ++i) {
    StateFromComplete(equation_, states_as_gradient_type[i], states[size_t(i)]);
  }
  ForEachComponent(
      [&](double& grad_x, double& grad_y, double uM, double u1, double u2,
          double u3, double u4) {
        std::array<double, 2> grad{0.0, 0.0};
        const std::array<double, 5> quantities{uM, u1, u2, u3, u4};
        ComputeGradients(grad, quantities, x);
        grad_x = grad[0];
        grad_y = grad[1];
      },
      gradient[0], gradient[1], states_as_gradient_type[0],
      states_as_gradient_type[1], states_as_gradient_type[2],
      states_as_gradient_type[3], states_as_gradient_type[4]);
}

template <typename GradientMethod> 
void MdGradients<GradientMethod>::LimitGradients(
    const std::array<StridedDataView<double, Rank>, Rank>& grad_u,
    StridedDataView<const double, Rank> u,
    const CutCellData<Rank>& geom, const Coordinates<Rank>& dx) const {
  ForEachIndex(u.Box(), [&](auto... is) {
    if (geom.volume_fractions(is...) > 0.0) {
      limiter_.LimitGradientsAtIndex(grad_u, u, geom, {is...}, dx);
    }
  });
}

template <typename GradientMethod> 
void MdGradients<GradientMethod>::ComputeGradients(
    span<double, 2> gradient, span<const double, 4> states,
    span<const Coordinates<Rank>, 4> x) {
  Eigen::Matrix<double, 3, 2> A;
  Eigen::Matrix<double, 3, 1> b;
  for (int i = 0; i < 3; ++i) {
    A.row(i) = (x[i + 1] - x[0]).transpose();
    b[i] = states[i + 1] - states[0];
  }
  Eigen::Vector2d grads = A.colPivHouseholderQr().solve(b);
  gradient[0] = grads[0];
  gradient[1] = grads[1];
}

template <typename GradientMethod> 
void MdGradients<GradientMethod>::ComputeGradients(
    span<double, 2> gradient, span<const double, 5> states,
    span<const Coordinates<Rank>, 5> x) {
  Eigen::Matrix<double, 4, 2> A;
  Eigen::Matrix<double, 4, 1> b;
  for (int i = 0; i < 4; ++i) {
    A.row(i) = (x[i + 1] - x[0]).transpose();
    b[i] = states[i + 1] - states[0];
  }
  Eigen::Vector2d grads = A.colPivHouseholderQr().solve(b);
  gradient[0] = grads[0];
  gradient[1] = grads[1];
}


} // namespace fub

#endif