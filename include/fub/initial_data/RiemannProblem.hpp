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

#ifndef FUB_INITIAL_DATA_RIEMANN_PROBLEM_HPP
#define FUB_INITIAL_DATA_RIEMANN_PROBLEM_HPP

#include "fub/Equation.hpp"
#include "fub/ForEach.hpp"
#include "fub/grid/AMReX/cutcell/PatchHierarchy.hpp"

#include <array>
#include <limits>

namespace fub {
namespace amrex {
namespace cutcell {

template <typename Eq, typename Geometry> struct RiemannProblem {
  using Equation = Eq;
  using Complete = fub::Complete<Equation>;
  using Conservative = fub::Conservative<Equation>;
  static constexpr int Rank = Eq::Rank();

  RiemannProblem(std::shared_ptr<PatchHierarchy> hier, const Eq& eq,
                 const Geometry& geom, const Complete& l, const Complete& r)
      : hierarchy{std::move(hier)}, equation{eq}, geometry_{geom}, left{l},
        right{r} {
    ForEachComponent<Complete>(
        [](double& x) { x = std::numeric_limits<double>::quiet_NaN(); },
        boundary);
  }

  void InitializeData(const View<Complete>& data, PatchHandle patch) {
    const ::amrex::Geometry& geom = hierarchy->GetGeometry(patch.level);
    const ::amrex::Box& box = patch.iterator->tilebox();
    const auto& factory = hierarchy->GetPatchLevel(patch.level).factory;
    const auto& flags = factory->getMultiEBCellFlagFab()[*patch.iterator];
    ::amrex::FabType type = flags.getType(box);
    if (type == ::amrex::FabType::covered) {
      ForEachVariable<Complete>(
          [](const auto& var) {
            span<double> span = var.Span();
            std::fill(span.begin(), span.end(),
                      std::numeric_limits<double>::quiet_NaN());
          },
          data);
      return;
    }

    CartesianCoordinates x = GetCartesianCoordinates(geom, box);
    if (type == ::amrex::FabType::regular) {
      ForEachIndex(Box<0>(data), [&](auto... is) {
        if (geometry_(x(is...)) < 0.0) {
          Store(data, left, {is...});
        } else {
          Store(data, right, {is...});
        }
      });
    } else {
      CutCellData<Rank> eb = hierarchy->GetCutCellData(patch, Direction::X);
      const PatchDataView<const ::amrex::EBCellFlag, Rank>& flags = eb.flags;
      ForEachIndex(Box<0>(data), [&](auto... is) {
        if (flags(is...).isCovered()) {
          Store(data, boundary, {is...});
        } else if (geometry_(x(is...)) < 0.0) {
          Store(data, left, {is...});
        } else {
          Store(data, right, {is...});
        }
      });
    }
  }

  std::shared_ptr<PatchHierarchy> hierarchy;
  Equation equation;
  Geometry geometry_;
  Complete left{equation};
  Complete right{equation};
  Complete boundary{equation};
};

template <typename Eq, typename Geom>
RiemannProblem(std::shared_ptr<PatchHierarchy>, const Eq&, const Geom&,
               nodeduce_t<const Complete<Eq>&>, nodeduce_t<const Complete<Eq>&>)
    ->RiemannProblem<Eq, Geom>;

} // namespace cutcell
} // namespace amrex
} // namespace fub

#endif
