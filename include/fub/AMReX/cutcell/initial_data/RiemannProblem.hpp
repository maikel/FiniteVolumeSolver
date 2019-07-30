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

#include "fub/AMReX/cutcell/PatchHierarchy.hpp"
#include "fub/Equation.hpp"
#include "fub/ForEach.hpp"

#include <AMReX.H>

#include <array>
#include <limits>

namespace fub::amrex::cutcell {

template <typename Eq, typename Geometry> struct RiemannProblem {
  using Equation = Eq;
  using Complete = fub::Complete<Equation>;
  using Conservative = fub::Conservative<Equation>;
  static constexpr int Rank = Eq::Rank();

  RiemannProblem(const Eq& eq, const Geometry& geom, const Complete& l,
                 const Complete& r);

  void InitializeData(::amrex::MultiFab& data, const ::amrex::Geometry& geom);

  Equation equation_;
  Geometry geometry_;
  Complete left_{equation_};
  Complete right_{equation_};
  Complete boundary_{equation_};
};

template <typename Eq, typename Geom>
RiemannProblem(const Eq&, const Geom&, nodeduce_t<const Complete<Eq>&>,
               nodeduce_t<const Complete<Eq>&>)
    ->RiemannProblem<Eq, Geom>;

template <typename Eq, typename Geometry>
RiemannProblem<Eq, Geometry>::RiemannProblem(const Eq& eq, const Geometry& geom,
                                             const Complete& l,
                                             const Complete& r)
    : equation_{eq}, geometry_{geom}, left_{l}, right_{r} {
  ForEachComponent([](double& x) { x = 0.0; }, boundary_);
}

template <typename Eq, typename Geometry>
void RiemannProblem<Eq, Geometry>::InitializeData(
    ::amrex::MultiFab& data, const ::amrex::Geometry& geom) {
  const ::amrex::EBFArrayBoxFactory& factory =
      dynamic_cast<const ::amrex::EBFArrayBoxFactory&>(data.Factory());
  const ::amrex::FabArray<::amrex::EBCellFlagFab>& flags =
      factory.getMultiEBCellFlagFab();
  ForEachFab(data, [&](const ::amrex::MFIter& mfi) {
    const ::amrex::FabType type = flags[mfi].getType();
    if (type == ::amrex::FabType::covered) {
      span<double> span(data[mfi].dataPtr(), data[mfi].size());
      std::fill(span.begin(), span.end(), 0.0);
      return;
    }
    View<Complete> states =
        MakeView<Complete>(data[mfi], equation_, mfi.tilebox());
    Eigen::Vector3d x = Eigen::Vector3d::Zero();
    if (type == ::amrex::FabType::regular) {
      ForEachIndex(Box<0>(states), [&](auto... is) {
        geom.CellCenter({int(is)...}, x.data());
        if (geometry_(x) < 0.0) {
          Store(states, left_, {is...});
        } else {
          Store(states, right_, {is...});
        }
      });
    } else {
      ForEachIndex(Box<0>(states), [&](auto... is) {
        geom.CellCenter({int(is)...}, x.data());
        if (flags[mfi]({int(is)...}).isCovered()) {
          Store(states, boundary_, {is...});
        } else if (geometry_(x) < 0.0) {
          Store(states, left_, {is...});
        } else {
          Store(states, right_, {is...});
        }
      });
    }
  });
}

} // namespace fub::amrex::cutcell

#endif
