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

#ifndef FUB_GRID_AMREX_BOUNDARY_CONDITION_HPP
#define FUB_GRID_AMREX_BOUNDARY_CONDITION_HPP

#include "fub/Direction.hpp"
#include "fub/Duration.hpp"

#include "fub/core/type_traits.hpp"

#include <AMReX_MultiFab.H>
#include <AMReX_PhysBCFunct.H>

namespace fub {
namespace amrex {

inline int GetSign(int side) { return (side == 0) - (side == 1); }

class GriddingAlgorithm;

struct BoundaryConditionBase {
  virtual ~BoundaryConditionBase() = default;
  virtual std::unique_ptr<BoundaryConditionBase> Clone() const = 0;
  virtual void FillBoundary(::amrex::MultiFab& mf,
                            const ::amrex::Geometry& geom, Duration time_point,
                            const GriddingAlgorithm& gridding) = 0;
};

template <typename BC>
struct BoundaryConditionWrapper : public BoundaryConditionBase {
  BoundaryConditionWrapper(const BC& bc) : boundary_condition_{bc} {}
  BoundaryConditionWrapper(BC&& bc) : boundary_condition_{std::move(bc)} {}

  std::unique_ptr<BoundaryConditionBase> Clone() const override;

  void FillBoundary(::amrex::MultiFab& mf, const ::amrex::Geometry& geom,
                    Duration time_point,
                    const GriddingAlgorithm& gridding) override;

  BC boundary_condition_;
};

class BoundaryCondition : public ::amrex::PhysBCFunctBase {
public:
  BoundaryCondition() = default;

  BoundaryCondition(const BoundaryCondition&);
  BoundaryCondition& operator=(const BoundaryCondition&);

  BoundaryCondition(BoundaryCondition&&) = default;
  BoundaryCondition& operator=(BoundaryCondition&&) = default;

  template <typename BC,
            typename = std::enable_if_t<!decays_to<BC, BoundaryCondition>()>>
  BoundaryCondition(BC&& bc);

  void FillBoundary(::amrex::MultiFab& mf, int, int, double time_point,
                    int) override;

  void FillBoundary(::amrex::MultiFab& mf, const ::amrex::Geometry& geom,
                    Duration timepoint, const GriddingAlgorithm& gridding);

  void FillBoundary(::amrex::MultiFab& mf, const ::amrex::Geometry& geom,
                    Duration timepoint, const GriddingAlgorithm& gridding, Direction dir);

  const GriddingAlgorithm* parent{};
  ::amrex::Geometry geometry{};

private:
  std::unique_ptr<BoundaryConditionBase> boundary_condition_{};
};

// Implementation

template <typename BC>
std::unique_ptr<BoundaryConditionBase>
BoundaryConditionWrapper<BC>::Clone() const {
  return std::make_unique<BoundaryConditionWrapper<BC>>(boundary_condition_);
}

template <typename BC>
void BoundaryConditionWrapper<BC>::FillBoundary(
    ::amrex::MultiFab& mf, const ::amrex::Geometry& geom, Duration time_point,
    const GriddingAlgorithm& gridding) {
  boundary_condition_.FillBoundary(mf, geom, time_point, gridding);
}

template <typename BC, typename>
BoundaryCondition::BoundaryCondition(BC&& bc)
    : boundary_condition_{
          std::make_unique<BoundaryConditionWrapper<std::decay_t<BC>>>(
              std::forward<BC>(bc))} {}

inline void BoundaryCondition::FillBoundary(::amrex::MultiFab& mf, int, int,
                                            double time_point, int) {
  if (boundary_condition_ && parent) {
    boundary_condition_->FillBoundary(mf, geometry, Duration(time_point),
                                      *parent);
  }
}

inline void BoundaryCondition::FillBoundary(::amrex::MultiFab& mf,
                                            const ::amrex::Geometry& geom,
                                            Duration timepoint,
                                            const GriddingAlgorithm& gridding) {
  if (boundary_condition_) {
    boundary_condition_->FillBoundary(mf, geom, timepoint, gridding);
  }
}

} // namespace amrex
} // namespace fub

#endif
