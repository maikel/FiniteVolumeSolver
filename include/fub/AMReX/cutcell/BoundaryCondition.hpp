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

#ifndef FUB_GRID_AMREX_CUTCELL_BOUNDARY_CONDITION_HPP
#define FUB_GRID_AMREX_CUTCELL_BOUNDARY_CONDITION_HPP

#include "fub/Direction.hpp"
#include "fub/Duration.hpp"

#include "fub/core/type_traits.hpp"

#include <AMReX_MultiFab.H>
#include <AMReX_PhysBCFunct.H>

namespace fub::amrex::cutcell {

class GriddingAlgorithm;

struct BoundaryConditionBase_ {
  virtual ~BoundaryConditionBase_() = default;
  virtual std::unique_ptr<BoundaryConditionBase_> Clone() const = 0;
  virtual void FillBoundary(::amrex::MultiFab& mf,
                            const ::amrex::Geometry& geom, Duration time_point,
                            const GriddingAlgorithm& gridding) = 0;
  virtual void FillBoundary(::amrex::MultiFab& mf,
                            const ::amrex::Geometry& geom, Duration time_point,
                            const GriddingAlgorithm& gridding,
                            Direction dir) = 0;
};

template <typename BC>
struct BoundaryConditionWrapper_ : public BoundaryConditionBase_ {
  BoundaryConditionWrapper_(const BC& bc) : boundary_condition_{bc} {}
  BoundaryConditionWrapper_(BC&& bc) : boundary_condition_{std::move(bc)} {}

  std::unique_ptr<BoundaryConditionBase_> Clone() const override;

  void FillBoundary(::amrex::MultiFab& mf, const ::amrex::Geometry& geom,
                    Duration time_point,
                    const GriddingAlgorithm& gridding) override;

  void FillBoundary(::amrex::MultiFab& mf, const ::amrex::Geometry& geom,
                    Duration time_point, const GriddingAlgorithm& gridding,
                    Direction dir) override;

  BC boundary_condition_;
};

class AnyBoundaryCondition : public ::amrex::PhysBCFunctBase {
public:
  AnyBoundaryCondition() = default;

  AnyBoundaryCondition(const AnyBoundaryCondition&);
  AnyBoundaryCondition& operator=(const AnyBoundaryCondition&);

  AnyBoundaryCondition(AnyBoundaryCondition&&) = default;
  AnyBoundaryCondition& operator=(AnyBoundaryCondition&&) = default;

  template <typename BC,
            typename = std::enable_if_t<!decays_to<BC, AnyBoundaryCondition>()>>
  AnyBoundaryCondition(BC&& bc);

  void FillBoundary(::amrex::MultiFab& mf, int, int, const ::amrex::IntVect&,
                    double time_point, int) override;

  void FillBoundary(::amrex::MultiFab& mf, const ::amrex::Geometry& geom,
                    Duration timepoint, const GriddingAlgorithm& gridding);

  void FillBoundary(::amrex::MultiFab& mf, const ::amrex::Geometry& geom,
                    Duration timepoint, const GriddingAlgorithm& gridding,
                    Direction dir);

  const GriddingAlgorithm* parent{};
  ::amrex::Geometry geometry{};

private:
  std::unique_ptr<BoundaryConditionBase_> boundary_condition_{};
};

// Implementation

template <typename BC>
std::unique_ptr<BoundaryConditionBase_>
BoundaryConditionWrapper_<BC>::Clone() const {
  return std::make_unique<BoundaryConditionWrapper_<BC>>(boundary_condition_);
}

template <typename BC>
void BoundaryConditionWrapper_<BC>::FillBoundary(
    ::amrex::MultiFab& mf, const ::amrex::Geometry& geom, Duration time_point,
    const GriddingAlgorithm& gridding) {
  boundary_condition_.FillBoundary(mf, geom, time_point, gridding);
}

template <typename BC>
void BoundaryConditionWrapper_<BC>::FillBoundary(
    ::amrex::MultiFab& mf, const ::amrex::Geometry& geom, Duration time_point,
    const GriddingAlgorithm& gridding, Direction dir) {
  boundary_condition_.FillBoundary(mf, geom, time_point, gridding, dir);
}

template <typename BC, typename>
AnyBoundaryCondition::AnyBoundaryCondition(BC&& bc)
    : boundary_condition_{
          std::make_unique<BoundaryConditionWrapper_<std::decay_t<BC>>>(
              std::forward<BC>(bc))} {}

inline void AnyBoundaryCondition::FillBoundary(::amrex::MultiFab& mf, int, int,
                                               const ::amrex::IntVect&,
                                               double time_point, int) {
  if (boundary_condition_ && parent) {
    boundary_condition_->FillBoundary(mf, geometry, Duration(time_point),
                                      *parent);
  }
}

inline void AnyBoundaryCondition::FillBoundary(
    ::amrex::MultiFab& mf, const ::amrex::Geometry& geom, Duration timepoint,
    const GriddingAlgorithm& gridding) {
  if (boundary_condition_) {
    boundary_condition_->FillBoundary(mf, geom, timepoint, gridding);
  }
}

inline void AnyBoundaryCondition::FillBoundary(
    ::amrex::MultiFab& mf, const ::amrex::Geometry& geom, Duration timepoint,
    const GriddingAlgorithm& gridding, Direction dir) {
  if (boundary_condition_) {
    boundary_condition_->FillBoundary(mf, geom, timepoint, gridding, dir);
  }
}

} // namespace fub::amrex::cutcell

#endif
