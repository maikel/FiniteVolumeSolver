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

#include "fub/Meta.hpp"

#include <memory>
#include <typeinfo>

namespace fub {

/// \ingroup GriddingAlgorithm
/// \defgroup BoundaryCondition Boundary Conditions
/// \brief This modules collects all components that fill the ghost layer of a
/// patch level.

namespace detail {
template <typename GriddingAlgorithm> struct BoundaryConditionBase {
  using DataReference = typename GridTraits<GriddingAlgorithm>::DataReference;

  virtual ~BoundaryConditionBase() = default;
  virtual std::unique_ptr<BoundaryConditionBase> Clone() const = 0;
  virtual void PreAdvanceHierarchy(const GriddingAlgorithm& grid) = 0;
  virtual void FillBoundary(DataReference mf, const GriddingAlgorithm& gridding,
                            int level) = 0;
  virtual void FillBoundary(DataReference mf, const GriddingAlgorithm& gridding,
                            int level, Direction dir) = 0;

  virtual const std::type_info& GetTypeInfo() const noexcept = 0;
  virtual void* GetPointer() noexcept = 0;
};
} // namespace detail

/// \ingroup BoundaryCondition PolymorphicValueType
///
/// \brief This is a polymorphic value type that wraps any BoundaryCondition
/// object.
template <typename GriddingAlgorithm> class AnyBoundaryCondition {
public:
  using DataReference = typename GridTraits<GriddingAlgorithm>::DataReference;

  /// @{
  /// \name Constructors

  /// \brief This constructs a method that does nothing on invocation.
  ///
  /// \throw Nothing.
  AnyBoundaryCondition() = default;

  /// \brief Stores any object which satisfies the BoundaryCondition concept.
  template <typename BC,
            typename = std::enable_if_t<!decays_to<BC, AnyBoundaryCondition>()>>
  AnyBoundaryCondition(BC&& bc);

  /// \brief Copies the implementation.
  AnyBoundaryCondition(const AnyBoundaryCondition& other);

  /// \brief Copies the implementation.
  AnyBoundaryCondition& operator=(const AnyBoundaryCondition& other);

  /// \brief Moves the `other` object without allocating and leaves an empty
  /// method.
  AnyBoundaryCondition(AnyBoundaryCondition&&) = default;

  /// \brief Moves the `other` object without allocating and leaves an empty
  /// method.
  AnyBoundaryCondition& operator=(AnyBoundaryCondition&&) = default;
  /// @}

  /// @{
  /// \name Actions

  /// \brief Perform some action on each time step if neccessary.
  ///
  /// This will do nothing if the boundary condition does not require to do
  /// anything.
  void PreAdvanceHierarchy(const GriddingAlgorithm& grid);

  /// \brief Fill the boundary layer of data.
  void FillBoundary(DataReference data, const GriddingAlgorithm& gridding,
                    int level);

  /// \brief Fill the boundary layer of data in direction `dir` only.
  void FillBoundary(DataReference data, const GriddingAlgorithm& gridding,
                    int level, Direction dir);
  /// @}

  /// \brief Try to cast this boundary condition to type BC.
  ///
  /// If the stored boundary condition is not of type BC then a nullptr is returned.
  template <typename BC>
  BC* Cast() {
    if (boundary_condition_) {
      const std::type_info& this_bc_type = boundary_condition_->GetTypeInfo();
      if (this_bc_type == typeid(BC)) {
        return static_cast<BC*>(boundary_condition_->GetPointer());
      }
    }
    return nullptr;
  }

private:
  std::unique_ptr<detail::BoundaryConditionBase<GriddingAlgorithm>>
      boundary_condition_{};
};

inline int GetSign(int side) { return (side == 0) - (side == 1); }

///////////////////////////////////////////////////////////////////////////////
//                                                              Implementation

namespace detail {
template <typename GriddingAlgorithm, typename BC>
struct BoundaryConditionWrapper
    : public BoundaryConditionBase<GriddingAlgorithm> {
  using DataReference = typename GridTraits<GriddingAlgorithm>::DataReference;

  BoundaryConditionWrapper(const BC& bc) : boundary_condition_{bc} {}
  BoundaryConditionWrapper(BC&& bc) : boundary_condition_{std::move(bc)} {}

  std::unique_ptr<BoundaryConditionBase<GriddingAlgorithm>>
  Clone() const override {
    return std::make_unique<BoundaryConditionWrapper<GriddingAlgorithm, BC>>(
        boundary_condition_);
  }

  void PreAdvanceHierarchy(const GriddingAlgorithm& grid) {
    if constexpr (is_detected<meta::PreAdvanceHierarchy, BC&,
                              const GriddingAlgorithm&>::value) {
      boundary_condition_.PreAdvanceHierarchy(grid);
    }
  }

  void FillBoundary(DataReference data, const GriddingAlgorithm& gridding,
                    int level) override {
    boundary_condition_.FillBoundary(data, gridding, level);
  }

  void FillBoundary(DataReference data, const GriddingAlgorithm& gridding,
                    int level, Direction dir) override {
    boundary_condition_.FillBoundary(data, gridding, level, dir);
  }

  const std::type_info& GetTypeInfo() const noexcept override {
    return typeid(BC);
  }

  void* GetPointer() noexcept override {
    return std::addressof(boundary_condition_);
  }

  BC boundary_condition_;
};
} // namespace detail

template <typename GriddingAlgorithm>
AnyBoundaryCondition<GriddingAlgorithm>::AnyBoundaryCondition(
    const AnyBoundaryCondition& other)
    : boundary_condition_(other.boundary_condition_
                              ? other.boundary_condition_->Clone()
                              : nullptr) {}

template <typename GriddingAlgorithm>
AnyBoundaryCondition<GriddingAlgorithm>&
AnyBoundaryCondition<GriddingAlgorithm>::
operator=(const AnyBoundaryCondition& other) {
  AnyBoundaryCondition tmp(other);
  return *this = std::move(tmp);
}

template <typename GriddingAlgorithm>
template <typename BC, typename>
AnyBoundaryCondition<GriddingAlgorithm>::AnyBoundaryCondition(BC&& bc)
    : boundary_condition_{std::make_unique<detail::BoundaryConditionWrapper<
          GriddingAlgorithm, std::decay_t<BC>>>(std::forward<BC>(bc))} {}

template <typename GriddingAlgorithm>
void AnyBoundaryCondition<GriddingAlgorithm>::PreAdvanceHierarchy(
    const GriddingAlgorithm& grid) {
  if (boundary_condition_) {
    boundary_condition_->PreAdvanceHierarchy(grid);
  }
}

template <typename GriddingAlgorithm>
void AnyBoundaryCondition<GriddingAlgorithm>::FillBoundary(
    DataReference data, const GriddingAlgorithm& gridding, int level) {
  if (boundary_condition_) {
    boundary_condition_->FillBoundary(data, gridding, level);
  }
}

template <typename GriddingAlgorithm>
void AnyBoundaryCondition<GriddingAlgorithm>::FillBoundary(
    DataReference data, const GriddingAlgorithm& gridding, int level,
    Direction dir) {
  if (boundary_condition_) {
    boundary_condition_->FillBoundary(data, gridding, level, dir);
  }
}

} // namespace fub

#endif
