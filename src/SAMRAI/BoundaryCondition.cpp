//
// Created by Patrick Denzler on 2019-07-30.
//

#include "fub/SAMRAI/BoundaryCondition.hpp"
#include "fub/SAMRAI/ViewFArrayBox.hpp"

namespace fub {
namespace samrai {

BoundaryCondition::BoundaryCondition(const BoundaryCondition& other)
    : geometry{other.geometry}, boundary_condition_{} {
  if (other.boundary_condition_) {
    boundary_condition_ = other.boundary_condition_->Clone();
  }
}

BoundaryCondition& BoundaryCondition::
operator=(const BoundaryCondition& other) {
  BoundaryCondition tmp{other};
  return (*this = std::move(tmp));
}

} // namespace amrex
} // namespace fub
