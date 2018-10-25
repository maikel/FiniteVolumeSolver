// Copyright (c) 2018 Maikel Nadolski
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

#ifndef FUB_SOLVER_EULER_IDEAL_GAS_SPLIT_INTEGRATOR_HPP
#define FUB_SOLVER_EULER_IDEAL_GAS_SPLIT_INTEGRATOR_HPP

#include "fub/core/span.hpp"
#include "fub/solver/BoundaryCondition.hpp"
#include "fub/solver/DimensionalSplitIntegrator.hpp"
#include "fub/solver/InitialCondition.hpp"
#include "fub/solver/euler/IdealGas.hpp"

#include "SAMRAI/hier/Patch.h"
#include "SAMRAI/pdat/CellData.h"
#include "SAMRAI/pdat/FaceData.h"

#include <array>
#include <memory>

namespace fub {
namespace euler {

class IdealGasSplitIntegrator : public fub::DimensionalSplitIntegrator {
public:
  using Variable = IdealGas::Variable;

  enum class FluxVariable : int { density, momentum, energy, size };

  static constexpr int variables_size = IdealGas::variables_size;

  static constexpr int flux_variables_size =
      static_cast<int>(FluxVariable::size);

  explicit IdealGasSplitIntegrator(const IdealGas& ideal_gas);

  const IdealGas& ideal_gas() const noexcept { return ideal_gas_; }

  void setInitialCondition(
      std::shared_ptr<const InitialCondition> initial_condition) {
    initial_condition_ = std::move(initial_condition);
  }

  void setBoundaryCondition(
      std::shared_ptr<const BoundaryCondition> boundary_condition) {
    boundary_condition_ = std::move(boundary_condition);
  }

  /////////////////////////////////////////////////////////////////////////////
  //                         Access SAMRAI data

  SAMRAI::pdat::CellData<double>&
  getScratchData(const SAMRAI::hier::Patch& patch, Variable var,
                 int dir) const {
    return *static_cast<SAMRAI::pdat::CellData<double>*>(
        patch.getPatchData(scratch(var, dir)).get());
  }

  SAMRAI::pdat::CellData<double>&
  getCurrentData(const SAMRAI::hier::Patch& patch, Variable var) const {
    return *static_cast<SAMRAI::pdat::CellData<double>*>(
        patch.getPatchData(current(var)).get());
  }

  SAMRAI::pdat::FaceData<double>& getFaceData(const SAMRAI::hier::Patch& patch,
                                              FluxVariable var) const {
    return *static_cast<SAMRAI::pdat::FaceData<double>*>(
        patch.getPatchData(face(var)).get());
  }

  span<const int, variables_size> current() const noexcept {
    return ideal_gas_.getDataIds();
  }

  int current(Variable var) const noexcept {
    const int val = static_cast<int>(var);
    FUB_ASSERT(0 <= val && val < variables_size);
    return current()[val];
  }

  span<const int> scratch() const noexcept {
    return span{scratch_}.first(ideal_gas_.getDim().getValue() *
                                variables_size);
  }

  span<const int, variables_size> scratch(int dir) const noexcept {
    FUB_ASSERT(0 <= dir && dir < 3);
    return {scratch_.data() + dir * variables_size, variables_size};
  }

  int scratch(Variable var, int dir) const noexcept {
    const int val = static_cast<int>(var);
    FUB_ASSERT(0 <= val && val < variables_size);
    return scratch(dir)[val];
  }

  span<const int, flux_variables_size> face() const noexcept { return face_; }

  int face(FluxVariable var) const noexcept {
    const int val = static_cast<int>(var);
    FUB_ASSERT(0 <= val && val < flux_variables_size);
    return face_[val];
  }

private:
  /////////////////////////////////////////////////////////////////////////////
  //          virtual overrides of DimensionalSplitIntergrator

  void initializePatchHierarchy(
      const std::shared_ptr<SAMRAI::hier::PatchHierarchy>& hierarchy,
      double time_point) const override;

  void doFillGhostLayer(
      const std::shared_ptr<SAMRAI::hier::PatchHierarchy>& hierarchy,
      double time_point, int dir) const override;

  double estimateStableTimeStepSize(const SAMRAI::hier::Patch& patch,
                                    double time_point) const override;

  void integratePatch(const SAMRAI::hier::Patch& patch, double time_point,
                      double time_step_size, int dir) const override;

  void coarsenFluxOnCoarseFineInteraces(
      const std::shared_ptr<SAMRAI::hier::PatchHierarchy>& hierarchy)
      const override;

  void fixPatchConservationOnCoarseFineBoundary(
      const SAMRAI::hier::Patch& patch,
      const SAMRAI::hier::BoundaryBox& boundary,
      double time_step_size) const override;

  void postIntegratePatchHierarchy(
      const std::shared_ptr<SAMRAI::hier::PatchHierarchy>& hierarchy,
      double time_point, double time_step_size, int dir) const override;

  IdealGas ideal_gas_;
  std::array<int, 3 * variables_size> scratch_;
  std::array<int, flux_variables_size> face_;
  std::shared_ptr<const InitialCondition> initial_condition_;
  std::shared_ptr<const BoundaryCondition> boundary_condition_;
};

} // namespace euler
} // namespace fub

#endif