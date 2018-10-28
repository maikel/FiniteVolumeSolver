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

#include "fub/solver/DimensionalSplitFluxMethod.hpp"
#include "fub/solver/DimensionalSplitTimeIntegrator.hpp"
#include "fub/solver/euler/IdealGas.hpp"

#include "SAMRAI/SAMRAI_config.h"

#include "SAMRAI/hier/Patch.h"
#include "SAMRAI/pdat/CellData.h"
#include "SAMRAI/pdat/FaceData.h"
#include "SAMRAI/xfer/RefineAlgorithm.h"

#include <array>
#include <memory>

namespace fub {
namespace euler {

class ForwardEulerTimeIntegrator : public fub::DimensionalSplitTimeIntegrator {
public:
  using Variable = IdealGas::Variable;
  using FluxState = IdealGas::FluxState;
  using ConsState = IdealGas::ConsState;
  using CompleteState = IdealGas::CompleteState;
  using FluxMethod = fub::DimensionalSplitFluxMethod<FluxState, CompleteState>;

  enum class FluxVariable : int { density, momentum, energy, size };
  enum class Scratch : int {};

  static constexpr int variables_size = static_cast<int>(Variable::size);
  static constexpr int flux_variables_size =
      static_cast<int>(FluxVariable::size);

  explicit ForwardEulerTimeIntegrator(const IdealGas& ideal_gas);

  const IdealGas& ideal_gas() const noexcept { return ideal_gas_; }

  /////////////////////////////////////////////////////////////////////////////
  /// \name Obtain Patch Data Id

  /// \brief Returns a patch data id which corresponds to a variable in the
  /// current context.
  int getPatchDataId(Variable variable) const {
    return ideal_gas_.getDataId(variable);
  }

  /// \brief Returns a patch data id which corresponds to a variable in a
  /// specified scratch context.
  int getPatchDataId(Variable variable, Scratch which_scratch) const {
    int dir = static_cast<int>(which_scratch);
    int var = static_cast<int>(variable);
    return scratch_[var + variables_size * dir];
  }

  /// \brief Returns a patch data id which corresponds to a flux variable.
  ///
  /// \note Flux variables exists only in one context.
  int getPatchDataId(FluxVariable variable) const {
    int var = static_cast<int>(variable);
    return face_[var];
  }

  /////////////////////////////////////////////////////////////////////////////
  /// \name Collect Patch Data 

  /// Returns a struct of CellDatas, one for each state variable.
  ///
  /// \param[in] patch The patch where CellData arrays live on.
  /// \param[in] Scratch A scratch value which specifies to which scratch
  ///                    context to take form.
  ///
  /// \return a CompleteState containing CellData to the specified scratch
  /// context.
  ///
  /// \throw Nothing.
  CompleteState getCompleteState(const SAMRAI::hier::Patch& patch,
                                 Scratch) const;

  /// Returns a struct of CellDatas, one for each state variable.
  ///
  /// \param[in] patch The patch where CellData arrays live on.
  ///
  /// \return a CompleteState containing CellData to the original context.
  ///
  /// \throw Nothing.
  CompleteState getCompleteState(const SAMRAI::hier::Patch& patch) const;

  /// Returns a struct of CellDatas, one for each conservative variable.
  ConsState getConsState(const SAMRAI::hier::Patch& patch, Scratch) const;

  ConsState getConsState(const SAMRAI::hier::Patch& patch) const;

  /// Returns a struct of FaceDatas, one for each flux variable.
  FluxState getFluxState(const SAMRAI::hier::Patch& patch) const;


  /// \name Virtual Overrides of DimensionalSplitTimeIntegrator

  void flagPatchDataIdsToAllocate(
      SAMRAI::hier::ComponentSelector& which_to_allocate) const override;

  double computeStableDtOnPatch(const SAMRAI::hier::Patch& patch,
                                double time_point) const override;

  void advanceTimeOnPatch(const SAMRAI::hier::Patch& patch, double time_point,
                          double time_step_size, Direction dir) const override;

  std::shared_ptr<SAMRAI::xfer::RefineAlgorithm>
  getFillGhostLayerRefineAlgorithm(Direction dir) const override;

  std::shared_ptr<SAMRAI::xfer::CoarsenAlgorithm>
  getFillGhostLayerCoarsenAlgorithm(Direction dir) const override;

  /////////////////////////////////////////////////////////////////////////////
  /// \name Class Member Variables

private:
  /// This is the object which stores patch data ids for state variables to some
  /// base context.
  IdealGas ideal_gas_;

  /// This array stores patch data ids in scratch contexts.
  /// We store variables in up to three scratch contexts, for each dimensional
  /// direction one scratch context. Scratch contexts have ghost cells only in
  /// one dimensional direction.
  std::array<int, SAMRAI_MAXIMUM_DIMENSION * variables_size> scratch_;

  /// This array stores patch data ids for the flux variables. This time
  /// integrator needs flux variables only once.
  std::array<int, flux_variables_size> face_;

  /// This array holds RefineAlgorithms which fill the ghost cell layer of
  /// patches. Data will be transfered from the base context to the a scratch
  /// context.
  std::array<std::shared_ptr<SAMRAI::xfer::RefineAlgorithm>,
             SAMRAI_MAXIMUM_DIMENSION>
      fill_ghost_layer_algorithm_;

  /// This is a polymorphic flux method object which we use to compute fluxes
  /// and to estimate a stable time step size.
  std::shared_ptr<const FluxMethod> flux_method_;
};

} // namespace euler
} // namespace fub

#endif