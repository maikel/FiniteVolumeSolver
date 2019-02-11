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

#ifndef FUB_SAMRAI_EULER_IDEAL_GAS_SPLIT_INTEGRATOR_HPP
#define FUB_SAMRAI_EULER_IDEAL_GAS_SPLIT_INTEGRATOR_HPP

#include "fub/SAMRAI/ideal_gas/IdealGasEquation.hpp"
#include "fub/SAMRAI/DimensionalSplitFluxMethod.hpp"
#include "fub/SAMRAI/DimensionalSplitTimeIntegrator.hpp"

#include "SAMRAI/SAMRAI_config.h"

#include "SAMRAI/hier/Patch.h"
#include "SAMRAI/pdat/CellData.h"
#include "SAMRAI/pdat/FaceData.h"
#include "SAMRAI/xfer/RefineAlgorithm.h"

#include <array>
#include <memory>

namespace fub {
namespace ideal_gas {
/// \ingroup Euler
class DimensionalSplitTimeIntegrator : public fub::DimensionalSplitTimeIntegrator {
public:
  using Variable = IdealGasEquation::Variable;
  using FluxPatchData = IdealGasEquation::FluxPatchData;
  using ConsPatchData = IdealGasEquation::ConsPatchData;
  using CompleteSpans = IdealGasEquation::CompleteSpans<double>;
  using FluxSpans = IdealGasEquation::FluxSpans<double>;
  using CompletePatchData = IdealGasEquation::CompletePatchData;
  using FluxMethod = fub::DimensionalSplitFluxMethod<FluxPatchData,
                                                     CompletePatchData>;

  enum class FluxVariable : int { density, momentum, energy, species, size };
  enum class Scratch : int { X, Y, Z };

  static constexpr int kVariablesSize = static_cast<int>(Variable::size);
  static constexpr int kFluxVariablesSize =
      static_cast<int>(FluxVariable::size);

  explicit DimensionalSplitTimeIntegrator(
      std::shared_ptr<const IdealGasEquation> ideal_gas);

  /// \brief Returns the equation of state object.
  const std::shared_ptr<const IdealGasEquation>& GetEquation() const noexcept {
    return ideal_gas_;
  }

  /////////////////////////////////////////////////////////////////////////////
  /// \name Obtain Patch Data Id

  /// \brief Returns a patch data id which corresponds to a variable in the
  /// current context.
  int GetPatchDataId(Variable variable) const {
    return ideal_gas_->GetPatchDataId(variable);
  }

  /// \brief Returns a patch data id which corresponds to a variable in a
  /// specified scratch context.
  int GetPatchDataId(Variable variable, Scratch which_scratch) const {
    int dir = static_cast<int>(which_scratch);
    int var = static_cast<int>(variable);
    return scratch_[var + kVariablesSize * dir];
  }

  /// \brief Returns a patch data id which corresponds to a flux variable.
  ///
  /// \param[in] variable  The flux variable which data id to return.
  ///
  /// \return An integer which is a patch data index.
  ///
  /// \require variable has to be a named value.
  ///
  /// \note Flux variables exists only in one context.
  int GetPatchDataId(FluxVariable variable) const {
    int var = static_cast<int>(variable);
    FUB_ASSERT(0 <= var && var < kFluxVariablesSize);
    return face_[var];
  }

  span<const int> GetScratch() const {
    return make_span(scratch_).first(kVariablesSize *
                                     ideal_gas_->GetDimension().getValue());
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
  CompletePatchData GetCompletePatchData(const SAMRAI::hier::Patch& patch,
                                          Scratch) const;

  /// Returns a struct of CellDatas, one for each state variable.
  ///
  /// \param[in] patch The patch where CellData arrays live on.
  ///
  /// \return a CompleteState containing CellData to the original context.
  ///
  /// \throw Nothing.
  CompletePatchData
  GetCompletePatchData(const SAMRAI::hier::Patch& patch) const;

  CompleteSpans GetCompleteSpans(const SAMRAI::hier::Patch& patch) const;
  CompleteSpans GetCompleteSpans(const SAMRAI::hier::Patch& patch, Scratch) const;

  /// Returns a struct of CellDatas, one for each conservative variable.
  ConsPatchData GetConsPatchData(const SAMRAI::hier::Patch& patch,
                                  Scratch) const;

  ConsPatchData GetConsPatchData(const SAMRAI::hier::Patch& patch) const;

  /// Returns a struct of FaceDatas, one for each flux variable.
  FluxPatchData GetFluxPatchData(const SAMRAI::hier::Patch& patch) const;

  FluxSpans GetFluxSpans(const SAMRAI::hier::Patch& patch, Direction face_normal) const;

  /// \name Virtual Overrides of DimensionalSplitTimeIntegrator

  void FlagPatchDataIdsToAllocate(
      SAMRAI::hier::ComponentSelector& which_to_allocate) const override;

  double ComputeStableDtOnPatch(const SAMRAI::hier::Patch& patch,
                                double time_point,
                                Direction dir) const override;

  void AdvanceTimeOnPatch(const SAMRAI::hier::Patch& patch, double time_point,
                          double time_step_size, Direction dir) const override;

  std::shared_ptr<SAMRAI::xfer::RefineAlgorithm>
  GetFillGhostLayerRefineAlgorithm(Direction dir) const override;

  std::shared_ptr<SAMRAI::xfer::CoarsenAlgorithm>
  GetFillGhostLayerCoarsenAlgorithm(Direction dir) const override;

  /////////////////////////////////////////////////////////////////////////////
  /// \name Class Member Variables

private:
  /// This is the object which stores patch data ids for state variables to some
  /// base context.
  std::shared_ptr<const IdealGasEquation> ideal_gas_;

  /// This array stores patch data ids in scratch contexts.
  /// We store variables in up to three scratch contexts, for each dimensional
  /// direction one scratch context. Scratch contexts have ghost cells only in
  /// one dimensional direction.
  std::array<int, SAMRAI_MAXIMUM_DIMENSION * kVariablesSize> scratch_;

  /// This array stores patch data ids for the flux variables. This time
  /// integrator needs flux variables only once.
  std::array<int, kFluxVariablesSize> face_;

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

} // namespace ideal_gas
} // namespace fub

#endif