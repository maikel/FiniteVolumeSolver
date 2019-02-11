// Copyright (c) 2018-2019 Maikel Nadolski
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

#ifndef FUB_SAMRAI_IDEAL_GAS_HLLE_GODUNOV_METHOD_HPP
#define FUB_SAMRAI_IDEAL_GAS_HLLE_GODUNOV_METHOD_HPP

#include "fub/SAMRAI/DimensionalSplitFluxMethod.hpp"
#include "fub/SAMRAI/ideal_gas/FlameMasterKinetics.hpp"

namespace fub {
namespace ideal_gas {
/// \ingroup IdealGas
/// This class is a first order flux method for the IdealGas problem which uses
/// the approximative HLL riemann problem solver with Einfeldts signals
/// velocities.
class MusclHancockMethod
    : public DimensionalSplitFluxMethod<IdealGasEquation::FluxPatchData,
                                        IdealGasEquation::CompletePatchData> {
public:
  using FluxPatchData = IdealGasEquation::FluxPatchData;
  using CompletePatchData = IdealGasEquation::CompletePatchData;
  using Variable = IdealGasEquation::Variable;

  explicit MusclHancockMethod(
      std::shared_ptr<const FlameMasterKinetics> equation)
      : equation_{std::move(equation)} {}

  /// \copydoc fub::DimensionalSplitFluxMethod::ComputeStableDtOnPatch
  double ComputeStableDtOnPatch(const CompletePatchData& state,
                                const SAMRAI::hier::Patch& patch,
                                Direction dir) const override;

  /// \copydoc fub::DimensionalSplitFluxMethod::ComputeFluxesOnPatch
  void ComputeFluxesOnPatch(const FluxPatchData& fluxes,
                            const CompletePatchData& state,
                            const SAMRAI::hier::Patch& patch, double dt,
                            Direction dir) const override;

  /// \copydoc fub::DimensionalSplitFluxMethod::GetStencilWidth
  ///
  /// \Return '2' in all directions.
  SAMRAI::hier::IntVector
  GetStencilWidth(const SAMRAI::tbox::Dimension& dim) const override;

  CompletePatchData
  GetHalfTimePatchData(const SAMRAI::hier::Patch& patch) const;

  SAMRAI::pdat::CellData<double>&
  GetHalfTimeCellData(const SAMRAI::hier::Patch& patch,
                      Variable variable) const;

private:
  std::shared_ptr<const FlameMasterKinetics> equation_;
};
} // namespace ideal_gas
} // namespace fub

#endif