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

/// \defgroup Euler Euler Equations
/// This module defines classes and methods to be used to concrete solver
/// implementations wrt the euler equations.

#ifndef FUB_SAMRAI_IDEAL_GAS_IDEAL_GAS_EQUATION_HPP
#define FUB_SAMRAI_IDEAL_GAS_IDEAL_GAS_EQUATION_HPP

#include "fub/core/assert.hpp"
#include "fub/core/mdspan.hpp"
#include "fub/core/span.hpp"

#include "fub/ideal_gas/IdealGasEquation.hpp"

#include "SAMRAI/hier/Patch.h"
#include "SAMRAI/hier/PatchHierarchy.h"
#include "SAMRAI/pdat/CellData.h"
#include "SAMRAI/pdat/FaceData.h"
#include "SAMRAI/tbox/Dimension.h"

#include <boost/optional.hpp>

#include <array>

namespace fub {
namespace ideal_gas {

struct IdealGasInvalidName : public std::runtime_error {
  IdealGasInvalidName()
      : runtime_error("Initialized IdealGas with an invalid name.") {}
};

span<double> MakeSpan(SAMRAI::pdat::CellData<double>& data) noexcept;

span<const double>
MakeSpan(const SAMRAI::pdat::CellData<double>& data) noexcept;

/// \ingroup IdealGas
///
/// \brief This class describes what data lives on SAMRAI patches and helps to
/// access their data arrays.
class IdealGasEquation : public IdealGasEquation {
public:
  template <typename T, typename V = T, typename S = V>
  using Complete = IdealGasEquation::CompleteState<T, V, T, T, T, T, S>;

  template <typename T, typename V = T, typename S = V>
  using Conservative = IdealGasEquation::ConservativeState<T, V, T, S>;

  template <typename T, typename V = T, typename S = V>
  using Primitive = IdealGasEquation::PrimitiveState<T, V, T, S>;

  using CompletePatchData = Complete<SAMRAI::pdat::CellData<double>&>;
  using ConsPatchData = Conservative<SAMRAI::pdat::CellData<double>&>;
  using FluxPatchData = Conservative<SAMRAI::pdat::FaceData<double>&>;
  using PrimPatchData = Primitive<SAMRAI::pdat::CellData<double>&>;

  template <typename T> using CompleteStateSpan = Complete<span<T>>;
  template <typename T> using ConsStateSpan = Conservative<span<T>>;
  template <typename T> using FluxStateSpan = Conservative<span<T>>;
  template <typename T> using PrimStateSpan = Primitive<span<T>>;

  template <typename T>
  using CompleteSpans = Complete<DynamicMdSpan<T, 3, layout_left>,
                                 DynamicMdSpan<T, 4, layout_left>>;
  template <typename T>
  using FluxSpans = Conservative<DynamicMdSpan<T, 3, layout_left>,
                                 DynamicMdSpan<T, 4, layout_left>>;

  static constexpr int kVariablesSize = static_cast<int>(Variable::size);

  IdealGasEquation(std::string name, SAMRAI::tbox::Dimension dim,
                   int n_species);

  virtual ~IdealGasEquation() = default;

  /// \brief Computes a complete state from a given conservative variables.
  virtual void FillFromCons(const CompletePatchData& q,
                            const ConsPatchData& u) const = 0;

  /// \brief Computes a complete state from a given primitive variables.
  virtual void FillFromPrim(const CompletePatchData& q,
                            const PrimPatchData& w) const = 0;

  // virtual void SetPressureIsentropic(const CompletePatchData& data,
  //                                    const SAMRAI::hier::Index& to,
  //                                    const SAMRAI::hier::Index& from) const =
  //                                    0;

  /////////////////////////////////////////////////////////////////////////////
  //                  SAMRAI related member functions

  /// \brief Returns the name which used as a prefix for SAMRAI's
  /// VariableDatabase.
  const std::string& GetName() const noexcept { return name_; }

  int GetNSpecies() const noexcept { return n_species_; }

  int GetNVariables() const noexcept { return kVariablesSize + n_species_ - 1; }

  /// \brief Returns the spatial dimension of the equations.
  const SAMRAI::tbox::Dimension& GetDimension() const noexcept {
    return dimension_;
  }

  /// \brief Returns the CellData handle for a given patch and variable.
  SAMRAI::pdat::CellData<double>& GetCellData(const SAMRAI::hier::Patch& patch,
                                              Variable var) const {
    return *static_cast<SAMRAI::pdat::CellData<double>*>(
        patch.getPatchData(GetPatchDataId(var)).get());
  }

  /// \brief Returns the array of all patch data ids registered by this problem
  /// class.
  span<const int, kVariablesSize> GetPatchDataIds() const noexcept {
    return patch_data_ids_;
  }

  /// \brief Map a variable value to a patch data id from SAMRAI.
  int GetPatchDataId(Variable var) const {
    int val = static_cast<int>(var);
    FUB_ASSERT(0 <= val && val < kVariablesSize);
    return patch_data_ids_[val];
  }

  //////////////////////////////////////////////////////////////////////////////
  // Access Patch Data on SAMRAI::hier::Patch

  CompletePatchData
  GetCompletePatchData(const SAMRAI::hier::Patch& patch) const;

  ConsPatchData GetConsPatchData(const SAMRAI::hier::Patch& patch) const;

  PrimPatchData GetPrimPatchData(const SAMRAI::hier::Patch& patch) const;

private:
  std::string name_;
  SAMRAI::tbox::Dimension dimension_;
  int n_species_;
  std::array<int, kVariablesSize> patch_data_ids_;
};

void CopyMassFractions(span<double> fractions,
                       const SAMRAI::pdat::CellData<double>& species,
                       const SAMRAI::pdat::CellIndex& index);

class IdealGasGrid {
public:
  using Variable = IdealGasEquation::Variable;

  IdealGasGrid(DynamicMdSpan<double, 2> buffer, double dx)
      : mdspan_(buffer), dx_{dx} {}

  double& operator()(Variable var, int cell) {
    return mdspan_(static_cast<int>(var), cell);
  }
  const double& operator()(Variable var, int cell) const {
    return mdspan_(static_cast<int>(var), cell);
  }

  DynamicMdSpan<double, 2> Mdspan() noexcept { return mdspan_; }
  DynamicMdSpan<const double, 2> Mdspan() const noexcept { return mdspan_; }

  DynamicExtents<2> Extents() const noexcept { return mdspan_.extents(); }
  int Extent(int dim) const { return mdspan_.extent(dim); }

  double Dx() const noexcept { return dx_; }

private:
  DynamicMdSpan<double, 2> mdspan_;
  double dx_;
};

boost::optional<IdealGasGrid>
GatherGrid(const SAMRAI::hier::PatchHierarchy& hier,
           const IdealGasEquation& equation, span<double> buffer, int root = 0);

IdealGasGrid AllGatherGrid(const SAMRAI::hier::PatchHierarchy& hier,
                           const IdealGasEquation& equation,
                           span<double> buffer, int root = 0);

} // namespace ideal_gas
} // namespace fub

#endif