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

#ifndef FUB_IDEAL_GAS_IDEAL_GAS_EQUATION_HPP
#define FUB_IDEAL_GAS_IDEAL_GAS_EQUATION_HPP

#include "fub/core/assert.hpp"
#include "fub/core/mdspan.hpp"
#include "fub/core/span.hpp"

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

template <typename T, typename VectorT, typename SpeciesT>
struct CompleteIdealGasState {
  T density;
  VectorT momentum;
  T energy;
  T pressure;
  T temperature;
  T speed_of_sound;
  SpeciesT species;
};

template <typename T, typename VectorT, typename SpeciesT>
struct ConservativeIdealGasState {
  T density;
  VectorT momentum;
  T energy;
  SpeciesT species;
};

template <typename T, typename VectorT, typename SpeciesT>
struct PrimitiveIdealGasState {
  T temperature;
  VectorT momentum;
  T pressure;
  SpeciesT species;
};

span<double> MakeSpan(SAMRAI::pdat::CellData<double>& data) noexcept;

span<const double>
MakeSpan(const SAMRAI::pdat::CellData<double>& data) noexcept;

/// \ingroup IdealGas
class IdealGasEquation {
public:
  template <typename T, typename V = T, typename S = T>
  using Complete = CompleteIdealGasState<T, V, S>;

  template <typename T, typename V = T, typename S = T>
  using Conservative = ConservativeIdealGasState<T, V, S>;

  template <typename T, typename V = T, typename S = T>
  using Primitive = PrimitiveIdealGasState<T, V, S>;

  using CompleteStatePatchData = Complete<SAMRAI::pdat::CellData<double>&>;
  using ConsStatePatchData = Conservative<SAMRAI::pdat::CellData<double>&>;
  using FluxStatePatchData = Conservative<SAMRAI::pdat::FaceData<double>&>;
  using PrimStatePatchData = Primitive<SAMRAI::pdat::CellData<double>&>;

  template <typename T> using CompleteStateSpan = Complete<span<T>>;
  template <typename T> using ConsStateSpan = Conservative<span<T>>;
  template <typename T> using FluxStateSpan = Conservative<span<T>>;
  template <typename T> using PrimStateSpan = Primitive<span<T>>;

  enum class Variable : int {
    density,
    momentum,
    energy,
    pressure,
    temperature,
    speed_of_sound,
    species,
    size ///< This is only a dummy value to get the number of variables
  };

  static constexpr int variables_size = static_cast<int>(Variable::size);

  IdealGasEquation(std::string name, SAMRAI::tbox::Dimension dim,
                   int n_species);

  virtual ~IdealGasEquation() = default;

  /// \brief Computes a complete state from a given conservative variables.
  virtual void FillFromCons(const CompleteStatePatchData& q,
                            const ConsStatePatchData& u) const = 0;

  /// \brief Computes a complete state from a given primitive variables.
  virtual void FillFromPrim(const CompleteStatePatchData& q,
                            const PrimStatePatchData& w) const = 0;

  // virtual void SetPressureIsentropic(const CompleteStatePatchData& data,
  //                                    const SAMRAI::hier::Index& to,
  //                                    const SAMRAI::hier::Index& from) const =
  //                                    0;

  /////////////////////////////////////////////////////////////////////////////
  //                  SAMRAI related member functions

  /// \brief Returns the name which used as a prefix for SAMRAI's
  /// VariableDatabase.
  const std::string& GetName() const noexcept { return name_; }

  int GetNSpecies() const noexcept { return n_species_; }

  int GetNVariables() const noexcept { return variables_size + n_species_ - 1; }

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
  span<const int, variables_size> GetPatchDataIds() const noexcept {
    return patch_data_ids_;
  }

  /// \brief Map a variable value to a patch data id from SAMRAI.
  int GetPatchDataId(Variable var) const {
    int val = static_cast<int>(var);
    FUB_ASSERT(0 <= val && val < variables_size);
    return patch_data_ids_[val];
  }

  CompleteStatePatchData
  GetCompleteStatePatchData(const SAMRAI::hier::Patch& patch) const;

  ConsStatePatchData
  GetConsStatePatchData(const SAMRAI::hier::Patch& patch) const;

  PrimStatePatchData
  GetPrimStatePatchData(const SAMRAI::hier::Patch& patch) const;

private:
  std::string name_;
  SAMRAI::tbox::Dimension dimension_;
  int n_species_;
  std::array<int, variables_size> patch_data_ids_;
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