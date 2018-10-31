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

#ifndef FUB_EULER_IDEAL_GAS_HPP
#define FUB_EULER_IDEAL_GAS_HPP

#include "fub/core/assert.hpp"
#include "fub/core/span.hpp"

#include "fub/solver/euler/Reactor.hpp"

#include "SAMRAI/hier/Patch.h"
#include "SAMRAI/pdat/CellData.h"
#include "SAMRAI/pdat/FaceData.h"
#include "SAMRAI/tbox/Dimension.h"

#include <array>
#include <vector>

namespace fub {
namespace euler {

struct IdealGasInvalidName : public std::runtime_error {
  IdealGasInvalidName()
      : runtime_error("Initialized IdealGas with an invalid name.") {}
};

template <typename T, typename VectorT = span<T>, typename SpeciesT = span<T>>
struct CompleteIdealGasState {
  T density;
  VectorT momentum;
  T energy;
  T pressure;
  T temperature;
  T speed_of_sound;
  SpeciesT species;
};

template <typename T, typename VectorT = span<T>, typename SpeciesT = span<T>>
struct ConservativeIdealGasState {
  T density;
  VectorT momentum;
  T energy;
  SpeciesT species;
};

class IdealGas {
public:
  template <typename T, typename V = span<T>, typename S = span<T>>
  using Complete = CompleteIdealGasState<T, V, S>;

  template <typename T, typename V = span<T>, typename S = span<T>>
  using Conservative = ConservativeIdealGasState<T, V, S>;

  using CompleteState =
      Complete<SAMRAI::pdat::CellData<double>&, SAMRAI::pdat::CellData<double>&,
               SAMRAI::pdat::CellData<double>&>;

  using ConsState = Conservative<SAMRAI::pdat::CellData<double>&,
                                 SAMRAI::pdat::CellData<double>&,
                                 SAMRAI::pdat::CellData<double>&>;

  using FluxState = Conservative<SAMRAI::pdat::FaceData<double>&,
                                 SAMRAI::pdat::FaceData<double>&,
                                 SAMRAI::pdat::FaceData<double>&>;

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

  IdealGas(const std::string& name, const SAMRAI::tbox::Dimension& dim);
  IdealGas(std::string&& name, const SAMRAI::tbox::Dimension& dim);

  /////////////////////////////////////////////////////////////////////////////
  //                  Problem related member functions

  /// \brief Returns a Reactor object which is used to make thermodynamic
  /// computations.
  /// @{
  const FlameMasterReactor& GetReactor() const noexcept { return reactor_; }
  FlameMasterReactor& GetReactor() noexcept { return reactor_; }
  /// @}

  /// \brief Computes a complete state from a given conservative state.
  ///
  /// \note This mutates the internal reactor state.
  void Reconstruct(const CompleteState& complete, const ConsState& cons);

  /// \brief Computes a complete state from a given conservative state.
  ///
  /// \note This mutates the internal reactor state.
  void AdvanceSourceTerm(const CompleteState& complete, double time_step_size);

  /////////////////////////////////////////////////////////////////////////////
  //                  SAMRAI related member functions

  /// \brief Returns the name which used as a prefix for SAMRAI's
  /// VariableDatabase.
  const std::string& name() const noexcept { return name_; }

  /// \brief Returns the spatial dimension of the equations.
  const SAMRAI::tbox::Dimension& getDim() const noexcept { return dimension_; }

  /// \brief Returns the CellData handle for a given patch and variable.
  SAMRAI::pdat::CellData<double>& getCellData(const SAMRAI::hier::Patch& patch,
                                              Variable var) const {
    return *static_cast<SAMRAI::pdat::CellData<double>*>(
        patch.getPatchData(getDataId(var)).get());
  }

  /// \brief Returns the array of all patch data ids registered by this problem
  /// class.
  span<const int, variables_size> getDataIds() const noexcept {
    return data_ids_;
  }

  /// \brief Map a variable value to a patch data id from SAMRAI.
  int getDataId(Variable var) const {
    int val = static_cast<int>(var);
    FUB_ASSERT(0 <= val && val < variables_size);
    return data_ids_[val];
  }

private:
  std::string name_;
  SAMRAI::tbox::Dimension dimension_;
  std::array<int, variables_size> data_ids_;
  FlameMasterReactor reactor_;
};

} // namespace euler
} // namespace fub

#endif