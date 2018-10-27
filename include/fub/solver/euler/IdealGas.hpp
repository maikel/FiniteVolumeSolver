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

#include "SAMRAI/hier/Patch.h"
#include "SAMRAI/pdat/CellData.h"
#include "SAMRAI/pdat/FaceData.h"
#include "SAMRAI/tbox/Dimension.h"

#include <array>

namespace fub {
namespace euler {

struct IdealGasInvalidName : public std::runtime_error {
  IdealGasInvalidName()
      : runtime_error("Initialized IdealGas with an invalid name.") {}
};

template <typename T> using Vector = std::array<T, SAMRAI_MAXIMUM_DIMENSION>;

template <typename T, typename VectorT = Vector<T>>
struct CompleteIdealGasState {
  T density;
  VectorT momentum;
  T energy;
  T pressure;
  T speed_of_sound;
};

template <typename T, typename VectorT = Vector<T>>
struct ConservativeIdealGasState {
  T density;
  VectorT momentum;
  T energy;
};

class IdealGas {
public:
  template <typename T, typename V = Vector<T>>
  using Complete = CompleteIdealGasState<T, V>;

  template <typename T, typename V = Vector<T>>
  using Conservative = ConservativeIdealGasState<T, V>;

  using CompleteState = Complete<SAMRAI::pdat::CellData<double>&,
                                 SAMRAI::pdat::CellData<double>&>;

  using ConsState = Conservative<SAMRAI::pdat::CellData<double>&,
                                 SAMRAI::pdat::CellData<double>&>;

  using FluxState = Conservative<SAMRAI::pdat::FaceData<double>&,
                                 SAMRAI::pdat::FaceData<double>&>;

  enum class Variable : int {
    density,
    momentum,
    energy,
    pressure,
    speed_of_sound,
    size // <-- This is only a dummy value to get the number of variables
  };

  static constexpr int variables_size = static_cast<int>(Variable::size);

  IdealGas(const std::string& name, const SAMRAI::tbox::Dimension& dim);
  IdealGas(std::string&& name, const SAMRAI::tbox::Dimension& dim);

  /////////////////////////////////////////////////////////////////////////////
  //                  Problem related member functions

  /// Returns the specific gas constant.
  static constexpr double getSpecificGasConstant() noexcept { return 1.0; }

  /// Returns the heat capacity ratio.
  static constexpr double getHeatCapacityRatio() noexcept { return 1.4; }

  /// Returns the heat capacity at constant volume.
  static constexpr double getCv() noexcept {
    // (i) R = cp - cv
    // (ii) gamma = cp / cv
    // R + cv = cp
    // gamma = (R + cv) / cv
    // cv gamma = R + cv
    // cv = R / (gamma - 1)
    constexpr double R = getSpecificGasConstant();
    constexpr double gamma = getHeatCapacityRatio();
    constexpr double cv = R / (gamma - 1.0);
    return cv;
  }

  /// Returns the heat capacity at constant pressure.
  static constexpr double getCp() noexcept {
    constexpr double gamma = getHeatCapacityRatio();
    constexpr double cv = getCv();
    return gamma * cv;
  }

  /// \brief Computes the pressure given a conservative state
  double computePressure(Conservative<double> q) const;

  /// \brief Computes the speed of sound given density and pressure.
  double computeSpeedOfSound(double density, double pressure) const;

  /////////////////////////////////////////////////////////////////////////////
  //                  SAMRAI related member functions

  const std::string& name() const noexcept { return name_; }

  const SAMRAI::tbox::Dimension& getDim() const noexcept { return dimension_; }

  SAMRAI::pdat::CellData<double>& getCellData(const SAMRAI::hier::Patch& patch,
                                              Variable var) const {
    return *static_cast<SAMRAI::pdat::CellData<double>*>(
        patch.getPatchData(getDataId(var)).get());
  }

  /// Returns the array of all patch data ids registered by this problem class.
  span<const int, variables_size> getDataIds() const noexcept {
    return data_ids_;
  }

  /// Map a variable value to a patch data id from SAMRAI.
  int getDataId(Variable var) const {
    int val = static_cast<int>(var);
    FUB_ASSERT(0 <= val && val < variables_size);
    return data_ids_[val];
  }

private:
  std::string name_;
  SAMRAI::tbox::Dimension dimension_;
  std::array<int, variables_size> data_ids_;
};

} // namespace euler
} // namespace fub

#endif