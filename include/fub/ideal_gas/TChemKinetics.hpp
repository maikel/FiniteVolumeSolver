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

#ifndef FUB_IDEAL_GAS_TCHEM_KINETICS_HPP
#define FUB_IDEAL_GAS_TCHEM_KINETICS_HPP

#include "fub/core/assert.hpp"
#include "fub/core/function_ref.hpp"
#include "fub/core/span.hpp"

#include "fub/ideal_gas/IdealGasKinetics.hpp"
#include "fub/ideal_gas/TChemReactor.hpp"

#include "SAMRAI/hier/Patch.h"
#include "SAMRAI/pdat/CellData.h"
#include "SAMRAI/pdat/FaceData.h"
#include "SAMRAI/tbox/Dimension.h"

#include <array>
#include <string>

namespace fub {
namespace ideal_gas {

struct TChemKineticsOptions {
  std::string chemfile;
  std::string thermofile;
};

/// \ingroup IdealGas
class TChemKinetics : public IdealGasKinetics {
public:
  TChemKinetics(std::string prefix, SAMRAI::tbox::Dimension dim,
                const TChemMechanism& mechanism) noexcept;

  /// \brief Computes a complete state from a given conservative state.
  ///
  /// \note This mutates the internal reactor state.
  void FillFromCons(const CompleteStatePatchData& complete,
                    const ConsStatePatchData& cons) const override;

  /// \brief Computes a complete state from a given primitive state.
  ///
  /// \note This mutates the internal reactor state.
  void FillFromPrim(const CompleteStatePatchData& complete,
                    const PrimStatePatchData& prim) const override;

  /// \brief Computes a complete state from a given conservative state.
  ///
  /// \note This mutates the internal reactor state.
  void AdvanceSourceTerm(const CompleteStatePatchData& complete,
                         double dt) const override;

  /// \brief Returns the underlying Reactor object.
  const TChemReactor& GetReactor() const { return reactor_; }

private:
  mutable TChemReactor reactor_;
};

} // namespace ideal_gas
} // namespace fub

#endif