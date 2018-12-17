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

#include "fub/ideal_gas/boundary_condition/IsentropicExpansionBoundary.hpp"
#include "fub/SAMRAI/utility.hpp"

#include "SAMRAI/pdat/CellIndex.h"

namespace fub {
namespace ideal_gas {
namespace {
double GetEnthalpy_(const FlameMasterReactor& reactor) {
  const double U = reactor.getInternalEnergy();
  const double rho = reactor.getDensity();
  const double V = 1.0 / rho;
  const double p = reactor.getPressure();
  return U + p * V;
}
} // namespace

// Expansion into a fixed-pressure plenum
void IsentropicExpansionBoundary::SetPhysicalBoundaryCondition(
    const SAMRAI::hier::Patch& patch, const SAMRAI::hier::Box& fill_box,
    double /* fill_time */, Direction dir, int side) const {
  using Scratch = HyperbolicTimeIntegrator::Scratch;
  IdealGasEquation::CompleteStatePatchData complete =
      integrator_->getCompleteState(patch, Scratch(dir));
  const int direction = static_cast<int>(dir);
  const int sign = side == 0 ? +1 : -1;
  FlameMasterReactor& reactor = kinetics_->GetReactor();
  std::vector<double> fractions(complete.species.getDepth());
  for (const SAMRAI::hier::Index& index : fill_box) {
    SAMRAI::pdat::CellIndex to(index);
    SAMRAI::pdat::CellIndex from(shift(index, Direction(direction), sign));
    const double u_before = complete.momentum(from) / complete.density(from);
    if (mass_fractions_.size() == 0) {
      CopyMassFractions(fractions, complete.species, from);
    } else {
      std::copy(mass_fractions_.begin(), mass_fractions_.end(),
                fractions.begin());
    }
    reactor.setMassFractions(fractions);
    reactor.setTemperature(complete.temperature(from));
    reactor.setPressure(complete.pressure(from));
    const double h_before = GetEnthalpy_(reactor);
    reactor.setPressureIsentropic(pressure_);
    const double h_after = GetEnthalpy_(reactor);
    const double dH = h_after - h_before;
    const double u2_after = 2 * std::abs(dH) + u_before * u_before;
    const double u_after =
        dH > 0.0 ? -std::sqrt(u2_after) : std::sqrt(u2_after);
    kinetics_->UpdateCellFromReactor(complete, to, reactor, u_after);
  }
}

} // namespace ideal_gas
} // namespace fub