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

#include "fub/ideal_gas/TChemKinetics.hpp"

#include <algorithm>
#include <numeric>

namespace fub {
namespace ideal_gas {
TChemKinetics::TChemKinetics(std::string prefix, SAMRAI::tbox::Dimension dim,
                             const TChemMechanism& mechanism) noexcept
    : IdealGasKinetics(prefix, dim, mechanism.getNSpecies()),
      reactor_(mechanism)
  {}

void TChemKinetics::FillFromCons(const CompleteStatePatchData& q,
                                 const ConsStatePatchData& cons) const {
  q.density.copy(cons.density);
  q.momentum.copy(cons.momentum);
  q.energy.copy(cons.energy);
  q.species.copy(cons.species);
  SAMRAI::hier::Box intersection(q.density.getDim());
  q.density.getGhostBox().intersect(cons.density.getGhostBox(), intersection);
  std::vector<double> Y(reactor_.GetNSpecies());
  for (const SAMRAI::hier::Index& index : intersection) {
    SAMRAI::pdat::CellIndex cell(index);
    CopyMassFractions(Y, q.species, cell);
    const double rho = q.density(cell);
    reactor_.SetMassFractions(Y);
    reactor_.SetDensity(rho);
    const double rhou = q.momentum(cell);
    const double rhoE = q.energy(cell);
    const double U = (rhoE - 0.5 * rhou * rhou / rho) / rho;
    reactor_.SetInternalEnergy(U);
    reactor_.UpdateThermoState();
    q.temperature(cell) = reactor_.GetTemperature();
    q.pressure(cell) = reactor_.GetPressure();
    q.speed_of_sound(cell) = reactor_.GetSpeedOfSound();
  }
}

void TChemKinetics::FillFromPrim(const CompleteStatePatchData& q,
                                 const PrimStatePatchData& prim) const {
  q.temperature.copy(prim.temperature);
  q.momentum.copy(prim.momentum);
  q.pressure.copy(prim.pressure);
  q.species.copy(prim.species);
  SAMRAI::hier::Box intersection(q.density.getDim());
  q.pressure.getGhostBox().intersect(prim.pressure.getGhostBox(), intersection);
  std::vector<double> X(reactor_.GetNSpecies());
  for (const SAMRAI::hier::Index& index : intersection) {
    SAMRAI::pdat::CellIndex cell(index);
    const double p = q.pressure(cell);
    CopyMassFractions(X, q.species, cell);
    reactor_.SetMoleFractions(X);
    reactor_.SetPressure(p);
    reactor_.SetTemperature(q.temperature(cell));
    reactor_.UpdateThermoState();
    double rho = reactor_.GetDensity();
    q.density(cell) = rho;
    const double rhou = q.momentum(cell);
    // internal energy = (total energy - kinetic energy) / density
    const double U = reactor_.GetInternalEnergy();
    const double rhoE = rho * U + 0.5 * rhou * rhou / rho;
    q.energy(cell) = rhoE;
    q.speed_of_sound(cell) = reactor_.GetSpeedOfSound();
    span<const double> Y = reactor_.GetMassFractions();
    for (int s = 0; s < Y.size(); ++s) {
      q.species(cell, s) = rho * Y[s];
    }
  }
}

void TChemKinetics::AdvanceSourceTerm(const CompleteStatePatchData& q,
                                      double time_step_size) const {
  const SAMRAI::hier::Box& box = q.density.getGhostBox();
  std::vector<double> rhoY(q.species.getDepth());
  for (const SAMRAI::hier::Index& index : box) {
    SAMRAI::pdat::CellIndex cell(index);
    CopyMassFractions(rhoY, q.species, cell);
    reactor_.SetMassFractions(rhoY);
    reactor_.SetPressure(q.pressure(cell));
    reactor_.SetTemperature(q.temperature(cell));
    reactor_.Advance(time_step_size);
    reactor_.UpdateThermoState();
    const double rho = reactor_.GetDensity();
    const double u = q.momentum(cell) / rho;
    q.density(cell) = rho;
    q.temperature(cell) = reactor_.GetTemperature();
    const double p = reactor_.GetPressure();
    q.pressure(cell) = p;
    q.energy(cell) = rho * reactor_.GetInternalEnergy() + 0.5 * rho * u * u;
    q.momentum(cell) = rho * u;
    q.speed_of_sound(cell) = reactor_.GetSpeedOfSound();
    span<const double> Y = reactor_.GetMassFractions();
    for (int s = 0; s < q.species.getDepth(); ++s) {
      q.species(cell, s) = rho * Y[s];
    }
  }
}

} // namespace ideal_gas
} // namespace fub