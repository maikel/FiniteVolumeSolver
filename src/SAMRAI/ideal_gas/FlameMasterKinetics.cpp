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

#include "fub/SAMRAI/ideal_gas/FlameMasterKinetics.hpp"

#include "fub/ideal_gas/mechanism/Burke2012.hpp"

#include "SAMRAI/hier/VariableDatabase.h"
#include "SAMRAI/pdat/CellVariable.h"

#include <algorithm>
#include <numeric>

namespace fub {
namespace ideal_gas {
using Variable = IdealGasEquation::Variable;

FlameMasterKinetics::FlameMasterKinetics(std::string name,
                                         SAMRAI::tbox::Dimension dim,
                                         FlameMasterReactor reactor)
    : IdealGasKinetics(name, dim, reactor.getNSpecies()), reactor_{std::move(
                                                              reactor)} {
  reactor_.setTemperature(300);
  span<const double> hs = reactor_.getEnthalpies();
  reference_enthalpies_.resize(hs.size());
  std::copy(hs.begin(), hs.end(), reference_enthalpies_.begin());
}

void FlameMasterKinetics::FillFromCons(const CompletePatchData& q,
                                       const ConsPatchData& cons) const {
  q.density.copy(cons.density);
  q.momentum.copy(cons.momentum);
  q.energy.copy(cons.energy);
  q.species.copy(cons.species);
  SAMRAI::hier::Box intersection(q.density.getDim());
  q.density.getGhostBox().intersect(cons.density.getGhostBox(), intersection);
  std::vector<double> fractions(q.species.getDepth());
  for (const SAMRAI::hier::Index& index : intersection) {
    SAMRAI::pdat::CellIndex cell(index);
    const double rho = q.density(cell);
    CopyMassFractions(fractions, q.species, cell);
    const double rhou = q.momentum(cell);
    const double rhoE = q.energy(cell);
    reactor_.setMassFractions(fractions);
    reactor_.setDensity(rho);
    // internal energy = (total energy - kinetic energy) / density
    const double e = (rhoE - 0.5 * rhou * rhou / rho) / rho;
    reactor_.setInternalEnergy(e);
    q.temperature(cell) = reactor_.getTemperature();
    q.pressure(cell) = reactor_.getPressure();
    const double gamma = reactor_.getCp() / reactor_.getCv();
    const double p = reactor_.getPressure();
    q.speed_of_sound(cell) = std::sqrt(gamma * p / rho);
  }
}

// void FlameMasterKinetics::FillFromCons(
//     const CompleteStateSpan<double>& q,
//     const ConsStateSpan<const double>& u) const {
//   std::copy(u.density.begin(), u.density.end(), q.density.begin());
//   std::copy(u.momentum.begin(), u.momentum.end(), q.momentum.begin());
//   std::copy(u.energy.begin(), u.energy.end(), q.energy.begin());
//   std::copy(u.species.begin(), u.species.end(), q.species.begin());
//   std::vector<double> masses(reactor_.getNSpecies());
//   const int size = q.density.size();
//   for (std::ptrdiff_t i = 0; i < size; ++i) {
//     const double rho = q.density[i];
//     const double rhou = q.momentum[i];
//     const double rhoE = q.energy[i];
//     for (int s = 0; s < reactor_.getNSpecies(); ++s) {
//       masses[s] = q.species[i + s * size];
//     }
//     reactor_.setMassFractions(masses);
//     reactor_.setDensity(rho);
//     // internal energy = (total energy - kinetic energy) / density
//     const double e = (rhoE - 0.5 * rhou * rhou / rho) / rho;
//     const double gamma = reactor_.getCp() / reactor_.getCv();
//     const double p = reactor_.getPressure();
//     reactor_.setInternalEnergy(e);
//     q.temperature[i] = reactor_.getTemperature();
//     q.pressure[i] = reactor_.getPressure();
//     q.speed_of_sound[i] = std::sqrt(gamma * p / rho);
//   }
// }

void FlameMasterKinetics::FillFromPrim(const CompletePatchData& q,
                                       const PrimPatchData& prim) const {
  q.temperature.copy(prim.temperature);
  q.momentum.copy(prim.momentum);
  q.pressure.copy(prim.pressure);
  q.species.copy(prim.species);
  SAMRAI::hier::Box intersection(q.density.getDim());
  q.density.getGhostBox().intersect(prim.pressure.getGhostBox(), intersection);
  std::vector<double> masses(q.species.getDepth());
  for (const SAMRAI::hier::Index& index : intersection) {
    SAMRAI::pdat::CellIndex cell(index);
    CopyMassFractions(masses, q.species, cell);
    reactor_.setMoleFractions(masses);
    reactor_.setTemperature(q.temperature(cell));
    reactor_.setPressure(q.pressure(cell));
    q.density(cell) = reactor_.getDensity();
    const double rho = reactor_.getDensity();
    q.momentum(cell) *= rho;
    const double rhou = q.momentum(cell);
    // internal energy = (total energy - kinetic energy) / density
    const double eps = reactor_.getInternalEnergy();
    const double rhoE = rho * eps + 0.5 * rhou * rhou / rho;
    q.energy(cell) = rhoE;
    const double gamma = reactor_.getCp() / reactor_.getCv();
    const double p = reactor_.getPressure();
    q.speed_of_sound(cell) = std::sqrt(gamma * p / rho);
    span<const double> fractions = reactor_.getMassFractions();
    for (int s = 0; s < fractions.size(); ++s) {
      q.species(cell, s) = rho * fractions[s];
    }
  }
}

// void FlameMasterKinetics::FillFromPrim(
//     const CompleteStateSpan<double>& q,
//     const PrimStateSpan<const double>& w) const {
//   FUB_ASSERT(q.temperature.size() == w.temperature.size());
//   std::copy(w.temperature.begin(), w.temperature.end(), q.temperature.begin());
//   std::copy(w.pressure.begin(), w.pressure.end(), q.pressure.begin());
//   std::copy(w.species.begin(), w.species.end(), q.species.begin());
//   std::copy(w.momentum.begin(), w.momentum.end(), q.momentum.begin());
//   std::vector<double> moles(reactor_.getNSpecies());
//   const int size = q.temperature.size();
//   for (std::ptrdiff_t i = 0; i < size; ++i) {
//     for (int s = 0; s < reactor_.getNSpecies(); ++s) {
//       moles[s] = q.species[i + s * size];
//     }
//     reactor_.setMoleFractions(moles);
//     reactor_.setTemperature(q.temperature[i]);
//     reactor_.setPressure(q.pressure[i]);
//     const double rho = reactor_.getDensity();
//     q.momentum[i] *= rho;
//     const double rhou = q.momentum[i];
//     // internal energy = (total energy - kinetic energy) / density
//     const double eps = reactor_.getInternalEnergy();
//     const double rhoE = rho * eps + 0.5 * rhou * rhou / rho;
//     const double gamma = reactor_.getCp() / reactor_.getCv();
//     const double p = reactor_.getPressure();
//     q.density[i] = rho;
//     q.energy[i] = rhoE;
//     q.speed_of_sound[i] = std::sqrt(gamma * p / rho);
//     span<const double> Y = reactor_.getMassFractions();
//     for (int s = 0; s < Y.size(); ++s) {
//       q.species[i + s * size] = rho * Y[s];
//     }
//   }
// }

namespace {
void Normalize_(span<double> x) {
  const double total = std::accumulate(x.begin(), x.end(), 0.0);
  std::transform(x.begin(), x.end(), x.begin(),
                 [total](double xi) { return xi / total; });
}
} // namespace

void FlameMasterKinetics::AdvanceSourceTerm(const CompletePatchData& q,
                                            double dt) const {
  const SAMRAI::hier::Box& box = q.density.getGhostBox();
  std::vector<double> fractions(q.species.getDepth());
  for (const SAMRAI::hier::Index& index : box) {
    SAMRAI::pdat::CellIndex cell(index);
    CopyMassFractions(fractions, q.species, cell);
    reactor_.setMassFractions(fractions);
    reactor_.setTemperature(q.temperature(cell));
    reactor_.setPressure(q.pressure(cell));
    reactor_.advance(dt);
    const double u = q.momentum(cell) / q.density(cell);
    const double rho = reactor_.getDensity();
    q.density(cell) = rho;
    q.temperature(cell) = reactor_.getTemperature();
    q.pressure(cell) = reactor_.getPressure();
    q.energy(cell) = rho * reactor_.getInternalEnergy() + 0.5 * rho * u * u;
    q.momentum(cell) = reactor_.getDensity() * u;
    const double gamma = reactor_.getCp() / reactor_.getCv();
    const double p = reactor_.getPressure();
    q.speed_of_sound(cell) = std::sqrt(gamma * p / rho);
    span<const double> Y = reactor_.getMassFractions();
    std::copy(Y.begin(), Y.end(), fractions.begin());
    Normalize_(fractions);
    for (int s = 0; s < q.species.getDepth(); ++s) {
      q.species(cell, s) = rho * fractions[s];
    }
  }
}

void FlameMasterKinetics::UpdateCellFromReactor(
    CompletePatchData data, const SAMRAI::pdat::CellIndex cell,
    const FlameMasterReactor& reactor, double velocity) const {
  const double rho = reactor.getDensity();
  data.density(cell) = rho;
  data.temperature(cell) = reactor.getTemperature();
  data.momentum(cell) = rho * velocity;
  data.energy(cell) =
      rho * reactor.getInternalEnergy() + 0.5 * rho * velocity * velocity;
  data.pressure(cell) = reactor.getPressure();
  const double gamma = reactor.getCp() / reactor.getCv();
  data.speed_of_sound(cell) = std::sqrt(gamma * reactor.getPressure() / rho);
  const int n_species = data.species.getDepth();
  span<const double> Y = reactor.getMassFractions();
  for (int s = 0; s < n_species; ++s) {
    data.species(cell, s) = rho * Y[s];
  }
}

void FlameMasterKinetics::UpdateStateFromReactor(
    span<double> buffer, const FlameMasterReactor& reactor,
    double velocity) const {
  const double rho = reactor.getDensity();
  buffer[int(Variable::density)] = rho;
  buffer[int(Variable::momentum)] = rho * velocity;
  buffer[int(Variable::temperature)] = reactor.getTemperature();
  buffer[int(Variable::energy)] =
      rho * reactor.getInternalEnergy() + 0.5 * rho * velocity * velocity;
  buffer[int(Variable::pressure)] = reactor.getPressure();
  const double gamma = reactor.getCp() / reactor.getCv();
  buffer[int(Variable::speed_of_sound)] =
      std::sqrt(gamma * reactor.getPressure() / rho);
  const int n_species = GetNSpecies();
  span<const double> Y = reactor.getMassFractions();
  std::transform(Y.begin(), Y.end(), &buffer[int(Variable::species)],
                 [rho](double Yi) { return rho * Yi; });
}

} // namespace ideal_gas
} // namespace fub
