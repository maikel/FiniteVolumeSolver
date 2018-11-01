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

#include "fub/solver/euler/IdealGas.hpp"
#include "fub/solver/euler/mechanism/Burke2012.hpp"

#include "SAMRAI/hier/VariableDatabase.h"
#include "SAMRAI/pdat/CellVariable.h"

#ifdef FUB_WITH_STD_STRING_VIEW
#include <string_view>
#endif

#include <numeric>

namespace fub {
namespace euler {
namespace {
using Variable = IdealGas::Variable;

#ifdef FUB_WITH_STD_STRING_VIEW
std::string getPrefixedCellVariableName_(const std::string& prefix,
                                         IdealGas::Variable var) {
  static constexpr std::array<std::string_view, IdealGas::variables_size> names{
      "/Density",    "/Momentum",     "/Energy", "/Pressure",
      "Temperature", "/SpeedOfSound", "Species"};
  std::string_view name = names[static_cast<int>(var)];
  std::string result;
  result.reserve(prefix.size() + name.size() + 1);
  result += prefix;
  result += name;
  return result;
}
#else
std::string getPrefixedCellVariableName_(const std::string& prefix,
                                         IdealGas::Variable var) {
  static constexpr std::array<const char*, IdealGas::variables_size> names{
      "/Density",    "/Momentum",     "/Energy", "/Pressure",
      "Temperature", "/SpeedOfSound", "Species"};
  const char* name = names[static_cast<int>(var)];
  std::string result;
  result.reserve(prefix.size() + std::strlen(name) + 1);
  result += prefix;
  result += name;
  return result;
}
#endif

int getDepth_(Variable var, int dim, int n_species) {
  int depth;
  switch (var) {
  case Variable::momentum:
    return dim;
  case Variable::species:
    return n_species - 1;
  default:
    return 1;
  }
}

int registerCellVariable_(const std::string& prefix, Variable var,
                          const SAMRAI::tbox::Dimension& dim, int n_species) {
  SAMRAI::hier::VariableDatabase& database =
      *SAMRAI::hier::VariableDatabase::getDatabase();
  // Construct a variable name using the given prefix
  std::string name = getPrefixedCellVariableName_(prefix, var);
  // Check if a variable with this name exists in the database
  std::shared_ptr<SAMRAI::hier::Variable> variable = database.getVariable(name);
  if (variable) {
    // This variable exists. Check if it is a CellVariable<double>.
    auto cell_variable =
        std::dynamic_pointer_cast<SAMRAI::pdat::CellVariable<double>>(variable);
    if (!cell_variable) {
      throw std::runtime_error(
          "A variable named '" + name +
          "' of differing type already exists in the variable database.");
    }
  } else {
    const int depth = getDepth_(var, dim.getValue(), n_species);
    // This variable does not exist yet, so we have to make a new one.
    variable =
        std::make_shared<SAMRAI::pdat::CellVariable<double>>(dim, name, depth);
  }
  // Register a (variable, "current")-pair w/o ghost cell width.
  // Return its patch data index.
  return database.registerVariableAndContext(
      variable, database.getContext("current"),
      SAMRAI::hier::IntVector::getZero(dim));
}

template <std::size_t N>
void registerAllVariables_(const std::string& prefix,
                           std::array<int, N>& data_ids,
                           const SAMRAI::tbox::Dimension& dim, int n_species) {
  const int size = N;
  for (int i = 0; i < size; ++i) {
    data_ids[i] =
        registerCellVariable_(prefix, IdealGas::Variable(i), dim, n_species);
  }
}
} // namespace

IdealGas::IdealGas(const std::string& name, const SAMRAI::tbox::Dimension& dim)
    : name_{name}, dimension_{dim}, reactor_{std::make_unique<Burke2012>()} {
  registerAllVariables_(name_, patch_data_ids_, dimension_,
                        reactor_.getNSpecies());
}

IdealGas::IdealGas(std::string&& name, const SAMRAI::tbox::Dimension& dim)
    : name_{std::move(name)},
      dimension_{dim}, reactor_{std::make_unique<Burke2012>()} {
  registerAllVariables_(name_, patch_data_ids_, dimension_,
                        reactor_.getNSpecies());
}

namespace {
void getMassFractions(span<double> fractions,
                      const SAMRAI::pdat::CellData<double>& species,
                      const SAMRAI::pdat::CellIndex& index) {
  FUB_ASSERT(fractions.size() == species.getDepth() + 1);
  const int max_depth = species.getDepth();
  for (int depth = 0; depth < max_depth; ++depth) {
    fractions[depth + 1] = species(index, depth);
  }
  const double total_fractions =
      std::accumulate(fractions.begin() + 1, fractions.end(), 0.0);
  fractions[0] = std::max(0.0, 1 - total_fractions);
}
} // namespace

void IdealGas::Reconstruct(const CompleteState& q, const ConsState& cons) {
  q.density.copy(cons.density);
  q.momentum.copy(cons.momentum);
  q.energy.copy(cons.energy);
  q.species.copy(cons.species);
  SAMRAI::hier::Box intersection(q.density.getDim());
  q.density.getGhostBox().intersect(cons.density.getGhostBox(), intersection);
  std::vector<double> fractions(q.species.getDepth() + 1);
  for (const SAMRAI::hier::Index& index : intersection) {
    SAMRAI::pdat::CellIndex cell(index);
    getMassFractions(fractions, q.species, cell);
    const double rho = q.density(cell);
    const double rhou = q.momentum(cell);
    const double rhoE = q.energy(cell);
    // internal energy = total energy - kinetic energy
    const double e = rhoE - 0.5 * rhou * rhou / rho;
    reactor_.setMassFractions(fractions);
    reactor_.setDensity(rho);
    reactor_.setInternalEnergy(e);
    q.temperature(cell) = reactor_.getTemperature();
    q.pressure(cell) = reactor_.getPressure();
    const double gamma = reactor_.getCp() / reactor_.getCv();
    const double p = reactor_.getPressure();
    q.speed_of_sound(cell) = std::sqrt(gamma * p / rho);
  }
}

void IdealGas::advanceSourceTerm(const CompleteState& q, double dt) {
  const SAMRAI::hier::Box& box = q.density.getGhostBox();
  std::vector<double> fractions(q.species.getDepth() + 1);
  for (const SAMRAI::hier::Index& index : box) {
    SAMRAI::pdat::CellIndex cell(index);
    getMassFractions(fractions, q.species, cell);
    reactor_.setMassFractions(fractions);
    reactor_.setTemperature(q.temperature(cell));
    reactor_.setPressure(q.pressure(cell));
    reactor_.advance(dt);
    const double u = q.momentum(cell) / q.density(cell);
    const double rho = reactor_.getDensity();
    q.density(cell) = rho;
    q.temperature(cell) = reactor_.getTemperature();
    q.pressure(cell) = reactor_.getPressure();
    q.energy(cell) = reactor_.getInternalEnergy() + 0.5 * rho * u * u;
    q.momentum(cell) = reactor_.getDensity() * u;
    const double gamma = reactor_.getCp() / reactor_.getCv();
    const double p = reactor_.getPressure();
    q.speed_of_sound(cell) = std::sqrt(gamma * p / rho);
  }
}

} // namespace euler
} // namespace fub