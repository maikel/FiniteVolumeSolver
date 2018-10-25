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

#include "SAMRAI/hier/VariableDatabase.h"
#include "SAMRAI/pdat/CellVariable.h"

#include <string_view>

namespace fub {
namespace euler {
namespace {
using Variable = IdealGas::Variable;

std::string getPrefixedCellVariableName_(const std::string& prefix,
                                         IdealGas::Variable var) {
  static constexpr std::array<std::string_view, IdealGas::variables_size> names{
      "/Density", "/Momentum", "/Energy", "/Pressure", "/SpeedOfSound"};
  std::string_view variable_name = names[static_cast<int>(var)];
  std::string result;
  result.reserve(prefix.size() + variable_name.size() + 1);
  result += prefix;
  result += variable_name;
  return result;
}

int registerCellVariable_(const std::string& prefix, Variable var,
                          const SAMRAI::tbox::Dimension& dim) {
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
    int depth = var == Variable::momentum ? dim.getValue() : 1;
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
                           const SAMRAI::tbox::Dimension& dim) {
  const int size = N;
  for (int i = 0; i < size; ++i) {
    data_ids[i] = registerCellVariable_(prefix, IdealGas::Variable(i), dim);
  }
}
} // namespace

IdealGas::IdealGas(const std::string& name, const SAMRAI::tbox::Dimension& dim)
    : name_{name}, dimension_{dim} {
  registerAllVariables_(name_, data_ids_, dimension_);
}

IdealGas::IdealGas(std::string&& name, const SAMRAI::tbox::Dimension& dim)
    : name_{std::move(name)}, dimension_{dim} {
  registerAllVariables_(name_, data_ids_, dimension_);
}

namespace {
double transform_reduce(Vector<double> x, Vector<double> y) {
  double total = 0.0;
  for (std::size_t i = 0; i < x.size(); ++i) {
    total += x[i] * y[i];
  }
  return total;
}
} // namespace

double IdealGas::computePressure(Conservative<double> q) const {
  const double kinetic_energy =
      0.5 * fub::euler::transform_reduce(q.momentum, q.momentum) / q.density;
  FUB_ASSERT(kinetic_energy <= q.energy);
  const double internal_energy = q.energy - kinetic_energy;
  const double gamma_1 = getHeatCapacityRatio() - 1.0;
  const double pressure = gamma_1 * internal_energy * q.density;
  return pressure;
}

double IdealGas::computeSpeedOfSound(double density, double pressure) const {
  return std::sqrt(getHeatCapacityRatio() * pressure / density);
}

IdealGas::Conservative<double>
IdealGas::computeFlux(Complete<double> state, int dir) const {
  const double velocity = state.momentum[dir] / state.density;
  IdealGas::Conservative<double> cons;
  cons.density = state.momentum[dir];
  cons.momentum = state.momentum;
  for (std::size_t i = 0; i < cons.momentum.size(); ++i) {
    cons.momentum[i] *= velocity;
  }
  cons.momentum[dir] += state.pressure;
  cons.energy = velocity * (state.energy + state.pressure);
  return cons;
}

} // namespace euler
} // namespace fub