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

#include "fub/ideal_gas/IdealGasEquation.hpp"

#include "SAMRAI/hier/VariableDatabase.h"
#include "SAMRAI/pdat/CellVariable.h"

#include <numeric>

namespace fub {
namespace ideal_gas {
namespace {
using Variable = IdealGasEquation::Variable;

std::string GetPrefixedCellVariableName_(const std::string& prefix,
                                         Variable var) {
  static constexpr std::array<const char*, IdealGasEquation::variables_size>
      names{"_Density",     "_Momentum",     "_Energy", "_Pressure",
            "_Temperature", "_SpeedOfSound", "_Species"};
  const char* name = names[static_cast<int>(var)];
  std::string result;
  result.reserve(prefix.size() + std::strlen(name) + 1);
  result += prefix;
  result += name;
  return result;
}

int GetDepth_(Variable var, int dim, int n_species) {
  int depth;
  switch (var) {
  case Variable::momentum:
    return dim;
  case Variable::species:
    return n_species;
  default:
    return 1;
  }
}

template <typename VariableType>
std::shared_ptr<VariableType>
GetVariableOr_(const std::string& name,
               const std::shared_ptr<VariableType>& alt) {
  std::shared_ptr<SAMRAI::hier::Variable> variable =
      SAMRAI::hier::VariableDatabase::getDatabase()->getVariable(name);
  if (variable) {
    auto var = std::dynamic_pointer_cast<VariableType>(variable);
    if (var && var->getDepth() == alt->getDepth()) {
      return var;
    }
    throw std::runtime_error(
        "A variable named '" + name +
        "' of differing type already exists in the variable database.");
  }
  return alt;
}

template <typename VariableType>
std::shared_ptr<VariableType> GetVariable_(const SAMRAI::tbox::Dimension& dim,
                                           const std::string& name, int depth) {
  return GetVariableOr_(name, std::make_shared<VariableType>(dim, name, depth));
}

int RegisterCellVariable_(const std::string& prefix, Variable var,
                          const SAMRAI::tbox::Dimension& dim, int n_species) {
  SAMRAI::hier::VariableDatabase& database =
      *SAMRAI::hier::VariableDatabase::getDatabase();
  // Construct a variable name using the given prefix
  std::string name = GetPrefixedCellVariableName_(prefix, var);
  // Check if a variable with this name exists in the database
  auto v = GetVariable_<SAMRAI::pdat::CellVariable<double>>(
      dim, name, GetDepth_(var, dim.getValue(), n_species));
  // Register a (variable, "current")-pair w/o ghost cell width.
  // Return its patch data index.
  return database.registerVariableAndContext(
      v, database.getContext("current"), SAMRAI::hier::IntVector::getZero(dim));
}

template <std::size_t N>
void RegisterAllVariables_(const std::string& prefix,
                           std::array<int, N>& data_ids,
                           const SAMRAI::tbox::Dimension& dim, int n_species) {
  const int size = N;
  for (int i = 0; i < size; ++i) {
    data_ids[i] = RegisterCellVariable_(prefix, IdealGasEquation::Variable(i),
                                        dim, n_species);
  }
}
} // namespace

IdealGasEquation::IdealGasEquation(std::string name,
                                   SAMRAI::tbox::Dimension dim, int n_species)
    : name_{std::move(name)}, dimension_{dim}, n_species_{n_species} {
  RegisterAllVariables_(name_, patch_data_ids_, dimension_, n_species);
}

void CopyMassFractions(span<double> fractions,
                       const SAMRAI::pdat::CellData<double>& species,
                       const SAMRAI::pdat::CellIndex& index) {
  FUB_ASSERT(fractions.size() == species.getDepth());
  const int max_depth = species.getDepth();
  for (int depth = 0; depth < max_depth; ++depth) {
    fractions[depth] = species(index, depth);
  }
}

IdealGasEquation::CompleteStatePatchData
IdealGasEquation::GetCompleteStatePatchData(
    const SAMRAI::hier::Patch& patch) const {
  // clang-format off
  SAMRAI::pdat::CellData<double>& density = GetCellData(patch, Variable::density);
  SAMRAI::pdat::CellData<double>& momentum = GetCellData(patch, Variable::momentum);
  SAMRAI::pdat::CellData<double>& energy = GetCellData(patch, Variable::energy);
  SAMRAI::pdat::CellData<double>& pressure = GetCellData(patch, Variable::pressure);
  SAMRAI::pdat::CellData<double>& temperature = GetCellData(patch, Variable::temperature);
  SAMRAI::pdat::CellData<double>& speed_of_sound = GetCellData(patch, Variable::speed_of_sound);
  SAMRAI::pdat::CellData<double>& species = GetCellData(patch, Variable::species);
  // clang-format on
  return {density,     momentum,       energy, pressure,
          temperature, speed_of_sound, species};
}

IdealGasEquation::ConsStatePatchData IdealGasEquation::GetConsStatePatchData(
    const SAMRAI::hier::Patch& patch) const {
  // clang-format off
  SAMRAI::pdat::CellData<double>& density = GetCellData(patch, Variable::density);
  SAMRAI::pdat::CellData<double>& momentum = GetCellData(patch, Variable::momentum);
  SAMRAI::pdat::CellData<double>& energy = GetCellData(patch, Variable::energy);
  SAMRAI::pdat::CellData<double>& species = GetCellData(patch, Variable::species);
  // clang-format on
  return {density, momentum, energy, species};
}

IdealGasEquation::PrimStatePatchData IdealGasEquation::GetPrimStatePatchData(
    const SAMRAI::hier::Patch& patch) const {
  // clang-format off
  SAMRAI::pdat::CellData<double>& temperature = GetCellData(patch, Variable::temperature);
  SAMRAI::pdat::CellData<double>& momentum = GetCellData(patch, Variable::momentum);
  SAMRAI::pdat::CellData<double>& pressure = GetCellData(patch, Variable::pressure);
  SAMRAI::pdat::CellData<double>& species = GetCellData(patch, Variable::species);
  // clang-format on
  return {temperature, momentum, pressure, species};
}

span<double> MakeSpan(SAMRAI::pdat::CellData<double>& data) noexcept {
  return span<double>(data.getPointer(), data.getArrayData().getOffset() *
                                             data.getArrayData().getDepth());
}

span<const double>
MakeSpan(const SAMRAI::pdat::CellData<double>& data) noexcept {
  return span<const double>(data.getPointer(),
                            data.getArrayData().getOffset() *
                                data.getArrayData().getDepth());
}

} // namespace ideal_gas
} // namespace fub