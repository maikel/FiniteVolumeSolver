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

#ifndef FUB_SAMRAI_REGISTER_VARIABLES_HPP
#define FUB_SAMRAI_REGISTER_VARIABLES_HPP

#include "fub/Equation.hpp"

#include <SAMRAI/pdat/CellVariable.h>
#include <SAMRAI/pdat/FaceVariable.h>
#include <SAMRAI/pdat/OuterfaceVariable.h>

#include <vector>

namespace fub {
namespace samrai {

struct DataDescription {
  std::vector<int> data_ids;
  std::vector<int> conservative;
};

std::string MakeVariableName(const std::string& prefix,
                             const std::string& name);

template <typename VariableType>
std::shared_ptr<VariableType> GetVariable(SAMRAI::tbox::Dimension dim,
                                          const std::string& name, int depth) {
  SAMRAI::hier::VariableDatabase* vardb =
      SAMRAI::hier::VariableDatabase::getDatabase();
  if (auto variable = vardb->getVariable(name)) {
    if (auto specialised = std::dynamic_pointer_cast<VariableType>(variable);
        specialised && specialised->getDepth() == depth &&
        specialised->getDim() == dim) {
      return specialised;
    }
    throw std::runtime_error(
        "A variable with name '" + name +
        "' was previously registered with differing depth or dimension.");
  }
  return std::make_shared<VariableType>(dim, name, depth);
}

template <typename State, typename VariableType, typename Equation>
void RegisterVariables(std::vector<int>& data_ids, const Equation& equation,
                       const SAMRAI::tbox::Dimension& dim,
                       const SAMRAI::hier::IntVector& ghost_layer_width,
                       const std::string& prefix,
                       const std::string& context_name) {
  constexpr auto names = State::Names();
  const auto sizes = equation.Shape(complete);
  SAMRAI::hier::VariableDatabase* vardb =
      SAMRAI::hier::VariableDatabase::getDatabase();
  std::shared_ptr<SAMRAI::hier::VariableContext> context =
      vardb->getContext(context_name);
  data_ids.reserve(boost::hana::length(names));
  boost::hana::for_each(boost::hana::zip(names, sizes), [&](auto xs) {
    const int depth = at_c<1>(xs);
    const char* name = at_c<0>(xs).c_str();
    const std::string variable_name = MakeVariableName(prefix, name);
    auto variable = GetVariable<VariableType>(dim, variable_name, depth);
    const int data_id =
        vardb->registerVariableAndContext(variable, context, ghost_layer_width);
    data_ids.push_back(data_id);
  });
}

/// This function registers all neccessary variables and contexts with SAMRAI.
/// The returned DataDescription shall be passed to solver classes.
///
/// \note It is currently assumed, that state and scratch variables are cell
/// centered and flux and coarse_fine_interfaces are face centered.
///
/// \param[in] equation The equation describes shape and size for each variable.
///
/// \param[in] prefix An optional prefix which will be prepend to all variable
/// names.
///
/// \return Returns the patch data ids for all registered variables.
template <typename Equation>
DataDescription RegisterVariables(const Equation& equation,
                                 std::string prefix = std::string()) {
  using Complete = typename Equation::Complete;
  using Conservative = typename Equation::Cons;

  DataDescription data_ids;

  const SAMRAI::tbox::Dimension dim(equation.Rank());

  const SAMRAI::hier::IntVector zero = SAMRAI::hier::IntVector::getZero(dim);
  RegisterVariables<Complete, SAMRAI::pdat::CellVariable<double>>(
      data_ids.data_ids, equation, dim, zero, prefix, "current");

  constexpr auto cons_names = Conservative::Names();
  constexpr auto complete_names = Complete::Names();

  data_ids.conservative.reserve(boost::hana::length(cons_names));
  boost::hana::for_each(cons_names, [&](auto name) {
    constexpr auto index =
        *boost::hana::index_if(complete_names, boost::hana::equal.to(name));
    data_ids.conservative.push_back(index);
  });

  return data_ids;
}

} // namespace samrai
} // namespace fub

#endif