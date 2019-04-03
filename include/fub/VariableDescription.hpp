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

#ifndef FUB_VARIABLE_DESCRIPTION_HPP
#define FUB_VARIABLE_DESCRIPTION_HPP

#include "fub/core/dynamic_extent.hpp"
#include "fub/core/type_traits.hpp"
#include "fub/ext/Eigen.hpp"
#include "fub/ext/hana.hpp"

namespace fub {

using Scalar = int_constant<1>;
using Vector2d = int_constant<2>;
using Vector3d = int_constant<3>;
using VectorXd = int_constant<dynamic_extent>;

template <typename T> struct RemoveTemplateParameter;

template <template <typename...> typename T, typename... Ts>
struct RemoveTemplateParameter<T<Ts...>> {
  using type = template_t<T>;
};

enum class StateType { Complete, Conservative };

template <StateType Type> using StateType_c =  std::integral_constant<StateType, Type>;

static constexpr StateType_c<StateType::Conservative> cons{};
static constexpr StateType_c<StateType::Complete> complete{};

template <typename ConsShape, typename CompleteShape = ConsShape>
struct VariableDescription {
  using ConsTemplate = typename RemoveTemplateParameter<ConsShape>::type;
  
  using CompleteTemplate =
      typename RemoveTemplateParameter<CompleteShape>::type;

  static constexpr auto ConservativeIsComplete() {
    return std::is_same<ConsShape, CompleteShape>{};
  }

  static constexpr auto Template(StateType_c<StateType::Conservative>) {
    return ConsTemplate{};
  }

  static constexpr auto Template(StateType_c<StateType::Complete>) {
    return CompleteTemplate{};
  }

  static constexpr auto Shape(StateType_c<StateType::Conservative>) {
    return boost::hana::members(ConsShape{});
  }

  static constexpr auto Shape(StateType_c<StateType::Complete>) {
    return boost::hana::members(CompleteShape{});
  }

  template <StateType Type>
  static constexpr auto StaticSize(StateType_c<Type> type) {
    constexpr auto shape = Shape(type);
    constexpr auto any_dynamic_size =
        boost::hana::any_of(shape, [](auto extents) {
          return extents == int_constant<dynamic_extent>{};
        });
    if constexpr (any_dynamic_size) {
      return int_constant<dynamic_extent>{};
    } else {
      auto sum = [](auto total, auto n) {
        return int_constant<total() + n()>{};
      };
      return boost::hana::fold_left(shape, sum);
    }
  }

  static constexpr auto ValueTypes(StateType_c<StateType::Conservative>) {
    return boost::hana::transform(boost::hana::members(ConsShape{}),
                                  [](auto) { return type_c<double>; });
  }
  static constexpr auto ValueTypes(StateType_c<StateType::Complete>) {
    return boost::hana::transform(boost::hana::members(CompleteShape{}),
                                  [](auto) { return type_c<double>; });
  }
};

} // namespace fub

#endif
