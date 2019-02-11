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

#ifndef FUB_EQUATION_HPP
#define FUB_EQUATION_HPP

#include "fub/StateFacade.hpp"

#include <boost/hana/at_key.hpp>
#include <boost/hana/map.hpp>
#include <boost/hana/string.hpp>

namespace fub {

template <typename Equation, typename String>
using SizeMemberFunction =
    decltype(std::declval<const Equation&>().Size(String{}));

template <typename Equation> struct EquationTraits {
  template <typename State, typename String>
  static constexpr auto Size(const Equation& equation,
                             boost::hana::basic_type<State>, String name) {
    namespace hana = boost::hana;
    constexpr auto names = State::names();
    constexpr auto value_types = State::value_types();
    constexpr auto name_to_value_type = hana::to<hana::map_tag>(
        hana::transform(hana::zip(names, value_types), [](auto zipped) {
          using namespace hana::literals;
          return hana::make_pair(zipped[0_c], zipped[1_c]);
        }));
    if constexpr (is_detected<SizeMemberFunction, Equation, String>{}) {
      return equation.Size(name);
    } else {
      return overloaded{
          [](const hana::type<double>&) {
            return std::integral_constant<int, 1>{};
          },
          [](const auto& x) {
            using T = typename std::decay_t<decltype(x)>::type;
            return std::integral_constant<int, std::tuple_size<T>::value>{};
          }}(name_to_value_type[name]);
    }
  }
};

template <typename State>
auto GetSizes(const typename State::Equation& equation) noexcept {
  constexpr auto names = State::names();
  using Equation = typename State::Equation;
  return Transform(
      [&](auto name) {
        return EquationTraits<Equation>::Size(equation,
                                              boost::hana::type_c<State>, name);
      },
      boost::hana::unpack(names, [](auto... names) {
        return boost::hana::make<typename State::hana_tag>(names...);
      }));
}

template <typename State> constexpr auto GetStaticSizes() noexcept {
  constexpr auto names = State::names();
  using Equation = typename State::Equation;
  return Transform(
      [&](auto name) {
        using SizeT = decltype(EquationTraits<Equation>::Size(
            std::declval<const Equation&>(), boost::hana::type_c<State>, name));
        if constexpr (std::is_integral_v<SizeT>) {
          return std::integral_constant<int, -1>{};
        } else {
          return SizeT{};
        }
      },
      boost::hana::unpack(names, [](auto... names) {
        return boost::hana::make<typename State::hana_tag>(names...);
      }));
}

template <typename SimdState> auto ToSingleState() noexcept {
  SimdState state;
  return Transform(
      [](auto x) {
        using T = std::decay_t<decltype(x)>;
        if constexpr (T::ColsAtCompileTime == 1) {
          return double{};
        } else if (T::ColsAtCompileTime == Eigen::Dynamic) {
          return std::vector<double>();
        } else {
          return std::array<double, T::ColsAtCompileTime>{};
        }
      },
      state);
};

template <typename State> auto ToStateArrayHelper() noexcept {
  State state;
  return Transform([](auto&&) { return std::vector<double>{}; }, state);
}

template <typename Equation>
using StateArray = decltype(ToStateArrayHelper<typename Equation::State>());

template <typename Equation>
using ConsArray = decltype(ToStateArrayHelper<typename Equation::Cons>());

template <typename Equation>
typename Equation::State MakeState(const Equation& eq) noexcept {
  using State = typename Equation::State;
  State state;
  const auto sizes = GetSizes<State>(eq);
  constexpr auto mem_accessors = State::accessors();
  constexpr auto size_accessors = std::decay_t<decltype(sizes)>::accessors();
  constexpr auto zipped = boost::hana::zip(mem_accessors, size_accessors);
  boost::hana::for_each(zipped, [&](auto accessor) {
    using namespace boost::hana::literals;
    boost::hana::second(accessor[0_c])(state).resize(
        kChunkSize, boost::hana::second(accessor[1_c])(sizes));
  });
  return state;
}

template <typename Equation>
typename Equation::Cons MakeCons(const Equation& eq) noexcept {
  using Cons = typename Equation::Cons;
  Cons state;
  const auto sizes = GetSizes<Cons>(eq);
  constexpr auto mem_accessors = Cons::accessors();
  constexpr auto size_accessors = std::decay_t<decltype(sizes)>::accessors();
  constexpr auto zipped = boost::hana::zip(mem_accessors, size_accessors);
  boost::hana::for_each(zipped, [&](auto accessor) {
    using namespace boost::hana::literals;
    boost::hana::second(accessor[0_c])(state).resize(
        kChunkSize, boost::hana::second(accessor[1_c])(sizes));
  });
  return state;
}

template <typename Equation>
StateArray<Equation> MakeStateArray(const Equation& eq, int n_cells) {
  const auto sizes = GetSizes<typename Equation::State>(eq);
  constexpr auto size_accessors = std::decay_t<decltype(sizes)>::accessors();
  return boost::hana::unpack(size_accessors, [&](auto... accessor) {
    return StateArray<Equation>{
        std::vector<double>(boost::hana::second(accessor)(sizes) * n_cells)...};
  });
}

template <typename Equation>
ConsArray<Equation> MakeConsArray(const Equation& eq, int n_cells) {
  const auto sizes = GetSizes<typename Equation::Cons>(eq);
  constexpr auto size_accessors = std::decay_t<decltype(sizes)>::accessors();
  return boost::hana::unpack(size_accessors, [&](auto... accessor) {
    return ConsArray<Equation>{
        std::vector<double>(boost::hana::second(accessor)(sizes) * n_cells)...};
  });
}

template <typename T, typename State, typename Layout = layout_left>
auto ToStateSpanHelper() noexcept {
  State state;
  return Transform(
      overloaded{[](const double&) {
                   return basic_mdspan<T, DynamicExtents<3>, Layout>{};
                 },
                 [](const auto&) {
                   return basic_mdspan<T, DynamicExtents<4>, Layout>{};
                 }},
      state);
}

template <typename T, typename Equation, typename Layout = layout_left>
using StateView =
    decltype(ToStateSpanHelper<T, typename Equation::State, Layout>());

template <typename T, typename Equation>
using ConsView = decltype(ToStateSpanHelper<T, typename Equation::Cons>());

template <typename State> auto Simdify(const State& state) {
  return Transform(
      overloaded{[](const double&) { return Scalar{}; },
                 [](const std::vector<double>& x) {
                   return ArrayXd(kChunkSize, x.size());
                 },
                 [](const auto& x) {
                   constexpr int array_size =
                       std::tuple_size<std::decay_t<decltype(x)>>::value;
                   return Array<double, array_size>{};
                 }},
      state);
}

template <typename State>
using SimdifyT = decltype(Simdify(std::declval<const State&>()));

template <typename State> auto AsCons(State&& state) {
  using StateT = std::decay_t<State>;
  using Equation = typename StateT::Equation;
  using ConsT = typename Equation::Cons;
  constexpr auto map = StateT::as_map();
  constexpr auto cnames = ConsT::names();
  return boost::hana::unpack(cnames, [&](auto... name) {
    return boost::hana::make<typename ConsT::hana_tag>(
        Ref(map[name](state))...);
  });
}

template <typename E, typename Function>
Function ForEachIndex(const layout_left::mapping<E>& mapping,
                      Function function) {
  if constexpr (E::rank() == 1) {
    for (int i = 0; i < mapping.extents().extent(0); ++i) {
      function(i);
    }
  } else if constexpr (E::rank() == 2) {
    for (int i = 0; i < mapping.extents().extent(1); ++i) {
      for (int j = 0; j < mapping.extents().extent(0); ++j) {
        function(j, i);
      }
    }
  } else if constexpr (E::rank() == 3) {
    for (int i = 0; i < mapping.extents().extent(2); ++i) {
      for (int j = 0; j < mapping.extents().extent(1); ++j) {
        for (int k = 0; k < mapping.extents().extent(0); ++k) {
          function(k, j, i);
        }
      }
    }
  }
  return function;
}

template <typename E, typename Function>
Function ForEachIndex(const layout_stride::mapping<E>& mapping,
                      Function function) {
  if constexpr (E::rank() == 1) {
    for (int i = 0; i < mapping.extents().extent(0); ++i) {
      function(i);
    }
  } else if constexpr (E::rank() == 2) {
    for (int i = 0; i < mapping.extents().extent(1); ++i) {
      for (int j = 0; j < mapping.extents().extent(0); ++j) {
        function(j, i);
      }
    }
  } else if constexpr (E::rank() == 3) {
    for (int i = 0; i < mapping.extents().extent(2); ++i) {
      for (int j = 0; j < mapping.extents().extent(1); ++j) {
        for (int k = 0; k < mapping.extents().extent(0); ++k) {
          function(k, j, i);
        }
      }
    }
  }
  return function;
}

template <typename State> auto ViewInnerRegion(const State& state, int gcw) {
  return Transform(
      [gcw](auto mdspan) {
        if constexpr (decltype(mdspan)::rank() == 2) {
          return subspan(mdspan, std::pair{gcw, mdspan.extent(0) - gcw}, all);
        }
        else if constexpr (decltype(mdspan)::rank() == 3) {
          return subspan(mdspan, std::pair{gcw, mdspan.extent(0) - gcw}, all,
                         all);
        } else {
          return subspan(mdspan, std::pair{gcw, mdspan.extent(0) - gcw}, all,
                         all, all);
        }
      },
      state);
}

} // namespace fub

#endif