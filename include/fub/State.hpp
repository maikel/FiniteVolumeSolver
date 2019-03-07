// Copyright (c) 2019 Maikel Nadolski
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

#ifndef FUB_STATE_HPP
#define FUB_STATE_HPP

#include "fub/VariableDescription.hpp"

namespace fub {
template <typename Base> struct StateBase : public Base {
  using BaseType = Base;
  using TemplateType = typename RemoveTemplateParameter<Base>::type;

  Base& base() noexcept { return *this; }
  const Base& base() const noexcept { return *this; }

  static constexpr auto Template() { return TemplateType{}; }

  static constexpr auto Names() {
    return decltype(boost::hana::keys(std::declval<Base>())){};
  }

  using MemberTypes = decltype(boost::hana::members(std::declval<Base>()));
  MemberTypes Members() const { return boost::hana::members(base()); }

  static constexpr auto Accessors() { return boost::hana::accessors<Base>(); }

  static constexpr auto PointersToMember() {
    return boost::hana::values(boost::hana::to_map(Accessors()));
  }

  static constexpr auto ValueTypes() {
    constexpr auto accessors = Accessors();
    return boost::hana::transform(accessors, [&](auto a) {
      using MemberType =
          decltype(boost::hana::second(a)(std::declval<StateBase<Base>&>()));
      using ValueType = remove_cvref_t<MemberType>;
      return type_c<ValueType>;
    });
  }

  template <typename... Xs> static constexpr auto Rebind(type<Xs>... types) {
    return Template()(types...);
  }
};

template <typename T>
constexpr auto ScalarComponentType(basic_type<T>, constant<1>) {
  return type_c<T>;
}

template <typename T>
constexpr auto ScalarComponentType(basic_type<T>, std::ptrdiff_t) {
  return type_c<Array<T, Eigen::Dynamic, 1>>;
}

template <typename T, int N, typename = std::enable_if_t<(N > 0)>>
constexpr auto ScalarComponentType(basic_type<T>, constant<N>) {
  return type_c<Array<T, N, 1>>;
}

template <typename Equation, StateType Type> constexpr auto ScalarBaseType() {
  constexpr StateType_c<Type> type{};
  constexpr auto value_types = Equation::ValueTypes(type);
  constexpr auto shape = Equation::Shape(type);
  constexpr auto member_types =
      boost::hana::transform(boost::hana::zip(value_types, shape), [](auto xs) {
        return ScalarComponentType(at_c<0>(xs), at_c<1>(xs));
      });
  constexpr auto template_ = Equation::Template(type);
  return boost::hana::unpack(
      member_types, [&](auto... types) { return template_(types...); });
}

template <typename Equation, StateType Type>
using ScalarBaseTypeT =
    typename decltype(ScalarBaseType<Equation, Type>())::type;

template <typename Derived, typename Base, bool IsDynamicSize = true>
struct StateVector : public StateBase<Base> {
  Derived operator+(const Derived& other) const {
    return boost::hana::unpack(
        boost::hana::zip(this->Members(), other.Members()),
        [](auto... xs) { return Derived{(at_c<0>(xs) + at_c<1>(xs))...}; });
  }

  Derived operator-(const Derived& other) const {
    return boost::hana::unpack(
        boost::hana::zip(this->Members(), other.Members()),
        [](auto... xs) { return Derived{(at_c<0>(xs) - at_c<1>(xs))...}; });
  }

  Derived operator*(const Derived& other) const {
    return boost::hana::unpack(
        boost::hana::zip(this->Members(), other.Members()),
        [](auto... xs) { return Derived{(at_c<0>(xs) * at_c<1>(xs))...}; });
  }

  StateVector() = default;

  template <typename... Args,
            typename = std::enable_if_t<
                StateBase<Base>::ValueTypes() ==
                boost::hana::make_tuple(type_c<remove_cvref_t<Args>>...)>>
  StateVector(Args&&... args) noexcept
      : StateBase<Base>{std::forward<Args>(args)...} {}
};

template <typename Equation>
struct Complete : StateVector<Complete<Equation>,
                              ScalarBaseTypeT<Equation, StateType::Complete>> {
  using EquationType = Equation;

  using Base = StateVector<Complete<Equation>,
                           ScalarBaseTypeT<Equation, StateType::Complete>>;

  using Base::Base;

  Complete(const Equation& eq [[maybe_unused]]) : Base{} {
    if constexpr (Equation::StaticSize(complete) == dynamic_extent) {
    }
  }
};

template <typename Equation>
struct Conservative
    : StateVector<Conservative<Equation>,
                  ScalarBaseTypeT<Equation, StateType::Conservative>> {
  using EquationType = Equation;

  using Base = StateVector<Conservative<Equation>,
                           ScalarBaseTypeT<Equation, StateType::Conservative>>;
  using Base::Base;

  Conservative(const Equation& eq [[maybe_unused]]) : Base{} {
    if constexpr (Equation::StaticSize(
                      std::integral_constant<StateType, cons>{}) ==
                  dynamic_extent) {
    }
  }
};

template <typename Equation>
Conservative<Equation> AsCons(const Complete<Equation>& complete) {
  constexpr auto names = Conservative<Equation>::Names();
  const auto mapped = boost::hana::to_map(complete);
  return boost::hana::unpack(names, [&](auto... name) {
    return Conservative<Equation>{mapped[name]...};
  });
}

template <typename T, int Width>
constexpr auto ArrayComponentType(basic_type<T>, constant<1>, constant<Width>) {
  return type_c<Array<T, 1, Width>>;
}

template <typename T, int Width>
constexpr auto ArrayComponentType(basic_type<T>, std::ptrdiff_t,
                                  constant<Width>) {
  return type_c<Array<T, Eigen::Dynamic, Width>>;
}

template <typename T, int N, int Width, typename = std::enable_if_t<(N > 0)>>
constexpr auto ArrayComponentType(basic_type<T>, constant<N>, constant<Width>) {
  return type_c<Array<T, N, Width>>;
}

template <int Width, typename Equation, StateType Type>
constexpr auto ArrayBaseType(constant<Width> width, basic_type<Equation>,
                             StateType_c<Type> type) {
  constexpr auto value_types = Equation::ValueTypes(type);
  constexpr auto shape = Equation::Shape(type);
  constexpr auto member_types = boost::hana::transform(
      boost::hana::zip(value_types, shape), [&](auto xs) {
        return ArrayComponentType(at_c<0>(xs), at_c<1>(xs), width);
      });
  constexpr auto template_ = Equation::Template(type);
  return boost::hana::unpack(
      member_types, [&](auto... types) { return template_(types...); });
}

template <int N, typename Equation, StateType Type>
using ArrayBaseTypeT = typename decltype(
    ArrayBaseType(constant<N>{}, type_c<Equation>, constant<Type>{}))::type;

template <typename Equation, int N = kDefaultChunkSize>
struct CompleteArray
    : StateVector<CompleteArray<Equation, N>,
                  ArrayBaseTypeT<N, Equation, StateType::Complete>> {
  using EquationType = Equation;
  using Base = StateVector<CompleteArray<Equation, N>,
                           ArrayBaseTypeT<N, Equation, StateType::Complete>>;
  using Base::Base;
  CompleteArray(const Equation& eq [[maybe_unused]]) : Base{} {
    if constexpr (Equation::StaticSize(
                      std::integral_constant<StateType, cons>{}) ==
                  dynamic_extent) {
    }
  }
};

template <typename Equation, int N = kDefaultChunkSize>
struct ConservativeArray
    : StateVector<ConservativeArray<Equation, N>,
                  ArrayBaseTypeT<N, Equation, StateType::Conservative>> {
  using EquationType = Equation;
  using Base =
      StateVector<ConservativeArray<Equation, N>,
                  ArrayBaseTypeT<N, Equation, StateType::Conservative>>;
  using Base::Base;
  ConservativeArray(const Equation& eq [[maybe_unused]]) : Base{} {
    if constexpr (Equation::StaticSize(
                      std::integral_constant<StateType, cons>{}) ==
                  dynamic_extent) {
    }
  }
};

template <typename StateType, typename F, typename... Ts>
void ForEachVariable(F function, Ts&&... states) {
  constexpr auto names = StateType::Names();
  constexpr auto maps = boost::hana::make_tuple(
      boost::hana::to_map(remove_cvref_t<Ts>::Accessors())...);
  boost::hana::for_each(names, [&](auto name) {
    boost::hana::unpack(boost::hana::zip(maps, boost::hana::make_tuple(
                                                   std::addressof(states)...)),
                        [&](auto... xs) {
                          return function(at_c<0>(xs)[name](*at_c<1>(xs))...);
                        });
  });
}

template <typename Layout, typename T, typename... Extents>
constexpr auto ViewComponent(constant<1>, T* pointer, Extents... extents) {
  return mdspan<T, sizeof...(Extents), Layout>(pointer, extents...);
}

template <typename Layout, int N, typename T, typename... Extents,
          typename = std::enable_if_t<(N > 0)>>
constexpr auto ViewComponent(constant<N>, T* pointer, Extents... extents) {
  return mdspan<T, sizeof...(Extents) + 1, Layout>(pointer, extents..., N);
}

template <typename Layout, typename T, typename... Extents>
constexpr auto ViewComponent(int n, T* pointer, Extents... extents) {
  return mdspan<T, sizeof...(Extents) + 1, Layout>(pointer, extents..., n);
}

template <typename Layout, typename T, int Rank>
constexpr auto ViewComponentType(basic_type<Layout>, constant<1>,
                                 basic_type<T*>, constant<Rank>) {
  return type_c<mdspan<T, Rank, Layout>>;
}

template <typename Layout, typename Depth, typename T, int Rank>
constexpr auto ViewComponentType(basic_type<Layout>, Depth, basic_type<T*>,
                                 constant<Rank>) {
  return type_c<mdspan<T, Rank + 1, Layout>>;
}

template <typename Equation, int Rank, typename Layout>
constexpr auto ViewBase(basic_type<Conservative<Equation>>, constant<Rank> rank,
                        basic_type<Layout> layout) {
  constexpr auto value_types = Equation::ValueTypes(cons);
  constexpr auto shape = Equation::Shape(cons);
  constexpr auto template_ = Equation::Template(cons);
  return boost::hana::unpack(
      boost::hana::zip(value_types, shape), [&](auto... xs) {
        auto view_type = [&](auto x) {
          constexpr auto pointer = boost::hana::traits::add_pointer(at_c<0>(x));
          constexpr auto depth = at_c<1>(x);
          return ViewComponentType(layout, depth, pointer, rank);
        };
        return template_(view_type(xs)...);
      });
};

template <typename Equation, int Rank, typename Layout>
constexpr auto ViewBase(basic_type<const Conservative<Equation>>,
                        constant<Rank> rank, basic_type<Layout> layout) {
  constexpr auto value_types = Equation::ValueTypes(cons);
  constexpr auto shape = Equation::Shape(cons);
  constexpr auto template_ = Equation::Template(cons);
  return boost::hana::unpack(
      boost::hana::zip(value_types, shape), [&](auto... xs) {
        auto view_type = [&](auto x) {
          constexpr auto pointer = boost::hana::traits::add_pointer(
              boost::hana::traits::add_const(at_c<0>(x)));
          constexpr auto depth = at_c<1>(x);
          return ViewComponentType(layout, depth, pointer, rank);
        };
        return template_(view_type(xs)...);
      });
};

template <typename Equation, int Rank, typename Layout>
constexpr auto ViewBase(basic_type<Complete<Equation>>, constant<Rank> rank,
                        basic_type<Layout> layout) {
  constexpr auto value_types = Equation::ValueTypes(complete);
  constexpr auto shape = Equation::Shape(complete);
  constexpr auto template_ = Equation::Template(complete);
  return boost::hana::unpack(
      boost::hana::zip(value_types, shape), [&](auto... xs) {
        auto view_type = [&](auto x) {
          constexpr auto pointer = boost::hana::traits::add_pointer(at_c<0>(x));
          constexpr auto depth = at_c<1>(x);
          return ViewComponentType(layout, depth, pointer, rank);
        };
        return template_(view_type(xs)...);
      });
};

template <typename Equation, int Rank, typename Layout>
constexpr auto ViewBase(basic_type<const Complete<Equation>>,
                        constant<Rank> rank, basic_type<Layout> layout) {
  constexpr auto value_types = Equation::ValueTypes(complete);
  constexpr auto shape = Equation::Shape(complete);
  constexpr auto template_ = Equation::Template(complete);
  return boost::hana::unpack(
      boost::hana::zip(value_types, shape), [&](auto... xs) {
        auto view_type = [&](auto x) {
          constexpr auto pointer = boost::hana::traits::add_pointer(
              boost::hana::traits::add_const(at_c<0>(x)));
          constexpr auto depth = at_c<1>(x);
          return ViewComponentType(layout, depth, pointer, rank);
        };
        return template_(view_type(xs)...);
      });
};

template <typename State, int Rank, typename Layout>
using ViewBaseT = typename decltype(
    ViewBase(type_c<State>, constant<Rank>(), type_c<Layout>))::type;

template <typename State, typename Layout = layout_left,
          int Rank = State::EquationType::Rank()>
struct View : StateBase<ViewBaseT<State, Rank, Layout>> {
  using EquationType = typename State::EquationType;
};

template <typename T> struct IsView : std::false_type {};

template <typename State, typename Layout, int Rank>
struct IsView<View<State, Layout, Rank>> : std::true_type {};

template <template <typename...> typename T, typename... Args>
StateBase<T<remove_cvref_t<Args>...>> MakeTemplate(template_t<T>,
                                                   Args&&... args) {
  return StateBase<T<remove_cvref_t<Args>...>>{std::forward<Args>(args)...};
}

template <typename Equation, typename Layout, int Rank>
View<Conservative<Equation>, Layout, Rank>
AsCons(View<Complete<Equation>, Layout, Rank> complete) {
  View<Conservative<Equation>, Layout, Rank> conservative;
  ForEachVariable<Conservative<Equation>>(
      [&](auto& cons, auto comp) { cons = comp; }, conservative, complete);
  return conservative;
}

template <typename Equation, typename Layout, int Rank>
View<const Conservative<Equation>, Layout, Rank>
AsCons(View<const Complete<Equation>, Layout, Rank> complete) {
  View<const Conservative<Equation>, Layout, Rank> conservative;
  ForEachVariable<Conservative<Equation>>(
      [&](auto& cons, auto comp) { cons = comp; }, conservative, complete);
  return conservative;
}

template <typename Equation, typename Layout, int Rank>
View<const Conservative<Equation>, Layout, Rank>
AsConst(View<Conservative<Equation>, Layout, Rank> view) {
  const auto members = view.Members();
  return boost::hana::unpack(members, [&](auto... mdspan) {
    return View<const Conservative<Equation>, Layout, Rank>{mdspan...};
  });
}

template <typename Equation, typename Layout, int Rank>
View<const Complete<Equation>, Layout, Rank>
AsConst(View<Complete<Equation>, Layout, Rank> view) {
  const auto members = view.Members();
  return boost::hana::unpack(members, [&](auto... mdspan) {
    return View<const Complete<Equation>, Layout, Rank>{mdspan...};
  });
}

template <typename Equation, typename Layout, int Rank>
struct View<Conservative<Equation>, Layout, Rank>
    : StateBase<ViewBaseT<Conservative<Equation>, Rank, Layout>> {
  using EquationType = Equation;
  using Base = StateBase<ViewBaseT<Conservative<Equation>, Rank, Layout>>;
  View() = default;
  View(const View&) = default;
  View(const View<Complete<Equation>, Layout>& complete)
      : View(AsCons(complete)) {}

  template <typename... Xs, typename = std::enable_if_t<(!IsView<Xs>() && ...)>>
  View(Xs&&... xs) : Base{std::forward<Xs>(xs)...} {}
};

template <typename Equation, typename Layout, int Rank>
struct View<const Conservative<Equation>, Layout, Rank>
    : StateBase<ViewBaseT<const Conservative<Equation>, Rank, Layout>> {
  using EquationType = Equation;
  using Base = StateBase<ViewBaseT<const Conservative<Equation>, Rank, Layout>>;
  using Base::Base;
  View() = default;
  View(const View&) = default;
  View(const View<Conservative<Equation>, Layout>& cons)
      : View(AsConst(cons)) {}
  View(const View<Complete<Equation>, Layout>& complete)
      : View(AsCons(complete)) {}
  View(const View<const Complete<Equation>, Layout>& complete)
      : View(AsCons(complete)) {}

  template <typename... Xs, typename = std::enable_if_t<(!IsView<Xs>() && ...)>>
  View(Xs&&... xs) : Base{std::forward<Xs>(xs)...} {}
};

template <typename Equation, typename Layout, int Rank>
struct View<const Complete<Equation>, Layout, Rank>
    : StateBase<ViewBaseT<const Complete<Equation>, Rank, Layout>> {
  using EquationType = Equation;
  using Base = StateBase<ViewBaseT<const Complete<Equation>, Rank, Layout>>;
  using Base::Base;
  View() = default;
  View(const View&) = default;
  View(const View<Complete<Equation>, Layout>& complete)
      : View(AsConst(complete)) {}

  template <typename... Xs,
            typename = std::enable_if_t<(!IsView<remove_cvref_t<Xs>>() && ...)>>
  View(Xs&&... xs) : Base{std::forward<Xs>(xs)...} {}
};

template <typename T> using StridedView = View<T, layout_stride>;

template <int N, typename State, typename Layout, int Rank>
auto Extents(View<State, Layout, Rank> view) {
  return at_c<N>(view.Members()).extents();
}

template <int N, typename State, typename Layout>
auto Mapping(View<State, Layout> view) {
  return at_c<N>(view.Members()).mapping();
}

template <typename D, typename T, typename Scalar,
          typename = std::enable_if_t<std::is_floating_point<Scalar>::value>>
constant<1> Components(const StateVector<D, T>&, const Scalar&) {
  return {};
}

template <typename D, typename Scalar, int M, int Options>
constant<M> Components(const Conservative<D>&,
                       const Eigen::Array<Scalar, 1, M, Options>&) {
  return {};
}

template <typename D, typename Scalar, int M, int Options>
constant<M> Components(const Complete<D>&,
                       const Eigen::Array<Scalar, 1, M, Options>&) {
  return {};
}

template <typename D, typename Scalar, int N, int M, int Options>
constant<M> Components(const ConservativeArray<D, N>&,
                       const Eigen::Array<Scalar, N, M, Options>&) {
  return {};
}

template <typename D, typename Scalar, int N, int M, int Options>
constant<M> Components(const CompleteArray<D, N>&,
                       const Eigen::Array<Scalar, N, M, Options>&) {
  return {};
}

template <typename S, typename L, int Rank, typename T, typename E>
auto Components(View<S, L, Rank>, basic_mdspan<T, E, L> span [[maybe_unused]]) {
  if constexpr (Rank == E::rank()) {
    return constant<1>{};
  } else {
    static_assert(Rank + 1 == E::rank());
    return span.extent(Rank);
  }
}

template <typename D, typename T, typename Scalar,
          typename = std::enable_if_t<std::is_floating_point<Scalar>::value>>
const Scalar& AtComponent(const StateVector<D, T>&, const Scalar& x, int) {
  return x;
}

template <typename D, typename T, typename Scalar,
          typename = std::enable_if_t<std::is_floating_point<Scalar>::value>>
Scalar& AtComponent(const StateVector<D, T>&, Scalar& x, int) {
  return x;
}

template <typename D, typename Scalar, int M, int Options>
const Scalar& AtComponent(const Conservative<D>&,
                          const Eigen::Array<Scalar, 1, M, Options>& x, int i) {
  return x[i];
}

template <typename D, typename Scalar, int M, int Options>
Scalar& AtComponent(const Conservative<D>&,
                    Eigen::Array<Scalar, 1, M, Options>& x, int i) {
  return x[i];
}

template <typename D, typename Scalar, int N, int Options>
Eigen::Array<Scalar, N, 1, Options>&
AtComponent(const ConservativeArray<D, N>&,
            Eigen::Array<Scalar, N, 1, Options>& x, int) {
  return x;
}

template <typename D, typename Scalar, int N, int Options>
const Eigen::Array<Scalar, N, 1, Options>&
AtComponent(const ConservativeArray<D, N>&,
            const Eigen::Array<Scalar, N, 1, Options>& x, int) {
  return x;
}

template <typename D, typename Scalar, int N, int M, int Options>
auto AtComponent(const ConservativeArray<D, N>&,
                 Eigen::Array<Scalar, N, M, Options>& x, int i) {
  return x.col(i);
}

template <typename D, typename Scalar, int N, int M, int Options>
Eigen::Array<Scalar, N, 1, Options>
AtComponent(const ConservativeArray<D, N>&,
            const Eigen::Array<Scalar, N, M, Options>& x, int i) {
  return x.col(i);
}

template <typename D, typename Scalar, int M, int Options>
const Scalar& AtComponent(const Complete<D>&,
                          const Eigen::Array<Scalar, 1, M, Options>& x, int i) {
  return x[i];
}

template <typename D, typename Scalar, int M, int Options>
Scalar& AtComponent(const Complete<D>&, Eigen::Array<Scalar, 1, M, Options>& x,
                    int i) {
  return x[i];
}

template <typename D, typename Scalar, int N, int Options>
Eigen::Array<Scalar, N, 1, Options>&
AtComponent(const CompleteArray<D, N>&, Eigen::Array<Scalar, N, 1, Options>& x,
            int) {
  return x;
}

template <typename D, typename Scalar, int N, int Options>
const Eigen::Array<Scalar, N, 1, Options>&
AtComponent(const CompleteArray<D, N>&,
            const Eigen::Array<Scalar, N, 1, Options>& x, int) {
  return x;
}

template <typename D, typename Scalar, int N, int M, int Options>
auto AtComponent(const CompleteArray<D, N>&,
                 Eigen::Array<Scalar, N, M, Options>& x, int i) {
  return x.col(i);
}

template <typename D, typename Scalar, int N, int M, int Options>
Eigen::Array<Scalar, N, 1, Options>
AtComponent(const CompleteArray<D, N>&,
            const Eigen::Array<Scalar, N, M, Options>& x, int i) {
  return x.col(i);
}

template <typename S, typename L, int Rank, typename T, typename E>
mdspan<T, Rank, L> AtComponent(View<S, L, Rank>, basic_mdspan<T, E, L> span,
                               int i [[maybe_unused]]) {
  if constexpr (Rank == E::rank()) {
    return span;
  } else {
    static_assert(Rank + 1 == E::rank());
    std::array<std::ptrdiff_t, E::rank()> index{};
    index[E::rank() - 1] = i;
    std::array<std::ptrdiff_t, E::rank() - 1> extents;
    std::array<std::ptrdiff_t, E::rank() - 1> strides;
    for (std::size_t r = 0; r < E::rank() - 1; ++r) {
      extents[r] = span.extent(r);
      strides[r] = span.stride(r);
    }
    if constexpr (std::is_same_v<L, layout_stride>) {
      layout_stride::mapping<dynamic_extents<Rank>> mapping{
          dynamic_extents<Rank>(extents), strides};
      return mdspan<T, Rank, L>(&span(index), mapping);
    } else {
      return mdspan<T, Rank, L>(&span(index), extents);
    }
  }
}

template <typename StateType, typename F, typename T, typename... Ts>
void ForEachComponent(F function, T&& state, Ts&&... states) {
  const auto pointers = boost::hana::make_tuple(std::addressof(states)...);
  ForEachVariable<StateType>(
      [&](auto&& var, auto&&... variables) {
        const auto n_components = Components(state, var);
        const auto pv =
            boost::hana::zip(pointers, boost::hana::make_tuple(&variables...));
        for (int i = 0; i < n_components; ++i) {
          boost::hana::unpack(pv, [&](auto... xs) {
            function(AtComponent(state, var, i),
                     AtComponent(*at_c<0>(xs), *at_c<1>(xs), i)...);
          });
        }
      },
      state, states...);
}

template <typename State, typename Layout>
void Load(State& state, View<const State, Layout> view,
          std::array<std::ptrdiff_t, State::EquationType::Rank()> index) {
  ForEachComponent<State>(
      [&](auto& component, auto mdspan) { component = mdspan(index); }, state,
      view);
}

template <typename State, typename Layout>
void Load(State& state, View<State, Layout> view,
          std::array<std::ptrdiff_t, State::EquationType::Rank()> index) {
  ForEachComponent<State>(
      [&](auto& component, auto mdspan) { component = mdspan(index); }, state,
      view);
}

template <typename Eq, typename Layout>
void Store(View<nodeduce_t<Conservative<Eq>>, Layout> view,
           const Conservative<Eq>& state,
           std::array<std::ptrdiff_t, Eq::Rank()> index) {
  ForEachComponent<Conservative<Eq>>(
      [&](auto mdspan, auto block) { Store(mdspan, block, index); }, view,
      state);
}

template <typename Eq, typename Layout>
void Store(View<nodeduce_t<Complete<Eq>>, Layout> view,
           const Complete<Eq>& state,
           std::array<std::ptrdiff_t, Eq::Rank()> index) {
  ForEachComponent<Complete<Eq>>(
      [&](auto mdspan, auto block) { Store(mdspan, block, index); }, view,
      state);
}

template <typename Eq, int N, typename Layout, int Rank>
void Load(ConservativeArray<Eq, N>& state,
          View<nodeduce_t<const Conservative<Eq>>, Layout, Rank> view,
          std::array<std::ptrdiff_t, Rank> index) {
  ForEachComponent<Conservative<Eq>>(
      [&](auto& component, auto mdspan) {
        component = std::apply(
            [&](auto... i) { return Load(constant<N>(), mdspan, i...); },
            index);
      },
      state, view);
}

template <typename Eq, int N, typename Layout, int Rank>
void Load(CompleteArray<Eq, N>& state,
          View<nodeduce_t<const Complete<Eq>>, Layout, Rank> view,
          std::array<std::ptrdiff_t, Rank> index) {
  ForEachComponent<Complete<Eq>>(
      [&](auto& component, auto mdspan) {
        component = Load(constant<N>(), mdspan, index);
      },
      state, view);
}

template <typename Eq, int N, typename Layout, int Rank>
void LoadN(CompleteArray<Eq, N>& state,
           const View<const Complete<Eq>, Layout, Rank>& view, int size,
           const std::array<std::ptrdiff_t, Rank>& pos) {
  ForEachComponent<Complete<Eq>>(
      [&](auto& s, auto v) { s = LoadN(constant<N>{}, v, size, pos); }, state,
      view);
}

template <typename Eq, int N, typename Layout, int Rank>
void LoadN(ConservativeArray<Eq, N>& state,
           const View<const Conservative<Eq>, Layout, Rank>& view, int size,
           const std::array<std::ptrdiff_t, Rank>& pos) {
  ForEachComponent<Conservative<Eq>>(
      [&](auto& s, auto v) { s = LoadN(constant<N>{}, v, size, pos); }, state,
      view);
}

template <typename Eq, int N, typename Layout, int Rank>
void Store(View<nodeduce_t<Conservative<Eq>>, Layout, Rank> view,
           const ConservativeArray<Eq, N>& state,
           std::array<std::ptrdiff_t, Rank> index) {
  ForEachComponent<Conservative<Eq>>(
      [&](auto mdspan, auto block) { Store(mdspan, block, index); }, view,
      state);
}

template <typename Eq, int N, typename Layout, int Rank>
void Store(View<nodeduce_t<Complete<Eq>>, Layout, Rank> view,
           const CompleteArray<Eq, N>& state,
           std::array<std::ptrdiff_t, Rank> index) {
  ForEachComponent<Complete<Eq>>(
      [&](auto mdspan, auto block) { Store(mdspan, block, index); }, view,
      state);
}

template <typename Eq, int N, typename Layout, int Rank>
void StoreN(const View<Complete<Eq>, Layout, Rank>& view,
            const CompleteArray<Eq, N>& state, int size,
            const std::array<std::ptrdiff_t, Rank>& pos) {
  ForEachComponent<Complete<Eq>>(
      [&](auto& s, auto v) { StoreN(v, size, s, pos[0]); }, state, view);
}

template <typename Eq, int N, typename Layout, int Rank>
void StoreN(const View<Conservative<Eq>, Layout, Rank>& view,
            const ConservativeArray<Eq, N>& state, int size,
            const std::array<std::ptrdiff_t, Rank>& pos) {
  ForEachComponent<Conservative<Eq>>(
      [&](auto& s, auto v) { StoreN(v, size, s, pos[0]); }, state, view);
}

template <typename State, typename Layout, int Rank,
          typename... SliceSpecifiers>
auto Subspan(const View<State, Layout, Rank>& view, SliceSpecifiers... slices) {
  static constexpr int R = SliceRank_<SliceSpecifiers...>;
  auto v = view.Members();
  return boost::hana::unpack(v, [&](auto... vs) {
    return View<State, layout_stride, R>{subspan(vs, slices...)...};
  });
}

template <typename F, typename View, typename... Views>
void ForEachRow(F function, View view, Views... views) {
  auto extents = Extents<0>(view);
  static constexpr int Rank = remove_cvref_t<decltype(extents)>::rank();
  if constexpr (Rank == 1) {
    function(view, views...);
  } else if constexpr (Rank == 2) {
    for (int j = 0; j < Extents<0>(view).extent(1); ++j) {
      function(Subspan(view, all, j), Subspan(views, all, j)...);
    }
  } else {
    static_assert(Rank == 3);
    for (int k = 0; k < Extents<0>(view).extent(2); ++k) {
      for (int j = 0; j < Extents<0>(view).extent(1); ++j) {
        function(Subspan(view, all, j, k), Subspan(views, all, j, k)...);
      }
    }
  }
}

template <Direction dir, typename T, typename E, typename L, typename A, typename SliceSpecifier>
auto Slice(basic_mdspan<T, E, L, A> span, SliceSpecifier slice) {
    if constexpr (E::rank() == 1) {
        static_assert(dir == Direction::X);
        return subspan(span, slice);
    } else if constexpr (E::rank() == 2) {
        if constexpr (dir == Direction::X) {
            return subspan(span, slice, all);
        } else {
            static_assert(dir == Direction::Y);
            return subspan(span, all, slice);
        }
    } else if constexpr (E::rank() == 3) {
        if constexpr (dir == Direction::X) {
            return subspan(span, slice, all, all);
        } else if constexpr (dir == Direction::Y) {
            return subspan(span, all, slice, all);
        } else {
            static_assert(dir == Direction::Z);
            return subspan(span, all, all, slice);
        }
    } else if constexpr (E::rank() == 4) {
        if constexpr (dir == Direction::X) {
            return subspan(span, slice, all, all, all);
        } else if constexpr (dir == Direction::Y) {
            return subspan(span, all, slice, all, all);
        } else {
            static_assert(dir == Direction::Z);
            return subspan(span, all, all, slice, all);
        }
    }
}


template <Direction dir, typename T, typename L, int Rank, typename SliceSpecifier>
auto Slice(View<T, L, Rank> view, SliceSpecifier slice) {
    View<T, layout_stride, Rank> sliced_view;
    ForEachVariable<remove_cvref_t<T>>([&](auto& s, auto v) {
        s = Slice<dir>(v, slice);
    }, sliced_view, view);
    return sliced_view;
}

} // namespace fub

#endif