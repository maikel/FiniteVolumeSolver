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

template <typename Equation, StateType Type>
constexpr auto ScalarBaseType() {
  constexpr constant<Type> type{};
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
using ScalarBaseTypeT = typename decltype(ScalarBaseType<Equation, Type>())::type;


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
struct Complete
    : StateVector<Complete<Equation>, ScalarBaseTypeT<Equation, StateType::Complete>> {
  using EquationType = Equation;

  using Base =
      StateVector<Complete<Equation>, ScalarBaseTypeT<Equation, StateType::Complete>>;

  using Base::Base;

  Complete(const Equation& eq [[maybe_unused]]) : Base{} {
    if constexpr (Equation::StaticSize(complete) == dynamic_extent) {
    }
  }
};

template <typename Equation>
struct Conservative
    : StateVector<Conservative<Equation>, ScalarBaseTypeT<Equation, StateType::Conservative>> {
  using EquationType = Equation;

  using Base =
      StateVector<Conservative<Equation>, ScalarBaseTypeT<Equation, StateType::Conservative>>;
  using Base::Base;

  Conservative(const Equation& eq [[maybe_unused]]) : Base{} {
    if constexpr (Equation::StaticSize(cons) == dynamic_extent) {
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
                             constant<Type> type) {
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


template <typename Equation, int N = kChunkSize>
struct CompleteArray : StateVector<CompleteArray<Equation, N>,
                                   ArrayBaseTypeT<N, Equation, StateType::Complete>> {
  using EquationType = Equation;
};

template <typename Equation, int N = kChunkSize>
struct ConsArray
    : StateVector<ConsArray<Equation, N>, ArrayBaseTypeT<N, Equation, StateType::Conservative>> {
  using EquationType = Equation;
};

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

template <typename State, typename Layout = layout_left>
struct View : StateBase<ViewBaseT<State, State::EquationType::Rank(), Layout>> {
  using EquationType = typename State::EquationType;
};

template <typename T> struct IsView : std::false_type {};

template <typename State, typename Layout>
struct IsView<View<State, Layout>> : std::true_type {};

template <template <typename...> typename T, typename... Args>
StateBase<T<remove_cvref_t<Args>...>> MakeTemplate(template_t<T>,
                                                   Args&&... args) {
  return StateBase<T<remove_cvref_t<Args>...>>{std::forward<Args>(args)...};
}

template <typename Equation, typename Layout>
View<Conservative<Equation>, Layout>
AsCons(View<Complete<Equation>, Layout> view) {
  const auto map = boost::hana::to_map(view);
  constexpr auto names = Conservative<Equation>::Names();
  return boost::hana::unpack(names, [&](auto... name) {
    return View<Conservative<Equation>, Layout>{map[name]...};
  });
}

template <typename Equation, typename Layout>
View<const Conservative<Equation>, Layout>
AsCons(View<const Complete<Equation>, Layout> view) {
  const auto map = boost::hana::to_map(view);
  constexpr auto names = Conservative<Equation>::Names();
  return boost::hana::unpack(names, [&](auto... name) {
    return View<const Conservative<Equation>, Layout>{map[name]...};
  });
}

template <typename Equation, typename Layout>
View<const Conservative<Equation>, Layout>
AsConst(View<Conservative<Equation>, Layout> view) {
  const auto members = view.Members();
  return boost::hana::unpack(members, [&](auto... mdspan) {
    return View<const Conservative<Equation>, Layout>{mdspan...};
  });
}

template <typename Equation, typename Layout>
View<const Complete<Equation>, Layout>
AsConst(View<Complete<Equation>, Layout> view) {
  const auto members = view.Members();
  return boost::hana::unpack(members, [&](auto... mdspan) {
    return View<const Complete<Equation>, Layout>{mdspan...};
  });
}

template <typename Equation, typename Layout>
struct View<Conservative<Equation>, Layout>
    : StateBase<ViewBaseT<Conservative<Equation>, Equation::Rank(), Layout>> {
  using EquationType = Equation;
  using Base =
      StateBase<ViewBaseT<Conservative<Equation>, Equation::Rank(), Layout>>;
  View() = default;
  View(const View&) = default;
  View(const View<Complete<Equation>, Layout>& complete)
      : View(AsCons(complete)) {}

  template <typename... Xs, typename = std::enable_if_t<(!IsView<Xs>() && ...)>>
  View(Xs&&... xs) : Base{std::forward<Xs>(xs)...} {}
};

template <typename Equation, typename Layout>
struct View<const Conservative<Equation>, Layout>
    : StateBase<
          ViewBaseT<const Conservative<Equation>, Equation::Rank(), Layout>> {
  using EquationType = Equation;
  using Base = StateBase<
      ViewBaseT<const Conservative<Equation>, Equation::Rank(), Layout>>;
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

template <typename Equation, typename Layout>
struct View<const Complete<Equation>, Layout>
    : StateBase<ViewBaseT<const Complete<Equation>, Equation::Rank(), Layout>> {
  using EquationType = Equation;
  using Base =
      StateBase<ViewBaseT<const Complete<Equation>, Equation::Rank(), Layout>>;
  using Base::Base;
  View() = default;
  View(const View&) = default;
  View(const View<Complete<Equation>, Layout>& complete)
      : View(AsConst(complete)) {}

  template <typename... Xs, typename = std::enable_if_t<(!IsView<Xs>() && ...)>>
  View(Xs&&... xs) : Base{std::forward<Xs>(xs)...} {}
};

template <typename T> using StridedView = View<T, layout_stride>;

template <int N, typename State, typename Layout>
auto Extents(View<State, Layout> view) {
  return at_c<N>(view.Members()).extents();
}

template <int N, typename State, typename Layout>
auto Mapping(View<State, Layout> view) {
  return at_c<N>(view.Members()).mapping();
}

template <typename StateType, typename F, typename... Ts>
void ForEachVariable(F function, Ts&&... states) {
  constexpr auto names = StateType::Names();
  constexpr auto maps =
      boost::hana::make_tuple(boost::hana::to_map(remove_cvref_t<Ts>::Accessors())...);
  boost::hana::for_each(names, [&](auto name) {
    boost::hana::unpack(boost::hana::zip(maps, boost::hana::make_tuple(
                                                   std::addressof(states)...)),
                        [&](auto... xs) {
                          return function(at_c<0>(xs)[name](*at_c<1>(xs))...);
                        });
  });
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
constant<M> Components(const ConsArray<D, N>&,
                       const Eigen::Array<Scalar, N, M, Options>&) {
  return {};
}

template <typename D, typename Scalar, int N, int M, int Options>
constant<M> Components(const CompleteArray<D, N>&,
                       const Eigen::Array<Scalar, N, M, Options>&) {
  return {};
}

template <typename S, typename L, typename T, typename E>
auto Components(View<S, L>, basic_mdspan<T, E, L> span [[maybe_unused]]) {
  if constexpr (S::EquationType::Rank() == E::rank()) {
    return constant<1>{};
  } else {
    static_assert(S::EquationType::Rank() + 1 == E::rank());
    return span.extent(S::EquationType::Rank());
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
AtComponent(const ConsArray<D, N>&, Eigen::Array<Scalar, N, 1, Options>& x,
            int) {
  return x;
}

template <typename D, typename Scalar, int N, int Options>
const Eigen::Array<Scalar, N, 1, Options>&
AtComponent(const ConsArray<D, N>&,
            const Eigen::Array<Scalar, N, 1, Options>& x, int) {
  return x;
}

template <typename D, typename Scalar, int N, int M, int Options>
auto AtComponent(const ConsArray<D, N>&, Eigen::Array<Scalar, N, M, Options>& x,
                 int i) {
  return x.col(i);
}

template <typename D, typename Scalar, int N, int M, int Options>
Eigen::Array<Scalar, N, 1, Options>
AtComponent(const ConsArray<D, N>&,
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

template <typename S, typename L, typename T, typename E>
mdspan<T, S::EquationType::Rank(), L>
AtComponent(View<S, L>, basic_mdspan<T, E, L> span, int i [[maybe_unused]]) {
  if constexpr (S::EquationType::Rank() == E::rank()) {
    return span;
  } else {
    static_assert(S::EquationType::Rank() + 1 == E::rank());
    std::array<std::ptrdiff_t, E::rank()> index{};
    index[E::rank() - 1] = i;
    std::array<std::ptrdiff_t, E::rank() - 1> extents;
    std::array<std::ptrdiff_t, E::rank() - 1> strides;
    for (std::size_t r = 0; r < E::rank() - 1; ++r) {
      extents[r] = span.extent(r);
      strides[r] = span.stride(r);
    }
    constexpr int Rank = E::rank() - 1;
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

// template <typename Eq, int N, typename Layout>
// void StoreArray(View<nodeduce_t<Conservative<Eq>>, Layout> view,
//                 const ConsArray<Eq, N>& state,
//                 std::array<std::ptrdiff_t, Eq::Rank()> index) {
//   ForEachComponent(
//       [&](auto mdspan, auto block) { Store(mdspan, block, index); }, view,
//       state);
// }


// template <typename Eq, int N, typename Layout>
// void LoadArray(ConsArray<Eq, N>& state,
//                View<nodeduce_t<const Conservative<Eq>>, Layout> view,
//                std::array<std::ptrdiff_t, Eq::Rank()> index) {
//   ForEachComponent(
//       [&](auto component, auto mdspan) {
//         component = std::apply(
//             [&](auto... i) { return Load(constant<N>(), mdspan, i...); },
//             index);
//       },
//       state, view);
// }

} // namespace fub

#endif