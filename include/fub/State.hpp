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
constexpr auto ScalarBaseType(basic_type<Equation>, constant<Type> type) {
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

template <typename Derived, typename Base>
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
};

template <typename Equation>
struct Complete
    : StateVector<Complete<Equation>, typename decltype(ScalarBaseType(
                                          type_c<Equation>, complete))::type> {
  using EquationType = Equation;
};

template <typename Equation>
struct Cons : StateVector<Cons<Equation>, typename decltype(ScalarBaseType(
                                              type_c<Equation>, cons))::type> {
  using EquationType = Equation;
};

template <typename Equation>
Cons<Equation> AsCons(const Complete<Equation>& complete) {
  constexpr auto names = Cons<Equation>::Names();
  const auto mapped = boost::hana::to_map(complete);
  return boost::hana::unpack(
      names, [&](auto... name) { return Cons<Equation>{mapped[name]...}; });
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

template <typename Equation, int N = kChunkSize>
struct CompleteArray : StateVector<CompleteArray<Equation, N>,
                                   ArrayBaseTypeT<N, Equation, complete>> {
  using EquationType = Equation;
};

template <typename Equation, int N = kChunkSize>
struct ConsArray
    : StateVector<ConsArray<Equation, N>, ArrayBaseTypeT<N, Equation, cons>> {
  using EquationType = Equation;
};

template <typename Layout, typename T, typename... Extents>
constexpr auto ViewComponent(constant<1>, T* pointer, Extents... extents) {
  return dynamic_mdspan<T, sizeof...(Extents), Layout>(pointer, extents...);
}

template <typename Layout, int N, typename T, typename... Extents,
          typename = std::enable_if_t<(N > 0)>>
constexpr auto ViewComponent(constant<N>, T* pointer, Extents... extents) {
  return dynamic_mdspan<T, sizeof...(Extents) + 1, Layout>(pointer, extents...,
                                                           N);
}

template <typename Layout, typename T, typename... Extents>
constexpr auto ViewComponent(int n, T* pointer, Extents... extents) {
  return dynamic_mdspan<T, sizeof...(Extents) + 1, Layout>(pointer, extents...,
                                                           n);
}

template <typename Layout, typename T, int Rank>
constexpr auto ViewComponentType(basic_type<Layout>, constant<1>,
                                 basic_type<T*>, constant<Rank>) {
  return type_c<dynamic_mdspan<T, Rank, Layout>>;
}

template <typename Layout, typename Depth, typename T, int Rank>
constexpr auto ViewComponentType(basic_type<Layout>, Depth, basic_type<T*>,
                                 constant<Rank>) {
  return type_c<dynamic_mdspan<T, Rank + 1, Layout>>;
}

template <typename Equation, int Rank, typename Layout>
constexpr auto ViewBase(basic_type<Cons<Equation>>, constant<Rank> rank,
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
constexpr auto ViewBase(basic_type<const Cons<Equation>>, constant<Rank> rank,
                        basic_type<Layout> layout) {
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
StateBase<T<remove_cvref_t<Args>...>> MakeTemplate(template_t<T> template_,
                                                   Args&&... args) {
  return StateBase<T<remove_cvref_t<Args>...>>{std::forward<Args>(args)...};
}

template <typename Equation, typename Layout>
View<Cons<Equation>, Layout> AsCons(View<Complete<Equation>, Layout> view) {
  const auto map = boost::hana::to_map(view);
  constexpr auto names = Cons<Equation>::Names();
  return boost::hana::unpack(names, [&](auto... name) {
    return View<Cons<Equation>, Layout>{map[name]...};
  });
}

template <typename Equation, typename Layout>
View<const Cons<Equation>, Layout>
AsCons(View<const Complete<Equation>, Layout> view) {
  const auto map = boost::hana::to_map(view);
  constexpr auto names = Cons<Equation>::Names();
  return boost::hana::unpack(names, [&](auto... name) {
    return View<const Cons<Equation>, Layout>{map[name]...};
  });
}

template <typename Equation, typename Layout>
View<const Cons<Equation>, Layout> AsConst(View<Cons<Equation>, Layout> view) {
  const auto members = view.Members();
  return boost::hana::unpack(members, [&](auto... mdspan) {
    return View<const Cons<Equation>, Layout>{mdspan...};
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
struct View<Cons<Equation>, Layout>
    : StateBase<ViewBaseT<Cons<Equation>, Equation::Rank(), Layout>> {
  using EquationType = Equation;
  using Base = StateBase<ViewBaseT<Cons<Equation>, Equation::Rank(), Layout>>;
  View() = default;
  View(const View&) = default;
  View(const View<Complete<Equation>, Layout>& complete)
      : View(AsCons(complete)) {}

  template <typename... Xs,
            typename = std::enable_if_t<
                (!IsView<Xs>() &&
                 ...)>> //,
                        // typename = std::enable_if_t<std::is_convertible<
                        //     typename Base::MemberTypes,
                        //     boost::hana::tuple<Xs...>>::value>>
  View(Xs&&... xs) : Base{std::forward<Xs>(xs)...} {}
};

template <typename Equation, typename Layout>
struct View<const Cons<Equation>, Layout>
    : StateBase<ViewBaseT<const Cons<Equation>, Equation::Rank(), Layout>> {
  using EquationType = Equation;
  using Base =
      StateBase<ViewBaseT<const Cons<Equation>, Equation::Rank(), Layout>>;
  using Base::Base;
  View() = default;
  View(const View&) = default;
  View(const View<Cons<Equation>, Layout>& cons) : View(AsConst(cons)) {}
  View(const View<Complete<Equation>, Layout>& complete)
      : View(AsCons(complete)) {}
  View(const View<const Complete<Equation>, Layout>& complete)
      : View(AsCons(complete)) {}

  template <typename... Xs,
            typename = std::enable_if_t<
                (!IsView<Xs>() &&
                 ...)>> //,
                        // typename = std::enable_if_t<std::is_convertible<
                        //     typename Base::MemberTypes,
                        //     boost::hana::tuple<Xs...>>::value>>
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

  template <typename... Xs,
            typename = std::enable_if_t<
                (!IsView<Xs>() &&
                 ...)>> //,
                        // typename = std::enable_if_t<std::is_convertible<
                        //     typename Base::MemberTypes,
                        //     boost::hana::tuple<Xs...>>::value>>
  View(Xs&&... xs) : Base{std::forward<Xs>(xs)...} {}
};

template <typename T> using StridedView = View<T, layout_stride>;

template <typename F, typename... Ts>
void ForEachVariable(F function, Ts&&... states) {
  constexpr auto pointers =
      boost::hana::zip(remove_cvref_t<Ts>::PointersToMember()...);
  boost::hana::for_each(pointers, [&](auto pointers_to_variable) {
    boost::hana::unpack(
        boost::hana::zip(pointers_to_variable,
                         boost::hana::make_tuple(std::addressof(states)...)),
        [&](auto... xs) { return function(at_c<0>(xs)(*at_c<1>(xs))...); });
  });
}

template <typename D, typename T, typename Scalar,
          typename = std::enable_if_t<std::is_floating_point<Scalar>::value>>
constant<1> Components(const StateVector<D, T>&, const Scalar&) {
  return {};
}

template <typename D, typename Scalar, int M, int Options>
constant<M> Components(const Cons<D>&,
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
auto Components(View<S, L>, basic_mdspan<T, E, L> mdspan) {
  if constexpr (S::EquationType::Rank() == E::rank()) {
    return constant<1>{};
  } else {
    static_assert(S::EquationType::Rank() + 1 == E::rank());
    return mdspan.extent(S::EquationType::Rank());
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
const Scalar& AtComponent(const Cons<D>&,
                          const Eigen::Array<Scalar, 1, M, Options>& x, int i) {
  return x[i];
}

template <typename D, typename Scalar, int M, int Options>
Scalar& AtComponent(const Cons<D>&, Eigen::Array<Scalar, 1, M, Options>& x,
                    int i) {
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
dynamic_mdspan<T, S::EquationType::Rank(), L>
AtComponent(View<S, L>, basic_mdspan<T, E, L> mdspan, int i) {
  if constexpr (S::EquationType::Rank() == E::rank()) {
    return mdspan;
  } else {
    static_assert(S::EquationType::Rank() + 1 == E::rank());
    std::array<std::ptrdiff_t, E::rank()> index{};
    index[E::rank() - 1] = i;
    std::array<std::ptrdiff_t, E::rank() - 1> extents;
    std::array<std::ptrdiff_t, E::rank() - 1> strides;
    for (std::size_t r = 0; r < E::rank() - 1; ++r) {
      extents[r] = mdspan.extent(r);
      strides[r] = mdspan.stride(r);
    }
    constexpr int Rank = E::rank() - 1;
    if constexpr (std::is_same_v<L, layout_stride>) {
      layout_stride::mapping<DynamicExtents<Rank>> mapping{
          DynamicExtents<Rank>(extents), strides};
      return dynamic_mdspan<T, Rank, L>(&mdspan(index), mapping);
    } else {
      return dynamic_mdspan<T, Rank, L>(&mdspan(index), extents);
    }
  }
}

template <typename F, typename T, typename... Ts>
void ForEachComponent(F function, T&& state, Ts&&... states) {
  const auto pointers = boost::hana::make_tuple(std::addressof(states)...);
  ForEachVariable(
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
  ForEachComponent(
      [&](auto& component, auto mdspan) { component = mdspan(index); }, state,
      view);
}

template <typename State, typename Layout>
void Load(State& state, View<State, Layout> view,
          std::array<std::ptrdiff_t, State::EquationType::Rank()> index) {
  ForEachComponent(
      [&](auto& component, auto mdspan) { component = mdspan(index); }, state,
      view);
}

template <typename Eq, int N, typename Layout>
void LoadArray(ConsArray<Eq, N>& state,
               View<nodeduce_t<const Cons<Eq>>, Layout> view,
               std::array<std::ptrdiff_t, Eq::Rank()> index) {
  ForEachComponent(
      [&](auto component, auto mdspan) {
        component = std::apply(
            [&](auto... i) { return Load(constant<N>(), mdspan, i...); },
            index);
      },
      state, view);
}

template <typename Eq, typename Layout>
void Store(View<nodeduce_t<Cons<Eq>>, Layout> view, const Cons<Eq>& state,
           std::array<std::ptrdiff_t, Eq::Rank()> index) {
  ForEachComponent(
      [&](auto mdspan, auto block) { Store(mdspan, block, index); }, view,
      state);
}

template <typename Eq, typename Layout>
void Store(View<nodeduce_t<Complete<Eq>>, Layout> view,
           const Complete<Eq>& state,
           std::array<std::ptrdiff_t, Eq::Rank()> index) {
  ForEachComponent(
      [&](auto mdspan, auto block) { Store(mdspan, block, index); }, view,
      state);
}

template <typename Eq, int N, typename Layout>
void StoreArray(View<nodeduce_t<Cons<Eq>>, Layout> view,
                const ConsArray<Eq, N>& state,
                std::array<std::ptrdiff_t, Eq::Rank()> index) {
  ForEachComponent(
      [&](auto mdspan, auto block) { Store(mdspan, block, index); }, view,
      state);
}

template <typename State, typename Layout>
auto Extents(View<State, Layout> view) {
  return at_c<0>(view.Members()).extents();
}

template <typename State, typename Layout>
auto Mapping(View<State, Layout> view) {
  return at_c<0>(view.Members()).mapping();
}

} // namespace fub

#endif