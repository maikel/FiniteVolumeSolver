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

#include "fub/PatchDataView.hpp"
#include "fub/core/tuple.hpp"
#include "fub/ext/Eigen.hpp"

#include <boost/mp11/algorithm.hpp>
#include <boost/mp11/tuple.hpp>

#include <numeric>

namespace fub {

template <typename T> struct StateTraits;

namespace meta {

template <typename Depths> struct Rank;

}

namespace detail {
template <std::size_t I, typename... Ts> constexpr auto ZipNested(Ts&&... ts) {
  return std::make_tuple(std::get<I>(ts)...);
}

template <std::size_t... Is, typename... Ts>
constexpr auto ZipHelper(std::index_sequence<Is...>, Ts&&... ts) {
  return std::make_tuple(ZipNested<Is>(std::forward<Ts>(ts)...)...);
}

template <std::size_t I, typename... Ts>
constexpr auto UnzipNested(const std::tuple<Ts...>& ts) {
  return std::apply([](const auto& t) { return std::get<I>(t); }, ts);
}

template <std::size_t... Is, typename... Ts>
constexpr auto UnzipHelper(std::index_sequence<Is...>,
                           const std::tuple<Ts...>& ts) {
  return std::make_tuple(UnzipNested<Is>(ts)...);
}
} // namespace detail

template <typename... Ts> constexpr auto Zip(Ts&&... ts) {
  using First =
      std::decay_t<boost::mp11::mp_front<boost::mp11::mp_list<Ts...>>>;
  return detail::ZipHelper(
      std::make_index_sequence<std::tuple_size<First>::value>(),
      std::forward<Ts>(ts)...);
}

template <typename... Ts>
constexpr auto Unzip(const std::tuple<Ts...>& zipped) {
  using First =
      std::decay_t<boost::mp11::mp_front<boost::mp11::mp_list<Ts...>>>;
  return detail::UnzipHelper(
      std::make_index_sequence<std::tuple_size<First>::value>(), zipped);
}

template <std::size_t I, typename State>
constexpr decltype(auto) get(State&& x) {
  using S = remove_cvref_t<State>;
  constexpr auto& pointers_to_member = S::Traits::pointers_to_member;
  return x.*std::get<I>(pointers_to_member);
}

template <typename F, typename... Ts>
void ForEachVariable(F function, Ts&&... states) {
  constexpr auto pointers =
      Zip(StateTraits<remove_cvref_t<Ts>>::pointers_to_member...);
  const auto xs = std::make_tuple(std::addressof(states)...);
  boost::mp11::tuple_for_each(pointers, [&function, &xs](const auto& ps) {
    boost::mp11::tuple_apply(
        [&function](auto... ys) {
          function((*std::get<0>(ys)).*std::get<1>(ys)...);
        },
        Zip(xs, ps));
  });
}

template <typename T>
struct NumVariables
    : std::integral_constant<
          std::size_t,
          std::tuple_size_v<std::decay_t<decltype(T::Traits::names)>>> {};

/// This type is used to tag scalar quantities
struct ScalarDepth : std::integral_constant<int, 1> {};

/// This type is used to tag quantities with a depth known at compile time.
template <int Depth> struct VectorDepth : std::integral_constant<int, Depth> {};

template <typename Depth> struct ToConcreteDepthImpl { using type = Depth; };

template <> struct ToConcreteDepthImpl<VectorDepth<-1>> { using type = int; };

template <typename Depth>
using ToConcreteDepth = typename ToConcreteDepthImpl<Depth>::type;

template <typename Depths>
using ToConcreteDepths = boost::mp11::mp_transform<ToConcreteDepth, Depths>;

namespace detail {
/// @{
/// This meta-function transforms a Depth into a value type used in a states or
/// state arrays.
template <typename T, int Width> struct DepthToStateValueTypeImpl;

/// The value-type for a single-cell state for a scalar quantity.
template <> struct DepthToStateValueTypeImpl<ScalarDepth, 1> {
  using type = double;
};

/// The value-type for a single-cell state for a scalar quantity.
template <int Width> struct DepthToStateValueTypeImpl<ScalarDepth, Width> {
  using type = Array<double, 1, Width>;
};

/// The value-type for a single-cell state for a scalar quantity.
template <int Depth, int Width>
struct DepthToStateValueTypeImpl<VectorDepth<Depth>, Width> {
  using type = Array<double, Depth, Width>;
};

/// For usage with boost::mp11
template <typename T>
using DepthToStateValueType = typename DepthToStateValueTypeImpl<T, 1>::type;
/// @}

template <typename T, typename Eq>
using DepthsT = typename StateTraits<T>::template Depths<Eq::Rank()>;

template <typename T, typename Eq> struct DepthsImpl {
  constexpr typename StateTraits<T>::template Depths<Eq::Rank()>
  operator()(const Eq&) const noexcept {
    return {};
  }
};
} // namespace detail

template <typename T> struct Type {};
inline constexpr struct DepthsFn {
  template <typename Eq, typename State>
  constexpr auto operator()(const Eq& eq,
                            [[maybe_unused]] Type<State> state) const
      noexcept(
          is_nothrow_tag_invocable<DepthsFn, const Eq&, Type<State>>::value ||
          !is_tag_invocable<DepthsFn, const Eq&, Type<State>>::value) {
    if constexpr (is_tag_invocable<DepthsFn, const Eq&, Type<State>>::value) {
      return fub::tag_invoke(*this, eq, std::move(state));
    } else if constexpr (is_detected<detail::DepthsT, State, Eq>::value) {
      detail::DepthsImpl<State, Eq> default_depths;
      return default_depths(eq);
    } else {
      static_assert(
          is_tag_invocable<DepthsFn, const Eq&, Type<State>>::value ||
          is_detected<detail::DepthsT, State, Eq>::value);
    }
  }
} Depths;

namespace meta {
template <typename T>
using Depths = decltype(
    ::fub::Depths(std::declval<typename T::Equation const&>(), Type<T>{}));
}

template <typename Depths>
using ScalarStateBase =
    boost::mp11::mp_transform<detail::DepthToStateValueType, Depths>;

template <typename Depths> struct ScalarState : ScalarStateBase<Depths> {
  using Traits = StateTraits<ScalarStateBase<Depths>>;
  using Base = ScalarStateBase<Depths>;

  using Base::Base;

  template <typename Equation>
  ScalarState(const Equation& eq) : Base{} {
    auto depths = ::fub::Depths(eq, Type<ScalarState>{});
    ForEachVariable(
        overloaded{
            [&](double& id, ScalarDepth) { id = 0.0; },
            [&](auto&& ids, auto depth) {
              if constexpr (std::is_same_v<std::decay_t<decltype(depth)>,
                                           int>) {
                ids = Array<double, 1, Eigen::Dynamic>::Zero(1, depth);
              } else {
                ids = Array<double, 1, decltype(depth)::value>::Zero();
              }
            },
        },
        *this, depths);
  }

  ScalarState& operator+=(const Base& other) {
    ForEachVariable([](auto&& that, auto&& other) { that += other; }, *this,
                    other);
    return *this;
  }

  ScalarState& operator-=(const Base& other) {
    ForEachVariable([](auto&& that, auto&& other) { that -= other; }, *this,
                    other);
    return *this;
  }

  ScalarState& operator*=(double lambda) {
    ForEachVariable([lambda](auto&& that, auto&& other) { that *= lambda; },
                    *this);
    return *this;
  }
};

template <typename Depths>
struct StateTraits<ScalarState<Depths>> : StateTraits<ScalarStateBase<Depths>> {
};

/// This type alias transforms state depths into a conservative state
/// associated with a specified equation.
template <typename Equation>
using ConservativeBase =
    boost::mp11::mp_transform<detail::DepthToStateValueType,
                              typename Equation::ConservativeDepths>;

/// This type has a constructor which takes an equation and might allocate any
/// dynamically sized member variable.
template <typename Eq> struct Conservative : ConservativeBase<Eq> {
  using Equation = Eq;
  using ValueType = double;
  using Traits = StateTraits<ConservativeBase<Equation>>;

  Conservative() = default;
  Conservative(const ConservativeBase<Eq>& x) : ConservativeBase<Eq>{x} {}
  Conservative& operator=(const ConservativeBase<Eq>& x) {
    static_cast<ConservativeBase<Eq>&>(*this) = x;
    return *this;
  }
  Conservative(const Equation& eq) : ConservativeBase<Eq>{} {
    auto depths = ::fub::Depths(eq, Type<Conservative<Equation>>{});
    ForEachVariable(
        overloaded{
            [&](double& id, ScalarDepth) { id = 0.0; },
            [&](auto&& ids, auto depth) {
              if constexpr (std::is_same_v<std::decay_t<decltype(depth)>,
                                           int>) {
                ids = Array<double, 1, Eigen::Dynamic>::Zero(1, depth);
              } else {
                ids = Array<double, 1, decltype(depth)::value>::Zero();
              }
            },
        },
        *this, depths);
  }

  Conservative& operator+=(const Conservative& other) {
    ForEachVariable([](auto&& that, auto&& other) { that += other; }, *this,
                    other);
    return *this;
  }
};

template <typename Eq>
struct StateTraits<Conservative<Eq>> : StateTraits<ConservativeBase<Eq>> {};

template <typename Eq>
struct Primitive : ScalarState<typename Eq::PrimitiveDepths> {
  using Base = ScalarState<typename Eq::PrimitiveDepths>;

  using Base::Base;
};

template <typename Eq>
struct StateTraits<Primitive<Eq>>
    : StateTraits<ScalarStateBase<typename Eq::PrimitiveDepths>> {};

template <typename Eq>
struct KineticState : ScalarState<typename Eq::KineticStateDepths> {
  using Base = ScalarState<typename Eq::KineticStateDepths>;

  using Base::Base;
};

template <typename Eq>
struct StateTraits<KineticState<Eq>>
    : StateTraits<ScalarStateBase<typename Eq::KineticStateDepths>> {};

template <typename Eq>
struct Characteristics : ScalarState<typename Eq::CharacteristicsDepths> {
  using Base = ScalarState<typename Eq::CharacteristicsDepths>;

  using Base::Base;
};

template <typename Eq>
struct StateTraits<Characteristics<Eq>>
    : StateTraits<ScalarStateBase<typename Eq::CharacteristicsDepths>> {};

template <typename State> auto StateToTuple(const State& x) {
  return boost::mp11::tuple_apply(
      [&x](auto... pointer) { return std::make_tuple((x.*pointer)...); },
      StateTraits<State>::pointers_to_member);
}

template <typename Equation>
using CompleteBase =
    boost::mp11::mp_transform<detail::DepthToStateValueType,
                              typename Equation::CompleteDepths>;

/// This type has a constructor which takes an equation and might allocate any
/// dynamically sized member variable.
template <typename Eq> struct Complete : CompleteBase<Eq> {
  using Equation = Eq;
  using ValueType = double;
  using Traits = StateTraits<CompleteBase<Equation>>;

  Complete() = default;
  Complete(const CompleteBase<Eq>& x) : CompleteBase<Eq>{x} {}
  Complete& operator=(const CompleteBase<Eq>& x) {
    static_cast<CompleteBase<Eq>&>(*this) = x;
  }
  Complete(const Equation& eq) : CompleteBase<Eq>{} {
    auto depths = ::fub::Depths(eq, Type<Complete<Equation>>{});
    ForEachVariable(
        overloaded{
            [&](double& id, ScalarDepth) { id = 0.0; },
            [&](auto&& ids, auto depth) {
              if constexpr (std::is_same_v<std::decay_t<decltype(depth)>,
                                           int>) {
                ids = Array<double, 1, Eigen::Dynamic>::Zero(1, depth);
              } else {
                ids = Array<double, 1, decltype(depth)::value>::Zero();
              }
            },
        },
        *this, depths);
  }

  Complete& operator+=(const Complete& other) {
    ForEachVariable([](auto&& that, auto&& other) { that += other; }, *this,
                    other);
    return *this;
  }
};

template <typename Eq>
struct StateTraits<Complete<Eq>> : StateTraits<CompleteBase<Eq>> {};

///// View type
namespace detail {
template <typename Depth, typename T, int Rank, typename Layout>
struct DepthToViewValueType {
  using type = PatchDataView<T, Rank + 1, Layout>;
};

template <typename T, int Rank, typename Layout>
struct DepthToViewValueType<ScalarDepth, T, Rank, Layout> {
  using type = PatchDataView<T, Rank, Layout>;
};

template <typename State, typename Layout, int Rank> struct ViewBaseImpl {
  template <typename Depth>
  using fn = typename DepthToViewValueType<Depth, double, Rank, Layout>::type;
  using type = boost::mp11::mp_transform<fn, meta::Depths<State>>;
};

template <typename State, typename Layout, int Rank>
struct ViewBaseImpl<const State, Layout, Rank> {
  template <typename Depth>
  using fn =
      typename DepthToViewValueType<Depth, const double, Rank, Layout>::type;
  using type = boost::mp11::mp_transform<fn, meta::Depths<State>>;
};
} // namespace detail

template <typename State, typename Layout, int Rank>
using ViewBase = typename detail::ViewBaseImpl<State, Layout, Rank>::type;

template <typename S, typename L = layout_left, int R = S::Equation::Rank()>
struct BasicView : ViewBase<S, L, R> {
  using State = S;
  using Layout = L;
  static constexpr int Rank = R;
  using ValueType =
      std::conditional_t<std::is_const<S>::value, const double, double>;
  using Equation = typename State::Equation;
  using Depths = typename Equation::ConservativeDepths;
  using Traits = StateTraits<ViewBase<S, L, R>>;

  static constexpr int rank() noexcept { return Rank; }

  BasicView() = default;

  BasicView(const ViewBase<S, L, R>& base) : ViewBase<S, L, R>{base} {}
  BasicView& operator=(const ViewBase<S, L, R>& base) {
    static_cast<ViewBase<S, L, R>&>(*this) = base;
  }
};

template <typename S, typename L, int R>
BasicView(const BasicView<S, L, R>&) -> BasicView<S, L, R>;

template <typename S, typename L, int R>
struct StateTraits<BasicView<S, L, R>> : StateTraits<ViewBase<S, L, R>> {};

template <typename State, typename Layout, int Rank>
BasicView<const State, Layout, Rank>
AsConst(const BasicView<State, Layout, Rank>& v) {
  BasicView<const State, Layout, Rank> w{};
  ForEachVariable([](auto& w, const auto& v) { w = v; }, w, v);
  return w;
}

template <typename Eq>
const Conservative<Eq>& AsCons(const Conservative<Eq>& x) {
  return x;
}

template <typename Eq> Conservative<Eq>& AsCons(Conservative<Eq>& x) {
  return x;
}

template <typename Eq>
const ConservativeBase<Eq>& AsCons(const Complete<Eq>& x) {
  return x;
}

template <typename Eq> ConservativeBase<Eq>& AsCons(Complete<Eq>& x) {
  return x;
}

template <typename State, typename L, int R>
auto AsCons(const BasicView<State, L, R>& view) {
  using T = std::remove_const_t<State>;
  using Equation = typename T::Equation;
  if constexpr (std::is_same<T, Conservative<Equation>>::value) {
    return view;
  } else {
    using Cons = std::conditional_t<std::is_const<State>::value,
                                    const Conservative<Equation>,
                                    Conservative<Equation>>;
    BasicView<Cons, L, R> cons(view);
    return cons;
  }
}

template <typename State, int Rank = State::Equation::Rank()>
using View = BasicView<State, layout_stride, Rank>;

template <typename T> struct IsView : std::false_type {};

template <typename State, typename Layout, int Rank>
struct IsView<BasicView<State, Layout, Rank>> : std::true_type {};

template <int N, typename State, typename Layout, int Rank>
dynamic_extents<static_cast<std::size_t>(Rank)>
Extents(const BasicView<State, Layout, Rank>& view) {
  static_assert(N >= 0, "N has to be greater than 0!");
  return get<static_cast<std::size_t>(N)>(view).Extents();
}

template <int N, typename State, typename Layout, int Rank>
IndexBox<Rank> Box(const BasicView<State, Layout, Rank>& view) {
  static_assert(N >= 0, "N has to be greater than 0!");
  return get<static_cast<std::size_t>(N)>(view).Box();
}

template <int N, typename State, typename Layout, int Rank>
typename Layout::template mapping<
    dynamic_extents<static_cast<std::size_t>(Rank)>>
Mapping(const BasicView<State, Layout, Rank>& view) {
  static_assert(N >= 0, "N has to be greater than 0!");
  return get<static_cast<std::size_t>(N)>(view).Mapping();
}

namespace detail {

template <typename State, typename Layout, int Rank, typename Eq>
struct DepthsImpl<BasicView<State, Layout, Rank>, Eq> {
  constexpr auto operator()(const Eq& eq) const noexcept {
    DepthsImpl<std::remove_const_t<State>, Eq> depths{};
    return depths(eq);
  }
};

template <typename T> struct GetNumberOfComponentsImpl {
  int_constant<1> operator()(double) const noexcept { return {}; }
  int_constant<1> operator()(int) const noexcept { return {}; }

  template <std::size_t N>
  int_constant<static_cast<int>(N)>
  operator()(const std::array<int, N>&) const noexcept {
    return {};
  }

  template <int N, int M, int O, int MR, int MC>
  int_constant<N>
  operator()(const Eigen::Array<double, N, M, O, MR, MC>&) const noexcept {
    return {};
  }

  template <int M, int O, int MR, int MC>
  int operator()(const Eigen::Array<double, Eigen::Dynamic, M, O, MR, MC>& x)
      const noexcept {
    return static_cast<int>(x.rows());
  }

  int operator()(const std::vector<int>& x) const noexcept {
    return static_cast<int>(x.size());
  }
};

template <typename S, typename Layout, int Rank>
struct GetNumberOfComponentsImpl<BasicView<S, Layout, Rank>> {
  int_constant<1>
  operator()(const PatchDataView<double, Rank, Layout>&) const noexcept {
    return {};
  }
  int operator()(
      const PatchDataView<double, Rank + 1, Layout>& pdv) const noexcept {
    return static_cast<int>(pdv.Extent(Rank));
  }
};

template <typename T>
static constexpr GetNumberOfComponentsImpl<T> GetNumberOfComponents{};

template <typename T> struct AtComponentImpl {
  template <typename FP>
  std::enable_if_t<std::is_floating_point<remove_cvref_t<FP>>::value, FP>
  operator()(FP&& x, int) const noexcept {
    return x;
  }

  template <int N>
  double& operator()(Eigen::Array<double, N, 1>& x, int n) const noexcept {
    return x[n];
  }

  template <int N>
  const double& operator()(const Eigen::Array<double, N, 1>& x,
                           int n) const noexcept {
    return x[n];
  }

  template <int N, int M, int O, int MR, int MC>
  auto operator()(Eigen::Array<double, N, M, O, MR, MC>& x,
                  int n) const noexcept {
    return x.row(n);
  }

  template <int N, int M, int O, int MR, int MC>
  auto operator()(const Eigen::Array<double, N, M, O, MR, MC>& x,
                  int n) const noexcept {
    return x.row(n);
  }
};

template <typename S, typename Layout, int Rank>
struct AtComponentImpl<BasicView<S, Layout, Rank>> {
  using ValueType = typename BasicView<S, Layout, Rank>::ValueType;

  const PatchDataView<ValueType, Rank, Layout>&
  operator()(const PatchDataView<ValueType, Rank, Layout>& x,
             int) const noexcept {
    return x;
  }

  PatchDataView<ValueType, Rank, Layout>
  operator()(const PatchDataView<ValueType, Rank + 1, Layout>& x,
             int n) const noexcept {
    constexpr std::size_t sRank = Rank;
    std::array<std::ptrdiff_t, sRank + 1> index{};
    index[sRank] = n;
    std::array<std::ptrdiff_t, sRank> extents;
    for (std::size_t r = 0; r < sRank; ++r) {
      extents[r] = x.Extent(r);
    }
    std::array<std::ptrdiff_t, sRank> strides;
    for (std::size_t r = 0; r < sRank; ++r) {
      strides[r] = x.Stride(r);
    }
    std::array<std::ptrdiff_t, sRank> origin;
    std::copy_n(x.Origin().begin(), Rank, origin.begin());
    if constexpr (std::is_same_v<Layout, layout_stride>) {
      layout_stride::mapping<dynamic_extents<sRank>> mapping{
          dynamic_extents<sRank>(extents), strides};
      mdspan<ValueType, sRank, Layout> mds(&x.MdSpan()(index), mapping);
      return PatchDataView<ValueType, Rank, Layout>(mds, origin);
    } else {
      mdspan<ValueType, sRank, Layout> mds(&x.MdSpan()(index), extents);
      return PatchDataView<ValueType, Rank, Layout>(mds, origin);
    }
  }
};

template <typename State> constexpr static AtComponentImpl<State> AtComponent{};
} // namespace detail

template <typename F, typename... Ts>
void ForEachComponent(F function, Ts&&... states) {
  ForEachVariable(
      [&](auto&&... variables) {
        const auto vs = std::make_tuple(std::addressof(variables)...);
        using T =
            remove_cvref_t<boost::mp11::mp_first<boost::mp11::mp_list<Ts...>>>;
        const int n_components =
            detail::GetNumberOfComponents<T>(*std::get<0>(vs));
        for (int i = 0; i < n_components; ++i) {
          function(detail::AtComponent<remove_cvref_t<Ts>>(variables, i)...);
        }
      },
      std::forward<Ts>(states)...);
}

template <typename State, typename Layout, int Rank>
void Load(State& state, const BasicView<const State, Layout, Rank>& view,
          const std::array<std::ptrdiff_t, State::Equation::Rank()>& index) {
  ForEachComponent(
      [&](double& component,
          const PatchDataView<const double, Rank, Layout>& pdv) {
        FUB_ASSERT(Contains(pdv.Box(), index));
        component = pdv(index);
      },
      state, view);
}

template <typename State, typename Layout, int Rank>
void Load(State& state, const BasicView<State, Layout, Rank>& view,
          const std::array<std::ptrdiff_t, State::Equation::Rank()>& index) {
  ForEachComponent(
      [&](double& component, const auto& pdv) {
        FUB_ASSERT(Contains(pdv.Box(), index));
        component = pdv(index);
      },
      state, view);
}

template <typename Eq, typename Layout>
void Store(const BasicView<Conservative<Eq>, Layout, Eq::Rank()>& view,
           const Conservative<Eq>& state,
           const std::array<std::ptrdiff_t, Eq::Rank()>& index) {
  ForEachComponent(
      [&](auto data, auto block) {
        FUB_ASSERT(Contains(data.Box(), index));
        Store(data, block, index);
      },
      view, state);
}

template <typename Eq, typename Layout>
void Store(const BasicView<Complete<Eq>, Layout, Eq::Rank()>& view,
           const Complete<Eq>& state,
           const std::array<std::ptrdiff_t, Eq::Rank()>& index) {
  ForEachComponent(
      [&](auto data, auto block) {
        FUB_ASSERT(Contains(data.Box(), index));
        Store(data, block, index);
      },
      view, state);
}

template <Direction dir, typename T, typename L, int Rank,
          typename SliceSpecifier>
auto Slice(const BasicView<T, L, Rank>& view, SliceSpecifier slice) {
  View<T, Rank> sliced_view;
  ForEachVariable([&](auto& s, auto v) { s = Slice<dir>(v, slice); },
                  sliced_view, view);
  return sliced_view;
}

template <typename State, int Rank, typename Layout>
View<State, Rank> Shrink(const BasicView<State, Layout, Rank>& view,
                         Direction dir, std::array<std::ptrdiff_t, 2> offsets) {
  IndexBox<Rank> box = Box<0>(view);
  IndexBox<Rank> shrinked_box = Shrink(box, dir, offsets);
  View<State, Rank> shrinked_view;
  using T = std::conditional_t<std::is_const_v<State>, const double, double>;
  ForEachVariable(
      overloaded{
          [&shrinked_box](PatchDataView<T, Rank, layout_stride>& dest,
                          PatchDataView<T, Rank, Layout> src) {
            dest = src.Subview(shrinked_box);
          },
          [&shrinked_box](PatchDataView<T, Rank + 1, layout_stride>& dest,
                          PatchDataView<T, Rank + 1, Layout> src) {
            IndexBox<Rank + 1> src_box = src.Box();
            IndexBox<Rank + 1> embed_shrinked_box = Embed<Rank + 1>(
                shrinked_box, {src_box.lower[Rank], src_box.upper[Rank]});
            dest = src.Subview(embed_shrinked_box);
          }},
      shrinked_view, view);
  return shrinked_view;
}

template <typename Equation> bool AnyNaN(const Complete<Equation>& state) {
  bool any_nan = false;
  ForEachComponent([&any_nan](double x) { any_nan |= std::isnan(x); }, state);
  return any_nan;
}

namespace detail {
template <typename T, typename Depth> struct DepthToPointerType {
  using type = std::pair<T*, std::ptrdiff_t>;
};

template <typename T> struct DepthToPointerType<T, ScalarDepth> {
  using type = T*;
};

template <typename State> struct ViewPointerBaseImpl {
  template <typename Depth>
  using fn = typename DepthToPointerType<double, Depth>::type;
  using type = boost::mp11::mp_transform<fn, meta::Depths<State>>;
};

template <typename State> struct ViewPointerBaseImpl<const State> {
  template <typename Depth>
  using fn = typename DepthToPointerType<const double, Depth>::type;
  using type = boost::mp11::mp_transform<fn, meta::Depths<State>>;
};
} // namespace detail

template <typename State>
using ViewPointerBase = typename detail::ViewPointerBaseImpl<State>::type;

template <typename State> struct ViewPointer : ViewPointerBase<State> {
  using Traits = StateTraits<ViewPointer<State>>;
};

template <typename State>
struct StateTraits<ViewPointer<State>> : StateTraits<ViewPointerBase<State>> {};

template <typename State>
ViewPointer(const ViewPointer<State>&) -> ViewPointer<State>;

template <typename State, typename Layout, int Rank>
ViewPointer<State> Begin(const BasicView<State, Layout, Rank>& view) {
  ViewPointer<State> pointer;
  ForEachVariable(
      overloaded{[&](double*& p, auto pdv) { p = pdv.Span().begin(); },
                 [&](const double*& p, auto pdv) { p = pdv.Span().begin(); },
                 [&](auto& p, auto pdv) {
                   p = std::pair{pdv.Span().begin(), pdv.Stride(Rank)};
                 }},
      pointer, view);
  return pointer;
}

template <typename State, typename Layout, int Rank>
ViewPointer<State> End(const BasicView<State, Layout, Rank>& view) {
  ViewPointer<State> pointer;
  ForEachVariable(
      overloaded{[&](double*& p, auto pdv) { p = pdv.Span().end(); },
                 [&](const double*& p, auto pdv) { p = pdv.Span().end(); },
                 [&](auto& p, auto pdv) {
                   p = std::pair{pdv.Span().end(), pdv.Stride(Rank)};
                 }},
      pointer, view);
  return pointer;
}

template <typename Eq>
void Load(Complete<Eq>& state,
          nodeduce_t<ViewPointer<const Complete<Eq>>> pointer) {
  ForEachVariable(
      overloaded{[](double& x, const double* p) { x = *p; },
                 [](auto& xs, std::pair<const double*, std::ptrdiff_t> ps) {
                   const double* p = ps.first;
                   for (std::size_t i = 0; i < xs.size(); ++i) {
                     xs[i] = *p;
                     p += ps.second;
                   }
                 }},
      state, pointer);
}

template <typename Eq>
void Load(Conservative<Eq>& state,
          nodeduce_t<ViewPointer<const Conservative<Eq>>> pointer) {
  ForEachVariable(
      overloaded{[](double& x, const double* p) { x = *p; },
                 [](auto& xs, std::pair<const double*, std::ptrdiff_t> ps) {
                   const double* p = ps.first;
                   for (std::size_t i = 0; i < xs.size(); ++i) {
                     xs[i] = *p;
                     p += ps.second;
                   }
                 }},
      state, pointer);
}

template <typename State>
void Advance(ViewPointer<State>& pointer, std::ptrdiff_t n) noexcept {
  ForEachVariable(overloaded{[n](double*& ptr) { ptr += n; },
                             [n](const double*& ptr) { ptr += n; },
                             [n](auto& ptr) { ptr.first += n; }},
                  pointer);
}

namespace detail {
template <typename T> struct DepthToIndexMappingTypeImpl;

/// The value-type for a single-cell state for a scalar quantity.
template <> struct DepthToIndexMappingTypeImpl<ScalarDepth> {
  using type = int;
};

/// The value-type for a single-cell state for a scalar quantity.
template <int Depth> struct DepthToIndexMappingTypeImpl<VectorDepth<Depth>> {
  using type = std::array<int, static_cast<std::size_t>(Depth)>;
};

/// The value-type for a single-cell state for a scalar quantity.
template <> struct DepthToIndexMappingTypeImpl<VectorDepth<dynamic_extent>> {
  using type = std::vector<int>;
};

template <typename T>
using DepthToIndexMappingType = typename DepthToIndexMappingTypeImpl<T>::type;
} // namespace detail

template <typename Equation>
using IndexMappingBase =
    boost::mp11::mp_transform<detail::DepthToIndexMappingType,
                              typename Equation::CompleteDepths>;

template <typename Equation> struct IndexMapping : IndexMappingBase<Equation> {
  IndexMapping() = default;
  IndexMapping(const Equation& equation) {
    auto depths = Depths(equation, Type<Complete<Equation>>{});
    int counter = 0;
    ForEachVariable(
        overloaded{
            [&](int& id, ScalarDepth) { id = counter++; },
            [&](auto&& ids, auto depth) {
              if constexpr (std::is_same_v<std::decay_t<decltype(depth)>,
                                           int>) {
                ids.resize(depth);
              }
              std::iota(ids.begin(), ids.end(), counter);
              counter += depth;
            },
        },
        *this, depths);
  }
};

template <typename Eq>
struct StateTraits<IndexMapping<Eq>> : StateTraits<IndexMappingBase<Eq>> {};

template <typename Equation>
void CopyFromBuffer(Complete<Equation>& state, span<const double> buffer) {
  int comp = 0;
  ForEachComponent(
      [&comp, &buffer](auto&& var) {
        var = buffer[comp];
        comp += 1;
      },
      state);
}

template <typename Equation>
void CopyFromBuffer(Conservative<Equation>& state, span<const double> buffer) {
  int comp = 0;
  ForEachComponent(
      [&comp, &buffer](auto&& var) {
        var = buffer[comp];
        comp += 1;
      },
      state);
}

template <typename Equation>
void CopyToBuffer(span<double> buffer, const Conservative<Equation>& state) {
  int comp = 0;
  ForEachComponent(
      [&comp, buffer](auto&& var) {
        buffer[comp] = var;
        comp += 1;
      },
      state);
}

} // namespace fub

#endif
