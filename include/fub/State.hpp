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

#include <Eigen/Dense>

#include <boost/mp11/algorithm.hpp>
#include <boost/mp11/tuple.hpp>

#include <tuple>

namespace fub {

template <typename T> struct StateTraits;

namespace detail {
template <std::size_t I, typename... Ts> constexpr auto ZipNested(Ts&&... ts) {
  return std::make_tuple(std::get<I>(ts)...);
}

template <std::size_t... Is, typename... Ts>
constexpr auto ZipHelper(std::index_sequence<Is...>, Ts&&... ts) {
  return std::make_tuple(ZipNested<Is>(std::forward<Ts>(ts)...)...);
}
} // namespace detail

template <typename... Ts> constexpr auto Zip(Ts&&... ts) {
  using First =
      std::decay_t<boost::mp11::mp_front<boost::mp11::mp_list<Ts...>>>;
  return detail::ZipHelper(
      std::make_index_sequence<std::tuple_size<First>::value>(),
      std::forward<Ts>(ts)...);
}

template <std::size_t I, typename State>
constexpr decltype(auto) get(State&& x) {
  using S = remove_cvref_t<State>;
  constexpr auto& pointers_to_member = S::Traits::pointers_to_member;
  return x.*std::get<I>(pointers_to_member);
}

template <typename StateType, typename F, typename... Ts>
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

template <typename State> auto ToTuple(const State& x) {
  return boost::mp11::tuple_apply(
      [&x](auto... pointer) { return std::make_tuple((x.*pointer)...); },
      StateTraits<State>::pointers_to_member);
}

/// This type is used to tag scalar quantities
struct ScalarDepth : std::integral_constant<int, 1> {};

/// This type is used to tag quantities with a depth known at compile time.
template <int Depth> struct VectorDepth : std::integral_constant<int, Depth> {};

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
  using type = Eigen::Array<double, 1, Width>;
};

/// The value-type for a single-cell state for a scalar quantity.
template <int Depth, int Width>
struct DepthToStateValueTypeImpl<VectorDepth<Depth>, Width> {
  using type = Eigen::Array<double, Depth, Width>;
};

/// For usage with boost::mp11
template <typename T>
using DepthToStateValueType = typename DepthToStateValueTypeImpl<T, 1>::type;
/// @}

/// This type alias transforms state depths into a conservative state associated
/// with a specified equation.
template <typename Equation>
using ConservativeBase =
    boost::mp11::mp_transform<DepthToStateValueType,
                              typename Equation::ConservativeDepths>;

template <typename Eq> struct Conservative;
template <typename Eq> struct Complete;

template <typename Equation>
void InitializeState(const Equation&, const Conservative<Equation>&) {}

template <typename Equation>
void InitializeState(const Equation&, const Complete<Equation>&) {}

/// This type has a constructor which takes an equation and might allocate any
/// dynamically sized member variable.
template <typename Eq> struct Conservative : ConservativeBase<Eq> {
  using Equation = Eq;
  using Depths = typename Equation::ConservativeDepths;

  Conservative() = default;
  Conservative(const ConservativeBase<Eq>& x) : ConservativeBase<Eq>{x} {}
  Conservative& operator=(const ConservativeBase<Eq>& x) {
    static_cast<ConservativeBase<Eq>&>(*this) = x;
  }
  Conservative(const Equation& eq) : ConservativeBase<Eq>{} {
    InitializeState(eq, *this);
  }
};

template <typename Eq>
struct StateTraits<Conservative<Eq>> : StateTraits<ConservativeBase<Eq>> {};

template <typename Equation>
using CompleteBase =
    boost::mp11::mp_transform<DepthToStateValueType,
                              typename Equation::CompleteDepths>;

/// This type has a constructor which takes an equation and might allocate any
/// dynamically sized member variable.
template <typename Eq> struct Complete : CompleteBase<Eq> {
  using Equation = Eq;
  using Depths = typename Equation::CompleteDepths;
  using Traits = StateTraits<CompleteBase<Equation>>;

  Complete() = default;
  Complete(const CompleteBase<Eq>& x) : CompleteBase<Eq>{x} {}
  Complete& operator=(const CompleteBase<Eq>& x) {
    static_cast<CompleteBase<Eq>&>(*this) = x;
  }
  Complete(const Equation& eq) : CompleteBase<Eq>{} {
    InitializeState(eq, *this);
  }
};

template <typename Eq>
struct StateTraits<Complete<Eq>> : StateTraits<CompleteBase<Eq>> {};

///// View type

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
  using type = boost::mp11::mp_transform<fn, typename State::Depths>;
};

template <typename State, typename Layout, int Rank>
struct ViewBaseImpl<const State, Layout, Rank> {
  template <typename Depth>
  using fn =
      typename DepthToViewValueType<Depth, const double, Rank, Layout>::type;
  using type = boost::mp11::mp_transform<fn, typename State::Depths>;
};

template <typename State, typename Layout, int Rank>
using ViewBase = typename ViewBaseImpl<State, Layout, Rank>::type;

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

  BasicView() = default;

  BasicView(const ViewBase<S, L, R>& base) : ViewBase<S, L, R>{base} {}
  BasicView& operator=(const ViewBase<S, L, R>& base) {
    static_cast<ViewBase<S, L, R>&>(*this) = base;
  }
};

template <typename State, typename Layout, int Rank>
BasicView<const State, Layout, Rank>
AsConst(const BasicView<State, Layout, Rank>& v) {
  BasicView<const State, Layout, Rank> w{};
  ForEachVariable<State>([](auto& w, const auto& v) { w = v; }, w, v);
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

template <typename S, typename L, int R>
struct StateTraits<BasicView<S, L, R>> : StateTraits<ViewBase<S, L, R>> {};

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

template <typename T, typename Eq> struct DepthsImpl;

template <typename Eq> struct DepthsImpl<Complete<Eq>, Eq> {
  constexpr typename Eq::CompleteDepths operator()(const Eq&) const noexcept {
    return {};
  }
};

template <typename Eq> struct DepthsImpl<Conservative<Eq>, Eq> {
  constexpr typename Eq::ConservativeDepths operator()(const Eq&) const
      noexcept {
    return {};
  }
};

template <typename State, typename Layout, int Rank, typename Eq>
struct DepthsImpl<BasicView<State, Layout, Rank>, Eq> {
  constexpr auto operator()(const Eq& eq) const noexcept {
    DepthsImpl<std::remove_const_t<State>, Eq> depths{};
    return depths(eq);
  }
};

template <typename T, typename Eq> auto Depths(const Eq& eq) {
  DepthsImpl<T, Eq> depths;
  return depths(eq);
}

template <typename T> struct GetNumberOfComponentsImpl {
  int_constant<1> operator()(double) const noexcept { return {}; }

  template <int N, int M>
  int_constant<N> operator()(const Eigen::Array<double, N, M>&) const noexcept {
    return {};
  }

  template <int M>
  int operator()(const Eigen::Array<double, Eigen::Dynamic, M>& x) const
      noexcept {
    return static_cast<int>(x.rows());
  }
};

template <typename S, typename Layout, int Rank>
struct GetNumberOfComponentsImpl<BasicView<S, Layout, Rank>> {
  int_constant<1> operator()(const PatchDataView<double, Rank, Layout>&) const
      noexcept {
    return {};
  }
  int operator()(const PatchDataView<double, Rank + 1, Layout>& pdv) const
      noexcept {
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
  const double& operator()(const Eigen::Array<double, N, 1>& x, int n) const
      noexcept {
    return x[n];
  }

  template <int N, int M>
  auto operator()(Eigen::Array<double, N, M>& x, int n) const noexcept {
    return x.row(n);
  }

  template <int N, int M>
  auto operator()(const Eigen::Array<double, N, M>& x, int n) const noexcept {
    return x.row(n);
  }
};

template <typename S, typename Layout, int Rank>
struct AtComponentImpl<BasicView<S, Layout, Rank>> {
  using ValueType = typename BasicView<S, Layout, Rank>::ValueType;

  const PatchDataView<ValueType, Rank, Layout>&
  operator()(const PatchDataView<ValueType, Rank, Layout>& x, int) const
      noexcept {
    return x;
  }

  PatchDataView<ValueType, Rank, Layout>
  operator()(const PatchDataView<ValueType, Rank + 1, Layout>& x, int n) const
      noexcept {
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

template <typename StateType, typename F, typename... Ts>
void ForEachComponent(F function, Ts&&... states) {
  ForEachVariable<StateType>(
      [&](auto&&... variables) {
        const auto vs = std::make_tuple(std::addressof(variables)...);
        using T =
            remove_cvref_t<boost::mp11::mp_first<boost::mp11::mp_list<Ts...>>>;
        const int n_components = GetNumberOfComponents<T>(*std::get<0>(vs));
        for (int i = 0; i < n_components; ++i) {
          function(AtComponent<remove_cvref_t<Ts>>(variables, i)...);
        }
      },
      std::forward<Ts>(states)...);
}

template <typename State, typename Layout, int Rank>
void Load(State& state, const BasicView<const State, Layout, Rank>& view,
          const std::array<std::ptrdiff_t, State::Equation::Rank()>& index) {
  ForEachComponent<State>(
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
  ForEachComponent<State>(
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
  ForEachComponent<Conservative<Eq>>(
      [&](auto data, auto block) {
        FUB_ASSERT(Contains(data.Box(), index));
        Store(data, block, index);
      }, view, state);
}

template <typename Eq, typename Layout>
void Store(const BasicView<Complete<Eq>, Layout, Eq::Rank()>& view,
           const Complete<Eq>& state,
           const std::array<std::ptrdiff_t, Eq::Rank()>& index) {
  ForEachComponent<Complete<Eq>>(
      [&](auto data, auto block) {
        FUB_ASSERT(Contains(data.Box(), index));
        Store(data, block, index);
  }, view, state);
}

template <typename State, typename Layout, int Rank,
          typename... SliceSpecifiers>
View<State, SliceRank_<SliceSpecifiers...>>
Subspan(const BasicView<State, Layout, Rank>& view, SliceSpecifiers... slices) {
  static constexpr int R = SliceRank_<SliceSpecifiers...>;
  View<State, R> subview{};
  ForEachVariable([&](auto& subcomp,
                      const auto& comp) { subcomp = subspan(comp, slices...); },
                  subview, view);
  return subview;
}

template <typename F, typename View, typename... Views>
void ForEachRow(F function, const View& view, const Views&... views) {
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

template <Direction dir, typename T, typename E, typename L, typename A,
          typename SliceSpecifier>
auto Slice(const basic_mdspan<T, E, L, A>& span, SliceSpecifier slice) {
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

template <Direction dir, typename T, typename L, int Rank,
          typename SliceSpecifier>
View<T, Rank> Slice(const BasicView<T, L, Rank>& view, SliceSpecifier slice) {
  View<T, Rank> sliced_view;
  ForEachVariable<remove_cvref_t<T>>(
      [&](auto& s, auto v) { s = Slice<dir>(v, slice); }, sliced_view, view);
  return sliced_view;
}

template <typename Equation> bool AnyNaN(const Complete<Equation>& state) {
  bool any_nan = false;
  ForEachComponent<Complete<Equation>>(
      [&any_nan](double x) { any_nan |= std::isnan(x); }, state);
  return any_nan;
}

} // namespace fub

#endif
