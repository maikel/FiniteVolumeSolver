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

#ifndef FUB_STATE_FACADE_HPP
#define FUB_STATE_FACADE_HPP

#include "fub/Direction.hpp"
#include "fub/ext/Eigen.hpp"

#include <boost/hana.hpp>
#include <boost/hana/core/to.hpp>
#include <boost/hana/define_struct.hpp>
#include <boost/hana/for_each.hpp>
#include <boost/hana/map.hpp>
#include <boost/hana/members.hpp>
#include <boost/hana/transform.hpp>
#include <boost/hana/zip.hpp>

namespace fub {
/// This is a tag which group states of the same state type together.
///
/// For example Equation::State and StateView<Equation::State, Equation> share
/// the same state tag because their member variables differ in types only.
template <template <typename...> class> struct StateTag {};

// Forward Decleration.
template <typename> struct StateFacade;

/// @{
/// This type trait returns true if the specified type T is a StateFacade.
template <typename T> struct IsStateFacade : std::false_type {};
template <typename T> struct IsStateFacade<StateFacade<T>> : std::true_type {};
/// @}

template <typename Eq> auto GetSizes(const Eq&);

inline auto Ref(double& x) { return std::ref(x); }

inline auto Ref(const double& x) { return x; }

template <typename T, std::size_t N> auto Ref(const std::array<T, N>& x) {
  return Eigen::Map<const Eigen::Array<T, N, 1>>(x.data());
}

template <typename T, std::size_t N> auto Ref(std::array<T, N>& x) {
  return Eigen::Map<Eigen::Array<T, N, 1>>(x.data());
}

template <typename T> auto Ref(const std::vector<T>& x) {
  return Eigen::Map<const Eigen::Array<T, Eigen::Dynamic, 1>>(x.data(),
                                                              x.size());
}

template <typename T> auto Ref(std::vector<T>& x) {
  return Eigen::Map<Eigen::Array<T, Eigen::Dynamic, 1>>(x.data(), x.size());
}

template <typename T, typename E, typename L, typename A>
auto Ref(basic_mdspan<T, E, L, A> mdspan) {
  return mdspan;
}

template <typename T> auto Ref(const T& x) { return T::Map(x.data()); }

template <typename T> auto Ref(T& x) { return T::Map(x.data()); }

template <typename Derived, typename Base, bool IsMdSpan>
struct StateFacade_MdSpan : Base {
  using Equation = typename Base::Equation;
  StateFacade_MdSpan() = default;
  constexpr StateFacade_MdSpan(const Base& base) : Base{base} {}
  constexpr StateFacade_MdSpan(Base&& base) : Base{std::move(base)} {}
  constexpr StateFacade_MdSpan(const Equation&) {}
};

template <typename Derived, typename Base>
struct StateFacade_MdSpan<Derived, Base, true> : Base {
  StateFacade_MdSpan() = default;
  constexpr StateFacade_MdSpan(const Base& base) : Base{base} {}

  template <typename... Is> auto operator()(Is... is) {
    return Transform([&](auto&& member) { return Ref(member(is...)); },
                     AsDerived());
  }

  template <typename... Is> auto operator()(Is... is) const {
    return Transform([&](auto&& member) { return Ref(member(is...)); },
                     AsDerived());
  }

  Derived& AsDerived() { return static_cast<Derived&>(*this); }
  const Derived& AsDerived() const {
    return static_cast<const Derived&>(*this);
  }
};

template <typename T> struct StateFacade;

template <typename T>
constexpr auto IsEigenType(boost::hana::basic_type<T> type) {
  return type == boost::hana::type_c<Scalar> ||
         type == boost::hana::type_c<ArrayXd> ||
         type == boost::hana::type_c<Array2d> ||
         type == boost::hana::type_c<Array3d>;
}

template <template <typename...> typename Data, typename... Ts>
constexpr auto
IsEigenFacade(boost::hana::basic_type<StateFacade<Data<Ts...>>>) {
  return (IsEigenType(boost::hana::type_c<Ts>) && ...);
}

template <typename T> struct IsMdSpan : std::false_type {};

template <template <typename...> typename C, typename... Ts>
struct IsMdSpan<C<Ts...>> : std::conjunction<is_mdspan<Ts>...> {};

template <typename Derived, typename Base, bool IsEigenType>
struct StateFacade_Eigen
    : StateFacade_MdSpan<Derived, Base, IsMdSpan<Base>::value> {
  using Base_ = StateFacade_MdSpan<Derived, Base, IsMdSpan<Base>::value>;
  using Base_::Base_;
};

template <typename Derived, typename Base>
struct StateFacade_Eigen<Derived, Base, true> : Base {
  using Equation = typename Base::Equation;

  StateFacade_Eigen() = default;
  constexpr StateFacade_Eigen(const Base& base) : Base{base} {}

  constexpr StateFacade_Eigen(const Equation&) {
    //   const auto sizes = GetSizes<StateFacade<Data<MemberTypes...>>>(eq);
    //   constexpr auto mem_accessors = accessors();
    //   constexpr auto size_accessors =
    //   std::decay_t<decltype(sizes)>::accessors(); constexpr auto zipped =
    //   boost::hana::zip(mem_accessors, size_accessors);
    //   boost::hana::for_each(zipped, [&](auto accessor) {
    //     using namespace boost::hana::literals;
    //     boost::hana::second(accessor[0_c])(*this).resize(
    //         kChunkSize, boost::hana::second(accessor[1_c])(sizes));
    //   });
    // }
  }
};

template <template <typename...> class Data, typename... MemberTypes>
struct StateFacade<Data<MemberTypes...>>
    : StateFacade_Eigen<StateFacade<Data<MemberTypes...>>, Data<MemberTypes...>,
                        IsEigenFacade(boost::hana::type_c<
                                      StateFacade<Data<MemberTypes...>>>)> {
  using hana_tag = StateTag<Data>;
  using Base = StateFacade_Eigen<
      StateFacade<Data<MemberTypes...>>, Data<MemberTypes...>,
      IsEigenFacade(boost::hana::type_c<StateFacade<Data<MemberTypes...>>>)>;
  using Equation = typename Data<MemberTypes...>::Equation;

  StateFacade() = default;
  StateFacade(const StateFacade&) = default;
  StateFacade& operator=(const StateFacade&) = default;
  StateFacade(StateFacade&&) = default;
  StateFacade& operator=(StateFacade&&) = default;
  ~StateFacade() = default;

  using Base::Base;

  template <
      typename... Ys,
      typename = std::enable_if_t<sizeof...(Ys) == sizeof...(MemberTypes)>,
      typename =
          std::enable_if_t<(std::is_convertible<Ys, MemberTypes>() && ...)>>
  constexpr StateFacade(Ys&&... ys)
      : Base(Data<MemberTypes...>{MemberTypes(std::forward<Ys>(ys))...}) {}

  template <typename... Ys>
  constexpr StateFacade(const StateFacade<Data<Ys...>>& other)
      : Base{boost::hana::unpack(
            boost::hana::accessors<Data<Ys...>>(), [&](auto&&... ys) {
              return Data<MemberTypes...>{(boost::hana::second(ys)(other))...};
            })} {}

  template <typename... Ys>
  // typename = std::enable_if_t<IsEigenType(boost::hana::type_c<Ys>...)>>
  constexpr StateFacade& operator=(const StateFacade<Data<Ys...>>& other) {
    constexpr auto this_accessors =
        boost::hana::accessors<Data<MemberTypes...>>();
    constexpr auto other_accessors = boost::hana::accessors<Data<Ys...>>();
    constexpr auto accessors =
        boost::hana::zip(this_accessors, other_accessors);
    boost::hana::for_each(accessors, [&](auto&& as) {
      using namespace boost::hana::literals;
      auto& this_member = boost::hana::second(as[0_c])(*this);
      using Member = std::decay_t<decltype(this_member)>;
      if constexpr (std::is_same_v<Member, ArrayXd>) {
        auto other_member = boost::hana::second(as[1_c])(other);
        FUB_ASSERT(this_member.cols() == other_member.cols());
        auto ov = other_member.reshaped();
        std::copy(ov.begin(), ov.end(), this_member.data());
      } else if constexpr (std::is_same_v<Member,
                                          std::reference_wrapper<double>>) {
        this_member.get() = boost::hana::second(as[1_c])(other);
      } else {
        this_member = boost::hana::second(as[1_c])(other);
      }
    });
    return *this;
  }

  static constexpr auto accessors() noexcept {
    return boost::hana::accessors<Data<MemberTypes...>>();
  }

  static constexpr auto size() noexcept {
    return boost::hana::length(accessors());
  }

  static constexpr auto names() noexcept {
    return boost::hana::keys(
        boost::hana::to<boost::hana::map_tag>(accessors()));
  }

  static constexpr auto value_types() noexcept {
    constexpr auto map = boost::hana::to<boost::hana::map_tag>(accessors());
    constexpr auto accessors = boost::hana::values(map);
    constexpr auto types = boost::hana::transform(accessors, [](auto accessor) {
      return boost::hana::type<remove_cvref_t<decltype(
          accessor(std::declval<Data<MemberTypes...>&>()))>>{};
    });
    return types;
  }

  static constexpr auto as_map() noexcept {
    constexpr auto as = accessors();
    return boost::hana::unpack(
        as, [](auto... a) { return boost::hana::make_map(a...); });
  }

  constexpr auto as_tuple() const noexcept {
    return boost::hana::transform(
        boost::hana::values(as_map()),
        [&](auto accessor) { return accessor(*this); });
  }
};

template <typename Fn, typename T>
constexpr auto Transform(Fn function, T&& state) {
  constexpr auto accessors = std::decay_t<T>::accessors();
  return boost::hana::unpack(accessors, [&](auto... accessors) {
    return boost::hana::make<typename std::decay_t<T>::hana_tag>(
        function(boost::hana::second(accessors)(state))...);
  });
}

template <typename Fn, typename T, typename S>
constexpr auto Transform(Fn function, T&& state1, S&& state2) {
  constexpr auto accessors1 = std::decay_t<T>::accessors();
  constexpr auto accessors2 = std::decay_t<S>::accessors();
  constexpr auto accessors = boost::hana::zip(accessors1, accessors2);
  return boost::hana::unpack(accessors, [&](auto... as) {
    using namespace boost::hana::literals;
    return boost::hana::make<typename std::decay_t<T>::hana_tag>(
        (function(boost::hana::second(as[0_c])(state1),
                  boost::hana::second(as[1_c])(state2)))...);
  });
}

template <typename T, typename S>
auto operator+(const StateFacade<T>& xs, const StateFacade<S>& ys) {
  return Transform([](auto&& x, auto&& y) { return x + y; }, xs, ys);
}

template <typename T, typename S>
auto operator-(const StateFacade<T>& xs, const StateFacade<S>& ys) {
  return Transform([](auto&& x, auto&& y) { return x - y; }, xs, ys);
}

template <typename Factor, typename T>
auto operator*(const Factor& x, const StateFacade<T>& ys) {
  return Transform([x](auto&& y) { return x * y; }, ys);
}

template <typename T, typename Factor>
auto operator*(const StateFacade<T>& ys, const Factor& x) {
  return Transform([x](auto&& y) { return y * x; }, ys);
}

template <typename T, typename S>
auto operator*(const StateFacade<T>& xs, const StateFacade<S>& ys) {
  return Transform([](auto&& x, auto&& y) { return x * y; }, xs, ys);
}

template <typename T, typename Factor>
auto operator/(const StateFacade<T>& ys, const Factor& x) {
  return Transform([x](auto&& y) { return y / x; }, ys);
}

template <typename T, typename S>
auto operator/(const StateFacade<T>& xs, const StateFacade<S>& ys) {
  return Transform([](auto&& x, auto&& y) { return x / y; }, xs, ys);
}

template <typename Chunk, typename MdSpans, typename... Indices>
void Load(StateFacade<Chunk>& dest, const StateFacade<MdSpans>& source,
          Indices... indices) {
  constexpr auto chunks = boost::hana::accessors<Chunk>();
  constexpr auto mdspans = boost::hana::accessors<MdSpans>();
  constexpr auto both = boost::hana::zip(chunks, mdspans);
  boost::hana::for_each(both, [&](auto zipped) {
    using namespace boost::hana::literals;
    auto chunks_acc = boost::hana::second(zipped[0_c]);
    auto mdspans_acc = boost::hana::second(zipped[1_c]);
    auto& chunk = chunks_acc(dest);
    auto& mdspan = mdspans_acc(source);
    if constexpr (std::decay_t<decltype(mdspan)>::rank() == 2) {
      FUB_ASSERT(mdspan.rank() == 2);
      FUB_ASSERT(chunk.cols() == mdspan.extent(1));
      const int cols = chunk.cols();
      for (int c = 0; c < cols; ++c) {
        chunk.col(c) = Load(subspan(mdspan, all, c), indices...);
      }
    } else {
      chunk = Load(mdspan, indices...);
    }
  });
}

template <std::ptrdiff_t N, typename Chunk, typename MdSpans>
void Load(span<StateFacade<Chunk>, N> dest, const StateFacade<MdSpans>& source,
          std::ptrdiff_t index) {
  for (StateFacade<Chunk>& state : dest) {
    Load(state, source, index++);
  }
}

template <typename T> auto NoAlias(T&& state) {
  return Transform([](auto&& x) { return x.noalias(); }, state);
}

template <typename Chunk, typename MdSpans, typename... Indices>
void LoadN(StateFacade<Chunk>& dest, int limited,
           const StateFacade<MdSpans>& source, Indices... indices) {
  constexpr auto chunks = boost::hana::accessors<Chunk>();
  constexpr auto mdspans = boost::hana::accessors<MdSpans>();
  constexpr auto both = boost::hana::zip(chunks, mdspans);
  boost::hana::for_each(both, [&](auto zipped) {
    using namespace boost::hana::literals;
    auto chunks_acc = boost::hana::second(zipped[0_c]);
    auto mdspans_acc = boost::hana::second(zipped[1_c]);
    auto& chunk = chunks_acc(dest);
    auto& mdspan = mdspans_acc(source);
    if constexpr (std::decay_t<decltype(mdspan)>::rank() == 2) {
      FUB_ASSERT(mdspan.rank() == 2);
      FUB_ASSERT(chunk.cols() == mdspan.extent(1));
      const int cols = chunk.cols();
      for (int c = 0; c < cols; ++c) {
        chunk.col(c) = LoadN(subspan(mdspan, all, c), limited, indices...);
      }
    } else {
      chunk = LoadN(mdspan, limited, indices...);
    }
  });
}

template <typename T> auto Extents(const StateFacade<T>& state) {
  constexpr auto members = boost::hana::accessors<T>();
  using namespace boost::hana::literals;
  return boost::hana::second(members[0_c])(state).extents();
}

template <typename T> auto Mapping(const StateFacade<T>& state) {
  constexpr auto members = boost::hana::accessors<T>();
  using namespace boost::hana::literals;
  return boost::hana::second(members[0_c])(state).mapping();
}

template <typename MdSpans, typename Chunk, typename... Indices>
void Store(StateFacade<MdSpans>& dest, const StateFacade<Chunk>& source,
           Indices... indices) {
  constexpr auto chunks = boost::hana::accessors<Chunk>();
  constexpr auto mdspans = boost::hana::accessors<MdSpans>();
  constexpr auto both = boost::hana::zip(chunks, mdspans);
  boost::hana::for_each(both, [&](auto zipped) {
    using namespace boost::hana::literals;
    auto chunks_acc = boost::hana::second(zipped[0_c]);
    auto mdspans_acc = boost::hana::second(zipped[1_c]);
    auto& chunk = chunks_acc(source);
    auto& mdspan = mdspans_acc(dest);
    if constexpr (std::decay_t<decltype(mdspan)>::rank() == 2) {
      FUB_ASSERT(chunk.cols() == mdspan.extent(1));
      const int cols = chunk.cols();
      for (int c = 0; c < cols; ++c) {
        Store(subspan(mdspan, all, c), chunk.col(c), indices...);
      }
    } else {
      Store(mdspan, chunk, indices...);
    }
  });
}

template <typename MdSpans, typename Chunk, typename... Indices>
void StoreN(StateFacade<MdSpans>& dest, int limited,
            const StateFacade<Chunk>& source, Indices... indices) {
  constexpr auto chunks = boost::hana::accessors<Chunk>();
  constexpr auto mdspans = boost::hana::accessors<MdSpans>();
  constexpr auto both = boost::hana::zip(chunks, mdspans);
  boost::hana::for_each(both, [&](auto zipped) {
    using namespace boost::hana::literals;
    auto chunks_acc = boost::hana::second(zipped[0_c]);
    auto mdspans_acc = boost::hana::second(zipped[1_c]);
    auto& chunk = chunks_acc(source);
    auto& mdspan = mdspans_acc(dest);
    if constexpr (std::decay_t<decltype(mdspan)>::rank() == 2) {
      FUB_ASSERT(chunk.cols() == mdspan.extent(1));
      const int cols = chunk.cols();
      for (int c = 0; c < cols; ++c) {
        StoreN(subspan(mdspan, all, c), limited, chunk.col(c), indices...);
      }
    } else {
      StoreN(mdspan, limited, chunk, indices...);
    }
  });
}

template <typename Extents, typename F, typename... Datas>
F ForEachRow(std::integral_constant<int, 0>, Extents extents, F function,
             const StateFacade<Datas>&... states) {
  if constexpr (Extents::rank() == 3) {
    constexpr auto view_subspan = [](const auto& state, std::ptrdiff_t y,
                                     std::ptrdiff_t z) {
      return Transform(
          [=](auto mdspan) {
            if constexpr (std::decay_t<decltype(mdspan)>::rank() == 3) {
              return subspan(mdspan, all, y, z);
            } else if constexpr (std::decay_t<decltype(mdspan)>::rank() == 4) {
              return subspan(mdspan, all, y, z, all);
            }
          },
          state);
    };
    for (std::ptrdiff_t z = 0; z < extents.extent(2); ++z) {
      for (std::ptrdiff_t y = 0; y < extents.extent(1); ++y) {
        function(view_subspan(states, y, z)...);
      }
    }
  } else if constexpr (Extents::rank() == 2) {
    constexpr auto view_subspan = [](const auto& state, std::ptrdiff_t y) {
      return Transform(
          [=](auto mdspan) {
            if constexpr (std::decay_t<decltype(mdspan)>::rank() == 2) {
              return subspan(mdspan, all, y);
            } else if constexpr (std::decay_t<decltype(mdspan)>::rank() == 3) {
              return subspan(mdspan, all, y, all);
            }
          },
          state);
    };
    for (std::ptrdiff_t y = 0; y < extents.extent(1); ++y) {
      function(view_subspan(states, y)...);
    }
  } else {
    function(states...);
  }
  return function;
}

template <typename Extents, typename F, typename... Datas>
F ForEachRow(std::integral_constant<int, 1>, Extents extents, F function,
             const StateFacade<Datas>&... states) {
  if constexpr (Extents::rank() == 3) {
    constexpr auto view_subspan = [](const auto& state, std::ptrdiff_t y,
                                     std::ptrdiff_t z) {
      return Transform(
          [=](auto mdspan) {
            if constexpr (std::decay_t<decltype(mdspan)>::rank() == 3) {
              return subspan(mdspan, y, all, z);
            } else if constexpr (std::decay_t<decltype(mdspan)>::rank() == 4) {
              return subspan(mdspan, y, all, z, all);
            }
          },
          state);
    };
    for (std::ptrdiff_t z = 0; z < extents.extent(2); ++z) {
      for (std::ptrdiff_t y = 0; y < extents.extent(1); ++y) {
        function(view_subspan(states, y, z)...);
      }
    }
  } else if constexpr (Extents::rank() == 2) {
    constexpr auto view_subspan = [](const auto& state, std::ptrdiff_t y) {
      return Transform(
          [=](auto mdspan) {
            if constexpr (std::decay_t<decltype(mdspan)>::rank() == 2) {
              return subspan(mdspan, y, all);
            } else if constexpr (std::decay_t<decltype(mdspan)>::rank() == 3) {
              return subspan(mdspan, y, all, all);
            }
          },
          state);
    };
    for (std::ptrdiff_t y = 0; y < extents.extent(1); ++y) {
      function(view_subspan(states, y)...);
    }
  }
  return function;
}

template <typename Extents, typename F, typename... Datas>
F ForEachRow(std::integral_constant<int, 2>, Extents extents, F function,
             const StateFacade<Datas>&... states) {
  if constexpr (Extents::rank() == 3) {
    constexpr auto view_subspan = [](const auto& state, std::ptrdiff_t y,
                                     std::ptrdiff_t z) {
      return Transform(
          [=](auto mdspan) {
            if constexpr (std::decay_t<decltype(mdspan)>::rank() == 3) {
              return subspan(mdspan, y, z, all);
            } else if constexpr (std::decay_t<decltype(mdspan)>::rank() == 4) {
              return subspan(mdspan, y, z, all, all);
            }
          },
          state);
    };
    for (std::ptrdiff_t z = 0; z < extents.extent(2); ++z) {
      for (std::ptrdiff_t y = 0; y < extents.extent(1); ++y) {
        function(view_subspan(states, y, z)...);
      }
    }
  }
  return function;
}

template <typename Extents, typename F, typename... Datas>
F ForEachRow(Extents extents, F function, const StateFacade<Datas>&... states) {
  return ForEachRow(std::integral_constant<int, 0>(), extents, function,
                    states...);
}

template <typename Extents, typename F, typename... Datas>
F ForEachRow(Direction dir, Extents extents, F function,
             const StateFacade<Datas>&... states) {

  if (dir == Direction::X) {
    return ForEachRow(std::integral_constant<int, 0>(), extents, function,
                      states...);
  } else if (dir == Direction::Y) {
    return ForEachRow(std::integral_constant<int, 1>(), extents, function,
                      states...);
  }
  FUB_ASSERT(dir == Direction::Z);
  return ForEachRow(std::integral_constant<int, 2>(), extents, function,
                    states...);
}

struct get_accessor_t {
  template <typename T>
  constexpr auto operator()(const StateFacade<T>&) const noexcept {
    return boost::hana::accessors<T>();
  }
};
inline constexpr get_accessor_t get_accessor;

template <typename F, typename... Ts>
F ForEachMember(F function, Ts&&... states) {
  auto broadcast_state = [](auto accessors, auto state) {
    return boost::hana::transform(accessors, [state](auto a) {
      return boost::hana::make_tuple(a, state);
    });
  };
  auto accessors =
      boost::hana::zip(broadcast_state(get_accessor(states), &states)...);
  boost::hana::for_each(accessors, [&function](auto as) {
    boost::hana::unpack(as, [&function](auto... a) {
      using namespace boost::hana::literals;
      function((boost::hana::second(a[0_c])(*a[1_c]))...);
    });
  });
  return function;
}

template <typename F, typename... Ts>
void ForEachVariable(F function, Ts&&... states) {
  ForEachMember([&](auto&&... members) { ForEachCol(function, members...); },
                states...);
}

} // namespace fub

namespace boost {
namespace hana {

template <template <typename...> class State>
struct make_impl<::fub::StateTag<State>> {
  template <typename... Args> static constexpr auto apply(Args&&... args) {
    return ::fub::StateFacade<State<::fub::remove_cvref_t<Args>...>>(
        std::forward<Args>(args)...);
  }
};

} // namespace hana
} // namespace boost

#endif