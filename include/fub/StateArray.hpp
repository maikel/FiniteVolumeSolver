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

#ifndef FUB_STATE_ARRAY_HPP
#define FUB_STATE_ARRAY_HPP

#include "fub/State.hpp"
#include "fub/ext/Eigen.hpp"
#include <boost/mp11/algorithm.hpp>

namespace fub {
namespace detail {
template <typename Depths, int Width> struct ArrayStateBaseImpl {
  template <typename T>
  using fn = typename DepthToStateValueTypeImpl<T, Width>::type;
  using type = boost::mp11::mp_transform<fn, Depths>;
};
} // namespace detail

template <typename Depths, int Width>
using ArrayStateBase = typename detail::ArrayStateBaseImpl<Depths, Width>::type;

template <typename Depths, int Width>
struct ArrayState : ArrayStateBase<Depths, Width> {
  using Traits = StateTraits<ArrayStateBase<Depths, Width>>;
  static constexpr int Rank = meta::Rank<Depths>::value;
  using Equation = typename Traits::template Equation<Rank>;
  using Base = ArrayStateBase<Depths, Width>;

  using Base::Base;

  ArrayState(const Equation& eq) : Base{} {
    auto depths = ::fub::Depths<ArrayState, Equation>(eq);
    ForEachVariable(
        overloaded{
            [&](Array1d& id, ScalarDepth) { id = Array1d::Zero(); },
            [&](auto&& ids, auto depth) {
              if constexpr (std::is_same_v<std::decay_t<decltype(depth)>,
                                           int>) {
                ids = ArrayXd::Zero(depth, kDefaultChunkSize);
              } else {
                ids = Array<double, decltype(depth)::value>::Zero();
              }
            },
        },
        *this, depths);
  }

  ArrayState& operator+=(const Base& other) {
    ForEachVariable([](auto&& that, auto&& other) { that += other; }, *this,
                    other);
    return *this;
  }

  ArrayState& operator-=(const Base& other) {
    ForEachVariable([](auto&& that, auto&& other) { that -= other; }, *this,
                    other);
    return *this;
  }

  ArrayState& operator*=(double lambda) {
    ForEachVariable([lambda](auto&& that) { that *= lambda; },
                    *this);
    return *this;
  }
};

template <typename Depths, int Width>
struct StateTraits<ArrayState<Depths, Width>>
    : StateTraits<ArrayStateBase<Depths, Width>> {};

template <typename Eq, int Width = kDefaultChunkSize>
struct PrimitiveArray : ArrayState<typename Eq::PrimitiveDepths, Width> {
  using Base = ArrayState<typename Eq::PrimitiveDepths, Width>;

  using Base::Base;
};

template <typename Eq, int Width>
struct StateTraits<PrimitiveArray<Eq, Width>>
    : StateTraits<ArrayState<typename Eq::PrimitiveDepths, Width>> {};

template <typename Eq, int Width = kDefaultChunkSize>
struct CharacteristicsArray
    : ArrayState<typename Eq::CharacteristicsDepths, Width> {
  using Base = ArrayState<typename Eq::CharacteristicsDepths, Width>;

  using Base::Base;
};

template <typename Eq, int Width>
struct StateTraits<CharacteristicsArray<Eq, Width>>
    : StateTraits<ArrayState<typename Eq::CharacteristicsDepths, Width>> {};

namespace detail {
template <typename Eq, int Width> struct ConservativeArrayBaseImpl {
  template <typename T>
  using fn = typename DepthToStateValueTypeImpl<T, Width>::type;
  using type = boost::mp11::mp_transform<fn, typename Eq::ConservativeDepths>;
};
} // namespace detail

template <typename Eq, int Width = kDefaultChunkSize>
using ConservativeArrayBase =
    typename detail::ConservativeArrayBaseImpl<Eq, Width>::type;

template <typename Eq, int Width = kDefaultChunkSize>
struct ConservativeArray : ConservativeArrayBase<Eq, Width> {
  using Equation = Eq;
  // using Depths = typename Equation::ConservativeDepths;
  using Traits = StateTraits<ConservativeArrayBase<Eq, Width>>;

  static constexpr int Size() { return Width; }

  ConservativeArray() = default;
  ConservativeArray(const Equation& eq) {
    auto depths = Depths<Conservative<Equation>>(eq);
    ForEachVariable(
        overloaded{
            [&](Array1d& id, ScalarDepth) { id = Array1d::Zero(); },
            [&](auto&& ids, auto depth) {
              if constexpr (std::is_same_v<std::decay_t<decltype(depth)>,
                                           int>) {
                ids = ArrayXd::Zero(depth, kDefaultChunkSize);
              } else {
                ids = Array<double, decltype(depth)::value>::Zero();
              }
            },
        },
        *this, depths);
  }
};

template <typename Eq, int Width>
struct StateTraits<ConservativeArray<Eq, Width>>
    : StateTraits<ConservativeArrayBase<Eq, Width>> {};

namespace detail {
template <typename Eq, int Width> struct CompleteArrayBaseImpl {
  template <typename T>
  using fn = typename DepthToStateValueTypeImpl<T, Width>::type;
  using type = boost::mp11::mp_transform<fn, typename Eq::CompleteDepths>;
};
} // namespace detail

template <typename Eq, int Width>
using CompleteArrayBase =
    typename detail::CompleteArrayBaseImpl<Eq, Width>::type;

template <typename Eq, int Width = kDefaultChunkSize>
struct CompleteArray : CompleteArrayBase<Eq, Width> {
  using Equation = Eq;
  // using Depths = typename Equation::CompleteDepths;
  using Traits = StateTraits<CompleteArrayBase<Eq, Width>>;

  static constexpr int Size() { return Width; }

  CompleteArray() = default;
  CompleteArray(const Equation& eq) {
    auto depths = Depths<Complete<Equation>>(eq);
    ForEachVariable(
        overloaded{
            [&](Array1d& id, ScalarDepth) { id = Array1d::Zero(); },
            [&](auto&& ids, auto depth) {
              if constexpr (std::is_same_v<std::decay_t<decltype(depth)>,
                                           int>) {
                ids = ArrayXd::Zero(depth, kDefaultChunkSize);
              } else {
                ids = Array<double, decltype(depth)::value>::Zero();
              }
            },
        },
        *this, depths);
  }
};

template <typename Eq, int Width>
struct StateTraits<CompleteArray<Eq, Width>>
    : StateTraits<CompleteArrayBase<Eq, Width>> {};

template <typename Eq, int N>
const ConservativeArray<Eq, N>& AsCons(const ConservativeArray<Eq, N>& x) {
  return x;
}

template <typename Eq, int N>
ConservativeArray<Eq, N>& AsCons(ConservativeArray<Eq, N>& x) {
  return x;
}

template <typename Eq, int N>
const ConservativeArrayBase<Eq, N>& AsCons(const CompleteArray<Eq, N>& x) {
  return x;
}

template <typename Eq, int N>
ConservativeArrayBase<Eq, N>& AsCons(CompleteArray<Eq, N>& x) {
  return x;
}

template <typename Eq, int N>
void Load(Conservative<Eq>& q, const ConservativeArray<Eq, N>& qs, int i) {
  ForEachComponent([&](auto& qi, auto qsi) { qi = qsi[i]; }, q, qs);
}

template <typename Eq, int N>
void Load(Complete<Eq>& q, const CompleteArray<Eq, N>& qs, int i) {
  ForEachComponent([&](auto& qi, auto qsi) { qi = qsi[i]; }, q, qs);
}

template <typename Eq, int N, typename Layout, std::size_t Rank>
void Load(ConservativeArray<Eq, N>& state,
          nodeduce_t<const BasicView<const Conservative<Eq>, Layout,
                                     static_cast<int>(Rank)>&>
              view,
          std::array<std::ptrdiff_t, Rank> index) {
  ForEachComponent(
      [&](auto& component, auto data) {
        component = std::apply(
            [&](auto... i) { return Load(int_constant<N>(), data, i...); },
            index);
      },
      state, view);
}

template <typename Eq, int N, typename Layout, int Rank>
void Load(
    ConservativeArray<Eq, N>& state,
    const BasicView<const Conservative<Eq>, Layout, Rank>& view,
    nodeduce_t<const std::array<std::ptrdiff_t, std::size_t(Rank)>&> index) {
  ForEachComponent(
      [&](auto&& component, auto data) {
        component = Load(int_constant<N>(), data, index);
      },
      state, view);
}

template <typename Eq>
void Load(ConservativeArray<Eq>& state,
          nodeduce_t<ViewPointer<const Conservative<Eq>>> pointer) {
  ForEachVariable(
      overloaded{
          [](Array1d& x, const double* p) { x = Eigen::Map<const Array1d>(p); },
          [](auto& xs, std::pair<const double*, std::ptrdiff_t> ps) {
            const double* p = ps.first;
            for (int i = 0; i < xs.rows(); ++i) {
              xs.row(i) = Eigen::Map<const Array1d>(p);
              p += ps.second;
            }
          }},
      state, pointer);
}

template <typename Eq, int N, typename Layout, int Rank>
void Load(
    CompleteArray<Eq, N>& state,
    const BasicView<const Complete<Eq>, Layout, Rank>& view,
    nodeduce_t<const std::array<std::ptrdiff_t, std::size_t(Rank)>&> index) {
  ForEachComponent(
      [&](auto&& component, auto data) {
        component = Load(int_constant<N>(), data, index);
      },
      state, view);
}

template <typename Eq>
void Load(CompleteArray<Eq>& state,
          nodeduce_t<ViewPointer<const Complete<Eq>>> pointer) {
  ForEachVariable(
      overloaded{
          [](Array1d& x, const double* p) { x = Eigen::Map<const Array1d>(p); },
          [](auto& xs, std::pair<const double*, std::ptrdiff_t> ps) {
            const double* p = ps.first;
            for (int i = 0; i < xs.rows(); ++i) {
              xs.row(i) = Eigen::Map<const Array1d>(p);
              p += ps.second;
            }
          }},
      state, pointer);
}

template <typename Eq, int N, typename Layout, int Rank>
void LoadN(
    CompleteArray<Eq, N>& state,
    const BasicView<const Complete<Eq>, Layout, Rank>& view, int size,
    nodeduce_t<const std::array<std::ptrdiff_t, std::size_t(Rank)>&> pos) {
  ForEachComponent(
      [&](auto&& s, auto v) { s = LoadN(int_constant<N>{}, v, size, pos); },
      state, view);
}

template <typename Eq, int N, typename Layout, int Rank>
void LoadN(
    ConservativeArray<Eq, N>& state,
    const BasicView<const Conservative<Eq>, Layout, Rank>& view, int size,
    nodeduce_t<const std::array<std::ptrdiff_t, std::size_t(Rank)>&> pos) {
  ForEachComponent(
      [&](auto& s, auto v) { s = LoadN(int_constant<N>{}, v, size, pos); },
      state, view);
}

template <typename Eq>
void LoadN(CompleteArray<Eq>& state,
           nodeduce_t<ViewPointer<const Complete<Eq>>> pointer, int n) {
  FUB_ASSERT(n <= kDefaultChunkSize);
  ForEachVariable(
      overloaded{[n](Array1d& x, const double* p) {
#ifdef __clang__
#pragma clang loop vectorize(disable)
#endif
                   for (int i = 0; i < n; ++i) {
                     x(i) = p[i];
                   }
                 },
                 [n](auto& xs, std::pair<const double*, std::ptrdiff_t> ps) {
                   const double* p = ps.first;
                   for (int i = 0; i < xs.rows(); ++i) {
#ifdef __clang__
#pragma clang loop vectorize(disable)
#endif
                     for (int j = 0; j < n; ++j) {
                       xs(i, j) = p[j];
                     }
                     p += ps.second;
                   }
                 }},
      state, pointer);
}

template <typename Eq>
void LoadN(ConservativeArray<Eq>& state,
           nodeduce_t<ViewPointer<const Conservative<Eq>>> pointer, int n) {
  FUB_ASSERT(n <= kDefaultChunkSize);
  ForEachVariable(
      overloaded{[n](Array1d& x, const double* p) {
#ifdef __clang__
#pragma clang loop vectorize(disable)
#endif
                   for (int i = 0; i < n; ++i) {
                     x(i) = p[i];
                   }
                 },
                 [n](auto& xs, std::pair<const double*, std::ptrdiff_t> ps) {
                   const double* p = ps.first;
                   for (int i = 0; i < xs.rows(); ++i) {
#ifdef __clang__
#pragma clang loop vectorize(disable)
#endif
                     for (int j = 0; j < n; ++j) {
                       xs(i, j) = p[j];
                     }
                     p += ps.second;
                   }
                 }},
      state, pointer);
}

template <typename Eq>
void Store(nodeduce_t<ViewPointer<Conservative<Eq>>> pointer,
           const ConservativeArray<Eq>& state) {
  ForEachVariable(overloaded{[](double* p, const Array1d& x) {
                               Eigen::Map<Array1d>{p} = x;
                             },
                             [](auto& ptr, const auto& x) {
                               for (int i = 0; i < x.rows(); ++i) {
                                 Eigen::Map<Array1d>(ptr.first) = x.row(i);
                                 ptr.first += ptr.second;
                               }
                             }},
                  pointer, state);
}

template <typename Eq>
void Store(nodeduce_t<ViewPointer<Complete<Eq>>> pointer,
           const CompleteArray<Eq>& state) {
  ForEachVariable(overloaded{[](double* p, const Array1d& x) {
                               Eigen::Map<Array1d>{p} = x;
                             },
                             [](auto& ptr, const auto& x) {
                               for (int i = 0; i < x.rows(); ++i) {
                                 Eigen::Map<Array1d>(ptr.first) = x.row(i);
                                 ptr.first += ptr.second;
                               }
                             }},
                  pointer, state);
}

template <typename Eq>
void StoreN(nodeduce_t<ViewPointer<Conservative<Eq>>> pointer,
            const ConservativeArray<Eq>& state, int n) {
  ForEachVariable(overloaded{[n](double* p, const Array1d& x) {
#ifdef __clang__
#pragma clang loop vectorize(disable)
#endif
                               for (int i = 0; i < n; ++i) {
                                 p[i] = x(i);
                               }
                             },
                             [n](auto& ptr, const auto& x) {
                               for (int i = 0; i < x.rows(); ++i) {
#ifdef __clang__
#pragma clang loop vectorize(disable)
#endif
                                 for (int j = 0; j < n; ++j) {
                                   ptr.first[i * ptr.second + j] = x(i, j);
                                 }
                               }
                             }},
                  pointer, state);
}

template <typename Eq>
void StoreN(nodeduce_t<ViewPointer<Complete<Eq>>> pointer,
            const CompleteArray<Eq>& state, int n) {
  ForEachVariable(overloaded{[n](double* p, const Array1d& x) {
#ifdef __clang__
#pragma clang loop vectorize(disable)
#endif
                               for (int i = 0; i < n; ++i) {
                                 p[i] = x(i);
                               }
                             },
                             [n](auto& ptr, const auto& x) {
                               for (int i = 0; i < x.rows(); ++i) {
#ifdef __clang__
#pragma clang loop vectorize(disable)
#endif
                                 for (int j = 0; j < n; ++j) {
                                   ptr.first[i * ptr.second + j] = x(i, j);
                                 }
                               }
                             }},
                  pointer, state);
}

template <typename Eq, int N, typename Layout, int Rank>
void Store(
    const BasicView<Conservative<Eq>, Layout, Rank>& view,
    const ConservativeArray<Eq, N>& state,
    nodeduce_t<const std::array<std::ptrdiff_t, std::size_t(Rank)>&> index) {
  ForEachComponent([&](auto data, auto block) { Store(data, block, index); },
                   view, state);
}

template <typename Eq, int N, typename Layout, int Rank>
void Store(
    const BasicView<Complete<Eq>, Layout, Rank>& view,
    const CompleteArray<Eq, N>& state,
    nodeduce_t<const std::array<std::ptrdiff_t, std::size_t(Rank)>&> index) {
  ForEachComponent([&](auto data, auto block) { Store(data, block, index); },
                   view, state);
}

template <typename Eq, int N, typename Layout, int Rank>
void StoreN(
    const BasicView<Complete<Eq>, Layout, Rank>& view,
    const CompleteArray<Eq, N>& state, int size,
    nodeduce_t<const std::array<std::ptrdiff_t, std::size_t(Rank)>&> pos) {
  ForEachComponent([&](auto& s, auto v) { StoreN(v, size, s, pos[0]); }, state,
                   view);
}

template <typename Eq, int N, typename Layout, int Rank>
void StoreN(
    const BasicView<Conservative<Eq>, Layout, Rank>& view,
    const ConservativeArray<Eq, N>& state, int size,
    nodeduce_t<const std::array<std::ptrdiff_t, std::size_t(Rank)>&> pos) {
  ForEachComponent([&](auto&& s, auto&& v) { StoreN(v, size, s, pos[0]); },
                   state, view);
}

} // namespace fub

#endif
