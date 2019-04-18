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

template <typename Eq, int Width> struct ConservativeArrayBaseImpl {
  template <typename T>
  using fn = typename DepthToStateValueTypeImpl<T, Width>::type;
  using type = boost::mp11::mp_transform<fn, typename Eq::ConservativeDepths>;
};

template <typename Eq, int Width>
using ConservativeArrayBase =
    typename ConservativeArrayBaseImpl<Eq, Width>::type;

template <typename Eq, int Width = kDefaultChunkSize>
struct ConservativeArray : ConservativeArrayBase<Eq, Width> {
  using Equation = Eq;
  using Depths = typename Equation::ConservativeDepths;
  using Traits = StateTraits<ConservativeArrayBase<Eq, Width>>;

  ConservativeArray() = default;
  ConservativeArray(const Equation&) {}
};

template <typename Eq, int Width> struct CompleteArrayBaseImpl {
  template <typename T>
  using fn = typename DepthToStateValueTypeImpl<T, Width>::type;
  using type = boost::mp11::mp_transform<fn, typename Eq::CompleteDepths>;
};

template <typename Eq, int Width>
using CompleteArrayBase = typename CompleteArrayBaseImpl<Eq, Width>::type;

template <typename Eq, int Width = kDefaultChunkSize>
struct CompleteArray : CompleteArrayBase<Eq, Width> {
  using Equation = Eq;
  using Depths = typename Equation::CompleteDepths;
  using Traits = StateTraits<CompleteArrayBase<Eq, Width>>;

  CompleteArray() = default;
  CompleteArray(const Equation&) {}
};

template <typename Eq, int N, typename Layout, std::size_t Rank>
void Load(ConservativeArray<Eq, N>& state,
          nodeduce_t<const BasicView<const Conservative<Eq>, Layout,
                                     static_cast<int>(Rank)>&>
              view,
          std::array<std::ptrdiff_t, Rank> index) {
  ForEachComponent<Conservative<Eq>>(
      [&](auto& component, auto data) {
        component = std::apply(
            [&](auto... i) { return Load(int_constant<N>(), data, i...); },
            index);
      },
      state, view);
}

template <typename Eq, int N, typename Layout, std::size_t Rank>
void Load(CompleteArray<Eq, N>& state,
          nodeduce_t<const BasicView<const Complete<Eq>, Layout,
                                     static_cast<int>(Rank)>&>
              view,
          std::array<std::ptrdiff_t, Rank> index) {
  ForEachComponent<Complete<Eq>>(
      [&](auto& component, auto data) {
        component = Load(int_constant<N>(), data, index);
      },
      state, view);
}

template <typename Eq, int N, typename Layout, std::size_t Rank>
void LoadN(CompleteArray<Eq, N>& state,
           nodeduce_t<const BasicView<const Complete<Eq>, Layout,
                                      static_cast<int>(Rank)>&>
               view,
           int size, const std::array<std::ptrdiff_t, Rank>& pos) {
  ForEachComponent<Complete<Eq>>(
      [&](auto& s, auto v) { s = LoadN(int_constant<N>{}, v, size, pos); },
      state, view);
}

template <typename Eq, int N, typename Layout, std::size_t Rank>
void LoadN(ConservativeArray<Eq, N>& state,
           nodeduce_t<const BasicView<const Conservative<Eq>, Layout,
                                      static_cast<int>(Rank)>&>
               view,
           int size, const std::array<std::ptrdiff_t, Rank>& pos) {
  ForEachComponent<Conservative<Eq>>(
      [&](auto& s, auto v) { s = LoadN(int_constant<N>{}, v, size, pos); },
      state, view);
}

template <typename Eq, int N, typename Layout, std::size_t Rank>
void Store(nodeduce_t<const BasicView<Conservative<Eq>, Layout,
                                      static_cast<int>(Rank)>&>
               view,
           const ConservativeArray<Eq, N>& state,
           std::array<std::ptrdiff_t, Rank> index) {
  ForEachComponent<Conservative<Eq>>(
      [&](auto data, auto block) { Store(data, block, index); }, view, state);
}

template <typename Eq, int N, typename Layout, std::size_t Rank>
void Store(
    nodeduce_t<const BasicView<Complete<Eq>, Layout, static_cast<int>(Rank)>&>
        view,
    const CompleteArray<Eq, N>& state, std::array<std::ptrdiff_t, Rank> index) {
  ForEachComponent<Complete<Eq>>(
      [&](auto data, auto block) { Store(data, block, index); }, view, state);
}

template <typename Eq, int N, typename Layout, std::size_t Rank>
void StoreN(
    nodeduce_t<const BasicView<Complete<Eq>, Layout, static_cast<int>(Rank)>&>
        view,
    const CompleteArray<Eq, N>& state, int size,
    const std::array<std::ptrdiff_t, Rank>& pos) {
  ForEachComponent<Complete<Eq>>(
      [&](auto& s, auto v) { StoreN(v, size, s, pos[0]); }, state, view);
}

template <typename Eq, int N, typename Layout, std::size_t Rank>
void StoreN(nodeduce_t<const BasicView<Conservative<Eq>, Layout,
                                       static_cast<int>(Rank)>&>
                view,
            const ConservativeArray<Eq, N>& state, int size,
            const std::array<std::ptrdiff_t, Rank>& pos) {
  ForEachComponent<Conservative<Eq>>(
      [&](auto& s, auto v) { StoreN(v, size, s, pos[0]); }, state, view);
}

} // namespace fub

#endif
