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

#ifndef FUB_TUPLE_HPP
#define FUB_TUPLE_HPP

#include <array>
#include <tuple>

namespace fub {

template <std::size_t... Is, typename... Ts>
[[nodiscard]] constexpr auto Take(std::index_sequence<Is...>,
                                  const std::tuple<Ts...>& t) {
  return std::tuple{std::get<Is>(t)...};
}

template <std::size_t N, typename... Ts>
[[nodiscard]] constexpr auto Take(const std::tuple<Ts...>& t) {
  return Take(std::make_index_sequence<N>(), t);
}

template <std::size_t N, std::size_t... Is, typename... Ts>
[[nodiscard]] constexpr auto Drop(std::index_sequence<Is...>,
                                  const std::tuple<Ts...>& t) {
  return std::tuple{std::get<N + Is>(t)...};
}

template <std::size_t N, typename... Ts>
[[nodiscard]] constexpr auto Drop(const std::tuple<Ts...>& t) {
  constexpr std::size_t size = sizeof...(Ts);
  return Drop<N>(std::make_index_sequence<size - N>(), t);
}

template <typename Tuple, typename Function>
[[nodiscard]] constexpr auto Transform(Tuple&& tuple, Function f) {
  return std::apply([&](auto&&... xs) { return std::tuple{f(xs)...}; }, tuple);
}

template <std::size_t... Is, typename T, typename... Ts>
[[nodiscard]] constexpr auto AsArray(std::index_sequence<Is...>,
                                     const std::tuple<T, Ts...>& t) {
  return std::array<T, sizeof...(Is)>{std::get<Is>(t)...};
}

template <typename T, typename... Ts>
[[nodiscard]] constexpr auto AsArray(const std::tuple<T, Ts...>& t) {
  return AsArray(std::make_index_sequence<sizeof...(Ts) + 1>(), t);
}

template <typename TupleLike>
[[nodiscard]] constexpr auto AsTuple(const TupleLike& t) {
  return std::apply([](const auto&... xs) { return std::tuple{xs...}; }, t);
}

} // namespace fub

#endif
