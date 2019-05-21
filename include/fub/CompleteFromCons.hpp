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

#ifndef FUB_COMPLETE_FROM_CONS_HPP
#define FUB_COMPLETE_FROM_CONS_HPP

#include "fub/State.hpp"
#include "fub/StateArray.hpp"

#include <cstring>
#include <type_traits>

namespace fub {
template <typename Eq, typename... Args>
using CompleteFromConsMemberFunction =
    decltype(std::declval<Eq>().CompleteFromCons(std::declval<Args>()...));

template <typename Equation>
void CompleteFromCons(Equation&& equation,
                      Complete<std::decay_t<Equation>>& complete,
                      const Conservative<std::decay_t<Equation>>& cons) {
  using Eq = std::decay_t<Equation>;
  if constexpr (is_detected<CompleteFromConsMemberFunction, Equation,
                            Complete<Eq>&, const Conservative<Eq>&>::value) {
    equation.CompleteFromCons(complete, cons);
  } else {
    static_assert(sizeof(Complete<Eq>) == sizeof(Conservative<Eq>));
    ForEachVariable([](auto& dest, const auto& src) { dest = src; },
                    AsCons(complete), cons);
  }
}

template <typename Equation>
void CompleteFromCons(Equation&& equation,
                      Complete<std::decay_t<Equation>>& complete,
                      const Complete<std::decay_t<Equation>>& cons) {
  using Eq = std::decay_t<Equation>;
  if constexpr (is_detected<CompleteFromConsMemberFunction, Equation,
                            Complete<Eq>&,
                            const ConservativeBase<Eq>&>::value) {
    equation.CompleteFromCons(complete, AsCons(cons));
  } else {
    static_assert(sizeof(Complete<Eq>) == sizeof(Conservative<Eq>));
    ForEachVariable([](auto& dest, const auto& src) { dest = src; },
                    AsCons(complete), AsCons(cons));
  }
}

template <typename Equation>
void CompleteFromCons(Equation&& equation,
                      CompleteArray<std::decay_t<Equation>>& complete,
                      const ConservativeArray<std::decay_t<Equation>>& cons) {
  using Eq = std::decay_t<Equation>;
  if constexpr (is_detected<CompleteFromConsMemberFunction, Equation,
                            CompleteArray<Eq>&,
                            const ConservativeArray<Eq>&>::value) {
    equation.CompleteFromCons(complete, cons);
  } else {
    static_assert(sizeof(Complete<Eq>) == sizeof(Conservative<Eq>));
    ForEachVariable([](auto& dest, const auto& src) { dest = src; },
                    AsCons(complete), cons);
  }
}

template <typename Equation>
void CompleteFromCons(Equation&& equation,
                      CompleteArray<std::decay_t<Equation>>& complete,
                      const CompleteArray<std::decay_t<Equation>>& cons) {
  using Eq = std::decay_t<Equation>;
  if constexpr (is_detected<CompleteFromConsMemberFunction, Equation,
                            CompleteArray<Eq>&,
                            const ConservativeArray<Eq>&>::value) {
    equation.CompleteFromCons(complete, AsCons(cons));
  } else {
    static_assert(sizeof(Complete<Eq>) == sizeof(Conservative<Eq>));
    ForEachVariable([](auto& dest, const auto& src) { dest = src; },
                    AsCons(complete), AsCons(cons));
  }
}

template <typename Equation>
void CompleteFromCons(
    Equation&& eq, const View<Complete<std::decay_t<Equation>>>& complete_view,
    const View<const Conservative<std::decay_t<Equation>>>& cons_view) {
  Complete<std::decay_t<Equation>> complete(eq);
  Conservative<std::decay_t<Equation>> cons(eq);
  ForEachIndex(Box<0>(complete_view), [&](auto... is) {
    Load(cons, cons_view, {is...});
    CompleteFromCons(eq, complete, cons);
    Store(complete_view, complete, {is...});
  });
}

} // namespace fub

#endif
