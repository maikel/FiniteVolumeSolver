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

template <typename Equation> struct CompleteFromConsImpl {
  static void apply(const Equation&, Complete<Equation>& complete,
                    const Conservative<Equation>& cons) {
    ForEachVariable<Conservative<Equation>>(
        [](auto& dest, const auto& src) { dest = src; }, complete, cons);
  }

  static void apply(const Equation&, Complete<Equation>&,
                    const Complete<Equation>&) {}
};

template <typename Equation>
void CompleteFromCons(Equation&& equation,
                      Complete<std::decay_t<Equation>>& complete,
                      const Conservative<std::decay_t<Equation>>& cons) {
  return CompleteFromConsImpl<std::decay_t<Equation>>::apply(equation, complete,
                                                             cons);
}

template <typename Equation>
void CompleteFromCons(Equation&& equation,
                      Complete<std::decay_t<Equation>>& complete,
                      const Complete<std::decay_t<Equation>>& cons) {
  return CompleteFromConsImpl<std::decay_t<Equation>>::apply(equation, complete,
                                                             cons);
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
