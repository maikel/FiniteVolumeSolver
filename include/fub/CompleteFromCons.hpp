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

#include "Execution.hpp"
#include "fub/State.hpp"
#include "fub/StateArray.hpp"
#include "fub/StateRow.hpp"

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

template <typename Equation> struct CompleteFromConsFn {
  Equation equation_;
  CompleteArray<Equation> complete_array_{equation_};
  ConservativeArray<Equation> cons_array_{equation_};
  Complete<Equation> complete_{equation_};
  Conservative<Equation> cons_{equation_};


  struct CompleteFromCons_Rows {
    CompleteFromConsFn<Equation>* this_;

    void operator()(const Row<Complete<Equation>>& complete_row,
                    const Row<const Conservative<Equation>>& cons_row) const {
      ViewPointer in = Begin(cons_row);
      ViewPointer end = End(cons_row);
      ViewPointer out = Begin(complete_row);
      Equation& equation = this_->equation_;
      CompleteArray<Equation>& complete = this_->complete_array_;
      ConservativeArray<Equation>& cons = this_->cons_array_;
      int n = static_cast<int>(get<0>(end) - get<0>(in));
      while (n >= kDefaultChunkSize) {
        Load(cons, in);
        ::fub::CompleteFromCons(equation, complete, cons);
        Store(out, complete);
        Advance(in, kDefaultChunkSize);
        Advance(out, kDefaultChunkSize);
        n = static_cast<int>(get<0>(end) - get<0>(in));
      }
      LoadN(cons, in, n);
      ::fub::CompleteFromCons(equation, complete, cons);
      StoreN(out, complete, n);
    }
  };

  void CompleteFromCons(execution::SimdTag, const View<Complete<Equation>>& complete_view,
                        const View<const Conservative<Equation>>& cons_view) {
    FUB_ASSERT(Box<0>(complete_view) == Box<0>(cons_view));
    ForEachRow(std::tuple{complete_view, cons_view},
               CompleteFromCons_Rows{this});
  }

  void CompleteFromCons(execution::SequentialTag,
                        const View<Complete<Equation>>& complete_view,
                        const View<const Conservative<Equation>>& cons_view) {
    FUB_ASSERT(Box<0>(complete_view) == Box<0>(cons_view));
    ForEachIndex(Box<0>(cons_view), [&](auto... is) {
      Load(cons_, cons_view, {is...});
      ::fub::CompleteFromCons(equation_, complete_, cons_);
      Store(complete_view, complete_, {is...});
    });
  }

  void CompleteFromCons(execution::OpenMpTag,
                        const View<Complete<Equation>>& complete_view,
                        const View<const Conservative<Equation>>& cons_view) {
    return CompleteFromCons(execution::seq, complete_view, cons_view);
  }

  void CompleteFromCons(execution::OpenMpSimdTag,
                        const View<Complete<Equation>>& complete_view,
                        const View<const Conservative<Equation>>& cons_view) {
    return CompleteFromCons(execution::simd, complete_view, cons_view);
  }
};

} // namespace fub

#endif
