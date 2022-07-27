#ifndef FUB_STATE_UTIL_HPP
#define FUB_STATE_UTIL_HPP

#include <type_traits>

namespace fub {
template <typename Equation, typename Complete, typename State>
void CompleteFromState(Equation& eq, Complete& dest,
                          const State& source) {
  if constexpr (std::is_same_v<State, Conservative<Equation>>) {
    CompleteFromCons(eq, dest, source);
  } else {
    CompleteFromPrim(eq, dest, source);
  }
}

template <typename Equation, typename Complete, typename State>
void StateFromComplete(Equation& eq, State& dest, const Complete& src) {
  if constexpr (std::is_same_v<State, Conservative<Equation>>) {
    dest = AsCons(src);
  } else {
    PrimFromComplete(eq, dest, src);
  }
}
}

#endif