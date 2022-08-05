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

template <typename Equation, typename State>
decltype(auto) GetComponent(const Equation& eq, const State& s, int ncomp) {
  auto depths = Depths(eq, Type<Complete<Equation>>{});
  int counter = ncomp;
  const double* pointer_to_component = nullptr;
  ForEachVariable(overloaded{
    [&](auto& component, ScalarDepth) {
      if (counter == 0) {
        pointer_to_component = &component;
      }
      counter -= 1;
    },
    [&](auto& variable, auto depth) {
      if (0 <= counter && counter < depth) {
        pointer_to_component = &variable[counter];
      }
      counter -= depth;
    }
  }, s, depths);
  FUB_ASSERT(pointer_to_component);
  return *pointer_to_component;
}

template <typename Equation, typename State>
decltype(auto) GetComponent(const Equation& eq, State& s, int ncomp) {
  auto depths = Depths(eq, Type<Complete<Equation>>{});
  int counter = ncomp;
  double* pointer_to_component = nullptr;
  ForEachVariable(overloaded{
    [&](auto& component, ScalarDepth) {
      if (counter == 0) {
        pointer_to_component = &component;
      }
      counter -= 1;
    },
    [&](auto& variable, auto depth) {
      if (0 <= counter && counter < depth) {
        pointer_to_component = &variable[counter];
      }
      counter -= depth;
    }
  }, s, depths);
  FUB_ASSERT(pointer_to_component);
  return *pointer_to_component;
}

template <typename Equation, typename State>
int GetNumberOfComponents(const Equation& eq, const State&) {
  auto depths = Depths(eq, Type<State>{});
  int counter = 0;
  ForEachVariable(
      overloaded{
          [&](ScalarDepth) { ++counter; },
          [&](auto depth) {
            counter += depth;
          },
      }, depths);
  return counter;
}
}

#endif