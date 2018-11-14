// Copyright (c) 2017-2018 Maikel Nadolski
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

#ifndef FUB_CORE_FUNCTION_REF
#define FUB_CORE_FUNCTION_REF

#include "fub/core/type_traits.hpp"

#include <cstdint>
#include <type_traits>
#include <utility>

namespace fub {
/// An efficient, type-erasing, non-owning reference to a callable. This is
/// intended for use as the type of a function parameter that is not used
/// after the function in question returns.
///
/// This class does not own the callable, so it is not in general safe to store
/// a function_ref.
template <typename Fn> class function_ref;

template <typename F, typename... Args>
using IsInvocableT_ = decltype(std::declval<F>()(std::declval<Args>()...));

template <typename F, typename R, typename... Args>
struct IsInvocable : is_detected_exact<R, IsInvocableT_, F, Args...> {};

/// An efficient, type-erasing, non-owning reference to a callable. This is
/// intended for use as the type of a function parameter that is not used
/// after the function in question returns.
///
/// This class does not own the callable, so it is not in general safe to store
/// a function_ref.
template <typename Ret, typename... Params> class function_ref<Ret(Params...)> {
  Ret (*callback)(std::intptr_t callable, Params... params);
  std::intptr_t callable;

  template <typename Callable>
  static Ret callback_fn(std::intptr_t callable, Params... params) {
    return (*reinterpret_cast<Callable*>(callable))(
        std::forward<Params>(params)...);
  }

public:
  template <
      typename Callable,
      typename = std::enable_if_t<IsInvocable<Callable, Ret, Params...>::value>>
  function_ref(Callable&& callable,
               typename std::enable_if<
                   !std::is_same<typename std::remove_reference<Callable>::type,
                                 function_ref>::value>::type* = nullptr)
      : callback(callback_fn<typename std::remove_reference<Callable>::type>),
        callable(reinterpret_cast<std::intptr_t>(&callable)) {}
  Ret operator()(Params... params) const {
    return callback(callable, std::forward<Params>(params)...);
  }
};
} // namespace fub

#endif