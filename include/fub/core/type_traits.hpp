// Copyright (c) 2017 Maikel Nadolski
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

/// @file type_traits.hpp This file adds basic type traits utilities which are
/// not yet implemented in all standard libraries.

#ifndef FUB_CORE_TYPE_TRAITS_HPP
#define FUB_CORE_TYPE_TRAITS_HPP

#include <type_traits>
#include <utility>

namespace fub {
#if defined(__cpp_lib_byte) && __cpp_lib_byte >= 201603
using std::byte;
#else
using byte = unsigned char;
#endif

////////////////////////////////////////////////////////////////////////////////
//                                                          [traits.is_detected]

/// \defgroup type-traits Type Traits
/// Type traits are used to enforce requirements on types with SFINAE.

template <class...> using void_t = void;

template <typename T> struct nodeduce { using type = T; };
template <typename T> using nodeduce_t = typename nodeduce<T>::type;

struct nonesuch {
  nonesuch() = delete;
  ~nonesuch() = delete;
  nonesuch(nonesuch const&) = delete;
  void operator=(nonesuch const&) = delete;
};

namespace detail {
template <class Default, class AlwaysVoid, template <class...> class Op,
          class... Args>
struct detector {
  using value_t = std::false_type;
  using type = Default;
};

template <class Default, template <class...> class Op, class... Args>
struct detector<Default, void_t<Op<Args...>>, Op, Args...> {
  using value_t = std::true_type;
  using type = Op<Args...>;
};

template <template <class...> class Op, class... Args>
using is_detected =
    typename detail::detector<nonesuch, void, Op, Args...>::value_t;

/// \ingroup type-traits
/// Returns the type of `Op<Args...>` or `nonesuch`
template <template <class...> class Op, class... Args>
using detected_t = typename detail::detector<nonesuch, void, Op, Args...>::type;

template <class Default, template <class...> class Op, class... Args>
using detected_or = detail::detector<Default, void, Op, Args...>;

template <class Default, template <class...> class Op, class... Args>
using detected_or_t = typename detected_or<Default, Op, Args...>::type;

template <class Expected, template <typename...> class Op, class... Args>
using is_detected_exact = std::is_same<Expected, detected_t<Op, Args...>>;
} // namespace detail

/// \ingroup type-traits
/// This is `std::true_type` if `Op<Args...>` is a valid SFINAE expression.
template <template <class...> class Op, class... Args>
struct is_detected : detail::is_detected<Op, Args...> {};

/// \ingroup type-traits
/// Returns the type of `Op<Args...>` or `nonesuch`
template <template <class...> class Op, class... Args>
using detected_t = detail::detected_t<Op, Args...>;

/// \ingroup type-traits
/// Returns the type of `Op<Args...>` or `Default`
template <class Default, template <class...> class Op, class... Args>
struct detected_or : detail::detected_or<Default, Op, Args...> {};

/// \ingroup type-traits
/// This is `std::true_type` if `Op<Args...>` is a valid SFINAE expression and
/// the return type is exactly `Expected`.
template <class Expected, template <typename...> class Op, class... Args>
struct is_detected_exact : detail::is_detected_exact<Expected, Op, Args...> {};

////////////////////////////////////////////////////////////////////////////////
//                                                          [traits.conjunction]
//                                                          [traits.disjunction]
//                                                          [traits.negation]
#if defined(__cpp_lib_logical_traits)
using std::conjunction;
using std::disjunction;
using std::negation;
#elif defined(__cpp_lib_experimental_logical_traits)
using std::experimental::conjunction;
using std::experimental::disjunction;
using std::experimental::negation;
#else
template <class...> struct conjunction : std::true_type {};
template <class B1> struct conjunction<B1> : B1 {};
template <class B1, class... Bn>
struct conjunction<B1, Bn...>
    : std::conditional_t<bool(B1::value), conjunction<Bn...>, B1> {};

template <class...> struct disjunction : std::false_type {};
template <class B1> struct disjunction<B1> : B1 {};
template <class B1, class... Bn>
struct disjunction<B1, Bn...>
    : std::conditional_t<bool(B1::value), B1, disjunction<Bn...>> {};

template <class Bool>
using negation = std::integral_constant<bool, !bool(Bool::value)>;
#endif

#if defined(__cpp_nontype_template_parameter_auto)
template <auto N> struct constant : std::integral_constant<decltype(N), N> {};
#endif

#if defined(__cpp_lib_integer_sequence)
using std::index_sequence;
using std::integer_sequence;
using std::make_index_sequence;
using std::make_integer_sequence;
#else
template <typename T, T... Ints> struct integer_sequence {
  using value_type = T;

  static constexpr std::size_t size() noexcept { return sizeof...(Ints); }
};

template <std::size_t... Is>
using index_sequence = integer_sequence<std::size_t, Is...>;

template <typename T, T N, T Last, T... Ints> struct MakeIntegerSequence_;

template <typename T, T N, T... Ints>
struct MakeIntegerSequence_<T, N, N, Ints...> {
  using type = integer_sequence<T, Ints...>;
};

template <typename T, T N, T L, T... Ints> struct MakeIntegerSequence_ {
  using type = typename MakeIntegerSequence_<T, N, L + 1, Ints..., L>::type;
};

template <typename T, T N>
using make_integer_sequence = typename MakeIntegerSequence_<T, N, T{0}>::type;

template <std::size_t N>
using make_index_sequence = make_integer_sequence<std::size_t, N>;
#endif

////////////////////////////////////////////////////////////////////////////////
//                                                               [meta.constant]

using index = std::ptrdiff_t;

template <bool Bool> using bool_constant = std::integral_constant<bool, Bool>;
template <index I> using index_constant = std::integral_constant<index, I>;
template <int I> using int_constant = std::integral_constant<int, I>;
template <std::size_t I>
using size_constant = std::integral_constant<std::size_t, I>;
template <int I> static constexpr int_constant<I> int_c{};
template <bool B> static constexpr bool_constant<B> bool_c{};
template <index I> static constexpr index_constant<I> index_c{};
template <std::size_t I> static constexpr size_constant<I> size_c{};

////////////////////////////////////////////////////////////////////////////////
//                                                         [traits.is_invocable]
//                                                        [traits.invoke_result]

#ifdef __cpp_lib_is_invocable
using std::invoke_result;
using std::invoke_result_t;
using std::is_invocable;
#else
/// \ingroup type-traits
/// This is `std::true_type` if `F` is a function type and can be invoked with
/// arguments of types `Args...`.
/// @{
template <typename F, typename... Args>
struct invoke_result : std::result_of<F(Args...)> {};

template <typename F, typename... Args>
using invoke_result_t = typename invoke_result<F, Args...>::type;
/// @}

/// \ingroup type-traits
/// This is `std::true_type` if a given object `f` of type `T` is callable by
/// `fub::invoke(f, args...)` for `Args... args`.
template <typename F, typename... Args>
struct is_invocable : is_detected<invoke_result_t, F, Args...> {};
#endif

////////////////////////////////////////////////////////////////////////////////
//                                                         [traits.remove_cvref]

/// \ingroup type-traits
/// This is equivalent to `std::remove_cv_t<std::remove_reference_t<T>>`.
template <typename T> struct remove_cvref {
  using type = std::remove_cv_t<std::remove_reference_t<T>>;
};
template <typename T> using remove_cvref_t = typename remove_cvref<T>::type;

////////////////////////////////////////////////////////////////////////////////
//                                                             [meta.value_type]

template <typename Var> struct value_type {
  using type = typename Var::value_type;
};
template <typename Var> using value_type_t = typename value_type<Var>::type;

////////////////////////////////////////////////////////////////////////////////
//                                               [traits.is_equality_comparable]
//                                       [traits.is_nothrow_equality_comparable]

template <typename T, typename S = T> struct is_equality_comparable_impl {
  template <typename L, typename R>
  using equality_t = decltype(std::declval<L>() == std::declval<R>());

  template <typename L, typename R>
  using inequality_t = decltype(std::declval<L>() != std::declval<R>());

  using type = conjunction<is_detected<equality_t, T, S>,
                           is_detected<inequality_t, T, S>>;
};

/// \ingroup type-traits
/// This is `std::true_type` if `T` and `S` are equality comparable.
template <typename T, typename S = T>
struct is_equality_comparable : is_equality_comparable_impl<T, S>::type {};

namespace detail {
template <typename T, typename S, bool IsComparable>
struct is_nothrow_equality_comparable_impl : bool_constant<false> {};

template <typename T, typename S>
struct is_nothrow_equality_comparable_impl<T, S, true> {
  static constexpr bool is_nothrow_equal =
      noexcept(std::declval<T>() == std::declval<S>());

  static constexpr bool is_nothrow_inequal =
      noexcept(std::declval<T>() != std::declval<S>());

  static constexpr bool value = is_nothrow_equal && is_nothrow_inequal;
};
} // namespace detail

/// \ingroup type-traits
/// This is `std::true_type` if `T` and `S` are nothrow equality comparable.
template <typename T, typename S = T> struct is_nothrow_equality_comparable {
  static constexpr bool value = detail::is_nothrow_equality_comparable_impl<
      T, S, is_equality_comparable<T, S>::value>::value;
};

////////////////////////////////////////////////////////////////////////////////
//                                                           [traits.is_regular]

/// \ingroup type-traits
/// This type trait checks if a specified type `T` fulfills the Regular concept.
template <typename T>
struct is_regular
    : conjunction<std::is_default_constructible<T>, std::is_copy_assignable<T>,
                  is_equality_comparable<T>> {};

// helper type for the visitor #4
#ifdef __cpp_deduction_guides
template <class... Ts> struct overloaded : Ts... { using Ts::operator()...; };
template <class... Ts> overloaded(Ts...)->overloaded<Ts...>;
#endif

} // namespace fub

#endif // !TYPE_TRAITS_HPP
