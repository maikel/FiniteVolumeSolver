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

#ifndef FUB_CORE_SPAN_HPP
#define FUB_CORE_SPAN_HPP

#include "fub/core/assert.hpp"
#include "fub/core/dynamic_extent.hpp"
#include "fub/core/type_traits.hpp"

#include <array>
#include <type_traits>

namespace fub {
/// \defgroup spans Spans
/// Types and function which are about non-owning ranges of contiguous memory.

namespace detail {
/// @{
/// Tests if the specified type `T` is a `std::array<U, N>` for some U and N.
template <typename T> struct is_std_array : std::false_type {};
template <typename T, std::size_t N>
struct is_std_array<std::array<T, N>> : std::true_type {};
/// @}

template <typename T> using data_t = decltype(std::declval<T>().data());
template <typename T> using element_type_ = std::remove_pointer_t<data_t<T>>;
template <typename T> using array_type_t = element_type_<T> (*)[];
} // namespace detail

template <typename T, std::ptrdiff_t N> class span;

///////////////////////////////////////////////////////////////////////////////
//                                                        [span.traits.is_span]

/// Returns true if the specified T is a `span<S, N>` for some type S and
/// integer N.
template <typename T> struct is_span : std::false_type {};

/// Returns true if the specified T is a `span<S, N>` for some type S and
/// integer N.
template <typename T, std::ptrdiff_t N>
struct is_span<span<T, N>> : std::true_type {};

/// Returns true if the specified T is a `span<S, N>` for some type S and
/// integer N.
template <typename T> static constexpr bool is_span_v = is_span<T>::value;

///////////////////////////////////////////////////////////////////////////////
//                                                                 [span.class]

/// \ingroup spans
/// A span is a view over a contiguous sequence of objects, the storage of which
/// is owned by some other object.
///
/// All member functions of span have constant time complexity.
///
/// `T` is required to be a complete object type that is not an abstract
/// class type.
///
/// If `N` is negative and not equal to dynamic_­extent, the program is
/// ill-formed.
template <typename T, std::ptrdiff_t N = dynamic_extent> class span {
public:
  static_assert(N > 0, "span's extent needs to be non-negative.");

  using element_type = T;
  using pointer = T*;
  using reference = T&;
  using value_type = std::remove_cv_t<T>;
  using index_type = std::ptrdiff_t;
  using difference_type = std::ptrdiff_t;
  using iterator = pointer;
  using const_iterator = std::add_const_t<T>*;
  using reverse_iterator = std::reverse_iterator<iterator>;
  using const_reverse_iterator = std::reverse_iterator<const_iterator>;

  static constexpr index_type extent = N;

  /// \name  Constructors, copy, and assignment [span.cons]

  /// Constructs a span from a `pointer` and `size` pair.
  ///
  /// This performs an assertion check in debug builds which will terminate the
  /// application if the specified size does not match the extent.
  ///
  /// \param[in] p  a valid address pointing to the first element of the span.
  /// \param[in] size  the size such that `[p, p + size)` is a valid range.
  ///
  /// \pre `size == extent`
  ///
  /// \post `data() == p`
  /// \post `size() == size`
  ///
  /// \throws Nothing.
  constexpr span(pointer p, [[maybe_unused]] index_type size) : pointer_{p} {
    FUB_ASSERT(size == extent);
  }

  /// Constructs a span from two pointers.
  ///
  /// This performs an assertion check in debug builds which will terminate the
  /// application if the specified size does not match the extent.
  ///
  /// \param[in] first  a valid address pointing to the first element of the
  /// span. \param[in] last  a pointer such that `[first, last)` is a valid
  /// range.
  ///
  /// \pre `last - first == extent`
  ///
  /// \post `data() == first`
  /// \post `size() == last - first`
  ///
  /// \throws Nothing.
  constexpr span(pointer first, pointer last) : span(first, last - first) {}

  /// Implicit conversion from a built-in C-style array.
  ///
  /// \param[in] arr  A native array which will be viewed on.
  ///
  /// \post `data() = &arr[0]`
  /// \post `size() == extent`
  ///
  /// \throws Nothing.
  constexpr span(element_type (&arr)[extent]) noexcept // NOLINT
      : pointer_{&arr[0]} {}

  /// Implicit conversion from const `std::array`s.
  ///
  /// \param[in] arr  A `std::array` which will be viewed on.
  ///
  /// \pre `std::is_const<element_type>`
  ///
  /// \post `data() = arr.data()`
  /// \post `size() == extent`
  ///
  /// \throws Nothing.
  template <typename ValueType, std::size_t Extent,
            typename = std::enable_if_t<Extent == extent>,
            typename = std::enable_if_t<std::is_convertible<
                const ValueType (*)[], element_type (*)[]>::value>>
  constexpr span(const std::array<ValueType, Extent>& arr) noexcept // NOLINT
      : pointer_{arr.data()} {}

  /// Implicit conversion from mutable `std::array`s.
  ///
  /// \param[in] arr  A `std::array` which will be viewed on.
  ///
  /// \pre `M == extent`
  ///
  /// \post `data() = arr.data()`
  /// \post `size() == extent`
  ///
  /// \throws Nothing.
  template <typename ValueType, std::size_t Extent,
            typename = std::enable_if_t<Extent == extent>,
            typename = std::enable_if_t<std::is_convertible<
                ValueType (*)[], element_type (*)[]>::value>>
  constexpr span(std::array<ValueType, Extent>& arr) noexcept // NOLINT
      : pointer_{arr.data()} {}

  /// Implicit conversion operator from a mutable container.
  ///
  /// \param[in] container  A container which owns the contiguous memory. The
  ///            container must outlive this span, otherwise the span is
  ///            dangling.
  ///
  /// \tparam Container  A container class which has the class member functions
  ///                    `Container::data()` and `Container::size()`.
  ///
  /// \pre `container.size() == extent`
  /// \pre `!is_span_v<Container>`
  ///
  /// \post  `data() == container.data()`
  /// \post `size() == container.size() == extent`
  ///
  /// \throws Nothing.
#ifdef DOXYGEN
  template <typename Container>
#else
  template <typename Container,
            typename = std::enable_if_t<!std::is_array<Container>::value>,
            typename = std::enable_if_t<
                !detail::is_std_array<std::remove_const_t<Container>>::value>,
            typename = std::enable_if_t<!is_span<Container>::value>,
            typename = std::enable_if_t<std::is_convertible<
                detected_t<detail::array_type_t, Container&>,
                element_type (*)[]>::value>>
#endif
  constexpr span(Container& container) // NOLINT
      : pointer_{container.data()} {
    FUB_ASSERT(container.size() == extent);
  }

  /// Implicit conversion operator from a constant container.
  ///
  /// \param[in] container  A container which owns the contiguous memory. The
  ///            container must outlive this span, otherwise the span is
  ///            dangling.
  ///
  /// \tparam Container  A container class which has the class member functions
  ///                    `Container::data()` and `Container::size()`.
  ///
  /// \pre `extent <= container.size()`
  /// \pre `!is_span_v<Container>`
  ///
  /// \post `data() == container.data()`
  /// \post `size() == container.size() == extent`
  ///
  /// \throws Nothing.
#ifdef DOXYGEN
  template <typename Container>
#else
  template <
      typename Container,
      typename = std::enable_if_t<!std::is_array<Container>::value>,
      typename = std::enable_if_t<!detail::is_std_array<Container>::value>,
      typename = std::enable_if_t<!is_span<Container>::value>,
      typename = std::enable_if_t<std::is_convertible<
          detected_t<detail::array_type_t, const Container&>,
          element_type (*)[]>::value>>
#endif
  constexpr span(const Container& container) // NOLINT
      : pointer_{container.data()} {
    FUB_ASSERT(container.size() == extent);
  }

  /// Implicit conversion from other span types.
  ///
  /// \tparam S  the element type of the other span.
  /// \tparam M  the static extents of the other span.
  ///
  /// \pre s.size() == extent
  /// \pre std::is_convertible<S(*)[], T(*)[]>::value
  ///
  /// \post `data() == s.data()`
  /// \post `size() == s.size() == extent`.
#ifdef DOXYGEN
  template <typename S, std::ptrdiff_t M>
#else
  template <
      typename S, std::ptrdiff_t M,
      typename = std::enable_if_t<std::is_convertible<
          typename span<S, M>::element_type (*)[], element_type (*)[]>::value>,
      typename = std::enable_if_t<extent == M>>
#endif
  constexpr span(const span<S, M>& s) noexcept // NOLINT
      : pointer_{s.data()} {
  }

  /// Defaulted copy constructor to trivially copy the class member variables.
  ///
  /// \throws Nothing.
  constexpr span(const span& s) noexcept = default;

  /// \name Subviews [span.sub]

  /// Returns a span of the first `Count`-many elements of this span.
  ///
  /// \tparam Count  Size of the returned subspan.
  ///
  /// \pre `0 <= Count && Count <= extent`.
  ///
  /// \throws Nothing.
  template <ptrdiff_t Count,
            typename = std::enable_if_t<(0 <= Count && Count <= extent)>>
  constexpr span<element_type, Count> first() const {
    return {pointer_, Count};
  }

  /// Returns a span of the first `count`-many elements of this span.
  ///
  /// \param[in] count  Size of the returned subspan.
  ///
  /// \pre `0 <= count && count <= extent`.
  ///
  /// \throws Nothing.
  constexpr span<element_type> first(index_type count) const {
    FUB_ASSERT(0 <= count && count <= size());
    return {pointer_, count};
  }

  /// Returns a span of the last `Count`-many elements of this span.
  ///
  /// \tparam Count  Size of the returned subspan.
  ///
  /// \pre `0 <= Count && Count <= extent`.
  ///
  /// \throws Nothing.
  template <ptrdiff_t Count,
            typename = std::enable_if_t<(0 <= Count && Count <= extent)>>
  constexpr span<element_type, Count> last() const {
    return {pointer_ + (extent - Count), Count};
  }

  /// Returns a span of the last `count`-many elements of this span.
  ///
  /// \param[in] count  Size of the returned subspan.
  ///
  /// \pre `0 <= count && count <= extent`.
  ///
  /// \throws Nothing.
  constexpr span<element_type> last(index_type count) const {
    FUB_ASSERT(0 <= count && count <= extent);
    return {pointer_ + (extent - count), count};
  }

  /// Returns a subspan viewing `Count` many elements from offset `Offset`.
  ///
  /// \tparam Offset  the offset of the returned subspan.
  /// \tparam Count  Size of the returned subspan.
  ///
  /// \pre `0 <= Count && Count <= extent`.
  ///
  /// \throws Nothing.
  template <
      ptrdiff_t Offset, ptrdiff_t Count = dynamic_extent,
      typename = std::enable_if_t<(0 <= Offset && Offset <= extent)>,
      typename = std::enable_if_t<(Count == dynamic_extent ||
                                   (0 <= Count && Count + Offset <= extent))>>
  constexpr auto subspan() const {
    if constexpr (Count == dynamic_extent) {
      return span<T, extent - Offset>{pointer_ + Offset, extent - Offset};
    } else {
      return span<T, Count>{pointer_ + Offset, Count};
    }
  }

  /// Returns a subspan viewing `count` many elements from offset `offset`.
  ///
  /// \param offset  the offset of the returned subspan.
  /// \param count  Size of the returned subspan.
  ///
  /// \pre `0 <= offset && offset <= extent`.
  /// \pre `count == dynamic_extent || 0 <= count && count + offset <= size()`.
  ///
  /// \throws Nothing.
  constexpr span<element_type>
  subspan(index_type offset, index_type count = dynamic_extent) const {
    FUB_ASSERT(0 <= offset && offset <= size());
    FUB_ASSERT(count == dynamic_extent ||
               (0 <= count && count + offset <= size()));
    return {pointer_ + offset,
            count == dynamic_extent ? size() - offset : count};
  }

  /// \name Observers [span.obs]

  /// Returns the number of elements in the span.
  ///
  /// \throws Nothing.
  constexpr index_type size() const noexcept { return extent; }

  /// Returns the number of bytes which are spanned by this span.
  ///
  /// \throws Nothing.
  constexpr index_type size_bytes() const noexcept {
    return sizeof(T) * extent;
  }

  /// Returns true if  `size() == 0`.
  ///
  /// \throws Nothing.
  constexpr bool empty() const noexcept { return false; }

  /// \name Element access [span.elem]

  /// Returns the underlying pointer.
  ///
  /// \throws Nothing.
  constexpr pointer data() const noexcept { return pointer_; }

  /// Accesses the n-th element of the spanned array.
  ///
  /// \throws Nothing.
  constexpr reference operator[](index_type n) const { return pointer_[n]; }

  /// Accesses the n-th element of the spanned array.
  ///
  /// \throws Nothing.
  constexpr reference operator()(index_type n) const { return pointer_[n]; }

  /// \name Iterator support [span.iterators]

  /// Returns an iterator pointing to the first element of the span.
  ///
  /// \throws Nothing.
  constexpr iterator begin() const noexcept { return pointer_; }

  /// Returns a const iterator pointing to the first element of the span.
  ///
  /// \throws Nothing.
  constexpr const_iterator cbegin() const noexcept { return pointer_; }

  /// Returns an iterator pointing one after the last element of the span.
  ///
  /// \throws Nothing.
  constexpr iterator end() const noexcept { return pointer_ + extent; }

  /// Returns a const iterator pointing one after the last element of the span.
  ///
  /// \throws Nothing.
  constexpr const_iterator cend() const noexcept { return pointer_ + extent; }

  /// Returns a reverse iterator pointing to the last element of the span.
  ///
  /// \throws Nothing.
  constexpr reverse_iterator rbegin() const noexcept {
    return std::make_reverse_iterator(end());
  }

  /// Returns a const reverse iterator pointing to the last element of the span.
  ///
  /// \throws Nothing.
  constexpr const_reverse_iterator crbegin() const noexcept {
    return std::make_reverse_iterator(cend());
  }

  /// Returns a reverse iterator pointing to the first element of the span.
  ///
  /// \throws Nothing.
  constexpr reverse_iterator rend() const noexcept {
    return std::make_reverse_iterator(begin());
  }

  /// Returns a const reverse iterator pointing to the first element of the
  /// span.
  ///
  /// \throws Nothing.
  constexpr const_reverse_iterator crend() const noexcept {
    return std::make_reverse_iterator(cbegin());
  }

private:
  pointer pointer_;
};

/// \ingroup spans
/// A span is a view over a contiguous sequence of objects, the storage of which
/// is owned by some other object.
///
/// All member functions of span have constant time complexity.
///
/// `T` is required to be a complete object type that is not an abstract
/// class type.
///
/// If `N` is negative and not equal to dynamic_­extent, the program is
/// ill-formed.
template <typename T> class span<T, 0> {
public:
  using element_type = T;
  using pointer = T*;
  using reference = T&;
  using value_type = std::remove_cv_t<T>;
  using index_type = std::ptrdiff_t;
  using difference_type = std::ptrdiff_t;
  using iterator = pointer;
  using const_iterator = std::add_const_t<T>*;
  using reverse_iterator = std::reverse_iterator<iterator>;
  using const_reverse_iterator = std::reverse_iterator<const_iterator>;

  static constexpr index_type extent = 0;

  /// \name  Constructors, copy, and assignment [span.cons]

  /// Constructs an empty span.
  ///
  /// \post `data() == nullptr`
  /// \post `size() == 0`
  ///
  /// \throws Nothing.
  constexpr span() noexcept = default;

  /// Constructs a span from a `pointer` + `size` pair.
  ///
  /// This performs an assertion check in debug builds which will terminate the
  /// application if the specified size does not match the extent.
  ///
  /// \pre `p` is a valid pointer.
  /// \pre `size == 0`
  ///
  /// \post `data() == p`
  /// \post `size() == 0`
  ///
  /// \throws Nothing.
  constexpr span(pointer p, [[maybe_unused]] index_type size) : pointer_{p} {
    FUB_ASSERT(size == 0);
  }

  /// Constructs a span from two `pointer`s.
  ///
  /// This performs an assertion check in debug builds which will terminate the
  /// application if the specified size does not match the extent.
  ///
  /// \pre `first` and `last` are valid pointers.
  /// \pre `first == last`
  ///
  /// \post `data() == first`
  /// \post `size() == last - first`
  ///
  /// \throws Nothing.
  constexpr span(pointer first, pointer last) : span(first, last - first) {}

  /// Implicit conversion operator from a mutable container.
  ///
  /// \pre `container.size() == extent`
  /// \pre `!is_span_v<Container>`
  /// \pre `container.size() == 0`
  ///
  /// \post  `data() == container.data()`
  /// \post `size() == container.size()`
  ///
  /// \throws Nothing.
#ifdef DOXYGEN
  template <typename Container>
#else
  template <typename Container,
            typename = std::enable_if_t<!std::is_array<Container>::value>,
            typename = std::enable_if_t<!is_span<Container>::value>,
            typename = std::enable_if_t<
                std::is_convertible<detected_t<detail::array_type_t, Container>,
                                    element_type (*)[]>::value>>
#endif
  constexpr span(Container& container) // NOLINT
      : pointer_{container.data()} {
    FUB_ASSERT(container.size() == extent);
  }

/// Implicit conversion operator from a constant container.
///
/// \pre `extent <= container.size()`
/// \pre `!is_span_v<Container>`
/// \pre `container.size() == 0`
///
/// \post `data() == container.data()`
/// \post `size() == container.size()`
///
/// \throws Nothing.
#ifdef DOXYGEN
  template <typename Container>
#else
  template <typename Container,
            typename = std::enable_if_t<!std::is_array<Container>::value>,
            typename = std::enable_if_t<!is_span<Container>::value>,
            typename = std::enable_if_t<
                std::is_convertible<detected_t<detail::array_type_t, Container>,
                                    element_type (*)[]>::value>>
#endif
  constexpr span(const Container& container) // NOLINT
      : pointer_{container.data()} {
    FUB_ASSERT(container.size() == extent);
  }

  /// Implicit conversion from other span types.
  ///
  /// \pre std::is_convertible<S(*)[], T(*)[]>::value
  ///
  /// \post data() == s.data()
  ///
  /// \throws Nothing.
#ifdef DOXYGEN
  template <typename S>
#else
  template <
      typename S,
      typename = std::enable_if_t<std::is_convertible<
          typename span<S, 0>::element_type (*)[], element_type (*)[]>::value>>
#endif
  constexpr span(const span<S, 0>& s) noexcept // NOLINT
      : pointer_{s.data()} {
  }

  /// Defaulted copy constructor to trivially copy the class member variables.
  ///
  /// \throws Nothing.
  constexpr span(const span& s) noexcept = default;

  /// \name Observers [span.obs]

  /// Returns the number of elements in the span.
  ///
  /// \throws Nothing.
  constexpr index_type size() const noexcept { return 0; }

  /// Returns the number of bytes which are spanned by this span.
  ///
  /// \throws Nothing.
  constexpr std::ptrdiff_t size_bytes() const noexcept { return 0; }

  /// Returns true if  `size() == 0`.
  ///
  /// \throws Nothing.
  constexpr bool empty() const noexcept { return true; }

  /// \name Element access [span.elem]

  /// Returns the underlying pointer.
  ///
  /// \throws Nothing.
  constexpr pointer data() const noexcept { return pointer_; }

  /// \name Iterator support [span.iterators]

  /// Returns an iterator pointing to the first element of the span.
  ///
  /// \throws Nothing.
  constexpr iterator begin() const noexcept { return data(); }

  /// Returns a const iterator pointing to the first element of the span.
  ///
  /// \throws Nothing.
  constexpr const_iterator cbegin() const noexcept { return data(); }

  /// Returns an iterator pointing one after the last element of the span.
  ///
  /// \throws Nothing.
  constexpr iterator end() const noexcept { return begin() + size(); }

  /// Returns a const iterator pointing one after the last element of the span.
  ///
  /// \throws Nothing.
  constexpr const_iterator cend() const noexcept { return begin() + size(); }

  /// Returns a reverse iterator pointing to the last element of the span.
  ///
  /// \throws Nothing.
  constexpr auto rbegin() const noexcept {
    return std::make_reverse_iterator(end());
  }

  /// Returns a const reverse iterator pointing to the last element of the span.
  ///
  /// \throws Nothing.
  constexpr auto crbegin() const noexcept {
    return std::make_reverse_iterator(cend());
  }

  /// Returns a reverse iterator pointing to the first element of the span.
  ///
  /// \throws Nothing.
  constexpr auto rend() const noexcept {
    return std::make_reverse_iterator(begin());
  }

  /// Returns a const reverse iterator pointing to the first element of the
  /// span.
  ///
  /// \throws Nothing.
  constexpr auto crend() const noexcept {
    return std::make_reverse_iterator(cbegin());
  }

private:
  pointer pointer_{nullptr};
};

/// \ingroup spans
/// span is a compact view over a contiguous range of data.
///
/// This class models something like
///
/// \code struct span { pointer ptr; std::ptrdiff_t length; }; \endcode
///
/// Where the `length` can be omitted for statically sized views.
template <typename T> class span<T, dynamic_extent> {
public:
  using element_type = T;
  using pointer = T*;
  using reference = T&;
  using value_type = std::remove_cv_t<T>;
  using index_type = std::ptrdiff_t;
  using difference_type = std::ptrdiff_t;
  using iterator = pointer;
  using const_iterator = std::add_const_t<T>*;
  using reverse_iterator = std::reverse_iterator<iterator>;
  using const_reverse_iterator = std::reverse_iterator<const_iterator>;

  static constexpr index_type extent = dynamic_extent;

  /////////////////////////////////////////////////////////////////////////////
  // Constructors

  /// Constructs an empty span of size 0.
  ///
  /// This initializes the underlying pointer to nullptr and the size to 0.
  ///
  /// \throws Nothing.
  ///
  /// \post data() == nullptr
  /// \post size() == 0
  constexpr span() = default;

  /// Constructs a span from a `pointer` + `size` pair.
  ///
  /// This performs an assertion check in debug builds which will terminate the
  /// application if the specified size does not match the extent.
  ///
  /// \throws Nothing.
  ///
  /// \post data() == p
  /// \post size() == size
  constexpr span(pointer p, index_type size) noexcept
      : pointer_{p}, size_{size} {}

  /// Constructs a span from two `pointer`s.
  ///
  /// This performs an assertion check in debug builds which will terminate the
  /// application if the specified size does not match the extent.
  ///
  /// \throws Nothing.
  ///
  /// \post data() == first
  /// \post size() == last - first
  constexpr span(pointer first, pointer last) noexcept
      : pointer_{first}, size_{last - first} {}

  /// Implicit conversion from a built-in C-style array.
  ///
  /// \throws Nothing.
  ///
  /// \post data() = pointer(&arr[0])
  /// \post size() == M
  template <std::size_t M>
  constexpr span(element_type (&arr)[M]) noexcept // NOLINT
      : pointer_{&arr[0]}, size_{static_cast<index_type>(M)} {}

  /// Implicit conversion from mutable Container-like types.
  ///
  /// \pre S(*)[] is convertible to element_type(*)[]
  ///
  /// \throws Nothing.
  ///
  /// \post data() = container.data()
  /// \post size() = container.size()
  template <
      typename Container,
      typename = std::enable_if_t<!std::is_array<Container>::value>,
      typename =
          std::enable_if_t<!is_span<std::remove_const_t<Container>>::value>,
      typename = std::enable_if_t<
          std::is_convertible<detected_t<detail::array_type_t, Container&>,
                              element_type (*)[]>::value>>
  constexpr span(Container& container) noexcept // NOLINT
      : pointer_{container.data()}, size_{static_cast<index_type>(
                                        container.size())} {}

  /// Implicit conversion from other span types.
  ///
  /// \pre size() <= M
  /// \pre std::is_convertible<S(*)[], T(*)[]>::value
  ///
  /// \post data() == s.data()
  template <
      typename S, index_type OtherExtents,
      typename = std::enable_if_t<std::is_convertible<S (*)[], T (*)[]>::value>>
  constexpr span(const span<S, OtherExtents>& s) noexcept // NOLINT
      : pointer_{s.data()}, size_{s.size()} {}

  constexpr span(const span& s) noexcept = default;

  /// \name Subviews [span.sub]

  /// Returns a span of the first `Count`-many elements of this span.
  ///
  /// \tparam Count  Size of the returned subspan.
  ///
  /// \pre `0 <= Count && Count <= extent`.
  ///
  /// \throws Nothing.
  template <ptrdiff_t Count, typename = std::enable_if_t<0 <= Count>>
  constexpr span<element_type, Count> first() const {
    FUB_ASSERT(Count <= size());
    return {pointer_, Count};
  }

  /// Returns a span of the first `count`-many elements of this span.
  ///
  /// \param[in] count  Size of the returned subspan.
  ///
  /// \pre `0 <= count && count <= extent`.
  ///
  /// \throws Nothing.
  constexpr span<element_type> first(index_type count) const {
    FUB_ASSERT(0 <= count && count <= size());
    return {pointer_, count};
  }

  /// Returns a span of the last `Count`-many elements of this span.
  ///
  /// \tparam Count  Size of the returned subspan.
  ///
  /// \pre `0 <= Count && Count <= extent`.
  ///
  /// \throws Nothing.
  template <ptrdiff_t Count, typename = std::enable_if_t<0 <= Count>>
  constexpr span<element_type, Count> last() const {
    FUB_ASSERT(Count <= size());
    return {pointer_ + (size() - Count), Count};
  }

  /// Returns a span of the last `count`-many elements of this span.
  ///
  /// \param[in] count  Size of the returned subspan.
  ///
  /// \pre `0 <= count && count <= extent`.
  ///
  /// \throws Nothing.
  constexpr span<element_type> last(index_type count) const {
    FUB_ASSERT(0 <= count && count <= size());
    return {pointer_ + (extent - count), count};
  }

  /// Returns a subspan viewing `Count` many elements from offset `Offset`.
  ///
  /// \tparam Offset  the offset of the returned subspan.
  /// \tparam Count  Size of the returned subspan.
  ///
  /// \pre `0 <= Count && Count <= extent`.
  ///
  /// \throws Nothing.
  template <ptrdiff_t Offset, ptrdiff_t Count = dynamic_extent,
            typename = std::enable_if_t<0 <= Offset>>
  constexpr span<element_type, Count> subspan() const {
    FUB_ASSERT(Offset <= size());
    FUB_ASSERT(Count == dynamic_extent ||
               (0 <= Count && Count + Offset <= size()));
    return {pointer_ + Offset,
            Count == dynamic_extent ? size() - Offset : Count};
  }

  /// Returns a subspan viewing `count` many elements from offset `offset`.
  ///
  /// \param offset  the offset of the returned subspan.
  /// \param count  Size of the returned subspan.
  ///
  /// \pre `0 <= offset && offset <= extent`.
  /// \pre `count == dynamic_extent || 0 <= count && count + offset <= size()`.
  ///
  /// \throws Nothing.
  constexpr span<element_type>
  subspan(index_type offset, index_type count = dynamic_extent) const {
    FUB_ASSERT(0 <= offset && offset <= size());
    FUB_ASSERT(count == dynamic_extent ||
               (0 <= count && count + offset <= size()));
    return {pointer_ + offset,
            count == dynamic_extent ? size() - offset : count};
  }

  /// \name Observers [span.obs]

  /// Returns the number of elements in the span.
  ///
  /// \throws Nothing.
  constexpr index_type size() const noexcept { return size_; }

  /// Returns the number of bytes which are spanned by this span.
  ///
  /// \throws Nothing.
  constexpr index_type size_bytes() const noexcept { return sizeof(T) * size_; }

  /// Returns true if  `size() == 0`.
  ///
  /// \throws Nothing.
  constexpr bool empty() const noexcept { return size_ == 0; }

  /// \name Element access [span.elem]

  /// Returns the underlying pointer.
  ///
  /// \throws Nothing.
  constexpr pointer data() const noexcept { return pointer_; }

  /// Accesses the n-th element of the spanned array.
  ///
  /// \throws Nothing.
  constexpr reference operator[](index_type n) const { return pointer_[n]; }

  /// Accesses the n-th element of the spanned array.
  ///
  /// \throws Nothing.
  constexpr reference operator()(index_type n) const { return pointer_[n]; }

  /// \name Iterator support [span.iterators]

  /// Returns an iterator pointing to the first element of the span.
  ///
  /// \throws Nothing.
  constexpr iterator begin() const noexcept { return pointer_; }

  /// Returns a const iterator pointing to the first element of the span.
  ///
  /// \throws Nothing.
  constexpr const_iterator cbegin() const noexcept { return pointer_; }

  /// Returns an iterator pointing one after the last element of the span.
  ///
  /// \throws Nothing.
  constexpr iterator end() const noexcept { return pointer_ + size_; }

  /// Returns a const iterator pointing one after the last element of the span.
  ///
  /// \throws Nothing.
  constexpr const_iterator cend() const noexcept { return pointer_ + size_; }

  /// Returns a reverse iterator pointing to the last element of the span.
  ///
  /// \throws Nothing.
  constexpr reverse_iterator rbegin() const noexcept {
    return std::make_reverse_iterator(end());
  }

  /// Returns a const reverse iterator pointing to the last element of the span.
  ///
  /// \throws Nothing.
  constexpr const_reverse_iterator crbegin() const noexcept {
    return std::make_reverse_iterator(cend());
  }

  /// Returns a reverse iterator pointing to the first element of the span.
  ///
  /// \throws Nothing.
  constexpr reverse_iterator rend() const noexcept {
    return std::make_reverse_iterator(begin());
  }

  /// Returns a const reverse iterator pointing to the first element of the
  /// span.
  ///
  /// \throws Nothing.
  constexpr const_reverse_iterator crend() const noexcept {
    return std::make_reverse_iterator(cbegin());
  }

private:
  pointer pointer_{nullptr};
  index_type size_{0};
};

template <class T, size_t N> span(T (&)[N])->span<T, N>;

template <class T, size_t N> span(std::array<T, N>&)->span<T, N>;

template <class T, size_t N> span(const std::array<T, N>&)->span<const T, N>;

template <class Container>
span(Container&)->span<typename Container::value_type>;

template <class Container>
span(const Container&)->span<const typename Container::value_type>;

} // namespace fub

#endif // !SPAN_HPP
