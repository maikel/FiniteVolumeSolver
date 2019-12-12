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

/// @file This file introduces `extents<E0, ..., En>`, a compact
/// multi-dimensional size type. Each integral extent `Ei` stands an upper bound
/// in dimension `i` and can either be a compile-time constant signed integral
/// value or `dynamic_extent`. Each compile-time sized extent does not take any
/// extra byte.

#ifndef FUB_CORE_MDSPAN_HPP
#define FUB_CORE_MDSPAN_HPP

#include "fub/core/assert.hpp"
#include "fub/core/dynamic_extent.hpp"
#include "fub/core/span.hpp"
#include "fub/core/type_traits.hpp"

#include <array>
#include <type_traits>

#include <boost/mp11/tuple.hpp>

namespace fub {
/// This is the storage type for the extents class and only takes storage for
/// dynamic extents.
///
/// This class is also responsible for accessing all extents.
template <std::size_t Rank, std::size_t DynamicRank,
          std::ptrdiff_t... StaticExtents>
class ExtentsStorage_;

/// This is the class template specialisation if all extents are statically
/// known.
///
/// In this case no extra storage is going to be used.
template <std::size_t Rank, std::ptrdiff_t... StaticExtents>
class ExtentsStorage_<Rank, 0, StaticExtents...> {
public:
  static_assert(Rank == sizeof...(StaticExtents),
                "Rank does not match sizeof...(StaticExtents)!");
  static_assert(
      conjunction<bool_constant<(StaticExtents != dynamic_extent)>...>::value,
      "Some static extent is equal to dynamic_extent!");

  // [mdspan.extents.obs], Observers of the domain multi-index space

  /// Returns sizeof...(StaticExtents)
  static constexpr std::size_t rank() noexcept { return Rank; }

  /// Returns 0.
  static constexpr std::size_t rank_dynamic() noexcept { return 0; }

  /// Returns the n-th Static Extent.
  static constexpr std::ptrdiff_t static_extent(std::size_t n) noexcept {
    const std::ptrdiff_t static_extents[Rank]{StaticExtents...};
    return static_extents[n];
  }

  /// Returns the N-th Extent.
  ///
  /// In this case where all extents are known at compile time, this returns
  /// `static_extent(n)`.
  constexpr std::ptrdiff_t extent(size_t n) const noexcept {
    return static_extent(n);
  }
};

template <std::size_t Rank, std::ptrdiff_t... StaticExtents>
class ExtentsStorage_<Rank, Rank, StaticExtents...> {
public:
  static_assert(Rank == sizeof...(StaticExtents),
                "Rank does not match sizeof...(StaticExtents)!");
  static_assert(
      conjunction<bool_constant<StaticExtents == dynamic_extent>...>::value,
      "Not all static extents are equal to dynamic_extent!");

  // [mdspan.extents.constructors]

  ExtentsStorage_() = default;

  template <typename... IndexType,
            typename = std::enable_if_t<conjunction<
                std::is_convertible<IndexType, std::ptrdiff_t>...>::value>,
            typename = std::enable_if_t<sizeof...(IndexType) == Rank>>
  constexpr explicit ExtentsStorage_(IndexType... extent) noexcept
      : m_dynamic_extents{static_cast<std::ptrdiff_t>(extent)...} {}

  constexpr explicit ExtentsStorage_(
      const std::array<std::ptrdiff_t, Rank>& extents) noexcept
      : m_dynamic_extents{extents} {}

  // [mdspan.extents.obs], Observers of the domain multi-index space

  /// Returns sizeof...(StaticExtents)
  static constexpr std::size_t rank() noexcept { return Rank; }

  /// Returns sizeof...(StaticExtents)
  static constexpr std::size_t rank_dynamic() noexcept { return Rank; }

  /// Returns dynamic_extent.
  static constexpr std::ptrdiff_t static_extent(std::size_t /* n */) noexcept {
    return dynamic_extent;
  }

  /// Returns the N-th Extent.
  constexpr std::ptrdiff_t extent(size_t n) const noexcept {
    return m_dynamic_extents[n];
  }

private:
  std::array<std::ptrdiff_t, Rank> m_dynamic_extents{};
};

template <std::size_t Rank, std::size_t RankDynamic,
          std::ptrdiff_t... StaticExtents>
class ExtentsStorage_ {
public:
  // [mdspan.extents.constructors]

  ExtentsStorage_() = default;

  template <typename... IndexType,
            typename = std::enable_if_t<conjunction<
                std::is_convertible<IndexType, std::ptrdiff_t>...>::value>,
            typename = std::enable_if_t<sizeof...(IndexType) == RankDynamic>>
  constexpr explicit ExtentsStorage_(IndexType... extent) noexcept
      : m_dynamic_extents{static_cast<std::ptrdiff_t>(extent)...} {}

  constexpr explicit ExtentsStorage_(
      const std::array<std::ptrdiff_t, RankDynamic>& extents) noexcept
      : m_dynamic_extents{extents} {}

  // [mdspan.extents.obs], Observers of the domain multi-index space

  /// Returns sizeof...(StaticExtents)
  static constexpr std::size_t rank() noexcept { return Rank; }

  /// Returns sizeof...(StaticExtents)
  static constexpr std::size_t rank_dynamic() noexcept { return RankDynamic; }

  /// Returns dynamic_extent.
  static constexpr std::ptrdiff_t static_extent(std::size_t n) noexcept {
    constexpr std::ptrdiff_t static_extents[Rank]{StaticExtents...};
    return static_extents[n];
  }

  /// Returns the N-th Extent.
  constexpr std::ptrdiff_t extent(size_t n) const noexcept {
    constexpr std::ptrdiff_t static_extents[Rank]{StaticExtents...};
    if (static_extents[n] != dynamic_extent) {
      return static_extents[n];
    }
    return m_dynamic_extents[find_dynamic_extent_index(n)];
  }

private:
  std::array<std::ptrdiff_t, RankDynamic> m_dynamic_extents{};

  static constexpr std::size_t
  find_dynamic_extent_index(std::size_t n) noexcept {
    constexpr std::ptrdiff_t static_extents[Rank]{StaticExtents...};
    std::size_t dynamic_extent_index = 0;
    for (std::size_t dim = 0; dim < n; ++dim) {
      if (static_extents[dim] == dynamic_extent) {
        ++dynamic_extent_index;
      }
    }
    return dynamic_extent_index;
  }
};

template <typename... IndexType>
constexpr std::size_t count_dynamic_extents(IndexType... extent) noexcept {
  constexpr std::size_t Rank = sizeof...(IndexType);
  std::ptrdiff_t extents[Rank]{static_cast<std::ptrdiff_t>(extent)...};
  std::size_t counter = 0;
  for (std::size_t dim = 0; dim < Rank; ++dim) {
    if (extents[dim] == dynamic_extent) {
      ++counter;
    }
  }
  return counter;
}

/// An extents object defines a multidimensional index space which is the
/// Cartesian product of integers extents `[0..N0) * [0..N1) * ...`
template <std::ptrdiff_t... StaticExtents>
class extents : private ExtentsStorage_<sizeof...(StaticExtents),
                                        count_dynamic_extents(StaticExtents...),
                                        StaticExtents...> {
private:
  using base_type = ExtentsStorage_<sizeof...(StaticExtents),
                                    count_dynamic_extents(StaticExtents...),
                                    StaticExtents...>;

public:
  // Type Definitions

  using index_type = std::ptrdiff_t;

  // [mdspan.extents.constructors]

  using base_type::base_type;

  // [mdspan.extents.static_observers]

  /// Returns sizeof...(StaticExtents)
  static constexpr std::size_t rank() noexcept { return base_type::rank(); }

  /// Returns sizeof...(StaticExtents)
  static constexpr std::size_t rank_dynamic() noexcept {
    return base_type::rank_dynamic();
  }

  /// Returns the `n`-th `StaticExtent`.
  static constexpr std::ptrdiff_t static_extent(std::size_t n) noexcept {
    return base_type::static_extent(n);
  }

  // [mdspan.extents.observers]

  /// Returns the `n`-th run-time extent.
  constexpr std::ptrdiff_t extent(size_t n) const noexcept {
    return base_type::extent(n);
  }
};

// [mdspan.extents.traits.is_extent]

/// \ingroup type-traits
/// This is true `std::true_type` iff `E` is `extents<Es...>`
/// for some `std::ptrdiff_t... Es`.
/// @{
template <typename E> struct is_extents : std::false_type {};
template <std::ptrdiff_t... StaticExtents>
struct is_extents<extents<StaticExtents...>> : std::true_type {};
template <typename E> static constexpr bool is_extents_v = is_extents<E>::value;
/// @}

/// Returns: `true` if `lhs.rank() == rhs.rank()` and `lhs.extents(r) ==
/// rhs.extents(r)` for all `r` in the range `[0, lhs.rank())`, or false
/// otherwise.
template <std::ptrdiff_t... StaticExtentsL, std::ptrdiff_t... StaticExtentsR>
constexpr bool operator==(const extents<StaticExtentsL...>& lhs,
                          const extents<StaticExtentsR...>& rhs) noexcept {
  if (lhs.rank() != rhs.rank()) {
    return false;
  }
  for (std::size_t r = 0; r < lhs.rank(); ++r) {
    if (lhs.extent(r) != rhs.extent(r)) {
      return false;
    }
  }
  return true;
}

template <std::ptrdiff_t... StaticExtentsL, std::ptrdiff_t... StaticExtentsR>
constexpr bool operator!=(const extents<StaticExtentsL...>& left,
                          const extents<StaticExtentsR...>& right) noexcept {
  return !(left == right);
}

template <std::ptrdiff_t... StaticExtents>
constexpr std::ptrdiff_t Size_(const extents<StaticExtents...> e) noexcept {
  const std::size_t rank = e.rank();
  std::ptrdiff_t size = 1;
  for (std::size_t r = 0; r < rank; ++r) {
    size *= e.extent(r);
  }
  return size;
}

/// This layout creates mappings which do row first indexing (as in Fortran).
///
/// It holds for all valid `i` and `j` the equation
/// \code mapping(i + 1, j) - mapping(1, j) == 1 \endcode
struct layout_left {

  /// This mapping does row first indexing (as in Fortran).
  ///
  /// It holds for all valid `i` and `j` the equation
  /// \code mapping(i + 1, j) - mapping(1, j) == 1 \endcode
  template <typename Extents> class mapping : private Extents {
  public:
    static_assert(is_extents_v<Extents>,
                  "Extents must satisfy the extents concept.");

    /// \name Constructors / Destructors
    /// @{
    constexpr mapping() = default;
    constexpr mapping(const mapping&) = default;
    constexpr mapping(mapping&&) = default;
    mapping& operator=(const mapping& other) noexcept = default;
    mapping& operator=(mapping&& other) noexcept = default;

    /// Implicit Conversion from Extents.
    ///
    /// \throws Nothing.
    ///
    /// \post extents() == extents.
    constexpr mapping(const Extents& extents) : Extents{extents} {} // NOLINT
    /// @}

    /// \name Observers
    /// @{
    constexpr const Extents& extents() const noexcept { return *this; }

    constexpr std::ptrdiff_t required_span_size() const noexcept {
      return Size_(extents());
    }

    /// \return Returns the codomain index for specified multi dimensional index
    /// coordinates.
    ///
    /// \effect Equivalent to offset in
    /// \code
    ///    index_type offset = 0;
    ///    for(size_t k = 0; k < Extents::rank(); ++k) {
    ///      offset += i[k] * stride(k);
    ///    }
    /// \endcode
    /// \throws Nothing.
    template <
        typename... IndexType,
        typename = std::enable_if_t<conjunction<
            std::is_convertible<IndexType, std::ptrdiff_t>...>::value>,
        typename = std::enable_if_t<(sizeof...(IndexType) == Extents::rank())>>
    constexpr std::ptrdiff_t operator()(IndexType... indices) const noexcept {
      return DoMapping_(make_index_sequence<Extents::rank()>(), indices...);
    }

    /// \return Returns the stride value in a specified dimension r.
    ///
    /// \throws Nothing.
    constexpr std::ptrdiff_t stride(std::size_t r) const noexcept {
      std::ptrdiff_t stride = 1;
      for (std::size_t dim = 0; dim < r; ++dim) {
        stride *= extents().extent(dim);
      }
      return stride;
    }

    static constexpr bool is_always_unique() noexcept { return true; }
    static constexpr bool is_always_contiguous() noexcept { return true; }
    static constexpr bool is_always_strided() noexcept { return true; }

    constexpr bool is_unique() const noexcept { return true; }
    constexpr bool is_contiguous() const noexcept { return true; }
    constexpr bool is_strided() const noexcept { return true; }

    template <class OtherExtents>
    constexpr bool operator==(const mapping<OtherExtents>& other) const
        noexcept {
      return extents() == other.extents();
    }

    template <class OtherExtents>
    constexpr bool operator!=(const mapping<OtherExtents>& other) const
        noexcept {
      return !(*this == other);
    }
    /// @}

  private:
    template <std::size_t... Is, typename... IndexType>
    constexpr std::ptrdiff_t DoMapping_(index_sequence<Is...>,
                                        IndexType... indices) const noexcept {
      const std::ptrdiff_t is[Extents::rank()]{
          static_cast<std::ptrdiff_t>(indices)...};
      const std::ptrdiff_t strides[Extents::rank()]{stride(Is)...};
      std::ptrdiff_t index = 0;
      for (std::size_t r = 0; r < Extents::rank(); ++r) {
        index += strides[r] * is[r];
      }
      return index;
    }
  };
};

struct layout_right {
  template <typename Extents> class mapping : private Extents {
  public:
    constexpr mapping() = default;
    constexpr mapping(const mapping&) = default;
    constexpr mapping(mapping&&) = default;
    mapping& operator=(const mapping&) noexcept = default;
    mapping& operator=(mapping&&) noexcept = default;

    constexpr mapping(const Extents& extents) : Extents{extents} {}

    constexpr const Extents& extents() const noexcept { return *this; }

    constexpr std::ptrdiff_t required_span_size() const noexcept {
      return Size_(extents());
    }

    template <
        typename... IndexType,
        typename = std::enable_if_t<conjunction<
            std::is_convertible<IndexType, std::ptrdiff_t>...>::value>,
        typename = std::enable_if_t<(sizeof...(IndexType) == Extents::rank())>>
    constexpr std::ptrdiff_t operator()(IndexType... indices) const noexcept {
      return DoMapping_(make_index_sequence<Extents::rank()>(), indices...);
    }

    constexpr std::ptrdiff_t stride(std::size_t r) const noexcept {
      std::ptrdiff_t stride = 1;
      for (std::size_t dim = r + 1; dim < Extents::rank(); ++dim) {
        stride *= extents().extent(dim);
      }
      return stride;
    }

    static constexpr bool is_always_unique() noexcept { return true; }
    static constexpr bool is_always_contiguous() noexcept { return true; }
    static constexpr bool is_always_strided() noexcept { return true; }

    constexpr bool is_unique() const noexcept { return true; }
    constexpr bool is_contiguous() const noexcept { return true; }
    constexpr bool is_strided() const noexcept { return true; }

    template <class OtherExtents>
    constexpr bool operator==(const mapping<OtherExtents>& other) const
        noexcept {
      return extents() == other.extents();
    }

    template <class OtherExtents>
    constexpr bool operator!=(const mapping<OtherExtents>& other) const
        noexcept {
      return !(*this == other);
    }

  private:
    template <std::size_t... Is, typename... IndexType>
    constexpr std::ptrdiff_t DoMapping_(index_sequence<Is...>,
                                        IndexType... indices) const noexcept {
      const std::ptrdiff_t is[Extents::rank()]{
          static_cast<std::ptrdiff_t>(indices)...};
      const std::ptrdiff_t strides[Extents::rank()]{stride(Is)...};
      std::ptrdiff_t index = 0;
      for (std::size_t r = 0; r < Extents::rank(); ++r) {
        index += strides[r] * is[r];
      }
      return index;
    }
  };
};

struct layout_stride {
  template <class Extents> class mapping : private Extents {
  public:
    // [mdspan.layout.stride.cons], layout_stride::mapping constructors
    constexpr mapping() = default;

    constexpr mapping(
        const Extents& e,
        const std::array<std::ptrdiff_t, Extents::rank()>& s) noexcept
        : Extents(e), strides_(s) {}

    std::array<std::ptrdiff_t, Extents::rank()>
    MakeStrides(const layout_left::mapping<Extents>& other) {
      std::array<std::ptrdiff_t, Extents::rank()> strides;
      for (std::size_t s = 0; s < Extents::rank(); ++s) {
        strides[s] = other.stride(s);
      }
      return strides;
    }

    constexpr mapping(const layout_left::mapping<Extents>& other)
        : mapping(other.extents(), MakeStrides(other)) {}

    // [mdspan.layout.stride.ops], layout_stride::mapping operations
    constexpr const Extents& extents() const noexcept { return *this; }

    constexpr const std::array<std::ptrdiff_t, Extents::rank()>& strides() const
        noexcept {
      return strides_;
    }

    constexpr std::ptrdiff_t required_span_size() const noexcept {
      std::ptrdiff_t max = extents().extent(0) * strides_[0];
      for (std::size_t r = 1; r < Extents::rank(); ++r) {
        max = std::max(extents().extent(r) * strides_[r], max);
      }
      return max;
    }

    template <
        typename... IndexType,
        typename = std::enable_if_t<conjunction<
            std::is_convertible<IndexType, std::ptrdiff_t>...>::value>,
        typename = std::enable_if_t<(sizeof...(IndexType) == Extents::rank())>>
    std::ptrdiff_t operator()(IndexType... is) const {
      const std::array<std::ptrdiff_t, Extents::rank()> idx{
          static_cast<std::ptrdiff_t>(is)...};
      std::ptrdiff_t offset = 0;
      for (std::size_t r = 0; r < Extents::rank(); ++r) {
        offset += idx[r] * strides_[r];
      }
      return offset;
    }

    static constexpr bool is_always_unique() noexcept { return true; }
    static constexpr bool is_always_contiguous() noexcept { return false; }
    static constexpr bool is_always_strided() noexcept { return true; }

    constexpr bool is_unique() const noexcept { return true; }
    constexpr bool is_contiguous() const noexcept { return false; }
    constexpr bool is_strided() const noexcept { return true; }

    std::ptrdiff_t stride(size_t rank) const noexcept { return strides_[rank]; }

    template <class OtherExtents>
    constexpr bool operator==(const mapping<OtherExtents>& other) const
        noexcept {
      return extents() == other.extents();
    }
    template <class OtherExtents>
    constexpr bool operator!=(const mapping<OtherExtents>& other) const
        noexcept {
      return !(*this == other);
    }

  private:
    std::array<std::ptrdiff_t, Extents::rank()> strides_;
  };
};

template <class ElementType> struct accessor_basic {
  using offset_policy = accessor_basic;
  using element_type = ElementType;
  using reference = ElementType&;
  using pointer = ElementType*;

  constexpr pointer offset(pointer p, ptrdiff_t i) const noexcept {
    return p + i;
  }

  constexpr reference access(pointer p, ptrdiff_t i) const noexcept {
    return *offset(p, i);
  }

  constexpr pointer decay(pointer p) const noexcept { return p; }
};

template <class ElementType, class Extents, class LayoutPolicy = layout_right,
          class AccessorPolicy = accessor_basic<ElementType>>
class basic_mdspan : AccessorPolicy, LayoutPolicy::template mapping<Extents> {
public:
  // Domain and codomain types
  using extents_type = Extents;
  using layout_type = LayoutPolicy;
  using accessor_type = AccessorPolicy;
  using mapping_type = typename layout_type::template mapping<extents_type>;
  using element_type = typename accessor_type::element_type;
  using value_type = std::remove_cv_t<element_type>;
  using index_type = std::ptrdiff_t;
  using difference_type = std::ptrdiff_t;
  using pointer = typename accessor_type::pointer;
  using reference = typename accessor_type::reference;

  // [mdspan.basic.cons], basic_mdspan constructors, assignment, and destructor
  constexpr basic_mdspan() noexcept = default;
  constexpr basic_mdspan(const basic_mdspan&) noexcept = default;
  constexpr basic_mdspan(basic_mdspan&&) noexcept = default;

  template <class... IndexType,
            typename = std::enable_if_t<conjunction<
                std::is_convertible<IndexType, std::ptrdiff_t>...>::value>,
            typename = std::enable_if_t<sizeof...(IndexType) ==
                                        Extents::rank_dynamic()>>
  explicit constexpr basic_mdspan(pointer p, IndexType... dynamic_extents)
      : basic_mdspan(p, Extents(dynamic_extents...)) {}

  template <class IndexType, std::size_t N>
  explicit constexpr basic_mdspan(
      pointer p, const std::array<IndexType, N>& dynamic_extents)
      : basic_mdspan(p, extents_type{dynamic_extents}) {}

  constexpr basic_mdspan(pointer p, const mapping_type& m)
      : accessor_type(), mapping_type(m), ptr_(p) {}

  constexpr basic_mdspan(pointer p, const mapping_type& m,
                         const accessor_type& a)
      : accessor_type(a), mapping_type(m), ptr_(p) {}

  template <class OtherElementType, class OtherExtents, class OtherLayoutPolicy,
            class OtherAccessorPolicy>
  constexpr basic_mdspan(
      const basic_mdspan<OtherElementType, OtherExtents, OtherLayoutPolicy,
                         OtherAccessorPolicy>& other)
      : accessor_type(), mapping_type(other.mapping()), ptr_(other.data()) {}

  ~basic_mdspan() = default;

  constexpr basic_mdspan& operator=(const basic_mdspan&) noexcept = default;
  constexpr basic_mdspan& operator=(basic_mdspan&&) noexcept = default;

  template <class OtherElementType, class OtherExtents, class OtherLayoutPolicy,
            class OtherAccessorPolicy>
  constexpr basic_mdspan& operator=(
      const basic_mdspan<OtherElementType, OtherExtents, OtherLayoutPolicy,
                         OtherAccessorPolicy>& other) noexcept {
    static_cast<mapping_type&>(*this) = other.mapping();
    ptr_ = other.data();
    return *this;
  }

  // [mdspan.basic.mapping], basic_mdspan mapping domain multi-index to access
  // codomain element
  constexpr reference operator[](index_type i) const noexcept {
    return this->operator()(i);
  }

  template <class... IndexType>
  constexpr reference operator()(IndexType... indices) const noexcept {
    const accessor_type acc = accessor();
    const mapping_type map = mapping();
    return acc.access(ptr_, map(indices...));
  }

  template <class IndexType, size_t N>
  constexpr reference operator()(const std::array<IndexType, N>& indices) const
      noexcept {
    const accessor_type acc = accessor();
    const mapping_type map = mapping();
    return acc.access(ptr_, boost::mp11::tuple_apply(map, indices));
  }

  // [mdspan.basic.domobs], basic_mdspan observers of the domain
  // multidimensional index space

  static constexpr int rank() noexcept { return Extents::rank(); }

  static constexpr int rank_dynamic() noexcept {
    return Extents::rank_dynamic();
  }

  static constexpr index_type static_extent(std::size_t rank) noexcept {
    return Extents::static_extent(rank);
  }

  constexpr Extents extents() const noexcept { return mapping().extents(); }

  constexpr index_type extent(std::size_t rank) const noexcept {
    return extents().extent(rank);
  }

  constexpr index_type size() const noexcept { return Size_(extents()); }

  constexpr index_type unique_size() const noexcept { return size(); }

  // [mdspan.basic.codomain], basic_mdspan observers of the codomain
  constexpr span<element_type> get_span() const noexcept {
    return {accessor().decay(ptr_), mapping().required_span_size()};
  }

  constexpr pointer data() const noexcept { return ptr_; }

  // [mdspan.basic.obs], basic_mdspan observers of the mapping
  static constexpr bool is_always_unique() noexcept {
    return mapping_type::is_always_unique();
  }
  static constexpr bool is_always_contiguous() noexcept {
    return mapping_type::is_always_contiguous();
  }
  static constexpr bool is_always_strided() noexcept {
    return mapping_type::is_always_strided();
  }

  constexpr mapping_type mapping() const noexcept { return *this; }
  constexpr accessor_type accessor() const noexcept { return *this; }
  constexpr bool is_unique() const noexcept { return mapping().is_unique(); }
  constexpr bool is_contiguous() const noexcept {
    return mapping().is_contiguous();
  }
  constexpr bool is_strided() const noexcept { return mapping().is_strided(); }

  constexpr index_type stride(std::size_t r) const {
    return mapping().stride(r);
  }

private:
  pointer ptr_{};
};

template <typename T, ptrdiff_t... Extents>
using static_mdspan = basic_mdspan<T, extents<Extents...>>;

template <typename T> struct dynamic_extents_;

template <typename I, I... Is>
struct dynamic_extents_<integer_sequence<I, Is...>> {
  using type = extents<((void)Is, fub::dynamic_extent)...>;
};

template <std::size_t Rank>
using dynamic_extents =
    typename dynamic_extents_<make_index_sequence<Rank>>::type;

// template <typename T, std::size_t Rank, typename Layout = layout_right>
// using DynamicMdSpan = basic_mdspan<T, dynamic_extents<Rank>, Layout>;

template <typename T, std::size_t Rank, typename Layout = layout_left>
using mdspan = basic_mdspan<T, dynamic_extents<Rank>, Layout>;

template <typename T, typename E, typename A = accessor_basic<T>>
using strided_mdspan = basic_mdspan<T, E, layout_stride, A>;

struct all_type {
  explicit all_type() = default;
};
static constexpr all_type all = all_type{};

/// \cond INTERNAL
template <typename T>
struct IsSliceRange_
    : disjunction<
          std::is_convertible<T, std::pair<std::ptrdiff_t, std::ptrdiff_t>>,
          std::is_convertible<T, all_type>> {};

template <typename I, typename... Is> constexpr I Accumulate_(I x0, Is... xi) {
  const I x[] = {x0, xi...};
  constexpr int size = 1 + static_cast<int>(sizeof...(Is));
  I total = I{};
  for (int i = 0; i < size; ++i) {
    total += x[i];
  }
  return total;
}

template <typename... SliceSpecifiers>
static constexpr std::size_t SliceRank_ =
    Accumulate_(static_cast<int>(IsSliceRange_<SliceSpecifiers>::value)...);
/// \endcond

// [mdspan.subspan], subspan creation
template <class ElementType, class Extents, class LayoutPolicy,
          class AccessorPolicy, class... SliceSpecifiers>
struct mdspan_subspan { // exposition only
  static constexpr std::size_t slice_rank = SliceRank_<SliceSpecifiers...>;
  using extents_t = dynamic_extents<slice_rank>;
  using layout_t = layout_stride;
  using type = basic_mdspan<ElementType, extents_t, layout_t,
                            typename AccessorPolicy::offset_policy>;
};

template <class ElementType, class Extents, class LayoutPolicy,
          class AccessorPolicy, class... SliceSpecifiers>
using mdspan_subspan_t = // exposition only
    typename mdspan_subspan<ElementType, Extents, LayoutPolicy, AccessorPolicy,
                            SliceSpecifiers...>::type;

constexpr std::ptrdiff_t
SliceExtent_(const std::pair<std::ptrdiff_t, std::ptrdiff_t>& p) {
  FUB_ASSERT(p.first <= p.second);
  FUB_ASSERT(0 <= p.first);
  return p.second - p.first;
}

/// \cond INTERNAL
constexpr std::ptrdiff_t SliceExtent_(all_type) { return dynamic_extent; }

template <typename T, typename = std::enable_if_t<!IsSliceRange_<T>::value>>
constexpr std::ptrdiff_t SliceExtent_(T) {
  return 0;
}

template <std::size_t Rank, typename Slice>
constexpr std::array<std::ptrdiff_t, Rank + 1>
SlicePushBack__(std::true_type, std::array<std::ptrdiff_t, Rank> array,
                Slice slice) {
  std::array<std::ptrdiff_t, Rank + 1> grown{};
  for (std::size_t i = 0; i < Rank; ++i) {
    grown[i] = array[i];
  }
  grown[Rank] = SliceExtent_(slice);
  return grown;
}

template <std::size_t Rank, typename Slice>
constexpr std::array<std::ptrdiff_t, Rank>
SlicePushBack__(std::false_type, std::array<std::ptrdiff_t, Rank> array,
                Slice /* slice */) {
  return array;
}

template <std::size_t Rank, typename Slice>
constexpr auto SlicePushBack_(std::array<std::ptrdiff_t, Rank> array,
                              Slice slice) {
  return SlicePushBack__(IsSliceRange_<Slice>(), array, slice);
}

template <std::size_t Rank>
constexpr auto MakeSliceExtents__(std::array<std::ptrdiff_t, Rank> result) {
  return result;
}

template <std::size_t Rank, typename Slice, typename... Rest>
constexpr auto MakeSliceExtents__(std::array<std::ptrdiff_t, Rank> prev,
                                  Slice slice, Rest... slices) {
  return MakeSliceExtents__(SlicePushBack_(prev, slice), slices...);
}

template <typename... Slices>
constexpr auto MakeSliceExtents_(Slices... slices) {
  return MakeSliceExtents__(std::array<std::ptrdiff_t, 0>{}, slices...);
}

template <class ElementType, class Extents, class LayoutPolicy,
          class AccessorPolicy, std::size_t... Is>
constexpr std::array<std::ptrdiff_t, sizeof...(Is)> MakeStrideArray_(
    const basic_mdspan<ElementType, Extents, LayoutPolicy, AccessorPolicy>& src,
    index_sequence<Is...>) {
  return {src.stride(Is)...};
}

template <typename T> struct Type_ {};

template <std::size_t N, std::size_t... Is>
constexpr std::array<std::size_t, N> AsStdArray__(const std::size_t (&array)[N],
                                                  index_sequence<Is...>) {
  return {array[Is]...};
}

template <std::size_t N>
constexpr std::array<std::size_t, N>
AsStdArray_(const std::size_t (&array)[N]) {
  return AsStdArray__(array, make_index_sequence<N>());
}

template <typename... Slices>
constexpr std::array<std::size_t, SliceRank_<Slices...>>
MakeRangeIndices__(Type_<Slices>... /* slices */) {
  constexpr std::size_t is_range[] = {IsSliceRange_<Slices>::value...};
  std::size_t ranges[SliceRank_<Slices...>]{};
  std::size_t ri = 0;
  for (std::size_t i = 0; i < sizeof...(Slices); ++i) {
    if (is_range[i]) {
      ranges[ri++] = i;
    }
  }
  return AsStdArray_(ranges);
}

template <typename... Slices, std::size_t... Is>
constexpr auto MakeRangeIndices_(index_sequence<Is...>,
                                 Type_<Slices>... /* slices */) {
  constexpr std::array<std::size_t, SliceRank_<Slices...>> indices{
      MakeRangeIndices__(Type_<Slices>{}...)};
  return index_sequence<indices[Is]...>{};
}

template <typename... Slices>
constexpr auto MakeRangeIndices_(Type_<Slices>... slices) {
  return MakeRangeIndices_(make_index_sequence<SliceRank_<Slices...>>{},
                           slices...);
}

template <typename... Slices> constexpr auto MakeRangeIndices(Slices...) {
  return MakeRangeIndices_(Type_<Slices>{}...);
}

template <class ElementType, class Extents, class LayoutPolicy,
          class AccessorPolicy, typename... Slices>
constexpr auto MakeSliceArray_(
    const basic_mdspan<ElementType, Extents, LayoutPolicy, AccessorPolicy>& src,
    Slices...) {
  constexpr auto range_indices = MakeRangeIndices_(Type_<Slices>{}...);
  return MakeStrideArray_(src, range_indices);
}

constexpr std::ptrdiff_t MakeOrigin__(std::ptrdiff_t n) { return n; }

constexpr std::ptrdiff_t
MakeOrigin__(std::pair<std::ptrdiff_t, std::ptrdiff_t> p) {
  return p.first;
}

constexpr std::ptrdiff_t MakeOrigin__(all_type /* all */) { return 0; }

template <typename... Slices>
constexpr std::array<std::ptrdiff_t, sizeof...(Slices)>
MakeOrigin_(Slices... slices) {
  return {MakeOrigin__(slices)...};
}

template <std::size_t... Is, typename F, std::size_t N>
constexpr auto Apply_(index_sequence<Is...>, F fn,
                      std::array<std::ptrdiff_t, N> indices) {
  return fn(indices[Is]...);
}

template <typename F, std::size_t N>
constexpr auto Apply_(F fn, std::array<std::ptrdiff_t, N> indices) {
  return Apply_(make_index_sequence<N>{}, fn, indices);
}

constexpr std::ptrdiff_t MapToLength_(all_type) { return -1; }

constexpr std::ptrdiff_t
MapToLength_(const std::pair<std::ptrdiff_t, std::ptrdiff_t>& rng) {
  return rng.second - rng.first;
}

template <typename T, typename = std::enable_if_t<!IsSliceRange_<T>()>>
constexpr std::ptrdiff_t MapToLength_(T) {
  return 0;
}

template <typename Extents, typename... SliceSpecifiers>
constexpr std::array<std::ptrdiff_t, SliceRank_<SliceSpecifiers...>>
MakeSliceExtents2_(const Extents& e, SliceSpecifiers... slices) {
  const std::array<std::ptrdiff_t, Extents::rank()> length{
      MapToLength_(slices)...};
  std::array<std::ptrdiff_t, Extents::rank()> dyn_length{};
  for (std::size_t i = 0; i < Extents::rank(); ++i) {
    if (length[i] < 0) {
      dyn_length[i] = e.extent(i);
    } else {
      dyn_length[i] = length[i];
    }
  }
  std::array<std::ptrdiff_t, SliceRank_<SliceSpecifiers...>> extents{};
  std::size_t slice = 0;
  for (std::size_t i = 0; i < Extents::rank(); ++i) {
    if (dyn_length[i] > 0) {
      extents[slice] = dyn_length[i];
      slice = slice + 1;
    }
  }
  return extents;
}
/// \endcond

template <class ElementType, class Extents, class LayoutPolicy,
          class AccessorPolicy, class... SliceSpecifiers>
mdspan_subspan_t<ElementType, Extents, LayoutPolicy, AccessorPolicy,
                 SliceSpecifiers...>
subspan(
    const basic_mdspan<ElementType, Extents, LayoutPolicy, AccessorPolicy>& src,
    SliceSpecifiers... slices) noexcept {
  const std::array<std::ptrdiff_t, Extents::rank()> origin =
      MakeOrigin_(slices...);
  const auto map = src.mapping();
  const auto acc = src.accessor();
  constexpr std::size_t slice_rank = SliceRank_<SliceSpecifiers...>;
  using SliceExtents = dynamic_extents<slice_rank>;
  const SliceExtents extents{MakeSliceExtents2_(src.extents(), slices...)};
  const std::array<std::ptrdiff_t, slice_rank> slice_array =
      MakeSliceArray_(src, slices...);
  const layout_stride::mapping<SliceExtents> slice_mapping(extents,
                                                           slice_array);
  return {acc.offset(src.data(), Apply_(map, origin)), slice_mapping};
}

template <typename> struct is_mdspan : std::false_type {};

template <typename T, typename E, typename L, typename A>
struct is_mdspan<basic_mdspan<T, E, L, A>> : std::true_type {};

template <typename T> static constexpr bool is_mdspan_v = is_mdspan<T>::value;

} // namespace fub

#endif
