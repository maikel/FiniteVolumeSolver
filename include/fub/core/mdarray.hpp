// Copyright (c) 2018-2020 Maikel Nadolski
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

#ifndef FUB_CORE_MDARRAY_HPP
#define FUB_CORE_MDARRAY_HPP

#include "fub/core/mdspan.hpp"
#include <vector>

namespace fub {

template <typename T, typename Allocator = std::allocator<T>>
struct std_vector_policy {
  using element_type = T;
  using container_type = std::vector<T, Allocator>;
  using allocator_type = typename container_type::allocator_type;
  using pointer = typename container_type::pointer;
  using const_pointer = typename container_type::const_pointer;
  using reference = typename container_type::reference;
  using const_reference = typename container_type::const_reference;
  using offset_policy = std_vector_policy<T, Allocator>;

  reference access(container_type& c, ptrdiff_t i);
  const_reference access(const container_type& c, ptrdiff_t i);

  container_type create(const T& val, std::size_t n,
                        const allocator_type& alloc = allocator_type()) {
    return container_type(n, val, alloc);
  }

  container_type create(std::size_t n,
                        const allocator_type& alloc = allocator_type()) {
    return container_type(n, alloc);
  }

  reference access(pointer p, std::ptrdiff_t i);
  const_reference access(const_pointer p, std::ptrdiff_t i);

  pointer offset(pointer p, std::ptrdiff_t i);
  const_pointer offset(const_pointer p, std::ptrdiff_t i);

  element_type* decay(pointer p);
  const element_type* decay(const_pointer p);
};

/// \brief This is a multi-dimensional array that is based on the interface of
/// mdspan.
template <typename ValueType, typename Extents, typename LayoutPolicy,
          typename ContainerPolicy>
class basic_mdarray {
public:
  using value_type = ValueType;
  using extents_type = Extents;
  using layout_policy = LayoutPolicy;
  using mapping_type = typename layout_policy::template mapping<Extents>;
  using container_policy = ContainerPolicy;
  using container_type = typename container_policy::container_type;
  using allocator_type = typename container_policy::allocator_type;
  using pointer = typename container_type::pointer;
  using const_pointer = typename container_type::const_pointer;
  using reference = typename container_type::reference;
  using const_reference = typename container_type::const_reference;
  using view_type =
      basic_mdspan<value_type, extents_type, layout_policy, container_policy>;
  using const_view_type = basic_mdspan<const value_type, extents_type,
                                       layout_policy, container_policy>;

  constexpr basic_mdarray() noexcept;

  constexpr basic_mdarray(const value_type& default_value,
                          const extents_type& extents,
                          const allocator_type& allocator = allocator_type());

  constexpr basic_mdarray(const extents_type& extents,
                          const allocator_type& allocator = allocator_type());

  constexpr basic_mdarray(const mapping_type& mapping,
                          const container_policy& cp = container_policy());

  constexpr view_type view() noexcept;
  constexpr const_view_type view() const noexcept;

  static constexpr int rank() noexcept { return extents_type::rank(); }
  static constexpr int dynamic_rank() noexcept {
    return extents_type::dynamic_rank();
  }
  static constexpr std::ptrdiff_t static_extent(int d) noexcept {
    return extents_type::static_extent(d);
  }

  constexpr extents_type extents() const noexcept;

  constexpr std::ptrdiff_t extent(int d) const noexcept;
  constexpr std::ptrdiff_t size(int d) const noexcept;

  constexpr reference operator()(std::array<std::ptrdiff_t, rank()> i);
  constexpr const_reference
  operator()(std::array<std::ptrdiff_t, rank()> i) const;

  template <typename... Is>
  constexpr std::enable_if_t<sizeof...(Is) == rank() &&
                                 (true && ... && std::is_integral_v<Is>),
                             reference>
  operator()(Is... is);

  template <typename... Is>
  constexpr std::enable_if_t<sizeof...(Is) == rank() &&
                                 (true && ... && std::is_integral_v<Is>),
                             const_reference>
  operator()(Is... is) const;

  constexpr const container_type& container() const noexcept;
  constexpr container_type& container() noexcept;

  constexpr std::size_t size() const noexcept;

private:
  [[no_unique_address]] mapping_type mapping_;
  [[no_unique_address]] container_policy container_policy_;
  container_type container_;
};

template <typename T, std::ptrdiff_t... Is>
using mdarray =
    basic_mdarray<T, extents<Is...>, layout_left, std_vector_policy<T>>;

template <typename T, std::ptrdiff_t Rank>
using dynamic_mdarray =
    basic_mdarray<T, dynamic_extents<Rank>, layout_left, std_vector_policy<T>>;

///////////////////////////////////////////////////////////////////////////////
//                                                              Implementation

template <typename T, typename Extents, typename LayoutPolicy,
          typename ContainerPolicy>
constexpr basic_mdarray<T, Extents, LayoutPolicy, ContainerPolicy>::
    basic_mdarray(const extents_type& extents, const allocator_type& allocator)
    : mapping_(extents), container_(ContainerPolicy().create(
                             mapping_.required_span_size(), allocator)) {}

template <typename T, typename Extents, typename LayoutPolicy,
          typename ContainerPolicy>
constexpr basic_mdarray<T, Extents, LayoutPolicy, ContainerPolicy>::
    basic_mdarray(const value_type& default_value, const extents_type& extents,
                  const allocator_type& allocator)
    : mapping_(extents),
      container_(ContainerPolicy().create(
          default_value, mapping_.required_span_size(), allocator)) {}

} // namespace fub

#endif