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

#include "fub/core/mdspan.hpp"

namespace fub {

template <typename ValueType, typename Extents, typename LayoutPolicy,
          typename ContainerPolicy>
class basic_mdarray {
public:
  using value_type = ValueType;
  using extents_type = Extents;
  using layout_policy = LayoutPolicy;
  using mapping_type = typename layout_policy::template mapping<Extents>;
  using container_policy = ContainerPolicy;
  using conatiner_type =
      typename container_policy::template container<value_type>;

  basic_mdarray(const value_type& default_value, const extents_type& extents,
                const allocator_type& allocator = allocator_type());

  basic_mdarray(const extents_type& extents,
                const allocator_type& allocator = allocator_type());

  view_type view() noexcept;
  const_view_type view() const noexcept;

  static constexpr std::ptrdiff_t rank() noexcept;

  reference operator()(std::array<std::ptrdiff_t, rank()> i);
  const_reference operator()(std::array<std::ptrdiff_t, rank()> i) const;

  template <typename... Is>
  std::enable_if_t<sizeof...(Is) == rank() &&
                       (true && ... && std::is_integral_v<Is>),
                   reference>
  operator()(Is... is);

  template <typename... Is>
  std::enable_if_t<sizeof...(Is) == rank() &&
                       (true && ... && std::is_integral_v<Is>),
                   const_reference>
  operator()(Is... is) const;

  std::ptrdiff_t extent(std::ptrdiff_t rank);

  const container_type& container() const noexcept;
  container_type& container() noexcept;

  std::size_t size() const noexcept;

private:
  container_type container_;
  [[no_unique_address]] mapping_type mapping_;
};

template <typename ValueType, typename... Is>
using mdarray =
    basic_mdarray<ValueType, extents<Is...>, layout_left, std_vector_policy>;

template <typename ValueType, std::ptrdiff_t Rank>
using dynamic_mdarray = basic_mdarray<ValueType, dynamic_extents<Rank>,
                                      layout_left, std_vector_policy>;

} // namespace fub