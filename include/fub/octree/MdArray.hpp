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

#ifndef FUB_P4EST_MDARRAY_HPP
#define FUB_P4EST_MDARRAY_HPP

namespace fub::p4est {

template <typename T, std::size_t Rank, typename Allocator = std::allocator<T>>
class MdArray {
public: 
  using container = std::vector<T, Allocator>;
  using value_type = typename container::value_type;
  using pointer = typename container::pointer;
  using const_pointer = typename container::const_pointer;
  using reference = typename container::reference;
  using const_reference = typename container::const_reference;

  MdArray();

  reference operator()(std::array<std::ptrdiff_t, Rank> i) noexcept;
  const_reference operator()(std::array<std::ptrdiff_t, Rank> i) const noexcept;

  template <typename... Is> reference operator()(Is... is) noexcept;
  template <typename... Is> const_reference operator()(Is... is) const noexcept;

  pointer data() noexcept;
  const_pointer data() const noexcept;
  const_pointer cdata() const noexcept;

private:
  std::vector<T, Allocator> container_;
  std::array<std::ptrdiff_t, Rank> extents_;
};

}

#endif