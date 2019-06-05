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

/// \file This file defines a small wrapper class which stores thread local
/// instances of a given object.

#ifndef FUB_EXT_OMP_HPP
#define FUB_EXT_OMP_HPP

#include "fub/core/assert.hpp"

#if defined(_OPENMP)
extern "C" {
#include <omp.h>
}
#endif

#include <vector>

namespace fub {

template <typename T, typename Allocator = std::allocator<T>> class OmpLocal {
public:
  explicit OmpLocal(const T& value, Allocator alloc = Allocator());

  T& Get() noexcept;
  const T& Get() const noexcept;

  T* operator->() noexcept;
  const T* operator->() const noexcept;

private:
  std::vector<T, Allocator> instances_;
};

#if defined(_OPENMP)
template <typename T, typename Allocator>
OmpLocal<T, Allocator>::OmpLocal(const T& value, Allocator alloc)
    : instances_(static_cast<std::size_t>(::omp_get_max_threads()), value,
                 alloc) {}

template <typename T, typename Allocator>
T& OmpLocal<T, Allocator>::Get() noexcept {
  const int thread_num = ::omp_get_thread_num();
  FUB_ASSERT(thread_num >= 0);
  const std::size_t index = static_cast<std::size_t>(thread_num);
  return instances_[index];
}

template <typename T, typename Allocator>
const T& OmpLocal<T, Allocator>::Get() const noexcept {
  const int thread_num = ::omp_get_thread_num();
  FUB_ASSERT(thread_num >= 0);
  const std::size_t index = static_cast<std::size_t>(thread_num);
  return instances_[index];
}
#else
template <typename T, typename Allocator>
OmpLocal<T, Allocator>::OmpLocal(const T& value, Allocator alloc)
    : instances_({value}, alloc) {}

template <typename T, typename Allocator>
T& OmpLocal<T, Allocator>::Get() noexcept {
  return instances_[0];
}

template <typename T, typename Allocator>
const T& OmpLocal<T, Allocator>::Get() const noexcept {
  return instances_[0];
}
#endif

template <typename T, typename Allocator>
T* OmpLocal<T, Allocator>::operator->() noexcept {
  return &Get();
}

template <typename T, typename Allocator>
const T* OmpLocal<T, Allocator>::operator->() const noexcept {
  return &Get();
}

} // namespace fub

#endif
