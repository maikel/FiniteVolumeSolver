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

/// \file This file defines several ForEachXXX functions to iterate through some
/// kind of range.

/// \defgroup ForEach For Each Loops
/// \brief This group contains all functions that help to iterate over some
/// range.

#ifndef FUB_FOR_EACH_HPP
#define FUB_FOR_EACH_HPP

#include "fub/State.hpp"
#include "fub/core/mdspan.hpp"

namespace fub {

/// \ingroup ForEach
/// @{
///
/// \brief Iterate through the multi-dimensional index space descibed by \p
/// mapping and invoke \p function for each such indices.
///
/// \param[in] mapping Describes how a multi dimensional index is mapped into a
/// linear space.
///
/// \param[in] function The callback object which is being invoked by this
/// function.
///
/// \return Returns the callback function obect.
///
/// The following example shows how to use this function.
///
/// ~~~~~~~~~~~~~{.cpp}
/// // Invoke a function in two one-dimensional index space
/// mdspan<double, 2> array;
/// ForEachIndex(array.mapping(), [&](int i, int j) {
///    std::cout << fmt::format("array({}, {}) = {}\n", i, j, array(i, j));
/// });
///
/// // Invoke a function in three one-dimensional index space
/// mdspan<double, 3> array;
/// ForEachIndex(array.mapping(), [](int i, int j, int k) {
///    std::cout << fmt::format("array({}, {}, {}) = \n", i, j, k, array(i, j,
///    k));
/// });
///
/// // Invoke a function where dimension is generic
/// mdspan<double, Rank> array;
/// ForEachIndex(array.mapping(), [](auto... is) {
///    std::array<int, Rank> index{is...};
///    std::cout << fmt::format("array({}) = \n", index, array(index));
/// });
/// ~~~~~~~~~~~~~
template <typename Extents, typename Function>
Function ForEachIndex(const layout_left::mapping<Extents>& mapping,
                      Function function) {
  static_assert(Extents::rank() == 1 || Extents::rank() == 2 ||
                Extents::rank() == 3 || Extents::rank() == 4);
  if constexpr (Extents::rank() == 1) {
    for (int i = 0; i < mapping.extents().extent(0); ++i) {
      function(i);
    }
  } else if constexpr (Extents::rank() == 2) {
    for (int i = 0; i < mapping.extents().extent(1); ++i) {
      for (int j = 0; j < mapping.extents().extent(0); ++j) {
        function(j, i);
      }
    }
  } else if constexpr (Extents::rank() == 3) {
    for (int i = 0; i < mapping.extents().extent(2); ++i) {
      for (int j = 0; j < mapping.extents().extent(1); ++j) {
        for (int k = 0; k < mapping.extents().extent(0); ++k) {
          function(k, j, i);
        }
      }
    }
  } else if constexpr (Extents::rank() == 4) {
    for (int c = 0; c < mapping.extents().extent(3); ++c) {
      for (int i = 0; i < mapping.extents().extent(2); ++i) {
        for (int j = 0; j < mapping.extents().extent(1); ++j) {
          for (int k = 0; k < mapping.extents().extent(0); ++k) {
            function(k, j, i, c);
          }
        }
      }
    }
  }
  return function;
}

template <typename Extents, typename Function>
Function ForEachIndex(const layout_stride::mapping<Extents>& mapping,
                      Function function) {
  static_assert(Extents::rank() == 1 || Extents::rank() == 2 ||
                Extents::rank() == 3 || Extents::rank() == 4);
  if constexpr (Extents::rank() == 1) {
    for (int i = 0; i < mapping.extents().extent(0); ++i) {
      function(i);
    }
  } else if constexpr (Extents::rank() == 2) {
    FUB_ASSERT(mapping.stride(0) < mapping.stride(1));
    for (int i = 0; i < mapping.extents().extent(1); ++i) {
      for (int j = 0; j < mapping.extents().extent(0); ++j) {
        function(j, i);
      }
    }
  } else if constexpr (Extents::rank() == 3) {
    FUB_ASSERT(mapping.stride(0) < mapping.stride(1));
    FUB_ASSERT(mapping.stride(1) < mapping.stride(2));
    for (int i = 0; i < mapping.extents().extent(2); ++i) {
      for (int j = 0; j < mapping.extents().extent(1); ++j) {
        for (int k = 0; k < mapping.extents().extent(0); ++k) {
          function(k, j, i);
        }
      }
    }
  } else if constexpr (Extents::rank() == 4) {
    for (int c = 0; c < mapping.extents().extent(3); ++c) {
      for (int i = 0; i < mapping.extents().extent(2); ++i) {
        for (int j = 0; j < mapping.extents().extent(1); ++j) {
          for (int k = 0; k < mapping.extents().extent(0); ++k) {
            function(k, j, i, c);
          }
        }
      }
    }
  }
  return function;
}
/// @}

template <int Rank, typename Function>
Function ForEachIndex(const IndexBox<Rank>& box, Function function) {
  static_assert(Rank == 1 || Rank == 2 || Rank == 3 || Rank == 4);
  if constexpr (Rank == 1) {
    for (std::ptrdiff_t i = box.lower[0]; i < box.upper[0]; ++i) {
      function(i);
    }
  } else if constexpr (Rank == 2) {
    for (std::ptrdiff_t i = box.lower[1]; i < box.upper[1]; ++i) {
      for (std::ptrdiff_t j = box.lower[0]; j < box.upper[0]; ++j) {
        function(j, i);
      }
    }
  } else if constexpr (Rank == 3) {
    for (std::ptrdiff_t i = box.lower[2]; i < box.upper[2]; ++i) {
      for (std::ptrdiff_t j = box.lower[1]; j < box.upper[1]; ++j) {
        for (std::ptrdiff_t k = box.lower[0]; k < box.upper[0]; ++k) {
          function(k, j, i);
        }
      }
    }
  } else if constexpr (Rank == 4) {
    for (std::ptrdiff_t c = box.lower[3]; c < box.upper[3]; ++c) {
      for (std::ptrdiff_t i = box.lower[2]; i < box.upper[2]; ++i) {
        for (std::ptrdiff_t j = box.lower[1]; j < box.upper[1]; ++j) {
          for (std::ptrdiff_t k = box.lower[0]; k < box.upper[0]; ++k) {
            function(k, j, i, c);
          }
        }
      }
    }
  }
  return function;
}

} // namespace fub

#endif