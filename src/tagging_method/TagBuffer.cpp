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

#include "fub/tagging_method/TagBuffer.hpp"

#include "fub/ForEach.hpp"

namespace fub {

void TagBuffer::TagCellsForRefinement(
    const PatchDataView<char, 1, layout_stride>& tags) {
  span<char> data = tags.Span();
  for (std::ptrdiff_t i = buffer_width_; i < data.size() - buffer_width_; ++i) {
    if (data[i] == '\1') {
      for (std::ptrdiff_t j = 0; j < buffer_width_; ++j) {
        data[i + j] = static_cast<char>('\2' - data[i + j]);
      }
      for (std::ptrdiff_t j = buffer_width_; j > 0; --j) {
        data[i - j - 1] = static_cast<char>('\2' - data[i - j - 1]);
      }
    }
  }
}

void TagBuffer::TagCellsForRefinement(
    const PatchDataView<char, 2, layout_stride>& tags) {
  mdspan<char, 2, layout_stride> data = tags.MdSpan();
  const std::ptrdiff_t nx = data.extent(0);
  const std::ptrdiff_t ny = data.extent(1);
  for (std::ptrdiff_t j = 0; j < ny; ++j) {
    for (std::ptrdiff_t i = 0; i < nx; ++i) {
      if (data(i, j) == '\1') {
        const std::ptrdiff_t k0 =
            std::max(i - buffer_width_, std::ptrdiff_t{0});
        const std::ptrdiff_t k1 = std::min(i + buffer_width_, nx);
        const std::ptrdiff_t l0 =
            std::max(j - buffer_width_, std::ptrdiff_t{0});
        const std::ptrdiff_t l1 = std::min(j + buffer_width_, ny);
        for (std::ptrdiff_t l = l0; l < l1; ++l) {
          for (std::ptrdiff_t k = k0; k < k1; ++k) {
            data(k, l) =
                std::max(data(k, l), static_cast<char>('\2' - data(k, l)));
          }
        }
      }
    }
  }
}

void TagBuffer::TagCellsForRefinement(
    const PatchDataView<char, 3, layout_stride>& tags) {
  const std::array<std::ptrdiff_t, 2> shift{buffer_width_, buffer_width_};
  const Direction X = Direction::X;
  const Direction Y = Direction::Y;
  const Direction Z = Direction::Z;
  IndexBox<3> box =
      Shrink(Shrink(Shrink(tags.Box(), X, shift), Y, shift), Z, shift);
  ForEachIndex(box, [&](std::ptrdiff_t i0, std::ptrdiff_t j0,
                        std::ptrdiff_t k0) {
    Index<3> lower{i0 - buffer_width_, j0 - buffer_width_, k0 - buffer_width_};
    Index<3> upper{i0 + buffer_width_, j0 + buffer_width_, k0 + buffer_width_};
    IndexBox<3> neighborhood{lower, upper};
    ForEachIndex(neighborhood, [&](std::ptrdiff_t i, std::ptrdiff_t j,
                                   std::ptrdiff_t k) {
      tags(i, j, k) =
          std::max(tags(i, j, k), static_cast<char>('\2' - tags(i, j, k)));
    });
  });
}

} // namespace fub