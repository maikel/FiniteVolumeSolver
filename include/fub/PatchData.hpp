// Copyright (c) 2020 Maikel Nadolski
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

#ifndef FUB_PATCH_DATA_HPP
#define FUB_PATCH_DATA_HPP

#include "fub/PatchDataView.hpp"

namespace fub {

template <typename T, int Rank> class PatchData {
private:
  dynamic_mdarray<T, Rank> array_;
  Index<Rank> origin_;

public:
  IndexBox<Rank> Box() const noexcept;

  PatchDataView<T, Rank, layout_left> View() noexcept;
  PatchDataView<const T, Rank, layout_left> View() const noexcept;
  PatchDataView<const T, Rank, layout_left> ConstView() const noexcept;
};

template <int Rank> struct CellIndex;
template <int Rank> struct FaceIndex;
template <int Rank> struct NodeIndex;

template <typename T, int Rank> class CellData {
private:
  PatchData<T, Rank> patch_data_;

public:
};

template <typename T, int Rank> class FaceData {
private:
  PatchData<T, Rank> patch_data_;

public:
};

template <typename T, int Rank> class NodeData {
private:
  PatchData<T, Rank> patch_data_;

public:
};

} // namespace fub

#endif