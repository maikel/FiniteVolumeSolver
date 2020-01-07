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

#include "fub/AMReX/ViewFArrayBox.hpp"

namespace fub {
namespace amrex {

std::array<std::ptrdiff_t, static_cast<std::size_t>(AMREX_SPACEDIM)>
AsArray(const ::amrex::IntVect& vec) {
  std::array<std::ptrdiff_t, static_cast<std::size_t>(AMREX_SPACEDIM)> array;
  for (std::size_t i = 0; i < static_cast<std::size_t>(AMREX_SPACEDIM); ++i) {
    array[i] = vec[static_cast<int>(i)];
  }
  return array;
}

std::array<::amrex::Box, 2>
GetCellsAndFacesInStencilRange(const ::amrex::Box& cell_tilebox,
                               const ::amrex::Box& face_validbox,
                               int stencil_width, Direction dir) {
  const int dir_v = static_cast<int>(dir);
  ::amrex::Box all_faces_in_tile = cell_tilebox;
  all_faces_in_tile.surroundingNodes(dir_v);
  all_faces_in_tile.grow(dir_v, -stencil_width);
  const ::amrex::Box face_tilebox = face_validbox & all_faces_in_tile;
  ::amrex::Box cells_in_stencil_range = face_tilebox;
  cells_in_stencil_range.enclosedCells();
  cells_in_stencil_range.grow(dir_v, stencil_width);
  return {cells_in_stencil_range, face_tilebox};
}

} // namespace amrex
} // namespace fub
