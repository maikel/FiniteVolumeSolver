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
template <>
::amrex::IntVect GetOptionOr(const ProgramOptions& map, const std::string& name,
                             const ::amrex::IntVect& value) {
  std::array<int, 3> iv{};
  std::copy_n(value.begin(), AMREX_SPACEDIM, iv.begin());
  iv = GetOptionOr(map, name, iv);
  return ::amrex::IntVect(AMREX_D_DECL(iv[0], iv[1], iv[2]));
}

template <>
::amrex::Box GetOptionOr(const ProgramOptions& map, const std::string& name,
                         const ::amrex::Box& value) {
  ProgramOptions box = GetOptions(map, name);
  ::amrex::IntVect lo = GetOptionOr(box, "lower", value.smallEnd());
  ::amrex::IntVect hi = GetOptionOr(box, "upper", value.bigEnd());
  return ::amrex::Box(lo, hi);
}

template <>
::amrex::RealBox GetOptionOr(const ProgramOptions& map, const std::string& name,
                             const ::amrex::RealBox& value) {
  ProgramOptions box = GetOptions(map, name);
  std::array<double, 3> lo{};
  std::array<double, 3> hi{};
  std::copy_n(value.lo(), AMREX_SPACEDIM, lo.begin());
  std::copy_n(value.hi(), AMREX_SPACEDIM, hi.begin());
  lo = GetOptionOr(box, "lower", lo);
  hi = GetOptionOr(box, "upper", hi);
  return ::amrex::RealBox(lo.data(), hi.data());
}

namespace amrex {

std::array<std::ptrdiff_t, static_cast<std::size_t>(AMREX_SPACEDIM)>
AsArray(const ::amrex::IntVect& vec) {
  std::array<std::ptrdiff_t, static_cast<std::size_t>(AMREX_SPACEDIM)> array{};
  for (std::size_t i = 0; i < static_cast<std::size_t>(AMREX_SPACEDIM); ++i) {
    array[i] = vec[static_cast<int>(i)];
  }
  return array;
}

std::array<::amrex::Box, 2>
GetCellsAndFacesInStencilRange(const ::amrex::Box& cell_validbox,
                               const ::amrex::Box& face_tilebox,
                               int stencil_width, Direction dir) {
  const int dir_v = static_cast<int>(dir);
  const ::amrex::Box cell_tilebox =
      cell_validbox &
      ::amrex::grow(::amrex::enclosedCells(face_tilebox), dir_v, stencil_width);
  const ::amrex::Box face_tilebox_in_range =
      face_tilebox &
      ::amrex::surroundingNodes(
          ::amrex::grow(cell_tilebox, dir_v, -stencil_width), dir_v);
  return {cell_tilebox, face_tilebox_in_range};
}

} // namespace amrex
} // namespace fub
