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

#ifndef FUB_AMREX_CUTCELL_INDEX_SPACE_HPP
#define FUB_AMREX_CUTCELL_INDEX_SPACE_HPP

#include <AMReX.H>

#ifdef AMREX_USE_EB
#include <AMReX_EB2.H>

#include <vector>

namespace fub {
namespace amrex {
namespace cutcell {

constexpr std::ptrdiff_t ipow(int base, int exponent) {
  std::ptrdiff_t prod{1};
  while (exponent > 0) {
    prod *= base;
    exponent -= 1;
  }
  return prod;
}

template <typename GShop>
std::vector<const ::amrex::EB2::IndexSpace*>
MakeIndexSpaces(GShop&& shop, const ::amrex::Geometry& coarse_geom,
                int n_level) {
  FUB_ASSERT(n_level > 0);
  std::vector<const ::amrex::EB2::IndexSpace*> index_spaces(
      static_cast<std::size_t>(n_level));
  ::amrex::Geometry geom = coarse_geom;
  for (int level = 0; level < n_level; ++level) {
    ::amrex::EB2::Build(shop, geom, level, level);
    index_spaces[static_cast<std::size_t>(level)] =
        &::amrex::EB2::IndexSpace::top();
    geom.refine({AMREX_D_DECL(2, 2, 2)});
  }
  return index_spaces;
}

template <typename GShop>
std::vector<const ::amrex::EB2::IndexSpace*>
MakeIndexSpaces(GShop&& shop, const ::amrex::Geometry& coarse_geom,
                int n_level, const ::amrex::IntVect& refine_ratio) {
  FUB_ASSERT(n_level > 0);
  std::vector<const ::amrex::EB2::IndexSpace*> index_spaces(
                                                            static_cast<std::size_t>(n_level));
  ::amrex::Geometry geom = coarse_geom;
  for (int level = 0; level < n_level; ++level) {
    ::amrex::EB2::Build(shop, geom, refine_ratio, level, level);
    index_spaces[static_cast<std::size_t>(level)] =
    &::amrex::EB2::IndexSpace::top();
    geom.refine(refine_ratio);
  }
  return index_spaces;
}

} // namespace cutcell
} // namespace amrex
} // namespace fub
#endif

#endif
