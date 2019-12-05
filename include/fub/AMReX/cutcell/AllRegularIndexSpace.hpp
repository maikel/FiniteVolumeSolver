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

#ifndef FUB_AMREX_CUTCELL_ALL_REGULAR_INDEX_SPACE_HPP
#define FUB_AMREX_CUTCELL_ALL_REGULAR_INDEX_SPACE_HPP

#include <AMReX_EB2.H>

namespace fub::amrex::cutcell {

class AllRegularIndexSpace final : public ::amrex::EB2::IndexSpace {
public:
  AllRegularIndexSpace(const ::amrex::Geometry& coarse_geometry,
                       const ::amrex::IntVect& refine_ratio,
                       int n_required_levels, int max_coarsening_level,
                       int n_grow = 4);

  const ::amrex::EB2::Level& getLevel(const ::amrex::Geometry& geom) const final;
  const ::amrex::Geometry& getGeometry(const ::amrex::Box& dom) const final;
  const ::amrex::Box& coarsestDomain() const final;

private:
  using LevelType = ::amrex::EB2::GShopLevel<
      ::amrex::EB2::GeometryShop<::amrex::EB2::AllRegularIF>>;
  std::vector<LevelType> levels_;
  std::vector<::amrex::Geometry> geometries_;
  std::vector<::amrex::Box> domains_;
  std::vector<int> n_grows_;
};

std::vector<const ::amrex::EB2::IndexSpace*> BuildRegularSpace(const ::amrex::Geometry& coarse_geometry,
                       const ::amrex::IntVect& refine_ratio,
                       int n_required_levels);

} // namespace fub::amrex

#endif