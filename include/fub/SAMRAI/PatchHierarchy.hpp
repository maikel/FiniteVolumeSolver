// Copyright (c) 2019 Maikel Nadolski
// Copyright (c) 2019 Patrick Denzler
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

#ifndef FUB_SAMRAI_PATCH_HIERARCHY_HPP
#define FUB_SAMRAI_PATCH_HIERARCHY_HPP

#include "fub/SAMRAI/RegisterVariables.hpp"
#include <SAMRAI/hier/PatchHierarchy.h>
#include <SAMRAI/geom/CartesianGridGeometry.h>

namespace fub::samrai {
  struct PatchHierarchyOptions {
    SAMRAI::hier::IntVector refine_ratio;
    int max_number_of_levels{1};
  };

  template <std::size_t Rank>
  struct CoordinateRange {
    std::array<double, Rank> lower;
    std::array<double, Rank> upper;
  };

  template <std::size_t Rank>
  std::shared_ptr<SAMRAI::geom::CartesianGridGeometry> MakeCartesianGridGeometry(const std::array<std::ptrdiff_t, Rank>& n_cells, const CoordinateRange<Rank>& coordinates, std::string prefix = std::string()) {
    const double* x_lo = coordinates.lower.data();
    const double* x_up = coordinates.upper.data();
    const SAMRAI::hier::Index idx_up = std::apply([](auto... is) {
      return SAMRAI::hier::Index((is - 1)...);
      }, n_cells);
    const SAMRAI::tbox::Dimension dim(Rank);
    const SAMRAI::hier::Box domain_box(SAMRAI::hier::Index::getZeroIndex(dim), idx_up, SAMRAI::hier::BlockId(0));
    SAMRAI::hier::BoxContainer domain{domain_box};
    return std::make_shared<SAMRAI::geom::CartesianGridGeometry>(prefix + "_Geometry", x_lo, x_up, domain);
  }

  class PatchHierarchy {
  public:
    /// \brief Constructs a PatchHierarchy object which is capable of holding data
    /// described by the secified data description on given geometry extents.
    template<typename Equation>
    PatchHierarchy(const Equation& eq, std::shared_ptr<SAMRAI::geom::CartesianGridGeometry> geom, PatchHierarchyOptions hier_opts, std::string prefix = std::string());

    /// \brief Return some additional patch hierarchy options.
    const PatchHierarchyOptions& GetOptions() const noexcept;

    const DataDescription& GetDataDescription() const noexcept;

    const std::shared_ptr<SAMRAI::hier::PatchHierarchy>& GetNative() const noexcept;

  private:
    std::shared_ptr<SAMRAI::hier::PatchHierarchy> hierarchy_;
    DataDescription data_desc_;
    PatchHierarchyOptions options_;
  };

  template <typename Equation>
  PatchHierarchy::PatchHierarchy(
      const Equation& eq,
      std::shared_ptr<SAMRAI::geom::CartesianGridGeometry> geom,
      PatchHierarchyOptions hier_opts, std::string prefix)
    : data_desc_{RegisterVariables(eq)}, options_{hier_opts} {
      hierarchy_ = std::make_shared<SAMRAI::hier::PatchHierarchy>(prefix + "_Hierarchy", geom);
      hierarchy_->setMaxNumberOfLevels(options_.max_number_of_levels);
      for(int i = 1; i < options_.max_number_of_levels; ++i) {
        hierarchy_->setRatioToCoarserLevel(options_.refine_ratio, i);
      }
  }


}

#endif
