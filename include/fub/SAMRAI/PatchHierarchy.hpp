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

#include "fub/ext/uuid.hpp"

#include <SAMRAI/geom/CartesianGridGeometry.h>
#include <SAMRAI/hier/PatchHierarchy.h>

namespace fub::samrai {

SAMRAI::hier::ComponentSelector SelectComponents(span<const int> data_ids);

struct PatchHierarchyOptions {
  SAMRAI::hier::IntVector refine_ratio;
  int max_number_of_levels{1};
};

template <std::size_t Rank> struct CoordinateRange {
  std::array<double, Rank> lower;
  std::array<double, Rank> upper;
};

template <std::size_t Rank>
std::shared_ptr<SAMRAI::geom::CartesianGridGeometry>
MakeCartesianGridGeometry(const std::array<std::ptrdiff_t, Rank>& n_cells,
                          const CoordinateRange<Rank>& coordinates) {
  const double* x_lo = coordinates.lower.data();
  const double* x_up = coordinates.upper.data();
  const SAMRAI::hier::Index idx_up = std::apply(
      [](auto... is) { return SAMRAI::hier::Index((is - 1)...); }, n_cells);
  const SAMRAI::tbox::Dimension dim(Rank);
  const SAMRAI::hier::Box domain_box(SAMRAI::hier::Index::getZeroIndex(dim),
                                     idx_up, SAMRAI::hier::BlockId(0));
  SAMRAI::hier::BoxContainer domain{domain_box};
  return std::make_shared<SAMRAI::geom::CartesianGridGeometry>(
      MakeUniqueName(), x_lo, x_up, domain);
}

class PatchHierarchy {
public:
  /// \brief Constructs a PatchHierarchy object which is capable of holding data
  /// described by the secified data description on given geometry extents.
  template <typename Equation>
  PatchHierarchy(const Equation& eq,
                 std::shared_ptr<SAMRAI::geom::CartesianGridGeometry> geom,
                 PatchHierarchyOptions hier_opts);

  PatchHierarchy(DataDescription dd,
                 std::shared_ptr<SAMRAI::geom::CartesianGridGeometry> geom,
                 PatchHierarchyOptions hier_opts);

  PatchHierarchy(const PatchHierarchy& ph);
  PatchHierarchy& operator=(const PatchHierarchy& ph) {
    PatchHierarchy tmp(ph);
    return (*this = std::move(tmp));
  }

  PatchHierarchy(PatchHierarchy&& ph) = default;
  PatchHierarchy& operator=(PatchHierarchy&& ph) = default;

  /// \brief Return some additional patch hierarchy options.
  [[nodiscard]] const PatchHierarchyOptions& GetOptions() const noexcept;

  [[nodiscard]] const DataDescription& GetDataDescription() const noexcept;

  [[nodiscard]] const std::shared_ptr<SAMRAI::hier::PatchHierarchy>& GetNative() const
      noexcept;

private:
  std::shared_ptr<SAMRAI::hier::PatchHierarchy> hierarchy_;
  DataDescription data_desc_;
  PatchHierarchyOptions options_;
};

template <typename Equation>
PatchHierarchy::PatchHierarchy(
    const Equation& eq,
    std::shared_ptr<SAMRAI::geom::CartesianGridGeometry> geom,
    PatchHierarchyOptions hier_opts)
    : PatchHierarchy(RegisterVariables(eq), std::move(geom),
                     std::move(hier_opts)) {}

} // namespace fub::samrai

#endif
