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

#include "fub/Duration.hpp"
#include "fub/ext/uuid.hpp"

#include <SAMRAI/geom/CartesianGridGeometry.h>
#include <SAMRAI/geom/CartesianPatchGeometry.h>
#include <SAMRAI/hier/PatchHierarchy.h>

namespace fub::samrai {

SAMRAI::hier::ComponentSelector
SelectComponents(const SAMRAI::hier::PatchDescriptor& desc);
SAMRAI::hier::ComponentSelector SelectComponents(span<const int> data_ids);

struct PatchHierarchyOptions {
  SAMRAI::hier::IntVector refine_ratio;
  int max_number_of_levels{1};
};

template <std::size_t Rank> struct CoordinateRange {
  std::array<double, Rank> lower;
  std::array<double, Rank> upper;
};

template <typename I, std::size_t Rank>
std::enable_if_t<std::is_integral_v<I>,
                 std::shared_ptr<SAMRAI::geom::CartesianGridGeometry>>
MakeCartesianGridGeometry(const std::array<I, Rank>& n_cells,
                          const CoordinateRange<Rank>& coordinates) {
  const double* x_lo = coordinates.lower.data();
  const double* x_up = coordinates.upper.data();
  const SAMRAI::hier::Index idx_up = std::apply(
      [](auto... is) { return SAMRAI::hier::Index(int(is - 1)...); }, n_cells);
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

  PatchHierarchy(PatchHierarchy&& ph) = default;
  PatchHierarchy& operator=(PatchHierarchy&& ph) = default;

  PatchHierarchy(const PatchHierarchy& ph);
  PatchHierarchy& operator=(const PatchHierarchy& ph) {
    PatchHierarchy tmp(ph);
    std::swap(*this, tmp);
    return *this;
  }

  /// \brief Return some additional patch hierarchy options.
  [[nodiscard]] const PatchHierarchyOptions& GetOptions() const noexcept;

  [[nodiscard]] int GetMaxNumberOfLevels() const noexcept;
  [[nodiscard]] int GetNumberOfLevels() const noexcept;

  [[nodiscard]] const DataDescription& GetDataDescription() const noexcept;

  [[nodiscard]] const std::shared_ptr<SAMRAI::hier::PatchHierarchy>&
  GetNative() const noexcept;

  [[nodiscard]] const SAMRAI::geom::CartesianGridGeometry&
  GetGeometry(int level) const noexcept;

  [[nodiscard]] span<const int> GetDataIds() const noexcept;

  [[nodiscard]] std::shared_ptr<SAMRAI::hier::PatchLevel>
  GetPatchLevel(int level) const;

  [[nodiscard]] std::ptrdiff_t GetCycles(int level = 0) const;
  [[nodiscard]] Duration GetTimePoint(int level = 0) const;

  void SetCycles(std::ptrdiff_t cycles, int level);
  void SetTimePoint(Duration time_point, int level);

private:
  std::shared_ptr<SAMRAI::hier::PatchHierarchy> hierarchy_;
  DataDescription data_desc_;
  PatchHierarchyOptions options_;
  std::vector<std::ptrdiff_t> cycles_;
  std::vector<Duration> time_points_;
};

template <typename Equation>
PatchHierarchy::PatchHierarchy(
    const Equation& eq,
    std::shared_ptr<SAMRAI::geom::CartesianGridGeometry> geom,
    PatchHierarchyOptions hier_opts)
    : PatchHierarchy(RegisterVariables(eq), std::move(geom),
                     std::move(hier_opts)) {}

template <typename... Is>
std::array<double, sizeof...(Is)>
GetCellCenter(const SAMRAI::geom::CartesianGridGeometry& geom, Is... is) {
  const double* dx = geom.getDx();
  const double* xlo = geom.getXLower();
  std::array<std::common_type_t<Is...>, sizeof...(Is)> index{
      std::common_type_t<Is...>(is)...};
  std::array<double, sizeof...(Is)> x;
  for (std::size_t d = 0; d < sizeof...(Is); ++d) {
    x[d] = xlo[d] + 0.5 * dx[d] + index[d] * dx[d];
  }
  return x;
}

} // namespace fub::samrai

#endif
