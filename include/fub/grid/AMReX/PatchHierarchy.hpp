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

#ifndef FUB_AMREX_PATCH_HIERARCHY_HPP
#define FUB_AMREX_PATCH_HIERARCHY_HPP

#include "fub/Duration.hpp"
#include "fub/Equation.hpp"
#include "fub/ext/Eigen.hpp"
#include "fub/grid/AMReX/CartesianGridGeometry.hpp"
#include "fub/grid/AMReX/PatchHandle.hpp"

#include <AMReX_FluxRegister.H>
#include <AMReX_Geometry.H>
#include <AMReX_MultiFab.H>
#include <AMReX_PlotFileUtil.H>

#include <fmt/format.h>

#include <vector>

namespace fub {
namespace amrex {

struct PatchHierarchyOptions {
  int max_number_of_levels{1};
};

struct PatchLevel {
  PatchLevel() = default;

  PatchLevel(const PatchLevel&) = delete;
  PatchLevel& operator=(const PatchLevel&) = delete;

  PatchLevel(PatchLevel&&) noexcept = default;
  PatchLevel& operator=(PatchLevel&&) noexcept = default;

  ~PatchLevel() = default;

  PatchLevel(int num, Duration tp, const ::amrex::BoxArray& ba,
             const ::amrex::DistributionMapping& dm, int n_components);

  PatchLevel(int num, Duration tp, const ::amrex::BoxArray& ba,
             const ::amrex::DistributionMapping& dm, int n_components,
             const ::amrex::FabFactory<::amrex::FArrayBox>& factory);

  int level_number{};
  Duration time_point{};
  std::ptrdiff_t cycles{};
  ::amrex::MultiFab data{};
};

struct DataDescription {
  int n_state_components;
  int first_cons_component;
  int n_cons_components;
};

class PatchHierarchy {
public:
  explicit PatchHierarchy(DataDescription description,
                          const CartesianGridGeometry& geometry,
                          const PatchHierarchyOptions& options);

  const ::amrex::Geometry& GetGeometry(int level) const noexcept {
    return patch_level_geometry_[static_cast<std::size_t>(level)];
  }
  const PatchHierarchyOptions& GetOptions() const noexcept { return options_; }

  const CartesianGridGeometry& GetGridGeometry() const noexcept {
    return grid_geometry_;
  }

  std::ptrdiff_t GetCycles(int level = 0) const {
    return patch_level_[static_cast<std::size_t>(level)].cycles;
  }

  Duration GetTimePoint(int level = 0) const {
    return patch_level_[static_cast<std::size_t>(level)].time_point;
  }

  int GetNumberOfLevels() const noexcept {
    return static_cast<int>(std::count_if(
        patch_level_.begin(), patch_level_.end(),
        [](const PatchLevel& level) { return !level.data.empty(); }));
  }

  int GetMaxNumberOfLevels() const noexcept {
    return static_cast<int>(patch_level_.size());
  }

  int GetRatioToCoarserLevel(int level) const noexcept {
    if (level) {
      return 2;
    }
    return 1;
  }

  PatchLevel& GetPatchLevel(int level) noexcept {
    return patch_level_[static_cast<std::size_t>(level)];
  }

  const PatchLevel& GetPatchLevel(int level) const noexcept {
    return patch_level_[static_cast<std::size_t>(level)];
  }

  const DataDescription& GetDataDescription() const noexcept {
    return description_;
  }

  template <typename Feedback>
  Feedback ForEachPatch(int level, Feedback feedback) {
    for (::amrex::MFIter mfi(GetPatchLevel(level).data); mfi.isValid(); ++mfi) {
      PatchHandle handle{level, &mfi};
      std::invoke(feedback, handle);
    }
    return feedback;
  }

private:
  DataDescription description_;
  CartesianGridGeometry grid_geometry_;
  PatchHierarchyOptions options_;
  std::vector<PatchLevel> patch_level_;
  std::vector<::amrex::Geometry> patch_level_geometry_;
};

template <typename Equation>
DataDescription MakeDataDescription(const Equation& equation) {
  const auto complete_shape = equation.Shape(complete);
  const int n_comp = boost::hana::fold_left(
      complete_shape, 0, [](int size, int ncomp) { return size + ncomp; });
  const auto cons_shape = equation.Shape(cons);
  const int n_cons_comp = boost::hana::fold_left(
      cons_shape, 0, [](int size, int ncomp) { return size + ncomp; });
  DataDescription desc;
  desc.n_state_components = n_comp;
  desc.first_cons_component = 0;
  desc.n_cons_components = n_cons_comp;
  return desc;
}

template <typename Equation>
void WritePlotFile(const std::string plotfilename,
                   const fub::amrex::PatchHierarchy& hier,
                   const Equation& equation) {
  const int nlevels = hier.GetNumberOfLevels();
  const double time_point = hier.GetTimePoint().count();
  FUB_ASSERT(nlevels >= 0);
  std::size_t size = static_cast<std::size_t>(nlevels);
  ::amrex::Vector<const ::amrex::MultiFab*> mf(size);
  ::amrex::Vector<::amrex::Geometry> geoms(size);
  ::amrex::Vector<int> level_steps(size);
  ::amrex::Vector<::amrex::IntVect> ref_ratio(size);
  for (std::size_t i = 0; i < size; ++i) {
    mf[i] = &hier.GetPatchLevel(static_cast<int>(i)).data;
    geoms[i] = hier.GetGeometry(static_cast<int>(i));
    level_steps[i] = static_cast<int>(hier.GetCycles(static_cast<int>(i)));
    ref_ratio[i] = hier.GetRatioToCoarserLevel(static_cast<int>(i)) *
                   ::amrex::IntVect::TheUnitVector();
  }
  constexpr auto names = Equation::Complete::Names();
  const auto shape = equation.Shape(complete);
  ::amrex::Vector<std::string> varnames;
  varnames.reserve(boost::hana::length(names));
  boost::hana::for_each(boost::hana::zip(names, shape), [&](auto xs) {
    const int ncomp = at_c<1>(xs);
    if (ncomp == 1) {
      varnames.push_back(at_c<0>(xs).c_str());
    } else {
      for (int i = 0; i < ncomp; ++i) {
        varnames.push_back(fmt::format("{}_{}", at_c<0>(xs).c_str(), i));
      }
    }
  });
  ::amrex::WriteMultiLevelPlotfile(plotfilename, nlevels, mf, varnames, geoms,
                                   time_point, level_steps, ref_ratio);
}

} // namespace amrex
} // namespace fub

#endif