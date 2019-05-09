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

#ifndef FUB_AMREX_CUTCELL_PATCH_HIERARCHY_HPP
#define FUB_AMREX_CUTCELL_PATCH_HIERARCHY_HPP

#include "fub/AMReX/CartesianGridGeometry.hpp"
#include "fub/AMReX/PatchHandle.hpp"
#include "fub/AMReX/PatchHierarchy.hpp"
#include "fub/CutCellData.hpp"
#include "fub/Duration.hpp"
#include "fub/Equation.hpp"
#include "fub/ext/Eigen.hpp"

#include <AMReX_EBFabFactory.H>
#include <AMReX_FluxRegister.H>
#include <AMReX_Geometry.H>
#include <AMReX_MultiCutFab.H>
#include <AMReX_MultiFab.H>
#include <AMReX_PlotFileUtil.H>

#include <fmt/format.h>

#include <functional>
#include <vector>

namespace fub {
namespace amrex {
namespace cutcell {

/// This class holds state data arrays for each refinement level of a patch
/// hierarchy.
///
/// This cut-cell version stores some embedded boundary informations in addition
/// to the normal patch level type.
struct PatchLevel : ::fub::amrex::PatchLevel {
  PatchLevel() = default;
  PatchLevel(const PatchLevel&);
  PatchLevel& operator=(const PatchLevel&);
  PatchLevel(PatchLevel&&) noexcept = default;
  PatchLevel& operator=(PatchLevel&&) noexcept = default;
  ~PatchLevel() = default;

  PatchLevel(int level, Duration tp, const ::amrex::BoxArray& ba,
             const ::amrex::DistributionMapping& dm, int n_components,
             std::shared_ptr<::amrex::EBFArrayBoxFactory> factory);

  using MultiCutFabs =
      std::array<std::unique_ptr<::amrex::MultiCutFab>, AMREX_SPACEDIM>;

  std::shared_ptr<::amrex::EBFArrayBoxFactory> factory;
  MultiCutFabs unshielded;
  MultiCutFabs shielded_left;
  MultiCutFabs shielded_right;
  MultiCutFabs doubly_shielded;
};

/// This class extents the normal hierarchy options with a pointer to an
/// embedded boundary index space for each possible refinement level.
struct PatchHierarchyOptions : public ::fub::amrex::PatchHierarchyOptions {
  std::vector<const ::amrex::EB2::IndexSpace*> index_spaces;
};

class PatchHierarchy {
public:
  using PatchHandle = ::fub::amrex::PatchHandle;

  PatchHierarchy(DataDescription description,
                 const CartesianGridGeometry& geometry,
                 const PatchHierarchyOptions& options);

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

  int GetRatioToCoarserLevel(int level, Direction dir) const noexcept;

  ::amrex::IntVect GetRatioToCoarserLevel(int level) const noexcept;

  const ::amrex::Geometry& GetGeometry(int level) const noexcept {
    FUB_ASSERT(0 <= level && level < GetMaxNumberOfLevels());
    const std::size_t level_num = static_cast<std::size_t>(level);
    return patch_level_geometry_[level_num];
  }

  PatchLevel& GetPatchLevel(int level) noexcept {
    FUB_ASSERT(0 <= level && level < GetMaxNumberOfLevels());
    const std::size_t level_num = static_cast<std::size_t>(level);
    return patch_level_[level_num];
  }

  const PatchLevel& GetPatchLevel(int level) const noexcept {
    FUB_ASSERT(0 <= level && level < GetMaxNumberOfLevels());
    const std::size_t level_num = static_cast<std::size_t>(level);
    return patch_level_[level_num];
  }

  const DataDescription& GetDataDescription() const noexcept {
    return description_;
  }

  const std::shared_ptr<::amrex::EBFArrayBoxFactory>&
  GetEmbeddedBoundary(int level) const noexcept {
    FUB_ASSERT(0 <= level && level < GetMaxNumberOfLevels());
    const std::size_t level_num = static_cast<std::size_t>(level);
    return patch_level_[level_num].factory;
  }

  template <typename Feedback>
  Feedback ForEachPatch(int level, Feedback feedback) const {
#ifdef _OPENMP
#pragma omp parallel firstprivate(feedback)
#endif
    {
      for (::amrex::MFIter mfi(GetPatchLevel(level).data); mfi.isValid();
           ++mfi) {
        PatchHandle handle{level, &mfi};
        feedback(handle);
      }
    }
    return feedback;
  }

  template <typename Feedback> double Minimum(int level, Feedback feedback) const {
    double global_min = std::numeric_limits<double>::infinity();
#ifdef _OPENMP
#pragma omp parallel reduction(min : global_min) firstprivate(feedback)
#endif
    {
      for (::amrex::MFIter mfi(GetPatchLevel(level).data); mfi.isValid();
           ++mfi) {
        PatchHandle handle{level, &mfi};
        const double local_min = feedback(handle);
        global_min = std::min(global_min, local_min);
      }
    }
    return global_min;
  }

  CutCellData<AMREX_SPACEDIM> GetCutCellData(PatchHandle patch,
                                             Direction dir) const;

private:
  DataDescription description_;
  CartesianGridGeometry grid_geometry_;
  PatchHierarchyOptions options_;
  std::vector<PatchLevel> patch_level_;
  std::vector<::amrex::Geometry> patch_level_geometry_;
};

template <typename Equation>
void WritePlotFile(const std::string plotfilename, const PatchHierarchy& hier,
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
  using Traits = StateTraits<Complete<Equation>>;
  constexpr auto names = Traits::names;
  const auto depths = Depths<Complete<Equation>>(equation);
  const std::size_t n_names =
      std::tuple_size<remove_cvref_t<decltype(names)>>::value;
  ::amrex::Vector<std::string> varnames;
  varnames.reserve(n_names);
  boost::mp11::tuple_for_each(Zip(names, ToTuple(depths)), [&](auto xs) {
    const int ncomp = std::get<1>(xs);
    if (ncomp == 1) {
      varnames.push_back(std::get<0>(xs));
    } else {
      for (int i = 0; i < ncomp; ++i) {
        varnames.push_back(fmt::format("{}_{}", std::get<0>(xs), i));
      }
    }
  });
  ::amrex::EB_WriteMultiLevelPlotfile(plotfilename, nlevels, mf, varnames,
                                      geoms, time_point, level_steps,
                                      ref_ratio);
}

void WriteCheckpointFile(const std::string checkpointname,
                         const PatchHierarchy& hier);

std::shared_ptr<PatchHierarchy>
ReadCheckpointFile(const std::string checkpointname, DataDescription desc,
                   const CartesianGridGeometry& geometry,
                   const PatchHierarchyOptions& options);

} // namespace cutcell
} // namespace amrex
} // namespace fub

#endif
