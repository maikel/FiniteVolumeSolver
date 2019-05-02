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
#include "fub/grid/AMReX/MultiFab.hpp"
#include "fub/grid/AMReX/PatchHandle.hpp"

#include <AMReX_FluxRegister.H>
#include <AMReX_Geometry.H>

#include <AMReX_ParallelDescriptor.H>
#include <AMReX_PlotFileUtil.H>
#include <AMReX_VisMF.H>

#include <fmt/format.h>

#include <vector>

namespace fub {
namespace amrex {

struct PatchHierarchyOptions {
  int max_number_of_levels{1};
  ::amrex::IntVect refine_ratio{AMREX_D_DECL(2, 2, 2)};
};

/// \brief The PatchLevel represents a distributed grid containing plain
/// simulation data without a ghost cell layer.
///
/// Copying a patch level object will deeply copy the data and creates a new
/// independent patch level. This includes making duplicate objects of box array
/// and distribution mapping and modifying the copy will not affect the original
/// patch level in any way.
struct PatchLevel {
  PatchLevel() = default;
  ~PatchLevel() noexcept = default;

  /// \brief Creates a independent copy of the patch level.
  PatchLevel(const PatchLevel& other);

  /// \brief Create a copy of the other patch level, deallocate old memory and
  /// allocate new memory for the copied data.
  PatchLevel& operator=(const PatchLevel& other);

  /// @{
  /// \brief Moves a patch level without any allocations happening.
  PatchLevel(PatchLevel&& other) noexcept = default;
  PatchLevel& operator=(PatchLevel&& other) = default;
  /// @}


  PatchLevel(int num, Duration tp, const ::amrex::BoxArray& ba,
             const ::amrex::DistributionMapping& dm, int n_components);

  PatchLevel(int num, Duration tp, const ::amrex::BoxArray& ba,
             const ::amrex::DistributionMapping& dm, int n_components,
             const ::amrex::FabFactory<::amrex::FArrayBox>& factory);

  int level_number{};
  Duration time_point{};
  std::ptrdiff_t cycles{};
  ::amrex::BoxArray box_array{};
  ::amrex::DistributionMapping distribution_mapping{};
  ::amrex::MultiFab data{};
};

/// The DataDescription class contains all information which is neccessary to
/// describe the complete and conservative state data of an equation.
struct DataDescription {
  int n_state_components;
  int first_cons_component;
  int n_cons_components;
  int dimension{AMREX_SPACEDIM};
};

class PatchHierarchy {
public:
  PatchHierarchy(DataDescription description,
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

  int GetRatioToCoarserLevel(int level, Direction dir) const noexcept;
  ::amrex::IntVect GetRatioToCoarserLevel(int level) const noexcept;

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
      feedback(handle);
    }
    return feedback;
  }

  template <typename Feedback> double Minimum(int level, Feedback feedback) {
    double global_min = std::numeric_limits<double>::infinity();
    for (::amrex::MFIter mfi(GetPatchLevel(level).data); mfi.isValid(); ++mfi) {
      PatchHandle handle{level, &mfi};
      const double local_min = feedback(handle);
      global_min = std::min(global_min, local_min);
    }
    return global_min;
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
  const auto complete_depths = Depths<Complete<Equation>>(equation);
  int n_comp = 0;
  ForEachVariable<Complete<Equation>>([&n_comp](int depth) { n_comp += depth; },
                                      complete_depths);

  const auto cons_depths = Depths<Conservative<Equation>>(equation);
  int n_cons_comp = 0;
  ForEachVariable<Conservative<Equation>>(
      [&n_cons_comp](int depth) { n_cons_comp += depth; }, cons_depths);

  DataDescription desc;
  desc.n_state_components = n_comp;
  desc.first_cons_component = 0;
  desc.n_cons_components = n_cons_comp;
  desc.dimension = Equation::Rank();
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
    ref_ratio[i] = hier.GetRatioToCoarserLevel(static_cast<int>(i));
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
  ::amrex::WriteMultiLevelPlotfile(plotfilename, nlevels, mf, varnames, geoms,
                                   time_point, level_steps, ref_ratio);
}

void WriteCheckpointFile(const std::string checkpointname,
                         const fub::amrex::PatchHierarchy& hier);

} // namespace amrex
} // namespace fub

#endif
