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

#include "fub/AMReX/CartesianGridGeometry.hpp"
#include "fub/Duration.hpp"
#include "fub/Equation.hpp"
#include "fub/Execution.hpp"
#include "fub/equations/IdealGasMix.hpp"
#include "fub/ext/Eigen.hpp"

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

  /// \brief Allocates arrays with specified box array and distribution mapping.
  ///
  /// \param num the refinement level number
  /// \param tp  the time point of the simulation
  /// \param ba  the box array of the distributed array
  /// \param dm  the distribution mapping for the array
  /// \param n_components the number of components of the array
  PatchLevel(int num, Duration tp, const ::amrex::BoxArray& ba,
             const ::amrex::DistributionMapping& dm, int n_components);

  /// Allocates arrays with specified box array and distribution mapping.
  ///
  /// \param num the refinement level number
  /// \param tp  the time point of the simulation
  /// \param ba  the box array of the distributed array
  /// \param dm  the distribution mapping for the array
  /// \param n_components the number of components of the array
  /// \param factory the FAB factory, important with EB.
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

template <typename Equation>
DataDescription MakeDataDescription(const Equation& equation);

/// The PatchHierarchy holds simulation data on multiple refinement levels. It
/// also holds a time stamp for each level.
class PatchHierarchy {
public:
  /// \brief Constructs a PatchHierarchy object which is capable of holding data
  /// described by the secified data description on given geometry extents.
  template <typename Equation>
  PatchHierarchy(const Equation& equation,
                 const CartesianGridGeometry& geometry,
                 const PatchHierarchyOptions& options);

  /// \brief Constructs a PatchHierarchy object which is capable of holding data
  /// described by the secified data description on given geometry extents.
  PatchHierarchy(DataDescription description,
                 const CartesianGridGeometry& geometry,
                 const PatchHierarchyOptions& options);

  [[nodiscard]] const DataDescription& GetDataDescription() const noexcept;

  /// \brief Return some additional patch hierarchy options.
  [[nodiscard]] const PatchHierarchyOptions& GetOptions() const noexcept;

  /// \brief Returns the Grid Geometry which was used to create the hierarchy
  /// with.
  [[nodiscard]] const CartesianGridGeometry& GetGridGeometry() const noexcept;

  [[nodiscard]] std::ptrdiff_t GetCycles(int level = 0) const;

  [[nodiscard]] Duration GetTimePoint(int level = 0) const;

  [[nodiscard]] int GetNumberOfLevels() const noexcept;

  [[nodiscard]] int GetMaxNumberOfLevels() const noexcept;

  [[nodiscard]] int GetRatioToCoarserLevel(int level, Direction dir) const noexcept;

  [[nodiscard]] ::amrex::IntVect GetRatioToCoarserLevel(int level) const noexcept;

  [[nodiscard]] PatchLevel& GetPatchLevel(int level);

  [[nodiscard]] const PatchLevel& GetPatchLevel(int level) const;

  /// \brief Returns a Geometry object for a specified level.
  ///
  /// \param[in] The refinement level number for this geometry obejct.
  [[nodiscard]] const ::amrex::Geometry& GetGeometry(int level) const;

  // Modifiers

  void PushBack(const PatchLevel& level);

  void PushBack(PatchLevel&& level);

  void PopBack();

  [[nodiscard]] span<const ::amrex::EB2::IndexSpace*> GetIndexSpaces();

private:
  DataDescription description_;
  CartesianGridGeometry grid_geometry_;
  PatchHierarchyOptions options_;
  std::vector<PatchLevel> patch_level_;
  std::vector<::amrex::Geometry> patch_level_geometry_;
  std::vector<const ::amrex::EB2::IndexSpace*> index_spaces_;
};

template <typename Equation>
DataDescription MakeDataDescription(const Equation& equation) {
  const auto complete_depths = Depths<Complete<Equation>>(equation);
  int n_comp = 0;
  ForEachVariable([&n_comp](int depth) { n_comp += depth; }, complete_depths);

  const auto cons_depths = Depths<Conservative<Equation>>(equation);
  int n_cons_comp = 0;
  ForEachVariable([&n_cons_comp](int depth) { n_cons_comp += depth; },
                  cons_depths);

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
  boost::mp11::tuple_for_each(Zip(names, StateToTuple(depths)), [&](auto xs) {
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

  template <int Rank>
  void WritePlotFile(const std::string plotfilename,
                     const fub::amrex::PatchHierarchy& hier,
                     const IdealGasMix<Rank>& equation) {
    using Equation = IdealGasMix<Rank>;
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
    boost::mp11::tuple_for_each(Zip(names, StateToTuple(depths)), [&](auto xs) {
      const int ncomp = std::get<1>(xs);
      if (ncomp == 1) {
        varnames.push_back(std::get<0>(xs));
      } else {
        span<const std::string> species = equation.GetReactor().GetSpeciesNames();
        for (int i = 0; i < ncomp; ++i) {
          if (std::get<0>(xs) == std::string{"Species"}) {
            varnames.push_back(species[i]);
          } else {
            varnames.push_back(fmt::format("{}_{}", std::get<0>(xs), i));
          }
        }
      }
    });
    ::amrex::WriteMultiLevelPlotfile(plotfilename, nlevels, mf, varnames, geoms,
                                     time_point, level_steps, ref_ratio);
  }

void WriteCheckpointFile(const std::string checkpointname,
                         const fub::amrex::PatchHierarchy& hier);

template <typename Equation>
PatchHierarchy::PatchHierarchy(const Equation& equation,
                               const CartesianGridGeometry& geometry,
                               const PatchHierarchyOptions& options)
    : PatchHierarchy(MakeDataDescription(equation), geometry, options) {}

} // namespace amrex
} // namespace fub

#endif
