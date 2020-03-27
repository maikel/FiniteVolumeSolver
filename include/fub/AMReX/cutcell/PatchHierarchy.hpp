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
#include <iosfwd>
#include <vector>

namespace fub::amrex::cutcell {

/// This class holds state data arrays for each refinement level of a patch
/// hierarchy.
///
/// This cut-cell version stores some embedded boundary informations in addition
/// to the normal patch level type.
struct PatchLevel : ::fub::amrex::PatchLevel {
  PatchLevel() = default;
  ~PatchLevel() = default;

  /// \brief Deeply copy a patch level on each MPI Rank and recompute shielded
  /// fractions.
  PatchLevel(const PatchLevel&);

  /// \brief Deeply copy a patch level on each MPI Rank and recompute shielded
  /// fractions.
  PatchLevel& operator=(const PatchLevel&);

  /// @{
  /// \brief Moves a patch level without any allocations happening.
  PatchLevel(PatchLevel&&) noexcept = default;
  PatchLevel& operator=(PatchLevel&&) noexcept = default;
  /// @}

  /// \brief Constructs a patch level and computes shielded fractions.
  ///
  /// \param[in] level  the refinement level number in a hierarchy
  /// \param[in] tp  the time point associated with this level
  /// \param[in] ba  the box array describing all boxes of this level
  /// \param[in] dm  the distribution mapping for all boxes of this level
  /// \param[in] factory  the boundary informations for cut cells.
  PatchLevel(int level, Duration tp, const ::amrex::BoxArray& ba,
             const ::amrex::DistributionMapping& dm, int n_components,
             std::shared_ptr<::amrex::EBFArrayBoxFactory> factory, int ngrow);

  using MultiCutFabs =
      std::array<std::unique_ptr<::amrex::MultiCutFab>, AMREX_SPACEDIM>;

  /// \brief This stores the EB factory for this refinement level.
  std::shared_ptr<::amrex::EBFArrayBoxFactory> factory;

  /// \brief Store unshielded face fractions for all faces which touch a cut
  /// cell.
  MultiCutFabs unshielded;

  /// \brief Store singly shielded from left face fractions for all faces which
  /// touch a cut cell.
  MultiCutFabs shielded_left;

  /// \brief Store singly shielded from right face fractions for all faces which
  /// touch a cut cell.
  MultiCutFabs shielded_right;

  /// \brief Store doubly shielded face fractions for all faces which touch a
  /// cut cell.
  MultiCutFabs doubly_shielded;
};

/// This class extents the normal hierarchy options with a pointer to an
/// embedded boundary index space for each possible refinement level.
struct PatchHierarchyOptions : public ::fub::amrex::PatchHierarchyOptions {
  PatchHierarchyOptions() = default;
  PatchHierarchyOptions(const ProgramOptions& options)
      : ::fub::amrex::PatchHierarchyOptions(options) {
    ngrow_eb_level_set =
        GetOptionOr(options, "ngrow_eb_level_set", ngrow_eb_level_set);
    cutcell_load_balance_weight =
        GetOptionOr(options, "cutcell_load_balance_weight", cutcell_load_balance_weight);
    remove_covered_grids =
        GetOptionOr(options, "remove_covered_grids", remove_covered_grids);
  }

  template <typename Log> void Print(Log& log) {
    ::fub::amrex::PatchHierarchyOptions::Print(log);
    BOOST_LOG(log) << " - ngrow_eb_level_set = " << ngrow_eb_level_set;
    BOOST_LOG(log) << " - cutcell_load_balance_weight = " << cutcell_load_balance_weight;
    BOOST_LOG(log) << " - remove_covered_grids = " << remove_covered_grids;
  }

  std::vector<const ::amrex::EB2::IndexSpace*> index_spaces{};
  int ngrow_eb_level_set{5};
  double cutcell_load_balance_weight{6.0};
  bool remove_covered_grids{true};
};

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

  /// @{
  /// \name Member Acccess

  const DataDescription& GetDataDescription() const noexcept;

  /// \brief Return some additional patch hierarchy options.
  const PatchHierarchyOptions& GetOptions() const noexcept;

  /// \brief Returns the Grid Geometry which was used to create the hierarchy
  /// with.
  const CartesianGridGeometry& GetGridGeometry() const noexcept;

  [[nodiscard]] const std::shared_ptr<CounterRegistry>&
  GetCounterRegistry() const noexcept;

  void SetCounterRegistry(std::shared_ptr<CounterRegistry> registry);
  /// @}

  /// @{
  /// \name Observers

  int GetNumberOfLevels() const noexcept;

  int GetMaxNumberOfLevels() const noexcept;
  /// @}

  /// @{
  /// \name Level specific Access

  /// \brief Returns the number of cycles done for a specific level.
  ///
  /// \param[in] level  the refinement level number
  std::ptrdiff_t GetCycles(int level = 0) const;

  /// \brief Returns the time point of a specified level number.
  ///
  /// \param[in] level  the refinement level number
  Duration GetTimePoint(int level = 0) const;

  /// \brief Returns the ratio from fine to coarse for specified fine level
  /// number and direction.
  ///
  /// \param[in] level  the fine refinement level number
  /// \param[in] dir  the direction
  int GetRatioToCoarserLevel(int level, Direction dir) const;

  ::amrex::IntVect GetRatioToCoarserLevel(int level) const;

  const ::amrex::Geometry& GetGeometry(int level) const;

  PatchLevel& GetPatchLevel(int level);

  const PatchLevel& GetPatchLevel(int level) const;

  const std::shared_ptr<::amrex::EBFArrayBoxFactory>&
  GetEmbeddedBoundary(int level) const;

  CutCellData<AMREX_SPACEDIM> GetCutCellData(int level,
                                             const ::amrex::MFIter& mfi) const;

  /// @}

private:
  DataDescription description_;
  CartesianGridGeometry grid_geometry_;
  PatchHierarchyOptions options_;
  std::vector<PatchLevel> patch_level_;
  std::vector<::amrex::Geometry> patch_level_geometry_;
  std::shared_ptr<CounterRegistry> registry_;
};

template <typename Equation>
PatchHierarchy::PatchHierarchy(const Equation& equation,
                               const CartesianGridGeometry& geometry,
                               const PatchHierarchyOptions& options)
    : PatchHierarchy(MakeDataDescription(equation), geometry, options) {}

template <typename Equation>
void WritePlotFile(const std::string& plotfilename, const PatchHierarchy& hier,
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
  ::amrex::EB_WriteMultiLevelPlotfile(plotfilename, nlevels, mf, varnames,
                                      geoms, time_point, level_steps,
                                      ref_ratio);
}

template <int Rank>
void WritePlotFile(const std::string& plotfilename, const PatchHierarchy& hier,
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
  ::amrex::EB_WriteMultiLevelPlotfile(plotfilename, nlevels, mf, varnames,
                                      geoms, time_point, level_steps,
                                      ref_ratio);
}

void WriteCheckpointFile(const std::string& checkpointname,
                         const PatchHierarchy& hier);

PatchHierarchy ReadCheckpointFile(const std::string& checkpointname,
                                  DataDescription desc,
                                  const CartesianGridGeometry& geometry,
                                  const PatchHierarchyOptions& options);

// void Write2Dfrom3D(std::ostream* out, const PatchHierarchy& hierarchy, const
// ::amrex::Box& finest_box,
//                 const IdealGasMix<3>& eq, fub::Duration time_point,
//                 std::ptrdiff_t cycle_number, MPI_Comm comm);

void WriteMatlabData(const std::string& name, const PatchHierarchy& hierarchy,
                     fub::Duration time_point, std::ptrdiff_t cycle_number,
                     MPI_Comm comm);

void Write2Dfrom3D(const std::string& name, const PatchHierarchy& hierarchy,
                   const ::amrex::Box& finest_box, const IdealGasMix<3>& eq,
                   fub::Duration time_point, std::ptrdiff_t cycle_number,
                   MPI_Comm comm);

std::vector<double> GatherStates(
    const PatchHierarchy& hierarchy,
    basic_mdspan<const double, extents<AMREX_SPACEDIM, dynamic_extent>> xs,
    MPI_Comm comm);

} // namespace fub::amrex::cutcell

#endif
