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

#include "fub/AMReX/cutcell/PatchHierarchy.hpp"
#include "fub/AMReX/ViewFArrayBox.hpp"

#include "fub/AMReX/FillCutCellData.hpp"

namespace fub {
namespace amrex {
namespace cutcell {
namespace {
using MultiCutFabs =
    std::array<std::unique_ptr<::amrex::MultiCutFab>, AMREX_SPACEDIM>;

MultiCutFabs MakeMultiCutFabs(const ::amrex::BoxArray& ba,
                              const ::amrex::DistributionMapping& dm,
                              const ::amrex::EBFArrayBoxFactory& factory) {
  MultiCutFabs mfabs;
  int dir = 0;
  for (std::unique_ptr<::amrex::MultiCutFab>& mf : mfabs) {
    ::amrex::IntVect unit = ::amrex::IntVect::TheDimensionVector(dir);
    mf = std::make_unique<::amrex::MultiCutFab>(
        ::amrex::convert(ba, unit), dm, 1, 4, factory.getMultiEBCellFlagFab());
    dir += 1;
  }
  return mfabs;
}
} // namespace

PatchLevel::PatchLevel(const PatchLevel& other)
    : ::fub::amrex::PatchLevel(other), factory(other.factory),
      unshielded(MakeMultiCutFabs(box_array, distribution_mapping, *factory)),
      shielded_left(
          MakeMultiCutFabs(box_array, distribution_mapping, *factory)),
      shielded_right(
          MakeMultiCutFabs(box_array, distribution_mapping, *factory)),
      doubly_shielded(
          MakeMultiCutFabs(box_array, distribution_mapping, *factory)) {
  const ::amrex::MultiFab& alphas = factory->getVolFrac();
  const ::amrex::MultiCutFab& normals = factory->getBndryNormal();
  const ::amrex::MultiCutFab& centeroids = factory->getBndryCent();
  const ::amrex::FabArray<::amrex::EBCellFlagFab>& flags =
      factory->getMultiEBCellFlagFab();
  for (std::size_t d = 0; d < static_cast<std::size_t>(AMREX_SPACEDIM); ++d) {
    const ::amrex::MultiCutFab& betas = *factory->getAreaFrac()[d];
    for (::amrex::MFIter mfi(box_array, distribution_mapping); mfi.isValid();
         ++mfi) {
      if (flags[mfi].getType(mfi.growntilebox(4)) ==
          ::amrex::FabType::singlevalued) {
        CutCellData<AMREX_SPACEDIM> cutcell_data;
        cutcell_data.flags = MakePatchDataView(flags[mfi], 0);
        cutcell_data.volume_fractions = MakePatchDataView(alphas[mfi], 0);
        cutcell_data.face_fractions = MakePatchDataView(betas[mfi], 0);
        cutcell_data.boundary_normals = MakePatchDataView(normals[mfi]);
        cutcell_data.boundary_centeroids = MakePatchDataView(centeroids[mfi]);
        PatchDataView<double, AMREX_SPACEDIM> us =
            MakePatchDataView((*unshielded[d])[mfi], 0);
        PatchDataView<double, AMREX_SPACEDIM> sL =
            MakePatchDataView((*shielded_left[d])[mfi], 0);
        PatchDataView<double, AMREX_SPACEDIM> sR =
            MakePatchDataView((*shielded_right[d])[mfi], 0);
        PatchDataView<double, AMREX_SPACEDIM> ds =
            MakePatchDataView((*doubly_shielded[d])[mfi], 0);
        FillCutCellData(us, sL, sR, ds, cutcell_data,
                        static_cast<Direction>(d));
      }
    }
  }
}

PatchLevel& PatchLevel::operator=(const PatchLevel& other) {
  PatchLevel tmp(other);
  return *this = std::move(tmp);
}

PatchLevel::PatchLevel(int level, Duration tp, const ::amrex::BoxArray& ba,
                       const ::amrex::DistributionMapping& dm, int n_components,
                       std::shared_ptr<::amrex::EBFArrayBoxFactory> f)
    : ::fub::amrex::PatchLevel(level, tp, ba, dm, n_components, *f),
      factory(std::move(f)), unshielded(MakeMultiCutFabs(ba, dm, *factory)),
      shielded_left(MakeMultiCutFabs(ba, dm, *factory)),
      shielded_right(MakeMultiCutFabs(ba, dm, *factory)),
      doubly_shielded(MakeMultiCutFabs(ba, dm, *factory)) {
  const ::amrex::MultiFab& alphas = factory->getVolFrac();
  const ::amrex::MultiCutFab& normals = factory->getBndryNormal();
  const ::amrex::MultiCutFab& centeroids = factory->getBndryCent();
  const ::amrex::FabArray<::amrex::EBCellFlagFab>& flags =
      factory->getMultiEBCellFlagFab();
  for (std::size_t d = 0; d < static_cast<std::size_t>(AMREX_SPACEDIM); ++d) {
    const ::amrex::MultiCutFab& betas = *factory->getAreaFrac()[d];
    for (::amrex::MFIter mfi(ba, dm); mfi.isValid(); ++mfi) {
      if (flags[mfi].getType(mfi.growntilebox(4)) ==
          ::amrex::FabType::singlevalued) {
        CutCellData<AMREX_SPACEDIM> cutcell_data;
        cutcell_data.flags = MakePatchDataView(flags[mfi], 0);
        cutcell_data.volume_fractions = MakePatchDataView(alphas[mfi], 0);
        cutcell_data.face_fractions = MakePatchDataView(betas[mfi], 0);
        cutcell_data.boundary_normals = MakePatchDataView(normals[mfi]);
        cutcell_data.boundary_centeroids = MakePatchDataView(centeroids[mfi]);
        PatchDataView<double, AMREX_SPACEDIM> us =
            MakePatchDataView((*unshielded[d])[mfi], 0);
        PatchDataView<double, AMREX_SPACEDIM> sL =
            MakePatchDataView((*shielded_left[d])[mfi], 0);
        PatchDataView<double, AMREX_SPACEDIM> sR =
            MakePatchDataView((*shielded_right[d])[mfi], 0);
        PatchDataView<double, AMREX_SPACEDIM> ds =
            MakePatchDataView((*doubly_shielded[d])[mfi], 0);
        FillCutCellData(us, sL, sR, ds, cutcell_data,
                        static_cast<Direction>(d));
      }
    }
  }
}

PatchHierarchy::PatchHierarchy(DataDescription desc,
                               const CartesianGridGeometry& geometry,
                               const PatchHierarchyOptions& options)
    : description_{std::move(desc)}, grid_geometry_{geometry},
      options_{options}, patch_level_{}, patch_level_geometry_{} {
  const std::size_t size =
      static_cast<std::size_t>(options.max_number_of_levels);
  patch_level_.resize(size);
  patch_level_geometry_.resize(size);
  ::amrex::IntVect lower{};
  ::amrex::IntVect upper{};
  for (int i = 0; i < AMREX_SPACEDIM; ++i) {
    upper[i] = geometry.cell_dimensions[static_cast<std::size_t>(i)] - 1;
  }
  ::amrex::Box level_box(lower, upper);
  for (::amrex::Geometry& geom : patch_level_geometry_) {
    geom = ::amrex::Geometry(level_box, &grid_geometry_.coordinates, -1,
                             grid_geometry_.periodicity.data());
    level_box.refine(options_.refine_ratio);
  }
}

// CutCellData<AMREX_SPACEDIM>
// PatchHierarchy::GetCutCellData(PatchHandle patch, Direction dir) const {
//   const std::size_t d = static_cast<std::size_t>(dir);
//   const ::amrex::MFIter& mfi = *patch.iterator;
//   CutCellData<AMREX_SPACEDIM> cutcell_data;
//   const PatchLevel& level = GetPatchLevel(patch.level);
//   cutcell_data.dir = dir;
//   cutcell_data.flags =
//       MakePatchDataView(level.factory->getMultiEBCellFlagFab()[mfi], 0);
//   cutcell_data.volume_fractions =
//       MakePatchDataView(level.factory->getVolFrac()[mfi], 0);
//   cutcell_data.face_fractions =
//       MakePatchDataView((*level.factory->getAreaFrac()[d])[mfi], 0);
//   cutcell_data.boundary_normals =
//       MakePatchDataView(level.factory->getBndryNormal()[mfi]);
//   cutcell_data.boundary_centeroids =
//       MakePatchDataView(level.factory->getBndryCent()[mfi]);
//   cutcell_data.unshielded_fractions =
//       MakePatchDataView((*level.unshielded[d])[mfi], 0);
//   cutcell_data.shielded_left_fractions =
//       MakePatchDataView((*level.shielded_left[d])[mfi], 0);
//   cutcell_data.shielded_right_fractions =
//       MakePatchDataView((*level.shielded_right[d])[mfi], 0);
//   cutcell_data.doubly_shielded_fractions =
//       MakePatchDataView((*level.doubly_shielded[d])[mfi], 0);
//   return cutcell_data;
// }

int PatchHierarchy::GetRatioToCoarserLevel(int level, Direction dir) const
    noexcept {
  if (level == 0) {
    return 1;
  }
  return options_.refine_ratio[static_cast<int>(dir)];
}

::amrex::IntVect PatchHierarchy::GetRatioToCoarserLevel(int level) const
    noexcept {
  if (level == 0) {
    return ::amrex::IntVect::TheUnitVector();
  }
  return options_.refine_ratio;
}

const ::amrex::Geometry& PatchHierarchy::GetGeometry(int level) const noexcept {
  return patch_level_geometry_[static_cast<std::size_t>(level)];
}

const PatchHierarchyOptions& PatchHierarchy::GetOptions() const noexcept {
  return options_;
}

const CartesianGridGeometry& PatchHierarchy::GetGridGeometry() const noexcept {
  return grid_geometry_;
}

std::ptrdiff_t PatchHierarchy::GetCycles(int level) const {
  return GetPatchLevel(level).cycles;
}

Duration PatchHierarchy::GetTimePoint(int level) const {
  return GetPatchLevel(level).time_point;
}

int PatchHierarchy::GetNumberOfLevels() const noexcept {
  return static_cast<int>(patch_level_.size());
}

int PatchHierarchy::GetMaxNumberOfLevels() const noexcept {
  return GetOptions().max_number_of_levels;
}

PatchLevel& PatchHierarchy::GetPatchLevel(int level) {
  return patch_level_.at(static_cast<std::size_t>(level));
}

const PatchLevel& PatchHierarchy::GetPatchLevel(int level) const {
  return patch_level_.at(static_cast<std::size_t>(level));
}

void WriteCheckpointFile(const std::string& checkpointname,
                         const PatchHierarchy& hier) {
  const int nlevels = hier.GetNumberOfLevels();
  ::amrex::PreBuildDirectorHierarchy(checkpointname, "Level_", nlevels, true);
  // write Header file
  if (::amrex::ParallelDescriptor::IOProcessor()) {
    const std::string header(checkpointname + "/Header");
    ::amrex::VisMF::IO_Buffer io_buffer(::amrex::VisMF::IO_Buffer_Size);
    std::ofstream hout;
    hout.rdbuf()->pubsetbuf(io_buffer.dataPtr(), io_buffer.size());
    hout.open(header.c_str(), std::ofstream::out | std::ofstream::trunc |
                                  std::ofstream::binary);
    if (!hout.good()) {
      ::amrex::FileOpenFailed(header);
    }

    hout.precision(17);

    // write out title line
    hout << "Checkpoint file for fub::cutcell::PatchHierarchy\n";

    // write out finest_level
    hout << nlevels - 1 << '\n';

    // write out array of istep
    for (int level = 0; level < nlevels; ++level) {
      hout << hier.GetCycles(level) << ' ';
    }
    hout << '\n';

    // write out array of t_new
    for (int level = 0; level < nlevels; ++level) {
      hout << hier.GetTimePoint(level).count() << ' ';
    }
    hout << '\n';

    // write the BoxArray at each level
    for (int level = 0; level < nlevels; ++level) {
      hier.GetPatchLevel(level).data.boxArray().writeOn(hout);
      hout << '\n';
    }
  }

  // write the MultiFab data
  for (int level = 0; level < nlevels; ++level) {
    ::amrex::VisMF::Write(hier.GetPatchLevel(level).data,
                          ::amrex::MultiFabFileFullPrefix(level, checkpointname,
                                                          "Level_", "data"));
  }
}

PatchHierarchy ReadCheckpointFile(const std::string& checkpointname,
                                  DataDescription desc,
                                  const CartesianGridGeometry& geometry,
                                  const PatchHierarchyOptions& options) {
  std::string File(checkpointname + "/Header");
  ::amrex::VisMF::IO_Buffer io_buffer(
      static_cast<std::size_t>(::amrex::VisMF::GetIOBufferSize()));
  ::amrex::Vector<char> fileCharPtr;
  ::amrex::ParallelDescriptor::ReadAndBcastFile(File, fileCharPtr);
  std::string fileCharPtrString(fileCharPtr.dataPtr());
  std::istringstream is(fileCharPtrString, std::istringstream::in);

  PatchHierarchy hierarchy(desc, geometry, options);

  // read in title line
  std::string line;
  std::getline(is, line);

  // read in finest_level
  std::getline(is, line);
  const int finest_level = std::stoi(line);
  std::vector<std::ptrdiff_t> cycles(
      static_cast<std::size_t>(finest_level + 1));
  std::vector<double> time_points(static_cast<std::size_t>(finest_level + 1));

  // read in array of istep
  std::getline(is, line);
  {
    std::istringstream lis(line);
    for (std::size_t i = 0; i <= static_cast<std::size_t>(finest_level); ++i) {
      lis >> cycles[i];
    }
  }

  // read in array of t_new
  std::getline(is, line);
  {
    std::istringstream lis(line);
    for (std::size_t i = 0; i <= static_cast<std::size_t>(finest_level); ++i) {
      lis >> time_points[i];
    }
  }

  for (int lev = 0; lev <= finest_level; ++lev) {
    // read in level 'lev' BoxArray from Header
    ::amrex::BoxArray ba;
    ba.readFrom(is);

    // create a distribution mapping
    ::amrex::DistributionMapping dm{ba, ::amrex::ParallelDescriptor::NProcs()};

    std::unique_ptr<::amrex::EBFArrayBoxFactory> eb_factory =
        ::amrex::makeEBFabFactory(hierarchy.GetGeometry(lev), ba, dm, {4, 4, 4},
                                  ::amrex::EBSupport::full);

    hierarchy.GetPatchLevel(lev) =
        PatchLevel(lev, Duration(time_points[static_cast<std::size_t>(lev)]),
                   ba, dm, desc.n_state_components, std::move(eb_factory));
    hierarchy.GetPatchLevel(lev).cycles = cycles[static_cast<std::size_t>(lev)];
  }

  // read in the MultiFab data
  for (int lev = 0; lev <= finest_level; ++lev) {
    ::amrex::VisMF::Read(
        hierarchy.GetPatchLevel(lev).data,
        ::amrex::MultiFabFileFullPrefix(lev, checkpointname, "Level_", "data"));
  }

  return hierarchy;
}

} // namespace cutcell
} // namespace amrex
} // namespace fub
