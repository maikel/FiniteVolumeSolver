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
#include "fub/AMReX/ForEachFab.hpp"
#include "fub/AMReX/ForEachIndex.hpp"
#include "fub/AMReX/ViewFArrayBox.hpp"
#include "fub/AMReX/output/DebugOutput.hpp"

#include "fub/AMReX/FillCutCellData.hpp"

#include "fub/equations/PerfectGas.hpp"

#include <boost/filesystem.hpp>

namespace fub {
namespace amrex {
namespace cutcell {
namespace {
using MultiCutFabs =
    std::array<std::shared_ptr<::amrex::MultiCutFab>, AMREX_SPACEDIM>;

MultiCutFabs MakeMultiCutFabs_(const ::amrex::BoxArray& ba,
                               const ::amrex::DistributionMapping& dm,
                               const ::amrex::EBFArrayBoxFactory& factory,
                               int ngrow) {
  MultiCutFabs mfabs;
  int dir = 0;
  const int ncomp = 2;
  for (std::shared_ptr<::amrex::MultiCutFab>& mf : mfabs) {
    ::amrex::IntVect unit = ::amrex::IntVect::TheDimensionVector(dir);
    mf = std::make_shared<::amrex::MultiCutFab>(
        ::amrex::convert(ba, unit), dm, ncomp, ngrow,
        factory.getMultiEBCellFlagFab());
    dir += 1;
  }
  return mfabs;
}
} // namespace

PatchLevel::PatchLevel(const PatchLevel& other) = default;

PatchLevel& PatchLevel::operator=(const PatchLevel& other) = default;

PatchLevel::PatchLevel(int level, Duration tp, const ::amrex::BoxArray& ba,
                       const ::amrex::DistributionMapping& dm, int n_components,
                       std::shared_ptr<::amrex::EBFArrayBoxFactory> f,
                       int ngrow)
    : ::fub::amrex::PatchLevel(level, tp, ba, dm, n_components, *f),
      factory(std::move(f)),
      unshielded(MakeMultiCutFabs_(ba, dm, *factory, ngrow)),
      shielded_left(MakeMultiCutFabs_(ba, dm, *factory, ngrow)),
      shielded_right(MakeMultiCutFabs_(ba, dm, *factory, ngrow)),
      doubly_shielded(MakeMultiCutFabs_(ba, dm, *factory, ngrow)) {
  const ::amrex::FabArray<::amrex::EBCellFlagFab>& flags =
      factory->getMultiEBCellFlagFab();
  for (std::size_t d = 0; d < AMREX_SPACEDIM; ++d) {
    unshielded[d]->setVal(0.0);
    shielded_left[d]->setVal(0.0);
    shielded_right[d]->setVal(0.0);
    doubly_shielded[d]->setVal(0.0);
  }
  for (std::size_t d = 0; d < AMREX_SPACEDIM; ++d) {
    const ::amrex::MultiCutFab& betas = *factory->getAreaFrac()[d];
    const ::amrex::BoxArray& ba = unshielded[d]->boxArray();
    const ::amrex::DistributionMapping& dm = unshielded[d]->DistributionMap();
    const int ngrow = unshielded[d]->nGrow();
    ForEachFab(execution::openmp, ba, dm, [&](const ::amrex::MFIter& mfi) {
      if (flags[mfi].getType() == ::amrex::FabType::singlevalued) {
        ::amrex::Box tilebox = mfi.growntilebox(ngrow);
        IndexBox<AMREX_SPACEDIM> face_box = AsIndexBox<AMREX_SPACEDIM>(tilebox);
        PatchDataView<const double, AMREX_SPACEDIM> beta =
            MakePatchDataView(betas[mfi], 0);
        StridedDataView<double, AMREX_SPACEDIM> us =
            MakePatchDataView((*unshielded[d])[mfi], 0).Subview(face_box);
        StridedDataView<double, AMREX_SPACEDIM> sL =
            MakePatchDataView((*shielded_left[d])[mfi], 0).Subview(face_box);
        StridedDataView<double, AMREX_SPACEDIM> sR =
            MakePatchDataView((*shielded_right[d])[mfi], 0).Subview(face_box);
        StridedDataView<double, AMREX_SPACEDIM> ds =
            MakePatchDataView((*doubly_shielded[d])[mfi], 0).Subview(face_box);
        StridedDataView<double, AMREX_SPACEDIM> us_rel =
            MakePatchDataView((*unshielded[d])[mfi], 1).Subview(face_box);
        StridedDataView<double, AMREX_SPACEDIM> sL_rel =
            MakePatchDataView((*shielded_left[d])[mfi], 1).Subview(face_box);
        StridedDataView<double, AMREX_SPACEDIM> sR_rel =
            MakePatchDataView((*shielded_right[d])[mfi], 1).Subview(face_box);
        StridedDataView<double, AMREX_SPACEDIM> ds_rel =
            MakePatchDataView((*doubly_shielded[d])[mfi], 1).Subview(face_box);
        FillCutCellData(us, sL, sR, ds, us_rel, sL_rel, sR_rel, ds_rel, beta,
                        Direction(d));
      }
    });
  }
}

PatchHierarchy::PatchHierarchy(DataDescription desc,
                               const CartesianGridGeometry& geometry,
                               const PatchHierarchyOptions& options)
    : description_{std::move(desc)},
      grid_geometry_{geometry}, options_{options}, patch_level_{},
      patch_level_geometry_{}, registry_{std::make_shared<CounterRegistry>()},
      debug_storage_{std::make_shared<DebugStorage>()} {
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

CutCellData<AMREX_SPACEDIM>
PatchHierarchy::GetCutCellData(int level_number,
                               const ::amrex::MFIter& mfi) const {
  CutCellData<AMREX_SPACEDIM> cutcell_data;
  const PatchLevel& level = GetPatchLevel(level_number);
  cutcell_data.volume_fractions =
      MakePatchDataView(level.factory->getVolFrac()[mfi], 0);
  cutcell_data.volume_centeroid =
      MakePatchDataView(level.factory->getCentroid()[mfi]);
  cutcell_data.boundary_normals =
      MakePatchDataView(level.factory->getBndryNormal()[mfi]);
  cutcell_data.boundary_centeroids =
      MakePatchDataView(level.factory->getBndryCent()[mfi]);
  for (std::size_t d = 0; d < AMREX_SPACEDIM; ++d) {
    cutcell_data.face_fractions[d] =
        MakePatchDataView((*level.factory->getAreaFrac()[d])[mfi], 0);
    cutcell_data.unshielded_fractions[d] =
        MakePatchDataView((*level.unshielded[d])[mfi], 0);
    cutcell_data.shielded_left_fractions[d] =
        MakePatchDataView((*level.shielded_left[d])[mfi], 0);
    cutcell_data.shielded_right_fractions[d] =
        MakePatchDataView((*level.shielded_right[d])[mfi], 0);
    cutcell_data.doubly_shielded_fractions[d] =
        MakePatchDataView((*level.doubly_shielded[d])[mfi], 0);
    cutcell_data.unshielded_fractions_rel[d] =
        MakePatchDataView((*level.unshielded[d])[mfi], 1);
    cutcell_data.shielded_left_fractions_rel[d] =
        MakePatchDataView((*level.shielded_left[d])[mfi], 1);
    cutcell_data.shielded_right_fractions_rel[d] =
        MakePatchDataView((*level.shielded_right[d])[mfi], 1);
    cutcell_data.doubly_shielded_fractions_rel[d] =
        MakePatchDataView((*level.doubly_shielded[d])[mfi], 1);
  }
  return cutcell_data;
}

const std::shared_ptr<::amrex::EBFArrayBoxFactory>&
PatchHierarchy::GetEmbeddedBoundary(int level) const {
  const std::size_t l = static_cast<std::size_t>(level);
  return patch_level_[l].factory;
}

int PatchHierarchy::GetRatioToCoarserLevel(int level, Direction dir) const {
  if (level == 0) {
    return 1;
  }
  return options_.refine_ratio[static_cast<int>(dir)];
}

::amrex::IntVect PatchHierarchy::GetRatioToCoarserLevel(int level) const {
  if (level == 0) {
    return ::amrex::IntVect::TheUnitVector();
  }
  return options_.refine_ratio;
}

const ::amrex::Geometry& PatchHierarchy::GetGeometry(int level) const {
  return patch_level_geometry_[static_cast<std::size_t>(level)];
}

const PatchHierarchyOptions& PatchHierarchy::GetOptions() const noexcept {
  return options_;
}

const CartesianGridGeometry& PatchHierarchy::GetGridGeometry() const noexcept {
  return grid_geometry_;
}

const std::shared_ptr<CounterRegistry>&
PatchHierarchy::GetCounterRegistry() const noexcept {
  return registry_;
}

const std::shared_ptr<DebugStorage>& PatchHierarchy::GetDebugStorage() const
    noexcept {
  return debug_storage_;
}

void PatchHierarchy::SetCounterRegistry(
    std::shared_ptr<CounterRegistry> registry) {
  registry_ = std::move(registry);
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

const DataDescription& PatchHierarchy::GetDataDescription() const noexcept {
  return description_;
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

    const int ngrow = options.ngrow_eb_level_set;

    std::unique_ptr<::amrex::EBFArrayBoxFactory> eb_factory =
        ::amrex::makeEBFabFactory(hierarchy.GetGeometry(lev), ba, dm,
                                  {ngrow, ngrow, ngrow},
                                  ::amrex::EBSupport::full);

    hierarchy.GetPatchLevel(lev) = PatchLevel(
        lev, Duration(time_points[static_cast<std::size_t>(lev)]), ba, dm,
        desc.n_state_components, std::move(eb_factory), ngrow - 1);
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
//
// namespace {
//
//::amrex::RealBox GetProbDomain_(const ::amrex::Geometry& geom,
//                                const ::amrex::Box& box) {
//  const double* dx = geom.CellSize();
//  double base[AMREX_SPACEDIM] = {};
//  geom.CellCenter({AMREX_D_DECL(0, 0, 0)}, base);
//  return ::amrex::RealBox(box, dx, base);
//}
//
//} // namespace

void WriteMatlabData(const std::string& name, const PatchHierarchy& hierarchy,
                     fub::Duration time_point, std::ptrdiff_t cycle_number,
                     MPI_Comm comm) {
  const std::size_t n_level =
      static_cast<std::size_t>(hierarchy.GetNumberOfLevels());
  std::vector<::amrex::FArrayBox> fabs{};
  fabs.reserve(n_level);
  for (std::size_t level = 0; level < n_level; ++level) {
    const int ilvl = static_cast<int>(level);
    const ::amrex::Geometry& level_geom = hierarchy.GetGeometry(ilvl);
    ::amrex::Box domain = level_geom.Domain();
    const ::amrex::MultiFab& level_data = hierarchy.GetPatchLevel(ilvl).data;
    ::amrex::FArrayBox local_fab(domain, level_data.nComp());
    local_fab.setVal(0.0);
    ForEachFab(level_data, [&](const ::amrex::MFIter& mfi) {
      const ::amrex::FArrayBox& patch_data = level_data[mfi];
      local_fab.copy(patch_data);
    });
    int rank = -1;
    ::MPI_Comm_rank(comm, &rank);
    if (rank == 0) {
      ::amrex::FArrayBox& fab = fabs.emplace_back(domain, level_data.nComp());
      fab.setVal(0.0);
      ::MPI_Reduce(local_fab.dataPtr(), fab.dataPtr(),
                   static_cast<int>(local_fab.size()), MPI_DOUBLE, MPI_SUM, 0,
                   comm);
      if (level > 0) {
        for (int comp = 1; comp < level_data.nComp(); ++comp) {
          for (int i = domain.smallEnd(0); i <= domain.bigEnd(0); ++i) {
            for (int j = domain.smallEnd(1); j <= domain.bigEnd(1); ++j) {
              ::amrex::IntVect fine_i{AMREX_D_DECL(i, j, domain.smallEnd(2))};
              ::amrex::IntVect coarse_i = fine_i;
              coarse_i.coarsen(hierarchy.GetRatioToCoarserLevel(ilvl));
              if (fab(fine_i, 0) == 0.0) {
                fab(fine_i, comp) = fabs[level - 1](coarse_i, comp);
              }
            }
          }
        }
        for (int i = domain.smallEnd(0); i <= domain.bigEnd(0); ++i) {
          for (int j = domain.smallEnd(1); j <= domain.bigEnd(1); ++j) {
            ::amrex::IntVect fine_i{AMREX_D_DECL(i, j, domain.smallEnd(2))};
            ::amrex::IntVect coarse_i = fine_i;
            coarse_i.coarsen(hierarchy.GetRatioToCoarserLevel(ilvl));
            if (fab(fine_i, 0) == 0.0) {
              fab(fine_i, 0) = fabs[level - 1](coarse_i, 0);
            }
          }
        }
      }
      if (level == n_level - 1) {
        boost::filesystem::path path(name);
        boost::filesystem::path dir = path.parent_path();
        boost::filesystem::create_directories(dir);
        {
          const ::amrex::Geometry& level_geom = hierarchy.GetGeometry(ilvl);
          std::ofstream out(name);
#if AMREX_SPACEDIM == 2
          out << fmt::format("size = ({}, {}, {})\n", domain.length(0),
                             domain.length(1), fab.nComp());
          out << fmt::format("dx = ({}, {})\n", level_geom.CellSize(0),
                             level_geom.CellSize(1));
          out << fmt::format("x0 = ({}, {})\n",
                             level_geom.CellCenter(domain.smallEnd(0), 0),
                             level_geom.CellCenter(domain.smallEnd(1), 1));
#else
          out << fmt::format("size = ({}, {}, {}, {})\n", domain.length(0),
                             domain.length(1), domain.length(2), fab.nComp());
          out << fmt::format("dx = ({}, {}, {})\n", level_geom.CellSize(0),
                             level_geom.CellSize(1), level_geom.CellSize(2));
          out << fmt::format("x0 = ({}, {}, {})\n",
                             level_geom.CellCenter(domain.smallEnd(0), 0),
                             level_geom.CellCenter(domain.smallEnd(1), 1),
                             level_geom.CellCenter(domain.smallEnd(2), 2));
#endif
          out << fmt::format("t = {}\n", time_point.count());
          out << fmt::format("cycle = {}\n", cycle_number);
          out << fmt::format("data_file = {}.bin\n", path.filename().string());
        }
        // Dump binary data
        std::ofstream bin(name + ".bin", std::ios::binary);
        char* pointer = static_cast<char*>(static_cast<void*>(fab.dataPtr()));
        bin.write(pointer, static_cast<std::streamsize>(fab.size()) *
                               static_cast<std::streamsize>(sizeof(double)));
      }
    } else {
      ::MPI_Reduce(local_fab.dataPtr(), nullptr,
                   static_cast<int>(local_fab.size()), MPI_DOUBLE, MPI_SUM, 0,
                   comm);
    }
  }
}

#if AMREX_SPACEDIM == 3
// void Write2Dfrom3D(std::string name, const PatchHierarchy& hierarchy,
//                   const ::amrex::Box& finest_box, const IdealGasMix<3>& eq,
//                   fub::Duration time_point, std::ptrdiff_t cycle_number,
//                   MPI_Comm comm) {
//  int rank = -1;
//  MPI_Comm_rank(comm, &rank);
//  if (rank == 0) {
//    boost::filesystem::path path(name);
//    boost::filesystem::path dir = path.parent_path();
//    boost::filesystem::create_directories(dir);
//    std::ofstream out(name);
//    if (!out) {
//      throw std::runtime_error("Could not open output file!");
//    }
//    Write2Dfrom3D(&out, hierarchy, finest_box, eq, time_point, cycle_number,
//                  comm);
//  } else {
//    Write2Dfrom3D(nullptr, hierarchy, finest_box, eq, time_point,
//    cycle_number,
//                  comm);
//  }
//}

void Write2Dfrom3D(const std::string& name, const PatchHierarchy& hierarchy,
                   const ::amrex::Box& finest_box,
                   const IdealGasMix<3>& /* eq */, fub::Duration time_point,
                   std::ptrdiff_t cycle_number, MPI_Comm comm) {
  int rank = -1;
  ::MPI_Comm_rank(comm, &rank);
  if (rank == 0) {
    boost::filesystem::path path(name);
    boost::filesystem::path dir = path.parent_path();
    boost::filesystem::create_directories(dir);
  }
  const std::size_t n_level =
      static_cast<std::size_t>(hierarchy.GetNumberOfLevels());
  std::vector<::amrex::Geometry> geoms{};
  std::vector<::amrex::MultiFab> data{};
  std::vector<::amrex::FArrayBox> fabs{};
  std::vector<::amrex::Box> boxes{};
  boxes.reserve(n_level);
  boxes.push_back(finest_box);
  ::amrex::Box box = finest_box;
  for (int level = static_cast<int>(n_level) - 1; level > 0; --level) {
    ::amrex::IntVect refine_ratio = hierarchy.GetRatioToCoarserLevel(level);
    box.coarsen(refine_ratio);
    boxes.push_back(box);
  }
  std::reverse(boxes.begin(), boxes.end());
  geoms.reserve(n_level);
  data.reserve(n_level);
  fabs.reserve(n_level);
  for (std::size_t level = 0; level < n_level; ++level) {
    const int ilvl = static_cast<int>(level);
    ::amrex::Box domain = boxes[level];
    //    ::amrex::RealBox probDomain = GetProbDomain_(level_geom, domain);
    //    ::amrex::Geometry& geom = geoms.emplace_back(domain, &probDomain);
    const ::amrex::MultiFab& level_data = hierarchy.GetPatchLevel(ilvl).data;
    ::amrex::FArrayBox local_fab(domain, level_data.nComp());
    local_fab.setVal(0.0);
    ForEachFab(level_data, [&](const ::amrex::MFIter& mfi) {
      const ::amrex::FArrayBox& patch_data = level_data[mfi];
      const ::amrex::Box box = mfi.tilebox() & domain;
      local_fab.copy(patch_data, box);
    });
    if (rank == 0) {
      ::amrex::FArrayBox& fab = fabs.emplace_back(domain, level_data.nComp());
      fab.setVal(0.0);
      ::MPI_Reduce(local_fab.dataPtr(), fab.dataPtr(),
                   static_cast<int>(local_fab.size()), MPI_DOUBLE, MPI_SUM, 0,
                   comm);
      if (level > 0) {
        for (int comp = 1; comp < level_data.nComp(); ++comp) {
          for (int k = domain.smallEnd(2); k <= domain.bigEnd(2); ++k) {
            for (int j = domain.smallEnd(1); j <= domain.bigEnd(1); ++j) {
              for (int i = domain.smallEnd(0); i <= domain.bigEnd(0); ++i) {
                ::amrex::IntVect fine_i{i, j, k};
                ::amrex::IntVect coarse_i = fine_i;
                coarse_i.coarsen(hierarchy.GetRatioToCoarserLevel(ilvl));
                if (fab(fine_i, 0) == 0.0) {
                  fab(fine_i, comp) = fabs[level - 1](coarse_i, comp);
                }
              }
            }
          }
        }
        for (int k = domain.smallEnd(2); k <= domain.bigEnd(2); ++k) {
          for (int j = domain.smallEnd(1); j <= domain.bigEnd(1); ++j) {
            for (int i = domain.smallEnd(0); i <= domain.bigEnd(0); ++i) {
              ::amrex::IntVect fine_i{i, j, k};
              ::amrex::IntVect coarse_i = fine_i;
              coarse_i.coarsen(hierarchy.GetRatioToCoarserLevel(ilvl));
              if (fab(fine_i, 0) == 0.0) {
                fab(fine_i, 0) = fabs[level - 1](coarse_i, 0);
              }
            }
          }
        }
      }
      if (level == n_level - 1) {
        boost::filesystem::path path(name);
        boost::filesystem::path dir = path.parent_path();
        boost::filesystem::create_directories(dir);
        // Write Header File
        //        {
        const ::amrex::Geometry& level_geom = hierarchy.GetGeometry(ilvl);
        WriteToHDF5(name, fab, level_geom, time_point, cycle_number);
      }
    } else {
      ::MPI_Reduce(local_fab.dataPtr(), nullptr,
                   static_cast<int>(local_fab.size()), MPI_DOUBLE, MPI_SUM, 0,
                   comm);
    }
  }
}
#endif

std::vector<double> GatherStates(
    const PatchHierarchy& hierarchy,
    basic_mdspan<const double, extents<AMREX_SPACEDIM, dynamic_extent>> xs,
    MPI_Comm comm) {
  const int nlevel = hierarchy.GetNumberOfLevels();
  const int finest_level = nlevel - 1;
  const int ncomp = hierarchy.GetDataDescription().n_state_components;
  std::vector<double> buffer(
      static_cast<std::size_t>(xs.extent(1) * ncomp * nlevel));
  mdspan<double, 3> states(buffer.data(), xs.extent(1), ncomp, nlevel);
  for (int level = 0; level < nlevel; ++level) {
    const ::amrex::MultiFab& level_data = hierarchy.GetPatchLevel(level).data;
    const ::amrex::Geometry& level_geom = hierarchy.GetGeometry(level);
    ForEachFab(level_data, [&](const ::amrex::MFIter& mfi) {
      ForEachIndex(mfi.tilebox(), [&](auto... is) {
        double lo[AMREX_SPACEDIM]{};
        double hi[AMREX_SPACEDIM]{};
        const ::amrex::IntVect iv{int(is)...};
        level_geom.LoNode(iv, lo);
        level_geom.HiNode(iv, hi);
        for (int k = 0; k < xs.extent(1); ++k) {
          bool is_in_range = true;
          for (int i = 0; i < AMREX_SPACEDIM; ++i) {
            is_in_range = is_in_range && lo[i] <= xs(i, k) && xs(i, k) < hi[i];
          }
          if (is_in_range) {
            for (int comp = 0; comp < level_data.nComp(); ++comp) {
              states(k, comp, level) = level_data[mfi](iv, comp);
            }
          }
        }
      });
    });
  }
  std::vector<double> global_buffer(buffer.size());
  ::MPI_Allreduce(buffer.data(), global_buffer.data(),
                  static_cast<int>(global_buffer.size()), MPI_DOUBLE, MPI_SUM,
                  comm);
  states = mdspan<double, 3>(global_buffer.data(), xs.extent(1), ncomp, nlevel);
  for (int level = 1; level < nlevel; ++level) {
    for (int comp = 0; comp < ncomp; ++comp) {
      for (int i = 0; i < xs.extent(1); ++i) {
        if (states(i, comp, level) == 0.0) {
          states(i, comp, level) = states(i, comp, level - 1);
        }
      }
    }
  }
  std::vector<double> result(&states(0, 0, finest_level),
                             &states(0, 0, finest_level) +
                                 xs.extent(1) * ncomp);
  return result;
}

} // namespace cutcell
} // namespace amrex
} // namespace fub
