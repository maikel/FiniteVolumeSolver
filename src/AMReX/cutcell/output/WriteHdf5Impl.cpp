// Copyright (c) 2021 Maikel Nadolski
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

#include "src/AMReX/cutcell/output/WriteHdf5Impl.hpp"
#include "fub/AMReX/ForEachFab.hpp"
#include "fub/AMReX/ForEachIndex.hpp"

namespace fub::amrex::cutcell {
namespace {
template <typename F>
void ForEachInterpolation(::amrex::FArrayBox& dest, const ::amrex::Box& dbox,
                          const ::amrex::FArrayBox& src,
                          const ::amrex::Box& /* sbox */,
                          const ::amrex::IntVect& ratio, ::amrex::SrcComp si,
                          ::amrex::DestComp di, ::amrex::NumComps nc,
                          F&& function) {
  for (int comp = 0; comp < nc.n; ++comp) {
    ForEachIndex(dbox, [&](auto... is) {
      ::amrex::IntVect fine_i{int(is)...};
      ::amrex::IntVect coarse_i = fine_i;
      coarse_i.coarsen(ratio);
      function(dest, fine_i, src, coarse_i, ::amrex::SrcComp(si.i + comp),
               ::amrex::DestComp(di.i + comp));
    });
  }
}

void ConstantInterpolation(::amrex::FArrayBox& dest, const ::amrex::Box& dbox,
                           const ::amrex::FArrayBox& src,
                           const ::amrex::Box& sbox,
                           const ::amrex::IntVect& ratio, ::amrex::SrcComp si,
                           ::amrex::DestComp di, ::amrex::NumComps nc) {
  ForEachInterpolation(
      dest, dbox, src, sbox, ratio, si, di, nc,
      [](::amrex::FArrayBox& fine, const ::amrex::IntVect& fine_i,
         const ::amrex::FArrayBox& coarse, const ::amrex::IntVect& coarse_i,
         ::amrex::SrcComp si, ::amrex::DestComp di) {
        fine(fine_i, di.i) = coarse(coarse_i, si.i);
      });
}

void ConstantInterpolation(::amrex::FArrayBox& dest,
                           const ::amrex::FArrayBox& src,
                           const ::amrex::IntVect& ratio, ::amrex::SrcComp si,
                           ::amrex::DestComp di, ::amrex::NumComps n) {
  const ::amrex::Box dbox = dest.box();
  const ::amrex::Box coarsened = ::amrex::coarsen(dbox, ratio);
  const ::amrex::Box sbox = src.box();
  const ::amrex::Box cbox = coarsened & sbox;
  const ::amrex::Box fbox = ::amrex::refine(cbox, ratio);
  ConstantInterpolation(dest, fbox, src, cbox, ratio, si, di, n);
}
}

void WriteHdf5UnRestricted(
    const std::string& name, const PatchHierarchy& hierarchy,
    span<const std::string> fields) {
  MPI_Comm comm = ::amrex::ParallelContext::CommunicatorAll();
  int rank = -1;
  ::MPI_Comm_rank(comm, &rank);
  if (rank == 0) {
    boost::filesystem::path path(name);
    boost::filesystem::path dir = boost::filesystem::absolute(
        path.parent_path(), boost::filesystem::current_path());
    if (!boost::filesystem::exists(dir)) {
      boost::filesystem::create_directories(dir);
    }
  }
  const std::size_t n_level =
      static_cast<std::size_t>(hierarchy.GetNumberOfLevels());
  std::vector<::amrex::FArrayBox> fabs{};
  fabs.reserve(n_level);
  for (std::size_t level = 0; level < n_level; ++level) {
    const int ilvl = static_cast<int>(level);
    const ::amrex::Geometry& level_geom = hierarchy.GetGeometry(ilvl);
    ::amrex::Box domain = level_geom.Domain();
    const ::amrex::MultiFab& level_data = hierarchy.GetPatchLevel(ilvl).data;
    const ::amrex::MultiFab& volume_fractions =
        hierarchy.GetEmbeddedBoundary(ilvl)->getVolFrac();
    const int ncomp = level_data.nComp();
    ::amrex::FArrayBox local_fab(domain, ncomp + 1);
    local_fab.setVal(0.0);
    // Initialize volume fraction data by using coarser level
    if (level > 0 && rank == 0) {
      FUB_ASSERT(fabs.size() == level);
      ConstantInterpolation(local_fab, fabs[level - 1],
                            hierarchy.GetRatioToCoarserLevel(ilvl),
                            ::amrex::SrcComp(ncomp), ::amrex::DestComp(ncomp),
                            ::amrex::NumComps(1));
    }
    ForEachFab(level_data, [&](const ::amrex::MFIter& mfi) {
      const ::amrex::FArrayBox& patch_data = level_data[mfi];
      const ::amrex::FArrayBox& alpha = volume_fractions[mfi];
      local_fab.copy(patch_data, 0, 0, ncomp);
      local_fab.copy(alpha, 0, ncomp, 1);
    });
    if (rank == 0) {
      ::amrex::FArrayBox& fab = fabs.emplace_back(domain, ncomp + 1);
      fab.setVal(0.0);
      ::MPI_Reduce(local_fab.dataPtr(), fab.dataPtr(),
                   static_cast<int>(local_fab.size()), MPI_DOUBLE, MPI_SUM, 0,
                   comm);
      if (level > 0) {
        for (int comp = 1; comp < ncomp; ++comp) {
          ForEachIndex(domain, [&](auto... is) {
            ::amrex::IntVect fine_i{int(is)...};
            ::amrex::IntVect coarse_i = fine_i;
            coarse_i.coarsen(hierarchy.GetRatioToCoarserLevel(ilvl));
            if (fab(fine_i, 0) == 0.0 && fab(fine_i, ncomp) > 0.0) {
              fab(fine_i, comp) = fabs[level - 1](coarse_i, comp);
            }
          });
        }
        ForEachIndex(domain, [&](auto... is) {
          ::amrex::IntVect fine_i{int(is)...};
          ::amrex::IntVect coarse_i = fine_i;
          coarse_i.coarsen(hierarchy.GetRatioToCoarserLevel(ilvl));
          if (fab(fine_i, 0) == 0.0 && fab(fine_i, ncomp) > 0.0) {
            fab(fine_i, 0) = fabs[level - 1](coarse_i, 0);
          }
        });
      }
      if (level == n_level - 1) {
        const fub::Duration time_point = hierarchy.GetTimePoint();
        const std::ptrdiff_t cycle_number = hierarchy.GetCycles();
        WriteToHDF5(name, fab, level_geom, time_point, cycle_number, fields);
      }
    } else {
      ::MPI_Reduce(local_fab.dataPtr(), nullptr,
                   static_cast<int>(local_fab.size()), MPI_DOUBLE, MPI_SUM, 0,
                   comm);
    }
  }
}

void WriteHdf5RestrictedToBox(
    const std::string& name, const PatchHierarchy& hierarchy,
    const ::amrex::Box& finest_box,
    span<const std::string> fields) {
  MPI_Comm comm = ::amrex::ParallelContext::CommunicatorAll();
  int rank = -1;
  ::MPI_Comm_rank(comm, &rank);
  if (rank == 0) {
    boost::filesystem::path path(name);
    boost::filesystem::path dir = boost::filesystem::absolute(
        path.parent_path(), boost::filesystem::current_path());
    if (!boost::filesystem::exists(dir)) {
      boost::filesystem::create_directories(dir);
    }
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
    const ::amrex::MultiFab& level_data = hierarchy.GetPatchLevel(ilvl).data;
    const ::amrex::MultiFab& volume_fractions =
        hierarchy.GetEmbeddedBoundary(ilvl)->getVolFrac();
    const int ncomp = level_data.nComp();
    ::amrex::FArrayBox local_fab(domain, ncomp + 1);
    local_fab.setVal(0.0);
    ForEachFab(level_data, [&](const ::amrex::MFIter& mfi) {
      const ::amrex::FArrayBox& patch_data = level_data[mfi];
      const ::amrex::FArrayBox& alpha = volume_fractions[mfi];
      const ::amrex::Box box = mfi.tilebox() & domain;
      local_fab.copy(patch_data, box, ::amrex::SrcComp(0), ::amrex::DestComp(0),
                     ::amrex::NumComps(ncomp));
      local_fab.copy(alpha, box, ::amrex::SrcComp(0), ::amrex::DestComp(ncomp),
                     ::amrex::NumComps(1));
    });
    if (rank == 0) {
      ::amrex::FArrayBox& fab = fabs.emplace_back(domain, ncomp + 1);
      fab.setVal(0.0);
      ::MPI_Reduce(local_fab.dataPtr(), fab.dataPtr(),
                   static_cast<int>(local_fab.size()), MPI_DOUBLE, MPI_SUM, 0,
                   comm);
      if (level > 0) {
        for (int comp = 1; comp < level_data.nComp(); ++comp) {
          ForEachIndex(domain, [&](auto... is) {
            ::amrex::IntVect fine_i{int(is)...};
            ::amrex::IntVect coarse_i = fine_i;
            coarse_i.coarsen(hierarchy.GetRatioToCoarserLevel(ilvl));
            if (fab(fine_i, 0) == 0.0 && fab(fine_i, ncomp) > 0.0) {
              fab(fine_i, comp) = fabs[level - 1](coarse_i, comp);
            }
          });
        }
        ForEachIndex(domain, [&](auto... is) {
          ::amrex::IntVect fine_i{int(is)...};
          ::amrex::IntVect coarse_i = fine_i;
          coarse_i.coarsen(hierarchy.GetRatioToCoarserLevel(ilvl));
          if (fab(fine_i, 0) == 0.0 && fab(fine_i, ncomp) > 0.0) {
            fab(fine_i, 0) = fabs[level - 1](coarse_i, 0);
          }
        });
      }
      if (level == n_level - 1) {
        const fub::Duration time_point = hierarchy.GetTimePoint();
        const std::ptrdiff_t cycle_number = hierarchy.GetCycles();
        const ::amrex::Geometry& level_geom = hierarchy.GetGeometry(ilvl);
        WriteToHDF5(name, fab, level_geom, time_point, cycle_number, fields);
      }
    } else {
      ::MPI_Reduce(local_fab.dataPtr(), nullptr,
                   static_cast<int>(local_fab.size()), MPI_DOUBLE, MPI_SUM, 0,
                   comm);
    }
  }
}

} // namespace cutcell