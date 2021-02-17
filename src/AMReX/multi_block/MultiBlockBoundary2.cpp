// Copyright (c) 2020 Maikel Nadolski
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

#include "fub/AMReX/multi_block/MultiBlockBoundary2.hpp"

#include "fub/AMReX/GriddingAlgorithm.hpp"
#include "fub/AMReX/cutcell/GriddingAlgorithm.hpp"
#include "fub/AMReX/multi_block/MultiBlockGriddingAlgorithm2.hpp"

#include "fub/AMReX/cutcell/boundary_condition/ReflectiveBoundary.hpp"

#include "fub/AMReX/ForEachFab.hpp"
#include "fub/AMReX/ForEachIndex.hpp"
#include "fub/AMReX/ViewFArrayBox.hpp"
#include "fub/ForEach.hpp"

#include "fub/equations/ideal_gas_mix/mechanism/Burke2012.hpp"

#include "fub/ext/Log.hpp"

#include <AMReX_EB2.H>

#include <boost/log/common.hpp>
#include <boost/log/trivial.hpp>

#include <mpi.h>

namespace fub::amrex {
namespace {
::amrex::IntVect ReduceDimension(const ::amrex::IntVect& idx, Direction dir) {
  ::amrex::IntVect ld = ::amrex::IntVect::TheZeroVector();
  ld[int(dir)] = idx[int(dir)];
  return ld;
}

::amrex::Box ReduceDimension(const ::amrex::Box& mirror,
                             Direction dir) noexcept {
  ::amrex::IntVect lower = ReduceDimension(mirror.smallEnd(), dir);
  ::amrex::IntVect upper = ReduceDimension(mirror.bigEnd(), dir);
  return ::amrex::Box{lower, upper};
}

double CellVolume(const ::amrex::Geometry& geom) {
  return AMREX_D_TERM(geom.CellSize(0), *geom.CellSize(1), *geom.CellSize(2));
}

::amrex::FArrayBox ReduceTotalVolume(const cutcell::PatchHierarchy& hierarchy,
                                     int level, const ::amrex::Box& mirror,
                                     Direction dir) {
  const ::amrex::Box ld_mirror = ReduceDimension(mirror, dir);
  const ::amrex::EBFArrayBoxFactory& factory =
      *hierarchy.GetEmbeddedBoundary(level);
  const ::amrex::MultiFab& alphas = factory.getVolFrac();
  ::amrex::FArrayBox local_reduced_volumes(ld_mirror);
  local_reduced_volumes.setVal(0.0);
  const double dx = CellVolume(hierarchy.GetGeometry(level));
  ForEachFab(alphas, [&](const ::amrex::MFIter& mfi) {
    const ::amrex::FArrayBox& alpha = alphas[mfi];
    IndexBox<AMREX_SPACEDIM> section =
        AsIndexBox<AMREX_SPACEDIM>(mirror & mfi.tilebox());
    ForEachIndex(section, [&](auto... is) {
      ::amrex::IntVect index{static_cast<int>(is)...};
      ::amrex::IntVect ld_index{AMREX_D_DECL(index[0], 0, 0)};
      const double frac = alpha(index);
      local_reduced_volumes(ld_index) += frac * dx;
    });
  });
  ::amrex::FArrayBox global_reduced_volumes(ld_mirror);
  ::MPI_Allreduce(
      local_reduced_volumes.dataPtr(), global_reduced_volumes.dataPtr(),
      int(local_reduced_volumes.size()), MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  return global_reduced_volumes;
}

void IntegrateConserativeStates(
    const PatchDataView<double, AMREX_SPACEDIM + 1>& dest,
    const PatchDataView<const double, AMREX_SPACEDIM>& total_volumes,
    const PatchDataView<const double, AMREX_SPACEDIM + 1>& src,
    const PatchDataView<const double, AMREX_SPACEDIM>& volumes,
    const ::amrex::Box& box, double cell_volume) {
  const int n_comps = static_cast<int>(dest.Extent(AMREX_SPACEDIM));
  for (int i = 0; i < n_comps; ++i) {
    ForEachIndex(AsIndexBox<AMREX_SPACEDIM>(box), [&](auto... is) {
      Index<AMREX_SPACEDIM> index{is...};
      Index<AMREX_SPACEDIM> lowdim_index{index[0]};
      if (volumes(index) > 0.0) {
        dest(AMREX_D_DECL(index[0], 0, 0), i) +=
            (volumes(index) * cell_volume / total_volumes(lowdim_index)) *
            src(is..., i);
      }
    });
  }
}

::amrex::FArrayBox AverageConservativeHierarchyStates(
    const ::amrex::MultiFab& datas, const cutcell::PatchHierarchy& hierarchy,
    int level, const ::amrex::Box& mirror, Direction dir) {
  const ::amrex::EBFArrayBoxFactory& eb = *hierarchy.GetEmbeddedBoundary(level);
  const ::amrex::MultiFab& volumes = eb.getVolFrac();
  const int n_comps = hierarchy.GetDataDescription().n_cons_components;
  const ::amrex::FArrayBox total_volume =
      ReduceTotalVolume(hierarchy, level, mirror, dir);
  ::amrex::FArrayBox fab(total_volume.box(), n_comps);
  fab.setVal(0.0);
  const double cell_volume = CellVolume(hierarchy.GetGeometry(level));
  PatchDataView<double, AMREX_SPACEDIM + 1> fabv = MakePatchDataView(fab);
  PatchDataView<const double, AMREX_SPACEDIM> tot_vol =
      MakePatchDataView(total_volume, 0);
  ForEachFab(datas, [&](const ::amrex::MFIter& mfi) {
    const ::amrex::FArrayBox& alpha = volumes[mfi];
    const ::amrex::FArrayBox& data = datas[mfi];
    PatchDataView<const double, AMREX_SPACEDIM + 1> src =
        MakePatchDataView(data);
    PatchDataView<const double, AMREX_SPACEDIM> vol =
        MakePatchDataView(alpha, 0);
    ::amrex::Box box = mfi.tilebox() & mirror;
    IntegrateConserativeStates(fabv, tot_vol, src, vol, box, cell_volume);
  });
  ::amrex::FArrayBox global_fab(fab.box(), fab.nComp());
  global_fab.setVal(0.0);
  ::MPI_Allreduce(fab.dataPtr(), global_fab.dataPtr(), int(fab.size()),
                  MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  return global_fab;
}

::amrex::FArrayBox
AverageConservativeHierarchyStates(const cutcell::PatchHierarchy& hierarchy,
                                   int level, const ::amrex::Box& mirror,
                                   Direction dir) {
  const ::amrex::MultiFab& datas = hierarchy.GetPatchLevel(level).data;
  return AverageConservativeHierarchyStates(datas, hierarchy, level, mirror,
                                            dir);
}

::amrex::FArrayBox
AverageConservativeHierarchyStates(const cutcell::IntegratorContext& context,
                                   int level, const ::amrex::Box& mirror,
                                   Direction dir) {
  const cutcell::PatchHierarchy& hierarchy = context.GetPatchHierarchy();
  const ::amrex::MultiFab& datas = context.GetScratch(level);
  return AverageConservativeHierarchyStates(datas, hierarchy, level, mirror,
                                            dir);
}

void InterpolateStates_Greater(span<double> dest, double dest_dx,
                               span<const double> src, double src_dx) {
  FUB_ASSERT(dest_dx > src_dx);
  std::ptrdiff_t src_index = 0LL;
  double last_remainder = 0.0;
  const double full_weight = src_dx / dest_dx;
  for (std::ptrdiff_t j = 0; j < dest.size(); ++j) {
    dest[j] = 0.0;
    double src_sum = 0.0;
    if (last_remainder > 0.0) {
      const double rest_dx = src_dx - last_remainder;
      const double rest_weight = rest_dx / dest_dx;
      dest[j] += rest_weight * src[src_index];
      src_index += 1;
      src_sum += rest_dx;
    }
    while (src_sum + src_dx <= dest_dx) {
      dest[j] += full_weight * src[src_index];
      src_sum += src_dx;
      src_index += 1;
    }
    const double remainder = dest_dx - src_sum;
    const double remainder_weight = remainder / dest_dx;
    FUB_ASSERT(remainder < src_dx);
    dest[j] += remainder_weight * src[src_index];
    last_remainder = remainder;
  }
}

void InterpolateStates_Less(span<double> dest, double dest_dx,
                            span<const double> src, double src_dx) {
  FUB_ASSERT(dest_dx < src_dx);
  for (std::ptrdiff_t i = 0; i < dest.size(); ++i) {
    const double x_lower = dest_dx * static_cast<double>(i);
    const double x_upper = dest_dx * static_cast<double>(i + 1);
    const std::ptrdiff_t j0 = static_cast<std::ptrdiff_t>(x_lower / src_dx);
    const std::ptrdiff_t j1 = static_cast<std::ptrdiff_t>(x_upper / src_dx);
    if (j0 == j1) {
      dest[i] = src[j0];
    } else {
      const double weight1 = ::fmod(x_upper, src_dx);
      const double weight0 = dest_dx - weight1;
      if (weight1 > 0.0 && j1 < src.size()) {
        dest[i] = weight0 / dest_dx * src[j0] + weight1 / dest_dx * src[j1];
      } else {
        dest[i] = src[j0];
      }
    }
  }
}

void InterpolateStates(::amrex::FArrayBox& dest, double dest_dx,
                       const ::amrex::FArrayBox& src, double src_dx, int side) {
  const int n_comps = dest.nComp();
  FUB_ASSERT(n_comps == src.nComp());
  if (dest_dx > src_dx) {
    // Coarsen Data
    for (int i = 0; i < n_comps; ++i) {
      span<double> destv = MakePatchDataView(dest, i).Span();
      span<const double> srcv = MakePatchDataView(src, i).Span();
      if (side == 1) {
        std::vector<double> rsrc(srcv.begin(), srcv.end());
        std::reverse(rsrc.begin(), rsrc.end());
        InterpolateStates_Greater(destv, dest_dx, rsrc, src_dx);
        std::reverse(destv.begin(), destv.end());
      } else {
        InterpolateStates_Greater(destv, dest_dx, srcv, src_dx);
      }
    }
  } else if (dest_dx < src_dx) {
    // Refine Data
    for (int i = 0; i < n_comps; ++i) {
      span<double> destv = MakePatchDataView(dest, i).Span();
      span<const double> srcv = MakePatchDataView(src, i).Span();
      if (side == 1) {
        std::vector<double> rsrc(srcv.begin(), srcv.end());
        std::reverse(rsrc.begin(), rsrc.end());
        InterpolateStates_Less(destv, dest_dx, rsrc, src_dx);
        std::reverse(destv.begin(), destv.end());
      } else {
        InterpolateStates_Less(destv, dest_dx, srcv, src_dx);
      }
    }
  } else {
    dest.copy(src);
  }
}

::amrex::FArrayBox AllgatherMirrorData(const ::amrex::MultiFab& mf,
                                       const ::amrex::Box& mirror,
                                       int first_comp, int n_comps) {
  ::amrex::FArrayBox local_fab(mirror, n_comps);
  local_fab.setVal(0.0);
  for (::amrex::MFIter mfi(mf); mfi.isValid(); ++mfi) {
    const ::amrex::FArrayBox& data = mf[mfi];
    const ::amrex::Box subbox = mfi.tilebox() & mirror;
    local_fab.plus(data, subbox, first_comp, 0, n_comps);
  }
  ::amrex::FArrayBox global_fab(local_fab.box(), local_fab.nComp());
  ::MPI_Allreduce(local_fab.dataPtr(), global_fab.dataPtr(),
                  int(local_fab.size()), MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  return global_fab;
}

::amrex::FArrayBox AllgatherMirrorData(const PatchHierarchy& hierarchy,
                                       int level, const ::amrex::Box& mirror) {
  const ::amrex::MultiFab& data = hierarchy.GetPatchLevel(level).data;
  const int n_comps = hierarchy.GetDataDescription().n_cons_components;
  return AllgatherMirrorData(data, mirror, 0, n_comps);
}

::amrex::FArrayBox AllgatherMirrorData(const IntegratorContext& context,
                                       int level, const ::amrex::Box& mirror) {
  const ::amrex::MultiFab& datas = context.GetScratch(level);
  const int n_comps =
      context.GetPatchHierarchy().GetDataDescription().n_cons_components;
  return AllgatherMirrorData(datas, mirror, 0, n_comps);
}

::amrex::Box MakeMirrorBox(const ::amrex::Box& box, int width, Direction dir,
                           int side, int level, const ::amrex::IntVect& ratio) {
  const int d = static_cast<int>(dir);
  ::amrex::Box result = box;
  for (int l = 0; l < level; ++l) {
    result.refine(ratio);
  }
  if (side == 1) {
    result.setSmall(d, result.bigEnd(d) - width + 1);
  } else {
    result.setBig(d, result.smallEnd(d) + width - 1);
  }
  return result;
}

::amrex::Box MakeGhostBox(const ::amrex::Box& box, int width, Direction dir,
                          int side) {
  const int d = static_cast<int>(dir);
  ::amrex::Box result = box;
  if (side == 1) {
    result.setSmall(d, box.bigEnd(d) + 1);
    result.setBig(d, box.bigEnd(d) + width);
  } else {
    result.setBig(d, box.smallEnd(d) - 1);
    result.setSmall(d, box.smallEnd(d) - width);
  }
  return result;
}

int Flip(int side) { return (side == 0) * 1 + (side != 0) * 0; }
//
// template <typename GriddingAlgorithm>
// int FindLevel(const ::amrex::Geometry& geom,
//              const GriddingAlgorithm& gridding) {
//  for (int level = 0; level <
//  gridding.GetPatchHierarchy().GetNumberOfLevels();
//       ++level) {
//    if (geom.Domain() ==
//        gridding.GetPatchHierarchy().GetGeometry(level).Domain()) {
//      return level;
//    }
//  }
//  return -1;
//}
} // namespace

AnyMultiBlockBoundary::AnyMultiBlockBoundary(const AnyMultiBlockBoundary& other)
    : impl_(other.impl_->Clone()), plenum_mirror_box_(other.plenum_mirror_box_),
      tube_mirror_box_(other.tube_mirror_box_),
      plenum_mirror_data_(std::make_unique<::amrex::FArrayBox>(
          other.plenum_mirror_data_->box(),
          other.plenum_mirror_data_->nComp())),
      tube_ghost_data_(std::make_unique<::amrex::FArrayBox>(
          other.tube_ghost_data_->box(), other.tube_ghost_data_->nComp())),
      tube_mirror_data_(std::make_unique<::amrex::FArrayBox>(
          other.tube_mirror_data_->box(), other.tube_mirror_data_->nComp())),
      plenum_ghost_data_(std::make_unique<::amrex::FArrayBox>(
          other.plenum_ghost_data_->box(), other.plenum_ghost_data_->nComp())),
      dir_(other.dir_),
      side_(other.side_), level_{other.level_}, gcw_{other.gcw_} {
  plenum_mirror_data_->copy(*other.plenum_mirror_data_);
  tube_ghost_data_->copy(*other.tube_ghost_data_);
  tube_mirror_data_->copy(*other.tube_mirror_data_);
  plenum_ghost_data_->copy(*other.plenum_ghost_data_);
}

void AnyMultiBlockBoundary::Initialize(
    const MultiBlockGriddingAlgorithm2& gridding,
    const BlockConnection& connection, int gcw, int level) {
  const std::ptrdiff_t pid = static_cast<std::ptrdiff_t>(connection.plenum.id);
  const cutcell::PatchHierarchy& plenum =
      gridding.GetPlena()[pid]->GetPatchHierarchy();
  const std::ptrdiff_t tid = static_cast<std::ptrdiff_t>(connection.tube.id);
  const PatchHierarchy& tube = gridding.GetTubes()[tid]->GetPatchHierarchy();
  const double plenum_dx =
      plenum.GetGeometry(level_).CellSize(static_cast<int>(dir_));
  const double tube_dx =
      tube.GetGeometry(level_).CellSize(static_cast<int>(dir_));
  // Allocate mirror data containing the plenum data (reduced to 1d)
  {
    const int mirror_width =
        static_cast<int>(std::ceil((gcw * tube_dx) / plenum_dx));
    plenum_mirror_box_ =
        MakeMirrorBox(connection.plenum.mirror_box, mirror_width, dir_, side_,
                      level_, plenum.GetRatioToCoarserLevel(level_));
    const int ncons = plenum.GetDataDescription().n_cons_components;
    const int ncomp = plenum.GetDataDescription().n_state_components;
    plenum_mirror_data_ = std::make_unique<::amrex::FArrayBox>(
        ReduceDimension(MakeMirrorBox(connection.plenum.mirror_box, gcw, dir_,
                                      side_, level_,
                                      plenum.GetRatioToCoarserLevel(level_)),
                        dir_),
        ncons);
    plenum_ghost_data_ = std::make_unique<::amrex::FArrayBox>(
        ReduceDimension(MakeGhostBox(plenum_mirror_box_, gcw, dir_, side_),
                        dir_),
        ncomp);
  }
  // Allocate mirror data for the tube data
  {
    const int mirror_width =
        static_cast<int>(std::ceil((gcw * plenum_dx) / tube_dx));
    tube_mirror_box_ =
        MakeMirrorBox(connection.tube.mirror_box, mirror_width, dir_,
                      Flip(side_), level_, tube.GetRatioToCoarserLevel(level));
    const int ncons = tube.GetDataDescription().n_cons_components;
    const int ncomp = tube.GetDataDescription().n_state_components;
    tube_mirror_data_ = std::make_unique<::amrex::FArrayBox>(
        MakeMirrorBox(connection.tube.mirror_box, gcw, dir_, Flip(side_),
                      level_, tube.GetRatioToCoarserLevel(level_)),
        ncons);
    tube_ghost_data_ = std::make_unique<::amrex::FArrayBox>(
        MakeGhostBox(tube_mirror_box_, gcw, dir_, Flip(side_)), ncomp);
  }
  // Fill internal mirror data
  ComputeBoundaryData(plenum, tube);
}

void AnyMultiBlockBoundary::PreAdvanceHierarchy(
    const MultiBlockGriddingAlgorithm2& grid) {
  const int plenum_id = connection_.plenum.id;
  const int tube_id = connection_.tube.id;
  const cutcell::PatchHierarchy& plenum =
      grid.GetPlena()[plenum_id]->GetPatchHierarchy();
  const PatchHierarchy& tube = grid.GetTubes()[tube_id]->GetPatchHierarchy();
  ComputeBoundaryData(plenum, tube);
}

template <int R> using Conservative = ::fub::Conservative<IdealGasMix<R>>;
template <int R> using Complete = ::fub::Complete<IdealGasMix<R>>;

void AnyMultiBlockBoundary::ComputeBoundaryData(
    const cutcell::PatchHierarchy& plenum, const PatchHierarchy& tube) {
  const int d = static_cast<int>(dir_);
  const double plenum_dx = plenum.GetGeometry(level_).CellSize(d);
  const double tube_dx = tube.GetGeometry(level_).CellSize(d);

  //////////////////////////////////////////////////////////////////////////////
  // Integrate plenum states over mirror volume and fill ghost cells of tube

  {
    // Integrate over mirror volume
    const ::amrex::FArrayBox plenum_data = AverageConservativeHierarchyStates(
        plenum, level_, plenum_mirror_box_, dir_);

    // Interpolate between grid cell sizes.
    InterpolateStates(*plenum_mirror_data_, tube_dx, plenum_data, plenum_dx,
                      side_);

    // Transform high dimensional states into low dimensional ones.
    // Store low dimensional states as the reference in the ghost cell region.
    impl_->FillTubeGhostLayer(*tube_ghost_data_, *plenum_mirror_data_);
  }

  {
    ::amrex::FArrayBox tube_data =
        AllgatherMirrorData(tube, level_, tube_mirror_box_);

    // Interpolate between grid cell sizes.
    InterpolateStates(*tube_mirror_data_, plenum_dx, tube_data, tube_dx,
                      Flip(side_));

    // Transform low dimensional states into high dimensional ones.
    impl_->FillPlenumGhostLayer(*plenum_ghost_data_, *tube_mirror_data_);
  }
}

void AnyMultiBlockBoundary::ComputeBoundaryData(
    const cutcell::IntegratorContext& plenum, const IntegratorContext& tube) {
  ComputeBoundaryDataForTube(plenum, tube);
  ComputeBoundaryDataForPlenum(plenum, tube);
}

/// \brief Integrate plenum states over mirror volume and fill ghost cells of
/// tube
void AnyMultiBlockBoundary::ComputeBoundaryDataForTube(
    const cutcell::IntegratorContext& plenum, const IntegratorContext& tube) {
  const int d = static_cast<int>(dir_);
  const double plenum_dx = plenum.GetGeometry(level_).CellSize(d);
  const double tube_dx = tube.GetGeometry(level_).CellSize(d);

  // Integrate over mirror volume
  const ::amrex::FArrayBox plenum_data = AverageConservativeHierarchyStates(
      plenum, level_, plenum_mirror_box_, dir_);

  // Interpolate between grid cell sizes.
  InterpolateStates(*plenum_mirror_data_, tube_dx, plenum_data, plenum_dx,
                    side_);

  // Transform high dimensional states into low dimensional ones.
  // Store low dimensional states as the reference in the ghost cell region.
  impl_->FillTubeGhostLayer(*tube_ghost_data_, *plenum_mirror_data_);
}

void AnyMultiBlockBoundary::ComputeBoundaryDataForPlenum(
    const cutcell::IntegratorContext& plenum, const IntegratorContext& tube) {
  const int d = static_cast<int>(dir_);
  const double plenum_dx = plenum.GetGeometry(level_).CellSize(d);
  const double tube_dx = tube.GetGeometry(level_).CellSize(d);

  ::amrex::FArrayBox tube_data =
      AllgatherMirrorData(tube, level_, tube_mirror_box_);

  // Interpolate between grid cell sizes.
  InterpolateStates(*tube_mirror_data_, plenum_dx, tube_data, tube_dx,
                    Flip(side_));

  // Transform low dimensional states into high dimensional ones.
  impl_->FillPlenumGhostLayer(*plenum_ghost_data_, *tube_mirror_data_);
}

void AnyMultiBlockBoundary::FillBoundary(::amrex::MultiFab& mf,
                                         const GriddingAlgorithm&, int level) {
  if (level == level_) {
    ::amrex::Box box = tube_ghost_data_->box();
    ForEachFab(execution::seq, mf, [&](const ::amrex::MFIter& mfi) {
      const ::amrex::Box b = mfi.growntilebox() & box;
      if (!b.isEmpty()) {
        for (int n = 0; n < mf.nComp(); ++n) {
          ForEachIndex(b, [&](auto... is) {
            const ::amrex::IntVect iv{int(is)...};
            if ((*tube_ghost_data_)(iv, 0) > 0.0) {
              mf[mfi](iv, n) = (*tube_ghost_data_)(iv, n);
            }
          });
        }
      }
    });
  }
}

void AnyMultiBlockBoundary::FillBoundary(::amrex::MultiFab& mf,
                                         const cutcell::GriddingAlgorithm&,
                                         int level) {
  if (level == level_) {
    ::amrex::Box ghost_box =
        MakeGhostBox(plenum_mirror_box_, gcw_, dir_, side_);
    ForEachFab(execution::openmp, mf, [&](const ::amrex::MFIter& mfi) {
      const ::amrex::Box box = mfi.growntilebox() & ghost_box;
      if (!box.isEmpty()) {
        ::amrex::FArrayBox& fab = mf[mfi];
        for (int c = 0; c < mf.nComp(); ++c) {
          ForEachIndex(box, [&](auto... is) {
            const ::amrex::IntVect dest{int(is)...};
            ::amrex::IntVect src = ReduceDimension(dest, dir_);
            fab(dest, c) = (*plenum_ghost_data_)(src, c);
          });
        }
      }
    });
  }
}

} // namespace fub::amrex
