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

#include "fub/AMReX/multi_block/MultiBlockBoundary.hpp"

#include "fub/AMReX/GriddingAlgorithm.hpp"
#include "fub/AMReX/cutcell/GriddingAlgorithm.hpp"
#include "fub/AMReX/multi_block/MultiBlockGriddingAlgorithm.hpp"

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

::amrex::FArrayBox
AverageConservativeHierarchyStates(const cutcell::PatchHierarchy& hierarchy,
                                   int level, const ::amrex::Box& mirror,
                                   Direction dir) {
  const ::amrex::EBFArrayBoxFactory& eb = *hierarchy.GetEmbeddedBoundary(level);
  const ::amrex::MultiFab& datas = hierarchy.GetPatchLevel(level).data;
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

void ReduceStateDimension(
    Complete<IdealGasMix<1>>& dest, IdealGasMix<1>& equation,
    const Conservative<IdealGasMix<AMREX_SPACEDIM>>& src) {
  dest.density = src.density;
  dest.momentum[0] = src.momentum[0];
  dest.species = src.species;
  dest.energy = src.energy;
  CompleteFromCons(equation, dest, AsCons(dest));
}

void EmbedState(Complete<IdealGasMix<AMREX_SPACEDIM>>& dest,
                IdealGasMix<AMREX_SPACEDIM>& equation,
                const Conservative<IdealGasMix<1>>& src) {
  dest.density = src.density;
  dest.momentum.setZero();
  dest.momentum[0] = src.momentum[0];
  dest.species = src.species;
  dest.energy = src.energy;
  CompleteFromCons(equation, dest, AsCons(dest));
}

::amrex::FArrayBox AllgatherMirrorData(const PatchHierarchy& hierarchy,
                                       int level, const ::amrex::Box& mirror) {
  const ::amrex::MultiFab& datas = hierarchy.GetPatchLevel(level).data;
  const int n_comps = hierarchy.GetDataDescription().n_cons_components;
  ::amrex::FArrayBox local_fab(mirror, n_comps);
  local_fab.setVal(0.0);
  for (::amrex::MFIter mfi(datas); mfi.isValid(); ++mfi) {
    const ::amrex::FArrayBox& data = datas[mfi];
    const ::amrex::Box subbox = mfi.tilebox() & mirror;
    local_fab.plus(data, subbox, 0, 0, n_comps);
  }
  ::amrex::FArrayBox global_fab(local_fab.box(), local_fab.nComp());
  ::MPI_Allreduce(local_fab.dataPtr(), global_fab.dataPtr(),
                  int(local_fab.size()), MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  return global_fab;
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

MultiBlockBoundary::MultiBlockBoundary(const MultiBlockBoundary& other)
    : log_(other.log_), valve_{other.valve_},
      plenum_equation_(other.plenum_equation_),
      tube_equation_(other.tube_equation_),
      plenum_mirror_box_(other.plenum_mirror_box_),
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

MultiBlockBoundary&
MultiBlockBoundary::operator=(const MultiBlockBoundary& other) {
  MultiBlockBoundary tmp(other);
  *this = std::move(tmp);
  return *this;
}

MultiBlockBoundary::MultiBlockBoundary(
    const MultiBlockGriddingAlgorithm& gridding,
    const BlockConnection& connection, int gcw,
    const FlameMasterReactor& reactor, int level)
    : MultiBlockBoundary("MultiBlockBoundary", gridding, connection, gcw,
                         reactor, level) {}

MultiBlockBoundary::MultiBlockBoundary(
    const std::string& channel, const MultiBlockGriddingAlgorithm& gridding,
    const BlockConnection& connection, int gcw,
    const FlameMasterReactor& reactor, int level)
    : MultiBlockBoundary(channel, gridding, connection, gcw, reactor, level,
                         nullptr) {}

MultiBlockBoundary::MultiBlockBoundary(
    const std::string& channel, const MultiBlockGriddingAlgorithm& gridding,
    const BlockConnection& connection, int gcw,
    const FlameMasterReactor& reactor, int level,
    std::shared_ptr<PressureValve> valve)
    : log_(boost::log::keywords::channel = channel,
           boost::log::keywords::severity = boost::log::trivial::debug),
      valve_{std::move(valve)}, plenum_equation_(reactor),
      tube_equation_(std::move(reactor)), dir_{connection.direction},
      side_{connection.side}, level_{level}, gcw_{gcw} {
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

template <int R> using Conservative = ::fub::Conservative<IdealGasMix<R>>;
template <int R> using Complete = ::fub::Complete<IdealGasMix<R>>;

void MultiBlockBoundary::ComputeBoundaryData(
    const cutcell::PatchHierarchy& plenum, const PatchHierarchy& tube) {
  const int d = static_cast<int>(dir_);
  const double plenum_dx = plenum.GetGeometry(level_).CellSize(d);
  const double tube_dx = tube.GetGeometry(level_).CellSize(d);

  BOOST_LOG_SCOPED_LOGGER_TAG(log_, "Time", plenum.GetTimePoint().count());

  //////////////////////////////////////////////////////////////////////////////
  // Integrate plenum states over mirror volume and fill ghost cells of tube
  // {{{
  // Integrate over mirror volume
  const ::amrex::FArrayBox plenum_data = AverageConservativeHierarchyStates(
      plenum, level_, plenum_mirror_box_, dir_);

  // Interpolate between grid cell sizes.
  InterpolateStates(*plenum_mirror_data_, tube_dx, plenum_data, plenum_dx,
                    side_);

  // Transform high dimensional states into low dimensional ones.
  // Store low dimensional states as the reference in the ghost cell region.
  {
    BasicView cons_states = MakeView<const Conservative<Plenum_Rank>>(
        *plenum_mirror_data_, plenum_equation_);
    BasicView complete_states =
        MakeView<Complete<Tube_Rank>>(*tube_ghost_data_, tube_equation_);
    Conservative<Plenum_Rank> cons(plenum_equation_);
    Complete<Tube_Rank> complete(tube_equation_);
    const std::ptrdiff_t i0 = plenum_mirror_data_->box().smallEnd(0);
    const std::ptrdiff_t j0 = tube_ghost_data_->box().smallEnd(0);
    ForEachIndex(Box<0>(complete_states), [&](std::ptrdiff_t j) {
      const std::ptrdiff_t k = j - j0;
      const std::ptrdiff_t i = i0 + k;
      Load(cons, cons_states, {i});
      if (cons.density > 0.0) {
        ReduceStateDimension(complete, tube_equation_, cons);
        Store(complete_states, complete, {j});
      }
    });

    const int ncomp = tube_ghost_data_->nComp();
    auto debug = boost::log::trivial::severity_level::debug;
    BOOST_LOG_SEV(log_, debug) << "Tube Ghost Data:";
    ForEachIndex(tube_ghost_data_->box(), [&](auto... is) {
      const ::amrex::IntVect iv{int(is)...};
      const double rho = (*tube_ghost_data_)(iv, 0);
      const double u = (*tube_ghost_data_)(iv, 1) / rho;
      const double p = (*tube_ghost_data_)(iv, ncomp - 5);
      BOOST_LOG_SEV(log_, debug)
          << fmt::format("xi = {}: {} kg/m3, {} m/s, {} Pa", iv[0], rho, u, p);
    });
  }
  // }}}

  // Interpolate between grid cell sizes.
  {
    ::amrex::FArrayBox tube_data =
        AllgatherMirrorData(tube, level_, tube_mirror_box_);
    InterpolateStates(*tube_mirror_data_, plenum_dx, tube_data, tube_dx,
                      Flip(side_));

    // Transform low dimensional states into high dimensional ones.
    BasicView cons_states =
        MakeView<Conservative<Tube_Rank>>(*tube_mirror_data_, tube_equation_);
    BasicView complete_states =
        MakeView<Complete<Plenum_Rank>>(*plenum_ghost_data_, plenum_equation_);
    Conservative<Tube_Rank> cons(tube_equation_);
    Complete<Plenum_Rank> complete(plenum_equation_);
    const std::ptrdiff_t i0 = Box<0>(cons_states).lower[0];
    const std::ptrdiff_t j0 = Box<0>(complete_states).lower[0];
    ForEachIndex(Box<0>(cons_states), [&](std::ptrdiff_t i) {
      const std::ptrdiff_t k = i - i0;
      const std::ptrdiff_t j = j0 + k;
      Load(cons, cons_states, {i});
      if (cons.density > 0.0) {
        EmbedState(complete, plenum_equation_, cons);
        Store(complete_states, complete, {j});
      }
    });

    const int ncomp = plenum_ghost_data_->nComp();
    BOOST_LOG(log_) << "Plenum Ghost Data:";
    ForEachIndex(plenum_ghost_data_->box(), [&](auto... is) {
      const ::amrex::IntVect iv{int(is)...};
      const double rho = (*plenum_ghost_data_)(iv, 0);
      const double u = (*plenum_ghost_data_)(iv, 1) / rho;
      const double p = (*plenum_ghost_data_)(iv, ncomp - 5);
      BOOST_LOG(log_) << fmt::format("{}: {} kg/m3, {} m/s, {} Pa", iv[0], rho,
                                     u, p);
    });
  }
}

const std::shared_ptr<PressureValve>& MultiBlockBoundary::GetValve() const
    noexcept {
  return valve_;
}

void MultiBlockBoundary::FillBoundary(::amrex::MultiFab& mf,
                                      const GriddingAlgorithm&, int) {
  if (valve_ && valve_->state == PressureValveState::closed) {
    //    ReflectiveBoundary closed{execution::openmp, tube_equation_, dir_,
    //      Flip(side_)};
    //    closed.FillBoundary(mf, geom, t, grid);
    return;
  }
  ::amrex::Box box = tube_ghost_data_->box();
  ForEachFab(execution::seq, mf, [&](const ::amrex::MFIter& mfi) {
    const ::amrex::Box b = mfi.growntilebox() & box;
    if (!b.isEmpty()) {
      for (int n = 0; n < mf.nComp(); ++n) {
        ForEachIndex(AsIndexBox<AMREX_SPACEDIM>(b), [&](auto... is) {
          const ::amrex::IntVect iv{int(is)...};
          if ((*tube_ghost_data_)(iv, 0) > 0.0) {
            mf[mfi](iv, n) = (*tube_ghost_data_)(iv, n);
          }
        });
      }

      if (valve_ && valve_->state == PressureValveState::open_fuel) {
        IdealGasMix<Tube_Rank>::Complete state(tube_equation_);
        auto view = MakeView<IdealGasMix<Tube_Rank>::Complete>(
            mf[mfi], tube_equation_, mfi.tilebox());
        std::vector<double> moles(static_cast<std::size_t>(
            tube_equation_.GetReactor().GetNSpecies()));
        ForEachIndex(Box<0>(view), [&](std::ptrdiff_t i) {
          if (2 < i && i < 5) {
            Load(state, view, {i});
            if (state.density > 0.0) {
              Array<double, 1, 1> velocity = state.momentum / state.density;
              tube_equation_.SetReactorStateFromComplete(state);
              span<const double> rmoles =
                  tube_equation_.GetReactor().GetMoleFractions();
              std::copy(rmoles.begin(), rmoles.end(), moles.begin());
              moles[Burke2012::sH2] = 2 * moles[Burke2012::sO2];
              tube_equation_.GetReactor().SetMoleFractions(moles);
              tube_equation_.CompleteFromReactor(state, velocity);
              Store(view, state, {i});
            }
          }
        });
      }
    }
  });
}

void MultiBlockBoundary::FillBoundary(::amrex::MultiFab& mf,
                                      const cutcell::GriddingAlgorithm& grid,
                                      int level) {
  ::amrex::Box ghost_box = MakeGhostBox(plenum_mirror_box_, gcw_, dir_, side_);
  if (valve_ && valve_->state == PressureValveState::closed) {
    const Eigen::Matrix<double, Plenum_Rank, 1> unit =
        UnitVector<Plenum_Rank>(dir_);
    const ::amrex::MultiFab& alphas =
        grid.GetPatchHierarchy().GetEmbeddedBoundary(level)->getVolFrac();
    OmpLocal<IdealGasMix<Plenum_Rank>::Complete> state(plenum_equation_);
    OmpLocal<IdealGasMix<Plenum_Rank>::Complete> reflected(plenum_equation_);
    ForEachFab(execution::openmp, mf, [&](const ::amrex::MFIter& mfi) {
      ::amrex::FArrayBox& fab = mf[mfi];
      const ::amrex::FArrayBox& alpha = alphas[mfi];
      ::amrex::Box box_to_fill = mfi.growntilebox() & ghost_box;
      if (!box_to_fill.isEmpty()) {
        auto states = MakeView<IdealGasMix<Plenum_Rank>::Complete>(
            fab, plenum_equation_, mfi.growntilebox());
        auto box = AsIndexBox<Plenum_Rank>(box_to_fill);
        ForEachIndex(box, [&](auto... is) {
          std::array<std::ptrdiff_t, Plenum_Rank> dest{is...};
          std::array<std::ptrdiff_t, Plenum_Rank> src =
              ReflectIndex(dest, box, dir_, side_);
          ::amrex::IntVect iv{
              AMREX_D_DECL(int(src[0]), int(src[1]), int(src[2]))};
          if (alpha(iv) > 0.0) {
            Load(*state, states, src);
            if (state->density > 0.0) {
              FUB_ASSERT(state->density > 0.0);
              Reflect(*reflected, *state, unit, plenum_equation_);
              Store(states, *reflected, dest);
            }
          }
        });
      }
    });
  } else {
    ForEachFab(execution::openmp, mf, [&](const ::amrex::MFIter& mfi) {
      const ::amrex::Box box = mfi.growntilebox() & ghost_box;
      if (!box.isEmpty()) {
        ::amrex::FArrayBox& fab = mf[mfi];
        for (int c = 0; c < mf.nComp(); ++c) {
          ForEachIndex(AsIndexBox<Plenum_Rank>(box), [&](auto... is) {
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
