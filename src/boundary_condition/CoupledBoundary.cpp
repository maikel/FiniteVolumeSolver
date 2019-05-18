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

#include "fub/boundary_condition/CoupledBoundary.hpp"
#include "fub/AMReX/ViewFArrayBox.hpp"
#include "fub/ForEach.hpp"

#include <AMReX_EB2.H>

#include <mpi.h>

namespace fub::amrex {
namespace {
::amrex::Box ReduceDimension(const ::amrex::Box& mirror) noexcept {
  ::amrex::IntVect lower{AMREX_D_DECL(mirror.smallEnd(0), 0, 0)};
  ::amrex::IntVect upper{AMREX_D_DECL(mirror.bigEnd(0), 0, 0)};
  return ::amrex::Box{lower, upper};
}

double CellVolume(const ::amrex::Geometry& geom) {
  return AMREX_D_TERM(geom.CellSize(0), *geom.CellSize(1), *geom.CellSize(2));
}

::amrex::FArrayBox ReduceTotalVolume(const cutcell::PatchHierarchy& hierarchy,
                                     int level, const ::amrex::Box& mirror) {
  const ::amrex::Box ld_mirror = ReduceDimension(mirror);
  const ::amrex::EBFArrayBoxFactory& factory =
      *hierarchy.GetEmbeddedBoundary(level);
  const ::amrex::MultiFab& alphas = factory.getVolFrac();
  ::amrex::FArrayBox local_reduced_volumes(ld_mirror);
  std::fill_n(local_reduced_volumes.dataPtr(),
              local_reduced_volumes.length()[0], 0.0);
  const double dx = CellVolume(hierarchy.GetGeometry(level));
  hierarchy.ForEachPatch(level, [&](PatchHandle patch) {
    const ::amrex::FArrayBox& alpha = alphas[*patch.iterator];
    IndexBox<AMREX_SPACEDIM> section =
        AsIndexBox<AMREX_SPACEDIM>(mirror & patch.iterator->tilebox());
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
  ForEachIndex(AsIndexBox<AMREX_SPACEDIM>(box), [&](auto... is) {
    Index<AMREX_SPACEDIM> index{is...};
    Index<AMREX_SPACEDIM> lowdim_index{index[0]};
    for (int i = 0; i < n_comps; ++i) {
      dest(AMREX_D_DECL(index[0], 0, 0), i) +=
          (volumes(index) * cell_volume / total_volumes(lowdim_index)) *
          src(is..., i);
    }
  });
}

::amrex::FArrayBox
AverageConservativeHierarchyStates(const cutcell::PatchHierarchy& hierarchy,
                                   int level, const ::amrex::Box& mirror) {
  const ::amrex::EBFArrayBoxFactory& eb = *hierarchy.GetEmbeddedBoundary(level);
  const ::amrex::MultiFab& datas = hierarchy.GetPatchLevel(level).data;
  const ::amrex::MultiFab& volumes = eb.getVolFrac();
  const int n_comps = hierarchy.GetDataDescription().n_cons_components;
  const ::amrex::FArrayBox total_volume =
      ReduceTotalVolume(hierarchy, level, mirror);
  ::amrex::FArrayBox fab(total_volume.box(), n_comps);
  fab.setVal(0.0);
  //  const std::ptrdiff_t linear =
  //      static_cast<std::ptrdiff_t>(fab.box().length()[0]);
  //  const std::ptrdiff_t size = linear * fab.nComp();
  //  std::fill_n(fab.dataPtr(), size, 0.0);
  const double cell_volume = CellVolume(hierarchy.GetGeometry(0));
  PatchDataView<double, AMREX_SPACEDIM + 1> fabv = MakePatchDataView(fab);
  PatchDataView<const double, AMREX_SPACEDIM> tot_vol =
      MakePatchDataView(total_volume, 0);
  hierarchy.ForEachPatch(level, [&](PatchHandle patch) {
    const ::amrex::FArrayBox& alpha = volumes[*patch.iterator];
    const ::amrex::FArrayBox& data = datas[*patch.iterator];
    PatchDataView<const double, AMREX_SPACEDIM + 1> src =
        MakePatchDataView(data);
    PatchDataView<const double, AMREX_SPACEDIM> vol =
        MakePatchDataView(alpha, 0);
    ::amrex::Box box = patch.iterator->tilebox() & mirror;
    IntegrateConserativeStates(fabv, tot_vol, src, vol, box, cell_volume);
  });
  ::amrex::FArrayBox global_fab(fab.box(), fab.nComp());
  global_fab.setVal(0.0);
  ::MPI_Allreduce(fab.dataPtr(), global_fab.dataPtr(), int(fab.size()), MPI_DOUBLE,
                  MPI_SUM, MPI_COMM_WORLD);
  return global_fab;
}

void InterpolateStates_Greater(span<double> dest, double dest_dx,
                               span<const double> src, double src_dx) {
  FUB_ASSERT(dest_dx > src_dx);
  const std::ptrdiff_t n = static_cast<std::ptrdiff_t>(dest_dx / src_dx);
  const double weight = src_dx / dest_dx;
  std::ptrdiff_t i0 = 0;
  double last_remainder = 0.0;
  for (std::ptrdiff_t j = 0; j < dest.size(); ++j) {
    dest[j] = 0.0;
    std::ptrdiff_t offset = 0.0;
    const double remainder = ::fmod(dest_dx - last_remainder, src_dx);
    if (last_remainder > 0.0) {
      dest[j] += (weight - last_remainder / dest_dx) * src[i0];
      offset = 1LL;
    }
    for (std::ptrdiff_t i = 0; i < n; ++i) {
      dest[j] += weight * src[i0 + offset + i];
    }
    if (i0 + offset + n < src.size()) {
      dest[j] += (remainder / dest_dx) * src[i0 + offset + n];
    }
    i0 += offset + n;
    last_remainder = remainder;
  }
}

void InterpolateStates_Less(span<double> dest, double dest_dx,
                            span<const double> src, double src_dx) {
  FUB_ASSERT(dest_dx < src_dx);
  for (std::ptrdiff_t i = 0; i < dest.size(); ++i) {
    const double x_lower = dest_dx * i;
    const double x_upper = dest_dx * (i + 1);
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
                       const ::amrex::FArrayBox& src, double src_dx) {
  const int n_comps = dest.nComp();
  FUB_ASSERT(n_comps == src.nComp());
  if (dest_dx > src_dx) {
    // Coarsen Data
    for (int i = 0; i < n_comps; ++i) {
      span<double> destv = MakePatchDataView(dest, i).Span();
      span<const double> srcv = MakePatchDataView(src, i).Span();
      InterpolateStates_Greater(destv, dest_dx, srcv, src_dx);
    }
  } else if (dest_dx < src_dx) {
    // Refine Data
    for (int i = 0; i < n_comps; ++i) {
      span<double> destv = MakePatchDataView(dest, i).Span();
      span<const double> srcv = MakePatchDataView(src, i).Span();
      InterpolateStates_Less(destv, dest_dx, srcv, src_dx);
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
  hierarchy.ForEachPatch(level, [&](PatchHandle patch) {
    const ::amrex::FArrayBox& data = datas[*patch.iterator];
    const ::amrex::Box subbox = patch.iterator->tilebox() & mirror;
    local_fab.plus(data, subbox, 0, 0, n_comps);
  });
  ::amrex::FArrayBox global_fab(local_fab.box(), local_fab.nComp());
  ::MPI_Allreduce(local_fab.dataPtr(), global_fab.dataPtr(), int(local_fab.size()),
                  MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  return global_fab;
}
} // namespace

CoupledBoundary::CoupledBoundary(const cutcell::PatchHierarchy& plenum,
                                 const PatchHierarchy& tube, int gcw,
                                 const FlameMasterReactor& reactor)
    : plenum_equation_(reactor), tube_equation_(std::move(reactor)) {
  const double plenum_dx = plenum.GetGeometry(0).CellSize(0);
  const double tube_dx = tube.GetGeometry(0).CellSize(0);
  // Allocate mirror data containing the plenum data (reduced to 1d)
  {
    const int mirror_width =
        static_cast<int>(std::ceil((gcw * tube_dx) / plenum_dx));
    ::amrex::Box domain = plenum.GetGeometry(0).Domain();
    ::amrex::IntVect lower = domain.smallEnd();
    ::amrex::IntVect upper = domain.bigEnd();
    upper[0] = mirror_width;
    plenum_mirror_box_ = ::amrex::Box{lower, upper};
    const int ncons = plenum.GetDataDescription().n_cons_components;
    const int ncomp = plenum.GetDataDescription().n_state_components;
    plenum_mirror_data_.resize(
        ::amrex::Box{{AMREX_D_DECL(domain.smallEnd(0), 0, 0)},
                     {AMREX_D_DECL(domain.smallEnd(0) + gcw - 1, 0, 0)}},
        ncons);
    plenum_ghost_data_.resize(
        ::amrex::Box{{AMREX_D_DECL(domain.smallEnd(0) - gcw, 0, 0)},
                     {AMREX_D_DECL(domain.smallEnd(0) - 1, 0, 0)}},
        ncomp);
  }
  // Allocate mirror data for the tube data
  {
    const int mirror_width =
        static_cast<int>(std::ceil((gcw * plenum_dx) / tube_dx));
    ::amrex::Box domain = tube.GetGeometry(0).Domain();
    ::amrex::IntVect lower = domain.smallEnd();
    ::amrex::IntVect upper = domain.bigEnd();
    lower[0] = upper[0] - mirror_width;
    tube_mirror_box_ = ::amrex::Box{lower, upper};
    const int ncons = tube.GetDataDescription().n_cons_components;
    const int ncomp = tube.GetDataDescription().n_state_components;
    tube_mirror_data_.resize(
        ::amrex::Box{{AMREX_D_DECL(domain.bigEnd(0) - gcw + 1, 0, 0)},
                     {AMREX_D_DECL(domain.bigEnd(0), 0, 0)}},
        ncons);
    tube_ghost_data_.resize(
        ::amrex::Box{{AMREX_D_DECL(domain.bigEnd(0) + 1, 0, 0)},
                     {AMREX_D_DECL(domain.bigEnd(0) + gcw, 0, 0)}},
        ncomp);
  }
  // Fill internal mirror data
  ComputeBoundaryData(plenum, tube);
}

void CoupledBoundary::ComputeBoundaryData(const cutcell::PatchHierarchy& plenum,
                                          const PatchHierarchy& tube) {
  const double plenum_dx = plenum.GetGeometry(0).CellSize(0);
  const double tube_dx = tube.GetGeometry(0).CellSize(0);
  {
    // Integrate plenum states.
    const ::amrex::FArrayBox plenum_data =
        AverageConservativeHierarchyStates(plenum, 0, plenum_mirror_box_);

    // Interpolate between grid cell sizes.
    InterpolateStates(plenum_mirror_data_, tube_dx, plenum_data, plenum_dx);

    // Transform high dimensional states into low dimensional ones.
    auto cons_states =
        MakeView<BasicView<Conservative<IdealGasMix<AMREX_SPACEDIM>>>>(
            MakePatchDataView(plenum_mirror_data_), plenum_equation_);
    auto complete_states = MakeView<BasicView<Complete<IdealGasMix<1>>>>(
        MakePatchDataView(tube_ghost_data_), tube_equation_);
    Conservative<IdealGasMix<AMREX_SPACEDIM>> cons(plenum_equation_);
    Complete<IdealGasMix<1>> complete(tube_equation_);
    const std::ptrdiff_t i0 = Box<0>(complete_states).lower[0];
    ForEachIndex(Box<0>(complete_states), [&](std::ptrdiff_t i) {
      Load(cons, cons_states, {i - i0});
      ReduceStateDimension(complete, tube_equation_, cons);
      Store(complete_states, complete, {i});
    });
  }

  {
    // Interpolate between grid cell sizes.
    ::amrex::FArrayBox tube_data =
        AllgatherMirrorData(tube, 0, tube_mirror_box_);
    InterpolateStates(tube_mirror_data_, plenum_dx, tube_data, tube_dx);

    // Transform low dimensional states into high dimensional ones.
    auto cons_states = MakeView<BasicView<Conservative<IdealGasMix<1>>>>(
        MakePatchDataView(tube_mirror_data_), tube_equation_);
    auto complete_states =
        MakeView<BasicView<Complete<IdealGasMix<AMREX_SPACEDIM>>>>(
            MakePatchDataView(plenum_ghost_data_), plenum_equation_);
    Conservative<IdealGasMix<1>> cons(tube_equation_);
    Complete<IdealGasMix<AMREX_SPACEDIM>> complete(plenum_equation_);
    const std::ptrdiff_t i0 = Box<0>(cons_states).lower[0];
    const std::ptrdiff_t j0 = Box<0>(complete_states).lower[0];
    ForEachIndex(Box<0>(cons_states), [&](std::ptrdiff_t i) {
      Load(cons, cons_states, {i});
      EmbedState(complete, plenum_equation_, cons);
      Store(complete_states, complete, {i + j0 - i0});
    });
  }
}

void CoupledBoundary::FillPlenumBoundary(
    const PatchDataView<double, AMREX_SPACEDIM + 1>& data,
    const cutcell::PatchHierarchy&, PatchHandle, Location loc, int, Duration) {
  if (loc.direction != 0 || loc.side != 0) {
    return;
  }

  auto dest = MakeView<BasicView<Complete<IdealGasMix<AMREX_SPACEDIM>>>>(
      data, plenum_equation_);

  auto src = MakeView<BasicView<Complete<IdealGasMix<AMREX_SPACEDIM>>>>(
      MakePatchDataView(plenum_ghost_data_), plenum_equation_);

  IndexBox<AMREX_SPACEDIM> fill_box = Box<0>(dest);
  fill_box.upper[0] = Box<0>(src).upper[0];
  IndexBox<1> ghost_data_box = AsIndexBox<1>(plenum_ghost_data_.box());

  Complete<IdealGasMix<AMREX_SPACEDIM>> complete(plenum_equation_);
  ForEachIndex(fill_box, [&](auto... is) {
    Index<AMREX_SPACEDIM> index{is...};
    if (Contains(ghost_data_box, Index<1>{index[0]})) {
      Load(complete, src, {index[0]});
      Store(dest, complete, index);
    }
  });
}

void CoupledBoundary::FillTubeBoundary(
    const PatchDataView<double, AMREX_SPACEDIM + 1>& data,
    const PatchHierarchy&, PatchHandle, Location loc, int fill_width,
    Duration) {
  if (loc.direction != 0 || loc.side != 1) {
    return;
  }

  auto dest =
      MakeView<BasicView<Complete<IdealGasMix<1>>>>(data, tube_equation_);

  auto src = MakeView<BasicView<Complete<IdealGasMix<1>>>>(
      MakePatchDataView(tube_ghost_data_), tube_equation_);

  IndexBox<1> fill_box = Box<0>(dest);
  fill_box.lower[0] = fill_box.upper[0] - fill_width;
  IndexBox<1> ghost_data_box = AsIndexBox<1>(tube_ghost_data_.box());

  Complete<IdealGasMix<1>> complete(tube_equation_);
  ForEachIndex(fill_box, [&](std::ptrdiff_t i) {
    if (Contains(ghost_data_box, Index<1>{i})) {
      Load(complete, src, {i});
      Store(dest, complete, {i});
    }
  });
}

} // namespace fub::amrex
