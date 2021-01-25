// Copyright (c) 2020 Maikel Nadolski
// Copyright (c) 2019 Stefan Vater
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

#include "fub/AMReX/solver/BK19Solver.hpp"

#include "fub/AMReX/ForEachFab.hpp"
#include "fub/AMReX/ForEachIndex.hpp"
#include "fub/AMReX/output/DebugOutput.hpp"
#include "fub/ForEach.hpp"
#include "fub/ext/Vc.hpp"

namespace fub::amrex {

BK19SolverOptions::BK19SolverOptions(const ProgramOptions& map) {
  do_initial_projection =
      GetOptionOr(map, "do_initial_projection", do_initial_projection);
  mlmg_tolerance_rel =
      GetOptionOr(map, "mlmg_tolerance_rel", mlmg_tolerance_rel);
  mlmg_tolerance_abs =
      GetOptionOr(map, "mlmg_tolerance_abs", mlmg_tolerance_abs);
  mlmg_max_iter = GetOptionOr(map, "mlmg_max_iter", mlmg_max_iter);
  mlmg_verbose = GetOptionOr(map, "mlmg_verbose", mlmg_verbose);
  bottom_tolerance_rel =
      GetOptionOr(map, "bottom_tolerance_rel", bottom_tolerance_rel);
  bottom_tolerance_abs =
      GetOptionOr(map, "bottom_tolerance_abs", bottom_tolerance_abs);
  bottom_max_iter = GetOptionOr(map, "bottom_max_iter", bottom_max_iter);
  bottom_verbose = GetOptionOr(map, "bottom_verbose", bottom_verbose);
  always_use_bnorm = GetOptionOr(map, "always_use_bnorm", always_use_bnorm);
  prefix = GetOptionOr(map, "prefix", prefix);
  output_between_steps =
      GetOptionOr(map, "output_between_steps", output_between_steps);
}

void BK19SolverOptions::Print(SeverityLogger& log) {
  BOOST_LOG(log) << fmt::format(" - do_initial_projection = {}",
                                do_initial_projection);
  BOOST_LOG(log) << fmt::format(" - mlmg_tolerance_rel = {}",
                                mlmg_tolerance_rel);
  BOOST_LOG(log) << fmt::format(" - mlmg_tolerance_abs = {}",
                                mlmg_tolerance_abs);
  BOOST_LOG(log) << fmt::format(" - mlmg_max_iter = {}", mlmg_max_iter);
  BOOST_LOG(log) << fmt::format(" - mlmg_verbose = {}", mlmg_verbose);
  BOOST_LOG(log) << fmt::format(" - bottom_tolerance_rel = {}",
                                bottom_tolerance_rel);
  BOOST_LOG(log) << fmt::format(" - bottom_tolerance_abs = {}",
                                bottom_tolerance_abs);
  BOOST_LOG(log) << fmt::format(" - bottom_max_iter = {}", bottom_max_iter);
  BOOST_LOG(log) << fmt::format(" - bottom_verbose = {}", bottom_verbose);
  BOOST_LOG(log) << fmt::format(" - always_use_bnorm = {}", always_use_bnorm);
  BOOST_LOG(log) << fmt::format(" - prefix = {}", prefix);
  BOOST_LOG(log) << fmt::format(" - output_between_steps = {}",
                                output_between_steps);
}

BK19PhysicalParameters::BK19PhysicalParameters(const ProgramOptions& map) {
  R_gas = GetOptionOr(map, "R_gas", R_gas);
  gamma = GetOptionOr(map, "gamma", gamma);
  c_p = GetOptionOr(map, "c_p", c_p);
  gravity = GetOptionOr(map, "gravity", gravity);
  f = GetOptionOr(map, "f", f);
  k_vect = GetOptionOr(map, "k_vect", k_vect);
  alpha_p = GetOptionOr(map, "alpha_p", alpha_p);
  Msq = GetOptionOr(map, "Msq", Msq);
}

void BK19PhysicalParameters::Print(SeverityLogger& log) {
  BOOST_LOG(log) << fmt::format(" - R_gas = {}", R_gas);
  BOOST_LOG(log) << fmt::format(" - gamma = {}", gamma);
  BOOST_LOG(log) << fmt::format(" - c_p = {}", c_p);
  BOOST_LOG(log) << fmt::format(" - gravity = {}", gravity);
  BOOST_LOG(log) << fmt::format(" - f = {}", f);
  BOOST_LOG(log) << fmt::format(" - k_vect = {}", k_vect);
  BOOST_LOG(log) << fmt::format(" - alpha_p = {}", alpha_p);
  BOOST_LOG(log) << fmt::format(" - Msq = {}", Msq);
}

namespace {

inline constexpr int no_ghosts = 0;
inline constexpr int one_ghost_cell_width = 1;

inline constexpr int one_component = 1;

template <typename... Pointers>
void AdvanceBy(std::ptrdiff_t n, Pointers&... ps) {
  ((ps += n), ...);
}

DebugSnapshot::ComponentNames GetCompleteVariableNames() {
  using Equation = CompressibleAdvection<2>;
  fub::CompressibleAdvection<2> equation{};
  using Traits = StateTraits<Complete<Equation>>;
  constexpr auto names = Traits::names;
  const auto depths = Depths(equation, Type<Complete<Equation>>{});
  const std::size_t n_names =
      std::tuple_size<remove_cvref_t<decltype(names)>>::value;
  DebugSnapshot::ComponentNames varnames;
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
  return varnames;
}

std::vector<::amrex::MultiFab>
CopyScratchFromContext(const IntegratorContext& context) {
  const int nlevel = context.GetPatchHierarchy().GetNumberOfLevels();
  std::vector<::amrex::MultiFab> scratch_vector;
  scratch_vector.reserve(nlevel);
  for (int level = 0; level < nlevel; ++level) {
    const ::amrex::MultiFab& scratch = context.GetScratch(level);
    ::amrex::MultiFab& scratch_aux = scratch_vector.emplace_back(
        scratch.boxArray(), scratch.DistributionMap(), scratch.nComp(),
        no_ghosts);
    scratch_aux.copy(scratch);
  }
  return scratch_vector;
}

void CopyScratchToContext(IntegratorContext& context,
                          span<const ::amrex::MultiFab> scratch) {
  for (int level = 0; level < scratch.size(); ++level) {
    ::amrex::MultiFab& scratch_out = context.GetScratch(level);
    scratch_out.copy(scratch[level]);
  }
}

void CopyPiToHierarchy(CompressibleAdvectionIntegratorContext& context,
                       span<const ::amrex::MultiFab> pi) {
  for (std::ptrdiff_t level = 0; level < pi.size(); ++level) {
    ::amrex::MultiFab& pi_out = context.GetPi(level);
    ::amrex::MultiFab::Copy(pi_out, pi[level], 0, 0, 1, ::amrex::IntVect(0));
  }
}

// Compute inplace dest = factor1 * (dest + factor2 * k x src)
void AddCrossProduct(std::array<span<double>, 2> dest,
                     std::array<span<const double>, 2> src, double factor1,
                     double factor2, span<const double, 3> k) {
  Vc::Vector<double> fac1(factor1);
  Vc::Vector<double> fac2(factor2);
  Vc::Vector<double> k2(k[2]);
  double* out_x = dest[0].begin();
  double* out_y = dest[1].begin();
  const double* in_x = src[0].begin();
  const double* in_y = src[1].begin();
  const double* end = dest[0].end();
  std::ptrdiff_t n = end - out_x;
  const auto size = static_cast<std::ptrdiff_t>(Vc::Vector<double>::size());
  while (n >= size) {
    Vc::Vector<double> du(out_x, Vc::Unaligned);
    Vc::Vector<double> dv(out_y, Vc::Unaligned);
    Vc::Vector<double> su(in_x, Vc::Unaligned);
    Vc::Vector<double> sv(in_y, Vc::Unaligned);
    Vc::Vector u_next = fac2 * (du - fac1 * k2 * sv);
    Vc::Vector v_next = fac2 * (dv + fac1 * k2 * su);
    u_next.store(out_x, Vc::Unaligned);
    v_next.store(out_y, Vc::Unaligned);
    AdvanceBy(size, out_x, out_y, in_x, in_y);
    n = end - out_x;
  }
  const auto mask = Vc::Vector<double>([](int i) { return i; }) < int(n);
  const Vc::Vector<double> du = mask_load(out_x, mask);
  const Vc::Vector<double> dv = mask_load(out_y, mask);
  const Vc::Vector<double> su = mask_load(in_x, mask);
  const Vc::Vector<double> sv = mask_load(in_y, mask);
  Vc::Vector<double> u_next = fac2 * (du - fac1 * k2 * sv);
  Vc::Vector<double> v_next = fac2 * (dv + fac1 * k2 * su);
  u_next.store(out_x, mask, Vc::Unaligned);
  v_next.store(out_y, mask, Vc::Unaligned);
}

// Compute inplace momentum = factor1 * (momentum + factor2 * k x momentum)
void AddCrossProduct(std::array<span<double>, 3> dest,
                     std::array<span<const double>, 3> src, double factor1,
                     double factor2, span<const double, 3> k) {
  Vc::Vector<double> fac1(factor1);
  Vc::Vector<double> fac2(factor2);
  Vc::Vector<double> k0(k[0]);
  Vc::Vector<double> k1(k[1]);
  Vc::Vector<double> k2(k[2]);
  double* out_x = dest[0].begin();
  double* out_y = dest[1].begin();
  double* out_z = dest[2].begin();
  const double* in_x = src[0].begin();
  const double* in_y = src[1].begin();
  const double* in_z = src[2].begin();
  const double* end = dest[0].end();
  std::ptrdiff_t n = end - out_x;
  const auto size = static_cast<std::ptrdiff_t>(Vc::Vector<double>::size());
  while (n >= size) {
    Vc::Vector<double> du(out_x, Vc::Unaligned);
    Vc::Vector<double> dv(out_y, Vc::Unaligned);
    Vc::Vector<double> dw(out_y, Vc::Unaligned);
    Vc::Vector<double> su(out_x, Vc::Unaligned);
    Vc::Vector<double> sv(out_y, Vc::Unaligned);
    Vc::Vector<double> sw(out_y, Vc::Unaligned);
    Vc::Vector u_next = fac2 * (du + fac1 * (k1 * sw - k2 * sv));
    Vc::Vector v_next = fac2 * (dv + fac1 * (k2 * su - k0 * sw));
    Vc::Vector w_next = fac2 * (dw + fac1 * (k0 * sv - k1 * su));
    u_next.store(out_x, Vc::Unaligned);
    v_next.store(out_y, Vc::Unaligned);
    w_next.store(out_z, Vc::Unaligned);
    AdvanceBy(size, out_x, out_y, out_z, in_x, in_y, in_z);
    n = end - out_x;
  }
  const auto mask = Vc::Vector<double>([](int i) { return i; }) < int(n);
  const Vc::Vector<double> du = mask_load(out_x, mask);
  const Vc::Vector<double> dv = mask_load(out_y, mask);
  const Vc::Vector<double> dw = mask_load(out_z, mask);
  const Vc::Vector<double> su = mask_load(in_x, mask);
  const Vc::Vector<double> sv = mask_load(in_y, mask);
  const Vc::Vector<double> sw = mask_load(in_z, mask);
  Vc::Vector u_next = fac2 * (du + fac1 * (k1 * sw - k2 * sv));
  Vc::Vector v_next = fac2 * (dv + fac1 * (k2 * su - k0 * sw));
  Vc::Vector w_next = fac2 * (dw + fac1 * (k0 * sv - k1 * su));
  u_next.store(out_x, Vc::Unaligned);
  v_next.store(out_y, Vc::Unaligned);
  w_next.store(out_z, Vc::Unaligned);
}

// Compute inplace dest = factor1 * (dest + factor2 * k x src)
template <typename I, std::size_t VelocityRank>
void AddCrossProduct(::amrex::MultiFab& dest,
                     const std::array<I, VelocityRank>& dest_index,
                     const ::amrex::MultiFab& src,
                     const std::array<I, VelocityRank>& src_index,
                     double factor1, double factor2, span<const double, 3> k) {
  ForEachFab(execution::openmp, dest, [&](const ::amrex::MFIter& mfi) {
    const ::amrex::Box box = mfi.growntilebox();
    if constexpr (VelocityRank == 2) {
      auto du = MakePatchDataView(dest[mfi], dest_index[0], box);
      auto dv = MakePatchDataView(dest[mfi], dest_index[1], box);
      auto su = MakePatchDataView(src[mfi], src_index[0], box);
      auto sv = MakePatchDataView(src[mfi], src_index[1], box);
      ForEachRow(std::tuple{du, dv, su, sv},
                 [factor1, factor2, k](span<double> du, span<double> dv,
                                       span<const double> su,
                                       span<const double> sv) {
                   AddCrossProduct(std::array<span<double>, 2>{du, dv},
                                   std::array<span<const double>, 2>{su, sv},
                                   factor1, factor2, k);
                 });
    } else if constexpr (VelocityRank == 3) {
      auto du = MakePatchDataView(dest[mfi], dest_index[0], box);
      auto dv = MakePatchDataView(dest[mfi], dest_index[1], box);
      auto dw = MakePatchDataView(dest[mfi], dest_index[2], box);
      auto su = MakePatchDataView(src[mfi], src_index[0], box);
      auto sv = MakePatchDataView(src[mfi], src_index[1], box);
      auto sw = MakePatchDataView(src[mfi], src_index[2], box);
      ForEachRow(
          std::tuple{du, dv, dw, su, sv, sw},
          [factor1, factor2, k](span<double> du, span<double> dv,
                                span<double> dw, span<const double> su,
                                span<const double> sv, span<const double> sw) {
            AddCrossProduct(std::array<span<double>, 3>{du, dv, dw},
                            std::array<span<const double>, 3>{su, sv, sw},
                            factor1, factor2, k);
          });
    }
  });
}
template <int Rank, int VelocityRank, std::size_t N>
void AddGravity(::amrex::MultiFab& dst, const std::array<int, N>& dest_index,
                const ::amrex::MultiFab& src,
                IndexMapping<CompressibleAdvection<Rank, VelocityRank>> index,
                const BK19PhysicalParameters& physical_parameters,
                Duration dt) {
  for (int i = 0; i < VelocityRank; ++i) {
    ::amrex::MultiFab::Saxpy(dst, dt.count() * physical_parameters.gravity[i],
                             src, index.density, dest_index[i], 1, dst.nGrow());
  }
}

template <int Rank, int VelocityRank>
void ApplyExplicitSourceTerms(
    BK19Solver<Rank, VelocityRank>& solver, Duration dt,
    const BK19PhysicalParameters& physical_parameters) {
  CompressibleAdvectionIntegratorContext& context =
      solver.GetAdvectionSolver().GetContext();
  const int nlevel = context.GetPatchHierarchy().GetNumberOfLevels();
  IndexMapping index = solver.GetEquation().GetIndexMapping();
  const double factor1 = -dt.count() * physical_parameters.f;
  const double factor2 = 1.0 / (1.0 + factor1 * factor1);
  for (int level = 0; level < nlevel; ++level) {
    ::amrex::MultiFab& scratch = context.GetScratch(level);
    AddCrossProduct(scratch, index.momentum, scratch, index.momentum, factor1,
                    factor2, physical_parameters.k_vect);
    AddGravity(scratch, index.momentum, scratch, index, physical_parameters,
               dt);
  }
}

std::vector<::amrex::MultiFab>
AllocateOnCells(const CompressibleAdvectionIntegratorContext& context,
                int n_comps, int n_grow) {
  const int nlevel = context.GetPatchHierarchy().GetNumberOfLevels();
  std::vector<::amrex::MultiFab> nodes;
  nodes.reserve(nlevel);
  for (int ilvl = 0; ilvl < nlevel; ++ilvl) {
    const ::amrex::MultiFab& scratch = context.GetScratch(ilvl);
    nodes.emplace_back(scratch.boxArray(), scratch.DistributionMap(), n_comps,
                       n_grow);
  }
  return nodes;
}

::amrex::MultiFab
AllocateOnCells(const CompressibleAdvectionIntegratorContext& context,
                int level, int n_comps, int n_grow) {
  const ::amrex::MultiFab& scratch = context.GetScratch(level);
  return ::amrex::MultiFab(scratch.boxArray(), scratch.DistributionMap(),
                           n_comps, n_grow);
}

std::vector<::amrex::MultiFab>
AllocateOnNodes(const CompressibleAdvectionIntegratorContext& context,
                int n_comps, int n_grow) {
  const int nlevel = context.GetPatchHierarchy().GetNumberOfLevels();
  std::vector<::amrex::MultiFab> nodes;
  nodes.reserve(nlevel);
  for (int ilvl = 0; ilvl < nlevel; ++ilvl) {
    const ::amrex::MultiFab& scratch = context.GetScratch(ilvl);
    const ::amrex::BoxArray grids =
        ::amrex::convert(scratch.boxArray(), ::amrex::IntVect::TheNodeVector());
    nodes.emplace_back(grids, scratch.DistributionMap(), n_comps, n_grow);
  }
  return nodes;
}

::amrex::MultiFab
AllocateOnNodes(const CompressibleAdvectionIntegratorContext& context,
                int level, int n_comps, int n_grow) {
  const ::amrex::MultiFab& scratch = context.GetScratch(level);
  const ::amrex::BoxArray grids =
      ::amrex::convert(scratch.boxArray(), ::amrex::IntVect::TheNodeVector());
  return ::amrex::MultiFab(grids, scratch.DistributionMap(), n_comps, n_grow);
}

std::vector<::amrex::MultiFab>
ZerosOnCells(const CompressibleAdvectionIntegratorContext& context, int n_comps,
             int n_grow) {
  std::vector<::amrex::MultiFab> cells =
      AllocateOnCells(context, n_comps, n_grow);
  for (::amrex::MultiFab& mf : cells) {
    mf.setVal(0.0);
  }
  return cells;
}

std::vector<::amrex::MultiFab>
ZerosOnNodes(const CompressibleAdvectionIntegratorContext& context, int n_comps,
             int n_grow) {
  std::vector<::amrex::MultiFab> nodes =
      AllocateOnNodes(context, n_comps, n_grow);
  for (::amrex::MultiFab& mf : nodes) {
    mf.setVal(0.0);
  }
  return nodes;
}

::amrex::Vector<::amrex::MultiFab*> ToPointers(span<::amrex::MultiFab> array) {
  ::amrex::Vector<::amrex::MultiFab*> pointers(array.size());
  for (int i = 0; i < array.size(); ++i) {
    pointers[i] = std::addressof(array[i]);
  }
  return pointers;
}

::amrex::Vector<const ::amrex::MultiFab*>
ToConstPointers(span<const ::amrex::MultiFab> array) {
  ::amrex::Vector<const ::amrex::MultiFab*> pointers(array.size());
  for (int i = 0; i < array.size(); ++i) {
    pointers[i] = std::addressof(array[i]);
  }
  return pointers;
}

::amrex::Vector<const ::amrex::Geometry*>
GetGeometries(const CompressibleAdvectionIntegratorContext& context) {
  const PatchHierarchy& hierarchy = context.GetPatchHierarchy();
  const int nlevel = hierarchy.GetNumberOfLevels();
  ::amrex::Vector<const ::amrex::Geometry*> pis(
      static_cast<std::size_t>(nlevel));
  for (int i = 0; i < nlevel; ++i) {
    pis[i] = std::addressof(hierarchy.GetGeometry(i));
  }
  return pis;
}

::amrex::Vector<const ::amrex::MultiFab*>
GetScratches(const CompressibleAdvectionIntegratorContext& context) {
  const PatchHierarchy& hierarchy = context.GetPatchHierarchy();
  const int nlevel = hierarchy.GetNumberOfLevels();
  ::amrex::Vector<const ::amrex::MultiFab*> scratches(
      static_cast<std::size_t>(nlevel));
  for (int i = 0; i < nlevel; ++i) {
    scratches[i] = std::addressof(context.GetScratch(i));
  }
  return scratches;
}

::amrex::Vector<const ::amrex::MultiFab*>
GetPvs(const CompressibleAdvectionIntegratorContext& context) {
  const PatchHierarchy& hierarchy = context.GetPatchHierarchy();
  const int nlevel = hierarchy.GetNumberOfLevels();
  ::amrex::Vector<const ::amrex::MultiFab*> pvs(
      static_cast<std::size_t>(nlevel));
  for (int i = 0; i < nlevel; ++i) {
    pvs[i] = std::addressof(context.GetAdvectiveFluxes(i).on_cells);
  }
  return pvs;
}

::amrex::Vector<const ::amrex::MultiFab*>
GetPis(const CompressibleAdvectionIntegratorContext& context) {
  const PatchHierarchy& hierarchy = context.GetPatchHierarchy();
  const int nlevel = hierarchy.GetNumberOfLevels();
  ::amrex::Vector<const ::amrex::MultiFab*> pis(
      static_cast<std::size_t>(nlevel));
  for (int i = 0; i < nlevel; ++i) {
    pis[i] = hierarchy.GetPatchLevel(i).nodes.get();
  }
  return pis;
}

::amrex::Vector<::amrex::MultiFab*>
GetPis(CompressibleAdvectionIntegratorContext& context) {
  PatchHierarchy& hierarchy = context.GetPatchHierarchy();
  const int nlevel = hierarchy.GetNumberOfLevels();
  ::amrex::Vector<::amrex::MultiFab*> pis(static_cast<std::size_t>(nlevel));
  for (int i = 0; i < nlevel; ++i) {
    pis[i] = hierarchy.GetPatchLevel(i).nodes.get();
  }
  return pis;
}

template <int Rank, int VelocityRank>
void ComputePvFromScratch(
    ::amrex::MultiFab& Pv, const ::amrex::MultiFab& scratch, int ngrow,
    const IndexMapping<CompressibleAdvection<Rank, VelocityRank>>& index) {
  for (std::size_t i = 0; i < index.momentum.size(); ++i) {
    const int dest_component = static_cast<int>(i);
    ::amrex::MultiFab::Copy(Pv, scratch, index.PTdensity, dest_component,
                            one_component, ngrow);
    ::amrex::MultiFab::Multiply(Pv, scratch, index.momentum[i], dest_component,
                                one_component, ngrow);
    ::amrex::MultiFab::Divide(Pv, scratch, index.density, dest_component,
                              one_component, ngrow);
  }
}

template <int Rank, int VelocityRank>
::amrex::MultiFab
ComputePvFromScratch(const BK19Solver<Rank, VelocityRank>& solver, int level,
                     int ngrow) {
  const CompressibleAdvectionIntegratorContext& context =
      solver.GetAdvectionSolver().GetContext();
  ::amrex::MultiFab Pv = AllocateOnCells(context, level, VelocityRank, ngrow);
  ComputePvFromScratch(Pv, context.GetScratch(level), ngrow,
                       solver.GetEquation().GetIndexMapping());
  return Pv;
}

template <int Rank, int VelocityRank>
std::vector<::amrex::MultiFab>
ComputePvFromScratch(const BK19Solver<Rank, VelocityRank>& solver, int ngrow) {
  const auto& context = solver.GetAdvectionSolver().GetContext();
  const int nlevel = context.GetPatchHierarchy().GetNumberOfLevels();
  std::vector<::amrex::MultiFab> Pvs;
  for (int level = 0; level < nlevel; ++level) {
    Pvs.push_back(ComputePvFromScratch(solver, level, ngrow));
  }
  return Pvs;
}

IndexBox<2> BoxFromStencil(const IndexBox<2>& box,
                           std::array<std::ptrdiff_t, 2> x_stencil,
                           std::array<std::ptrdiff_t, 2> y_stencil) {
  return Shrink(Shrink(box, Direction::X, x_stencil), Direction::Y, y_stencil);
}

void AverageCellToFace(::amrex::MultiFab& mf_faces, int face_component,
                       const ::amrex::MultiFab& mf_cells, int cell_component,
                       Direction dir) {
  if constexpr (AMREX_SPACEDIM == 2) {
    FUB_ASSERT(dir == Direction::X || dir == Direction::Y);
    if (dir == Direction::X) {
      ForEachFab(mf_cells, [&](const ::amrex::MFIter& mfi) {
        auto cells = MakePatchDataView(mf_cells[mfi], cell_component);
        auto all_faces = MakePatchDataView(mf_faces[mfi], face_component);
        // We are cautious and only compute the average for faces which exist
        // and are in range for the cells which exist on our grid
        auto cell_box_for_stencil = BoxFromStencil(cells.Box(), {1, 0}, {1, 1});
        auto face_box_for_stencil = Shrink(cell_box_for_stencil, dir, {0, 1});
        auto face_box = Intersect(all_faces.Box(), face_box_for_stencil);
        auto faces = all_faces.Subview(face_box);
        ForEachIndex(faces.Box(), [&](std::ptrdiff_t i, std::ptrdiff_t j) {
          // clang-format off
          faces(i, j) = 1.0 * cells(i - 1, j - 1) + 1.0 * cells(i, j - 1) +
                        2.0 * cells(i - 1,     j) + 2.0 * cells(i,     j) +
                        1.0 * cells(i - 1, j + 1) + 1.0 * cells(i, j + 1);
          faces(i, j) *= 0.125;
          // clang-format on
        });
      });
    } else {
      ForEachFab(mf_cells, [&](const ::amrex::MFIter& mfi) {
        auto cells = MakePatchDataView(mf_cells[mfi], cell_component);
        auto all_faces = MakePatchDataView(mf_faces[mfi], face_component);
        // We are cautious and only compute the average for faces which exist
        // and are in range for the cells which exist on our grid
        auto cell_box_for_stencil = BoxFromStencil(cells.Box(), {1, 1}, {1, 0});
        auto face_box_for_stencil = Shrink(cell_box_for_stencil, dir, {0, 1});
        auto face_box = Intersect(all_faces.Box(), face_box_for_stencil);
        auto faces = all_faces.Subview(face_box);
        ForEachIndex(faces.Box(), [&](std::ptrdiff_t i, std::ptrdiff_t j) {
          // clang-format off
          faces(i, j) = 1.0 * cells(i - 1,     j) + 2.0 * cells(i,     j) + 1.0 * cells(i + 1,     j) +
                        1.0 * cells(i - 1, j - 1) + 2.0 * cells(i, j - 1) + 1.0 * cells(i + 1, j - 1);
          faces(i, j) *= 0.125;
          // clang-format on
        });
      });
    }
  }
}

void AverageCellToNode(::amrex::MultiFab& mf_nodes, int node_component,
                       const ::amrex::MultiFab& mf_cells, int cell_component) {
  if constexpr (AMREX_SPACEDIM == 2) {
    ForEachFab(mf_cells, [&](const ::amrex::MFIter& mfi) {
      auto cells = MakePatchDataView(mf_cells[mfi], cell_component);
      auto nodes = MakePatchDataView(mf_nodes[mfi], node_component);
      ForEachIndex(nodes.Box(), [&](std::ptrdiff_t i, std::ptrdiff_t j) {
        // clang-format off
        nodes(i, j) = 0.25 * cells(i - 1, j - 1) + 0.25 * cells(i, j - 1) +
                      0.25 * cells(i - 1,     j) + 0.25 * cells(i,     j);
        // clang-format on
      });
    });
  }
}

template <int Rank, int VelocityRank>
void ComputeKCrossMomentum(
    ::amrex::MultiFab& result, const ::amrex::MultiFab& scratch,
    const std::array<double, 3>& k,
    const IndexMapping<CompressibleAdvection<Rank, VelocityRank>>& index);

template <int Rank, int VelocityRank>
std::vector<::amrex::MultiFab>
ComputeEulerBackwardRHSAndSetAlphaForLinearOperator(
    BK19Solver<Rank, VelocityRank>& solver, Duration dt,
    const BK19PhysicalParameters& physical_parameters,
    DebugSnapshotProxy dbg_sn = DebugSnapshotProxy()) {
  CompressibleAdvectionIntegratorContext& context =
      solver.GetAdvectionSolver().GetContext();
  const IndexMapping index = solver.GetEquation().GetIndexMapping();
  const int nlevel = context.GetPatchHierarchy().GetNumberOfLevels();
  // first compute the divergence term
  // vector field needs one ghost cell width to compute divergence!
  std::vector<::amrex::MultiFab> rhs_hierarchy =
      ZerosOnNodes(context, one_component, no_ghosts);
  std::vector<::amrex::MultiFab> UV =
      ComputePvFromScratch(solver, one_ghost_cell_width);
  solver.GetLinearOperator()->compDivergence(ToPointers(rhs_hierarchy),
                                             ToPointers(UV));
  dbg_sn.SaveData(ToConstPointers(rhs_hierarchy), "rhs",
                  GetGeometries(context));
  // Now add the diagonal term to rhs if alpha_p > 0.0
  // We also set the alpha coefficients for the linear operator in this case
  if (physical_parameters.alpha_p > 0.0) {
    const double factor1 = -physical_parameters.alpha_p *
                           physical_parameters.Msq /
                           (physical_parameters.gamma - 1.0) / dt.count();
    FUB_ASSERT(factor1 < 0.0);
    std::vector<::amrex::MultiFab> diagonal_term_on_nodes =
        AllocateOnNodes(context, one_component, no_ghosts);
    for (int ilvl = 0; ilvl < nlevel; ++ilvl) {
      const std::size_t level = static_cast<std::size_t>(ilvl);
      const ::amrex::MultiFab& scratch = context.GetScratch(ilvl);
      ::amrex::MultiFab diagonal_term_on_cells =
          AllocateOnCells(context, ilvl, one_component, one_ghost_cell_width);
      ForEachFab(
          execution::openmp, diagonal_term_on_cells,
          [&](const ::amrex::MFIter& mfi) {
            ::amrex::Box tilebox = mfi.growntilebox();
            auto diagfac =
                MakePatchDataView(diagonal_term_on_cells[mfi], 0, tilebox);
            auto PTdensity =
                MakePatchDataView(scratch[mfi], index.PTdensity, tilebox);
            ForEachRow(std::tuple{diagfac, PTdensity},
                       [physical_parameters, factor1,
                        dt](span<double> dfac, span<const double> PTdens) {
                         for (std::ptrdiff_t i = 0; i < dfac.size(); ++i) {
                           dfac[i] = factor1 *
                                     std::pow(PTdens[i],
                                              2.0 - physical_parameters.gamma);
                         }
                       });
          });
      AverageCellToNode(diagonal_term_on_nodes[level], 0,
                        diagonal_term_on_cells, 0);
      // This copies the termns into a local array within the linear operator
      solver.GetLinearOperator()->setAlpha(ilvl, diagonal_term_on_nodes[level]);
    }
    dbg_sn.SaveData(ToConstPointers(diagonal_term_on_nodes), "alpha",
                    GetGeometries(context));
    for (int ilvl = 0; ilvl < nlevel; ++ilvl) {
      // Now we weigth the diagonal term with our old solution and add this to
      // the RHS.
      const std::size_t level = static_cast<std::size_t>(ilvl);
      const ::amrex::MultiFab& pi = context.GetPi(ilvl);
      ::amrex::MultiFab::Multiply(diagonal_term_on_nodes[level], pi, 0, 0,
                                  one_component, no_ghosts);
      ::amrex::MultiFab::Add(rhs_hierarchy[level],
                             diagonal_term_on_nodes[level], 0, 0, one_component,
                             no_ghosts);
    }
  }
  dbg_sn.SaveData(ToConstPointers(rhs_hierarchy), "rhs_plus_alpha",
                  GetGeometries(context));
  return rhs_hierarchy;
}

template <int Rank, int VelocityRank>
std::vector<::amrex::MultiFab>
SetSigmaForLinearOperator(BK19Solver<Rank, VelocityRank>& solver, Duration dt,
                          const BK19PhysicalParameters& physical_parameters,
                          DebugSnapshotProxy dbg_sn = DebugSnapshotProxy()) {
  CompressibleAdvectionIntegratorContext& context =
      solver.GetAdvectionSolver().GetContext();
  const int nlevel = context.GetPatchHierarchy().GetNumberOfLevels();
  auto index = solver.GetEquation().GetIndexMapping();
  std::vector<::amrex::MultiFab> sigma =
      AllocateOnCells(context, one_component, no_ghosts);
  for (int level = 0; level < nlevel; ++level) {
    const ::amrex::MultiFab& scratch = context.GetScratch(level);
    sigma[level].setVal(
        physical_parameters.c_p * dt.count() /
        (1.0 + std::pow(dt.count() * physical_parameters.f, 2)));
    ::amrex::MultiFab::Multiply(sigma[level], scratch, index.PTdensity, 0,
                                one_component, sigma[level].nGrow());
    ::amrex::MultiFab::Divide(sigma[level], scratch, index.PTinverse, 0,
                              one_component, sigma[level].nGrow());
    // set weights in linear operator
    solver.GetLinearOperator()->setSigma(level, sigma[level]);
  }
  dbg_sn.SaveData(ToConstPointers(sigma), "sigma", GetGeometries(context));
  return sigma;
}

template <int Rank, int VelocityRank>
void SetSigmaCrossForLinearOperator(
    BK19Solver<Rank, VelocityRank>& solver, span<const ::amrex::MultiFab> sigma,
    Duration dt, const BK19PhysicalParameters& physical_parameters,
    DebugSnapshotProxy dbg_sn = DebugSnapshotProxy()) {
  CompressibleAdvectionIntegratorContext& context =
      solver.GetAdvectionSolver().GetContext();
  const int nlevel = context.GetPatchHierarchy().GetNumberOfLevels();
  std::vector<::amrex::MultiFab> sigmacross = AllocateOnCells(
      context, AMREX_SPACEDIM * (AMREX_SPACEDIM - 1), no_ghosts);
  for (int level = 0; level < nlevel; ++level) {
    for (int i = 0; i < AMREX_SPACEDIM * (AMREX_SPACEDIM - 1); ++i) {
      ::amrex::MultiFab::Copy(sigmacross[level], sigma[level], 0, i,
                              one_component, no_ghosts);
      const double facsign = std::pow(-1.0, i);
      sigmacross[level].mult(facsign * dt.count() * physical_parameters.f, i,
                             1);
    }
    solver.GetLinearOperator()->setSigmaCross(level, sigmacross[level]);
  }
  dbg_sn.SaveData(ToConstPointers(sigmacross),
                  DebugSnapshot::ComponentNames{"sigmac0", "sigmac1"},
                  GetGeometries(context));
}

template <int Rank, int VelocityRank>
void ApplyDivergenceCorrectionOnScratch(
    BK19Solver<Rank, VelocityRank>& solver, Duration dt,
    span<::amrex::MultiFab> UV_correction,
    const BK19PhysicalParameters& physical_parameters,
    DebugSnapshotProxy dbg_sn = DebugSnapshotProxy()) {
  CompressibleAdvectionIntegratorContext& context =
      solver.GetAdvectionSolver().GetContext();
  const int nlevel = context.GetPatchHierarchy().GetNumberOfLevels();
  auto index = solver.GetEquation().GetIndexMapping();
  for (int level = 0; level < nlevel; ++level) {
    std::array<int, VelocityRank> UV_index{};
    std::iota(UV_index.begin(), UV_index.end(), 0); // fill with 0, 1, 2

    const double factor1 = -dt.count() * physical_parameters.f;
    const double factor2 = 1.0;
    AddCrossProduct(UV_correction[level], UV_index, UV_correction[level],
                    UV_index, factor1, factor2, physical_parameters.k_vect);

    ::amrex::MultiFab& scratch = context.GetScratch(level);
    for (std::size_t i = 0; i < index.momentum.size(); ++i) {
      const int UV_component = static_cast<int>(i);
      ::amrex::MultiFab::Multiply(UV_correction[level], scratch,
                                  index.PTinverse, UV_component, one_component,
                                  no_ghosts);
      // UV_correction is now a momentum correction. Thus add it.
      ::amrex::MultiFab::Add(scratch, UV_correction[level], UV_component,
                             index.momentum[i], one_component, no_ghosts);
    }
  }
  dbg_sn.SaveData(
      ToConstPointers(UV_correction),
      DebugSnapshot::ComponentNames{"Momentum_corr0", "Momentum_corr1"},
      GetGeometries(context));
}

template <int Rank, int VelocityRank>
void ApplyPiCorrectionOnScratch(
    BK19Solver<Rank, VelocityRank>& solver, span<::amrex::MultiFab> div,
    Duration dt, DebugSnapshotProxy dbg_sn = DebugSnapshotProxy()) {
  CompressibleAdvectionIntegratorContext& context =
      solver.GetAdvectionSolver().GetContext();
  dbg_sn.SaveData(ToConstPointers(div), "div_Pv", GetGeometries(context));
  const int nlevel = context.GetPatchHierarchy().GetNumberOfLevels();
  const BK19PhysicalParameters& phys_param = solver.GetPhysicalParameters();
  auto index = solver.GetEquation().GetIndexMapping();
  for (int level = 0; level < nlevel; ++level) {
    const ::amrex::MultiFab& scratch = context.GetScratch(level);
    ::amrex::MultiFab dpidP_cells =
        AllocateOnCells(context, level, one_component, one_ghost_cell_width);
    ForEachFab(execution::openmp, dpidP_cells, [&](const ::amrex::MFIter& mfi) {
      ::amrex::Box tilebox = mfi.growntilebox();
      auto dpidP = MakePatchDataView(dpidP_cells[mfi], 0, tilebox);
      auto PTdensity =
          MakePatchDataView(scratch[mfi], index.PTdensity, tilebox);
      ForEachRow(std::tuple{dpidP, PTdensity},
                 [phys_param, dt](span<double> dpi, span<const double> PTdens) {
                   for (std::ptrdiff_t i = 0; i < dpi.size(); ++i) {
                     dpi[i] = -dt.count() * (phys_param.gamma - 1.0) /
                              phys_param.Msq *
                              std::pow(PTdens[i], phys_param.gamma - 2.0);
                   }
                 });
    });
    ::amrex::MultiFab dpidP_nodes =
        AllocateOnNodes(context, level, one_component, no_ghosts);
    AverageCellToNode(dpidP_nodes, 0, dpidP_cells, 0);
    ::amrex::MultiFab::Multiply(div[level], dpidP_nodes, 0, 0, one_component,
                                no_ghosts);
    ::amrex::MultiFab& pi = context.GetPi(level);
    ::amrex::MultiFab::Add(pi, div[level], 0, 0, one_component, no_ghosts);
  }
  dbg_sn.SaveData(ToConstPointers(div), "Pi_correction",
                  GetGeometries(context));
}

namespace {
/// \brief Recompute the state variable rho_chi_fast from the other variables on
/// the scratch of the specified context
///
/// The formula to recompute rho chi fast is
///
///           rho_chi_fast = rho / (PT - S0) * dS0_dy
template <int Rank, int VelocityRank>
void SnychronizeRhoChiFastWithIntegratorContext(
    CompressibleAdvectionIntegratorContext& context, 
    IndexMapping<CompressibleAdvection<Rank, VelocityRank>> index, Duration dt,
    span<const ::amrex::MultiFab> S0s) {
  const int n_levels = context.GetPatchHierarchy().GetNumberOfLevels();
  for (int level = 0; level < n_levels; ++level) {
    const ::amrex::Geometry& geom = context.GetGeometry(level);
    const double dy = geom.CellSize(1);
    ::amrex::MultiFab& scratch = context.GetScratch(level);
    const ::amrex::MultiFab& S0c = S0s[level];
    ForEachFab(S0c, [&](const ::amrex::MFIter& mfi) {
      const ::amrex::FArrayBox& S0 = S0c[mfi];
      ::amrex::FArrayBox& cells = scratch[mfi];
      ForEachIndex(cells.box(), [&](auto... is) {
        const ::amrex::IntVect iv(int(is)...);                      // iv + {0, 0}
        const ::amrex::IntVect upper = Shift(iv, Direction::Y, 1);  // iv + {0,  1}
        const ::amrex::IntVect lower = Shift(iv, Direction::Y, -1); // iv + {0, -1}
        const double dSdy = 0.5 * (S0c(upper) - S0c(lower)) / dy;
        // Sol.rhoX[i1] = (Sol.rho[i1] * (Sol.rho[i1] / Sol.rhoY[i1] - S0c)) + dt * (-v * dSdy) * Sol.rho[i1]
        const double rho = cells(iv, index.density);
        const double v = cells(iv, index.momentum[1]) / rho;
        const double PTdensity = cells(iv, index.PTdensity);
        const double rho_chi_fast = rho * (rho / PTdensity) + dt.count() * (-v * dSdy) * rho;
        cells(iv, index.rho_chi_fast) = rho_chi_fast;
      });
    });
  }
}
} // namespace

template <int Rank, int VelocityRank>
void DoEulerForward(BK19Solver<Rank, VelocityRank>& solver, Duration dt,
                    DebugSnapshotProxy dbg_sn = DebugSnapshotProxy()) {
  CompressibleAdvectionIntegratorContext& context =
      solver.GetAdvectionSolver().GetContext();

  auto index = solver.GetEquation().GetIndexMapping();

  std::vector<::amrex::MultiFab> UV =
      ComputePvFromScratch(solver, one_ghost_cell_width);

  std::vector<::amrex::MultiFab> div =
      ZerosOnNodes(context, one_component, no_ghosts);
  solver.GetLinearOperator()->compDivergence(ToPointers(div), ToPointers(UV));

  const std::shared_ptr<::amrex::MLNodeHelmholtz>& linear_operator =
      solver.GetLinearOperator();

  const BK19PhysicalParameters& physical_parameters =
      solver.GetPhysicalParameters();

  // construct sigma as in DoEulerBackward
  SetSigmaForLinearOperator(solver, dt, physical_parameters, dbg_sn);

  // Add explicit source terms to momentum on scratch
  std::vector<::amrex::MultiFab> source_terms =
      ZerosOnCells(context, AMREX_SPACEDIM, no_ghosts);
  std::array<int, AMREX_SPACEDIM> source_terms_index{AMREX_D_DECL(0, 1, 2)};

  // this computes: -sigma Grad(pi)
  linear_operator->getFluxes(ToPointers(source_terms), GetPis(context));

  SnychronizeRhoChiFastWithStates(context, index, dt, S0s);
  const int nlevel = context.GetPatchHierarchy().GetNumberOfLevels();
  const double factor1 = -dt.count() * physical_parameters.f;
  const double factor2 = 1.0;
  for (int ilvl = 0; ilvl < nlevel; ++ilvl) {
    const std::size_t level = static_cast<std::size_t>(ilvl);
    const ::amrex::MultiFab& scratch = context.GetScratch(ilvl);
    // Add coriolis terms to source_terms multifab
    AddCrossProduct(source_terms[level], source_terms_index, scratch,
                    index.momentum, factor1, factor2,
                    physical_parameters.k_vect);
    // Add gravity to source_terms
    AddGravity(source_terms[level], source_terms_index, scratch, index,
               physical_parameters, dt);
    // Add all source_terms to current momentum on scratch
    ::amrex::MultiFab::Add(context.GetScratch(ilvl), source_terms[level], 0,
                           index.momentum[0], index.momentum.size(), no_ghosts);
  }

  // compute update for pi (compressible case), (equation (16) in [BK19])

  if (physical_parameters.alpha_p > 0.0) {
    ApplyPiCorrectionOnScratch(solver, div, dt, dbg_sn);
  }
}

template <int Rank, int VelocityRank>
std::vector<::amrex::MultiFab>
DoEulerBackward(BK19Solver<Rank, VelocityRank>& solver, Duration dt,
                const BK19PhysicalParameters& physical_parameters,
                DebugSnapshotProxy dbg_sn = DebugSnapshotProxy()) {
  // If we play with a non-trivial coriolis term we need an explicit source
  // term integration step here
  ApplyExplicitSourceTerms(solver, dt, physical_parameters);
  solver.FillAllGhostLayers();

  const std::vector<::amrex::MultiFab> rhs =
      ComputeEulerBackwardRHSAndSetAlphaForLinearOperator(
          solver, dt, physical_parameters, dbg_sn);

  const std::vector<::amrex::MultiFab> sigma =
      SetSigmaForLinearOperator(solver, dt, physical_parameters, dbg_sn);

  SetSigmaCrossForLinearOperator(solver, sigma, dt, physical_parameters,
                                 dbg_sn);

  const CompressibleAdvectionIntegratorContext& context =
      solver.GetAdvectionSolver().GetContext();

  // solve elliptic equation for pi
  // weiqun recommends to use one ghost cell width to prevent internal copies
  std::vector<::amrex::MultiFab> pi =
      ZerosOnNodes(context, one_component, one_ghost_cell_width);

  std::shared_ptr linear_operator = solver.GetLinearOperator();
  const BK19SolverOptions& options = solver.GetSolverOptions();
  {
    Timer _ = context.GetCounterRegistry()->get_timer(
        "BK19LevelIntegrator::EulerBackward::solve");
    ::amrex::MLMG nodal_solver(*linear_operator);
    nodal_solver.setMaxIter(options.mlmg_max_iter);
    nodal_solver.setVerbose(options.mlmg_verbose);
    nodal_solver.setBottomVerbose(options.bottom_verbose);
    nodal_solver.setBottomMaxIter(options.bottom_max_iter);
    nodal_solver.setBottomToleranceAbs(options.bottom_tolerance_abs);
    nodal_solver.setBottomTolerance(options.bottom_tolerance_rel);
    nodal_solver.setAlwaysUseBNorm(options.always_use_bnorm);
    nodal_solver.solve(ToPointers(pi), ToConstPointers(rhs),
                       options.mlmg_tolerance_rel, options.mlmg_tolerance_abs);
  }
  dbg_sn.SaveData(ToConstPointers(pi), "solution", GetGeometries(context));

  // compute momentum correction
  std::vector<::amrex::MultiFab> UV_correction =
      ZerosOnCells(context, AMREX_SPACEDIM, no_ghosts);

  // this computes: -sigma Grad(pi)
  linear_operator->getFluxes(ToPointers(UV_correction), ToPointers(pi));

  // Given this solution we apply the correction for the momentum on the
  // current scratch.
  // This correction also takes the coriolis force into consideration.
  ApplyDivergenceCorrectionOnScratch(solver, dt, UV_correction,
                                     physical_parameters, dbg_sn);

  return pi;
}

} // namespace

///////////////////////////////////////////////////////////////////////////////
//                                                             IMPLEMENTATION
///////////////////////////////////////////////////////////////////////////////

template <int Rank, int VelocityRank>
BK19Solver<Rank, VelocityRank>::BK19Solver(
    const Equation& equation, AdvectionSolver advection,
    std::shared_ptr<::amrex::MLNodeHelmholtz> linop,
    const BK19PhysicalParameters& physical_parameters,
    const BK19SolverOptions& options)
    : equation_{equation},
      physical_parameters_{physical_parameters}, options_{options},
      advection_{std::move(advection)}, lin_op_{std::move(linop)} {}

template <int Rank, int VelocityRank>
void BK19Solver<Rank, VelocityRank>::RecomputeAdvectiveFluxes() {
  // Get all neccessary objects
  IndexMapping index = equation_.GetIndexMapping();
  const int nlevel =
      GetGriddingAlgorithm()->GetPatchHierarchy().GetNumberOfLevels();
  for (int level = 0; level < nlevel; ++level) {
    IntegratorContext& context = advection_.GetContext();
    const ::amrex::MultiFab& scratch = context.GetScratch(level);
    ::amrex::MultiFab& Pv_cells = context.GetAdvectiveFluxes(level).on_cells;
    std::array<::amrex::MultiFab, AMREX_SPACEDIM>& Pv_faces =
        context.GetAdvectiveFluxes(level).on_faces;

    // Here we assume that boundary conditions are satisified on scratch
    ComputePvFromScratch(Pv_cells, scratch, Pv_cells.nGrow(), index);

    // Average Pv_i for each velocity direction
    // TODO Something is off here if Rank != Velocity Rank
    constexpr int face_component = 0;
    for (std::size_t dir = 0; dir < std::size_t(Rank); ++dir) {
      const int cell_component = static_cast<int>(dir);
      AverageCellToFace(Pv_faces[dir], face_component, Pv_cells, cell_component,
                        Direction(dir));
    }
  }
}

template <int Rank, int VelocityRank>
void BK19Solver<Rank, VelocityRank>::DoInitialProjection() {
  std::shared_ptr counters = advection_.GetContext().GetCounterRegistry();
  Timer _ = counters->get_timer("BK19Solver::DoInitialProjection");
  BK19PhysicalParameters phys_param_aux{physical_parameters_};
  phys_param_aux.alpha_p = 0.0;
  phys_param_aux.f = 0.0;
  FillAllGhostLayers();
  DoEulerBackward(*this, Duration(1.0), phys_param_aux);
}

template <int Rank, int VelocityRank>
void BK19Solver<Rank, VelocityRank>::FillAllGhostLayers() {
  const int nlevel =
      GetGriddingAlgorithm()->GetPatchHierarchy().GetNumberOfLevels();
  IntegratorContext& context = advection_.GetContext();
  context.FillGhostLayerSingleLevel(0);
  for (int level = 1; level < nlevel; ++level) {
    context.FillGhostLayerTwoLevels(level, level - 1);
  }
}

template <int Rank, int VelocityRank>
Result<void, TimeStepTooLarge>
BK19Solver<Rank, VelocityRank>::AdvanceHierarchy(Duration dt) {
  std::shared_ptr<CounterRegistry> counters = advection_.GetCounterRegistry();
  Timer measure_everything =
      counters->get_timer("BK19LevelIntegrator::AdvanceLevelNonRecursively");
  AdvectionSolver& advection = GetAdvectionSolver();
  CompressibleAdvectionIntegratorContext& context = advection.GetContext();

  const int current_cycle = context.GetCycles(0);
  const Duration current_timepoint = context.GetTimePoint(0);

  const Duration half_dt = 0.5 * dt;

  // Save data on current time level for later use
  std::vector<::amrex::MultiFab> scratch_hierarchy =
      CopyScratchFromContext(context);

  // Get objects to debug this algorithm at various stages

  DebugStorage& debug = *context.GetPatchHierarchy().GetDebugStorage();

  DebugSnapshotProxy dbgPreStep = debug.AddSnapshot("BK19_pre-step");
  dbgPreStep.SaveData(GetScratches(context), GetCompleteVariableNames(),
                      GetGeometries(context));
  dbgPreStep.SaveData(GetPis(std::as_const(context)), "pi",
                      GetGeometries(context));

  // 1) Compute current Pv and interpolate to face centered quantity
  //    Current Pv is given by: Pv = PTdensity * momentum / density
  RecomputeAdvectiveFluxes();

  // 2) Do the advection with the face-centered Pv
  Result<void, TimeStepTooLarge> result = boost::outcome_v2::success();
  {
    Timer _ = counters->get_timer("BK19Solver::Advection_1");
    result = advection.AdvanceHierarchy(half_dt);

    DebugSnapshotProxy dbgAdvect = debug.AddSnapshot("BK19_advect");
    dbgAdvect.SaveData(GetScratches(context), GetCompleteVariableNames(),
                       GetGeometries(context));
    dbgAdvect.SaveData(GetPis(std::as_const(context)), "pi",
                       GetGeometries(context));
  }
  if (!result) {
    return result;
  }

  // 3) Do the first euler backward integration step for the source term
  {
    Timer _ = counters->get_timer("BK19Solver::EulerBackward_1");
    FillAllGhostLayers();

    DebugSnapshotProxy dbgBackward1 = debug.AddSnapshot("BK19_advect-backward");

    std::vector<::amrex::MultiFab> pi_aux =
        DoEulerBackward(*this, half_dt, physical_parameters_, dbgBackward1);

    // NOTE: the following update of pi in the pseudo-incompressible case is
    // not present in BK19, but a further development in the work of Ray Chow
    if (physical_parameters_.alpha_p == 0) {
      CopyPiToHierarchy(context, pi_aux);
    }

    dbgBackward1.SaveData(GetScratches(context), GetCompleteVariableNames(),
                          GetGeometries(context));
    dbgBackward1.SaveData(GetPis(std::as_const(context)), "pi",
                          GetGeometries(context));

    // 4) Recompute Pv at half time
    FillAllGhostLayers();
    RecomputeAdvectiveFluxes();

    dbgBackward1.SaveData(GetPvs(context),
                          DebugSnapshot::ComponentNames{"Pu", "Pv"},
                          GetGeometries(context));
  }

  // 5) Explicit Euler with old scratch data
  //   - We need a current pi_n here. What is the initial one?

  {
    // Copy data from old time level back to scratch
    Timer _ = counters->get_timer("BK19Solver::RestoreInitialScratchData");
    CopyScratchToContext(context, scratch_hierarchy);
    FillAllGhostLayers();
  }

  {
    Timer _ = counters->get_timer("BK19Solver::EulerForward");
    DebugSnapshotProxy dbgAdvBF =
        debug.AddSnapshot("BK19_advect-backward-forward");

    DoEulerForward(*this, half_dt, dbgAdvBF);

    dbgAdvBF.SaveData(GetScratches(context), GetCompleteVariableNames(),
                      GetGeometries(context));
    dbgAdvBF.SaveData(GetPvs(context),
                      DebugSnapshot::ComponentNames{"Pu", "Pv"},
                      GetGeometries(context));
    dbgAdvBF.SaveData(GetPis(std::as_const(context)), "pi",
                      GetGeometries(context));
  }

  // 6) Do the second advection step with half-time Pv and full time step
  //   - Currently, scratch contains the result of euler forward step,
  //     which started at the old time level.
  {
    Timer _ = counters->get_timer("BK19Solver::Advection_2");
    context.SetCycles(current_cycle, 0);
    context.SetTimePoint(current_timepoint, 0);
    FillAllGhostLayers();
    result = advection.AdvanceHierarchy(dt);

    DebugSnapshotProxy dbgAdvBFA =
        debug.AddSnapshot("BK19_advect-backward-forward-advect");
    dbgAdvBFA.SaveData(GetPvs(context),
                       DebugSnapshot::ComponentNames{"Pu", "Pv"},
                       GetGeometries(context));
    dbgAdvBFA.SaveData(GetScratches(context), GetCompleteVariableNames(),
                       GetGeometries(context));
    dbgAdvBFA.SaveData(GetPis(std::as_const(context)), "pi",
                       GetGeometries(context));
  }
  if (!result) {
    return result;
  }

  // 7) Do the second euler backward integration step for the source term
  {
    Timer _ = counters->get_timer("BK19Solver::EulerBackward_2");

    DebugSnapshotProxy dbgAdvBFAB =
        debug.AddSnapshot("BK19_advect-backward-forward-advect-backward");
    FillAllGhostLayers();
    std::vector<::amrex::MultiFab> pi_new =
        DoEulerBackward(*this, half_dt, physical_parameters_, dbgAdvBFAB);

    // Copy pi_n+1 to pi_n
    CopyPiToHierarchy(context, pi_new);

    for (int level = 0; level < context.GetPatchHierarchy().GetNumberOfLevels();
         ++level) {
      context.CopyScratchToData(level);
      context.SetCycles(current_cycle + 1, level);
      context.SetTimePoint(current_timepoint + dt, level);
    }

    dbgAdvBFAB.SaveData(GetScratches(context), GetCompleteVariableNames(),
                        GetGeometries(context));
    dbgAdvBFAB.SaveData(GetPis(std::as_const(context)), "pi",
                        GetGeometries(context));
  }

  return boost::outcome_v2::success();
}

template class BK19Solver<2, 2>;

} // namespace fub::amrex