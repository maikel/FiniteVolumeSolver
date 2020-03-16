// Copyright (c) 2019 Maikel Nadolski
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

#include "fub/AMReX/bk19/BK19LevelIntegrator.hpp"
#include "fub/AMReX/ForEachFab.hpp"
#include "fub/AMReX/ForEachIndex.hpp"
#include "fub/AMReX/MLMG/MLNodeHelmDualCstVel.hpp"
#include "fub/AMReX/bk19/BK19IntegratorContext.hpp"
#include "fub/AMReX/output/DebugOutput.hpp"
#include <AMReX_MLMG.H>

namespace fub::amrex {

DebugSnapshot::ComponentNames GetCompleteVariableNames() {
  using Equation = CompressibleAdvection<2>;
  fub::CompressibleAdvection<2> equation{};
  using Traits = StateTraits<Complete<Equation>>;
  constexpr auto names = Traits::names;
  const auto depths = Depths<Complete<Equation>>(equation);
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

void WriteBK19Plotfile::operator()(const GriddingAlgorithm& grid) const {
  const fub::amrex::PatchHierarchy& hier = grid.GetPatchHierarchy();
  std::string name = fmt::format("{}/plt{:09}", plotfilename, grid.GetCycles());
  const int nlevels = hier.GetNumberOfLevels();
  const double time_point = hier.GetTimePoint().count();
  FUB_ASSERT(nlevels >= 0);
  std::size_t size = static_cast<std::size_t>(nlevels);
  ::amrex::Vector<const ::amrex::MultiFab*> mf(size);
  ::amrex::Vector<const ::amrex::MultiFab*> mfnodes(size);
  ::amrex::Vector<::amrex::Geometry> geoms(size);
  ::amrex::Vector<int> level_steps(size);
  ::amrex::Vector<::amrex::IntVect> ref_ratio(size);
  for (std::size_t i = 0; i < size; ++i) {
    mf[i] = &hier.GetPatchLevel(static_cast<int>(i)).data;
    geoms[i] = hier.GetGeometry(static_cast<int>(i));
    level_steps[i] = static_cast<int>(hier.GetCycles(static_cast<int>(i)));
    ref_ratio[i] = hier.GetRatioToCoarserLevel(static_cast<int>(i));
  }

  std::vector<std::string> varnames = GetCompleteVariableNames();
  ::amrex::Vector<std::string> vnames(varnames.begin(), varnames.end());

  ::amrex::Vector<std::string> rfs{"raw_fields"};
  ::amrex::WriteMultiLevelPlotfile(name, nlevels, mf, vnames, geoms, time_point,
                                   level_steps, ref_ratio, "HyperCLaw-V1.1",
                                   "Level_", "Cell", rfs);

  // write nodal raw fields
  for (int level = 0; level < int(size); ++level) {
    if (hier.GetPatchLevel(level).nodes) {
      const ::amrex::MultiFab& nodes = *hier.GetPatchLevel(level).nodes;
      WriteRawField(name, "pi", nodes, level);
    }
  }
}

void WriteRawField(const std::string& path, const std::string& name,
                   const ::amrex::MultiFab& data, int level) {
  ::amrex::VisMF::SetHeaderVersion(::amrex::VisMF::Header::Version_v1);
  const std::string raw_pltname = fmt::format("{}/raw_fields", path);
  const std::string level_prefix = "Level_";
  const std::string full_prefix =
      ::amrex::MultiFabFileFullPrefix(level, raw_pltname, level_prefix, name);
  ::amrex::VisMF::Write(data, full_prefix);
}

using ::amrex::MFIter;
using ::amrex::MultiFab;

// We give names to some magic zeros and ones.
inline constexpr int no_ghosts = 0;
inline constexpr int one_ghost_cell_width = 1;
inline constexpr int one_component = 1;

inline constexpr int Rank = AMREX_SPACEDIM;

using Equation = CompressibleAdvection<Rank>;

// For now the implementation assumes Rank == 2
static_assert(Rank == 2);

namespace {
IndexBox<2> BoxFromStencil_(const IndexBox<2>& box,
                            std::array<std::ptrdiff_t, 2> x_stencil,
                            std::array<std::ptrdiff_t, 2> y_stencil) {
  return Shrink(Shrink(box, Direction::X, x_stencil), Direction::Y, y_stencil);
}

// Apply the cell to face average to a specified cell_component on mf_cells and
// write its result into face_component of mf_faces.
void AverageCellToFace_(MultiFab& mf_faces, int face_component,
                        const MultiFab& mf_cells, int cell_component,
                        Direction dir) {
  if constexpr (Rank == 2) {
    FUB_ASSERT(dir == Direction::X || dir == Direction::Y);
    if (dir == Direction::X) {
      ForEachFab(mf_cells, [&](const MFIter& mfi) {
        auto cells = MakePatchDataView(mf_cells[mfi], cell_component);
        auto all_faces = MakePatchDataView(mf_faces[mfi], face_component);
        // We are cautious and only compute the average for faces which exist
        // and are in range for the cells which exist on our grid
        auto cell_box_for_stencil =
            BoxFromStencil_(cells.Box(), {1, 0}, {1, 1});
        auto face_box_for_stencil = Shrink(cell_box_for_stencil, dir, {0, 1});
        auto face_box = Intersect(all_faces.Box(), face_box_for_stencil);
        auto faces = all_faces.Subview(face_box);
        ForEachIndex(faces.Box(), [&](int i, int j) {
          // clang-format off
          faces(i, j) = 1.0 * cells(i - 1, j - 1) + 1.0 * cells(i, j - 1) +
                        2.0 * cells(i - 1,     j) + 2.0 * cells(i,     j) +
                        1.0 * cells(i - 1, j + 1) + 1.0 * cells(i, j + 1);
          faces(i, j) *= 0.125;
          // clang-format on
        });
      });
    } else {
      ForEachFab(mf_cells, [&](const MFIter& mfi) {
        auto cells = MakePatchDataView(mf_cells[mfi], cell_component);
        auto all_faces = MakePatchDataView(mf_faces[mfi], face_component);
        // We are cautious and only compute the average for faces which exist
        // and are in range for the cells which exist on our grid
        auto cell_box_for_stencil =
            BoxFromStencil_(cells.Box(), {1, 1}, {1, 0});
        auto face_box_for_stencil = Shrink(cell_box_for_stencil, dir, {0, 1});
        auto face_box = Intersect(all_faces.Box(), face_box_for_stencil);
        auto faces = all_faces.Subview(face_box);
        ForEachIndex(faces.Box(), [&](int i, int j) {
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

// Apply the cell to node average to a specified cell_component on mf_cells and
// write its result into node_component of mf_nodes.
void AverageCellToNode_(MultiFab& mf_nodes, int node_component,
                        const MultiFab& mf_cells, int cell_component) {
  if constexpr (Rank == 2) {
    ForEachFab(mf_cells, [&](const MFIter& mfi) {
      auto cells = MakePatchDataView(mf_cells[mfi], cell_component);
      auto nodes = MakePatchDataView(mf_nodes[mfi], node_component);
      // We are cautious and only compute the average for nodes which exist
      // and are in range for the cells which exist on our grid
      ForEachIndex(nodes.Box(), [&](int i, int j) {
        // clang-format off
        nodes(i, j) = 0.25 * cells(i - 1, j - 1) + 0.25 * cells(i, j - 1) +
                      0.25 * cells(i - 1,     j) + 0.25 * cells(i,     j);
        // clang-format on
      });
    });
  }
}

void ComputePvFromScratch_(const IndexMapping<Equation>& index, MultiFab& dest,
                           const MultiFab& scratch,
                           const ::amrex::Periodicity& periodicity) {
  // Shall be: Pv[i] = PTdensity * v[i]
  // Compute Pv_i for each velocity direction
  for (std::size_t i = 0; i < index.momentum.size(); ++i) {
    const int dest_component = static_cast<int>(i);
    MultiFab::Copy(dest, scratch, index.PTdensity, dest_component,
                   one_component, no_ghosts);
    MultiFab::Multiply(dest, scratch, index.velocity[i], dest_component,
                       one_component, no_ghosts);
  }
  dest.FillBoundary(periodicity);
}
} // namespace

void RecomputeAdvectiveFluxes(const IndexMapping<Equation>& index,
                              std::array<MultiFab, Rank>& Pv_faces,
                              MultiFab& Pv_cells, const MultiFab& scratch,
                              const ::amrex::Periodicity& periodicity) {
  ComputePvFromScratch_(index, Pv_cells, scratch, periodicity);
  // Average Pv_i for each velocity direction
  constexpr int face_component = 0;
  for (std::size_t dir = 0; dir < index.velocity.size(); ++dir) {
    const int cell_component = static_cast<int>(dir);
    AverageCellToFace_(Pv_faces[dir], face_component, Pv_cells, cell_component,
                       Direction(dir));
  }
}

namespace {
Result<void, TimeStepTooLarge>
Advect_(BK19LevelIntegrator::AdvectionSolver& advection, int level, Duration dt,
        std::pair<int, int> subcycle) {
  return advection.AdvanceLevelNonRecursively(level, dt, subcycle);
}

void RecoverVelocityFromMomentum_(MultiFab& scratch,
                                  const IndexMapping<Equation>& index) {
  // MultiFab::Copy(dest, src, src_comp, dest_comp, n_comp, n_grow);
  // MultiFab::Divide(dest, src, src_comp, dest_comp, n_comp, n_grow);
  for (std::size_t i = 0; i < index.momentum.size(); ++i) {
    MultiFab::Copy(scratch, scratch, index.momentum[i], index.velocity[i],
                   one_component, no_ghosts);
    MultiFab::Divide(scratch, scratch, index.density, index.velocity[i],
                     one_component, no_ghosts);
  }
}

::amrex::MultiFab DoEulerBackward_(const IndexMapping<Equation>& index,
                                   ::amrex::MLNodeHelmDualCstVel& lin_op,
                                   const BK19PhysicalParameters& phys_param,
                                   const BK19LevelIntegratorOptions& options,
                                   BK19IntegratorContext& context, int level,
                                   Duration dt, DebugSnapshotProxy& dbg_sn) {

  MultiFab& scratch = context.GetScratch(level);
  const ::amrex::Geometry& geom = context.GetGeometry(level);
  const ::amrex::Periodicity periodicity = geom.periodicity();
  ::amrex::DistributionMapping distribution_map = scratch.DistributionMap();
  ::amrex::BoxArray on_cells = scratch.boxArray();
  ::amrex::BoxArray on_nodes = on_cells;
  on_nodes.surroundingNodes();

  // compute RHS for elliptic solve (equation (28) in [BK19] divided by -dt)
  // first compute diagonal part for compressibility
  MultiFab diagfac_cells(on_cells, distribution_map, one_component,
                         scratch.nGrow());
  MultiFab diagfac_nodes(on_nodes, distribution_map, one_component, no_ghosts);
  ForEachFab(diagfac_cells, [&](const MFIter& mfi) {
    span<double> diagfac = MakePatchDataView(diagfac_cells[mfi], 0).Span();
    span<double> PTdensity =
        MakePatchDataView(scratch[mfi], index.PTdensity).Span();
    FUB_ASSERT(diagfac.size() == PTdensity.size());
    for (std::ptrdiff_t i = 0; i < diagfac.size(); ++i) {
      diagfac[i] = -phys_param.alpha_p * phys_param.Msq / (phys_param.gamma - 1.0) *
                   std::pow(PTdensity[i], 2.0 - phys_param.gamma) / dt.count();
    }
  });
  AverageCellToNode_(diagfac_nodes, 0, diagfac_cells, 0);

  // get pi from scratch
  const MultiFab& pi_old = context.GetPi(level);
  MultiFab diagcomp(on_nodes, distribution_map, one_component, no_ghosts);
  MultiFab::Copy(diagcomp, pi_old, 0, 0, one_component, no_ghosts);
  MultiFab::Multiply(diagcomp, diagfac_nodes, 0, 0, one_component, no_ghosts);

  // then compute divergence term
  // This assumes f = 0 and N = 0
  //
  // vector field needs one ghost cell width to compute divergence!
  MultiFab UV(on_cells, distribution_map, index.momentum.size(),
              one_ghost_cell_width);
  ComputePvFromScratch_(index, UV, scratch, periodicity);

  MultiFab rhs(on_nodes, distribution_map, one_component, no_ghosts);
  rhs.setVal(0.0);
  lin_op.compDivergence({&rhs}, {&UV});
  dbg_sn.SaveData(rhs, "rhs", geom);

  // take sum of diagonal term and divergence
  MultiFab::Add(rhs, diagcomp, 0, 0, one_component, no_ghosts);

  // This is to check if RHS sums up to zero in psinc case:
  //  ::amrex::Box node_domain = geom.Domain();
  //  node_domain.surroundingNodes();
  //  node_domain.setBig(0, node_domain.bigEnd(0) - 1);
  //  node_domain.setBig(1, node_domain.bigEnd(1) - 1);
  //  double rhs_sum = 0.0;
  //  ForEachFab(rhs, [&](const ::amrex::MFIter& mfi) {
  //    const ::amrex::Box subbox = mfi.validbox() & node_domain;
  //    rhs_sum += rhs[mfi].sum(subbox, 0);
  //  });
  //  if (rhs_sum > 1e-6) {
  //    throw std::runtime_error("Fehler!");
  //  }

  // Construct sigma / weight of Laplacian (equation (27) in [BK19]
  // divided by -dt)
  MultiFab sigma(on_cells, distribution_map, one_component, no_ghosts);
  sigma.setVal(phys_param.c_p * dt.count());
  MultiFab::Multiply(sigma, scratch, index.PTdensity, 0, one_component,
                     sigma.nGrow());
  MultiFab::Divide(sigma, scratch, index.PTinverse, 0, one_component,
                   sigma.nGrow());
  dbg_sn.SaveData(sigma, "sigma", geom);

  // set weights in linear operator
  lin_op.setSigma(level, sigma);
  // We only set the diagonal term for alpha_p > 0, since this also tells then
  // linear operator that it is non-singular and the RHS must not summ up to
  // zero.
  if (phys_param.alpha_p > 0.0) {
    lin_op.setAlpha(level, diagfac_nodes);
  }

  // solve elliptic equation for pi
  MultiFab pi(on_nodes, distribution_map, one_component, no_ghosts);
  pi.setVal(0.0);

  ::amrex::MLMG nodal_solver(lin_op);
  nodal_solver.setMaxIter(options.mlmg_max_iter);
  nodal_solver.setVerbose(options.mlmg_verbose);
  nodal_solver.setBottomVerbose(options.bottom_verbose);
  nodal_solver.setBottomMaxIter(options.bottom_max_iter);
  nodal_solver.setBottomToleranceAbs(options.bottom_tolerance_abs);
  nodal_solver.setBottomTolerance(options.bottom_tolerance_rel);
  nodal_solver.setAlwaysUseBNorm(options.always_use_bnorm);

  nodal_solver.solve({&pi}, {&rhs}, options.mlmg_tolerance_rel,
                     options.mlmg_tolerance_abs);
  dbg_sn.SaveData(pi, "pi", geom);

  // compute momentum correction
  MultiFab UV_correction(on_cells, distribution_map, index.momentum.size(),
                         no_ghosts);
  UV_correction.setVal(0.0);
  // this computes: -sigma Grad(pi)
  lin_op.getFluxes({&UV_correction}, {&pi});

  for (std::size_t i = 0; i < index.momentum.size(); ++i) {
    const int UV_component = static_cast<int>(i);
    MultiFab::Multiply(UV_correction, scratch, index.PTinverse, UV_component,
                       one_component, no_ghosts);
    // UV_correction is now a momentum correction. Thus add it.
    MultiFab::Add(scratch, UV_correction, UV_component, index.momentum[i],
                  one_component, no_ghosts);
  }
  dbg_sn.SaveData(
      UV_correction,
      DebugSnapshot::ComponentNames{"Momentum_corr0", "Momentum_corr1"}, geom);

  RecoverVelocityFromMomentum_(scratch, index);

  return pi;
}

void DoEulerForward_(const IndexMapping<Equation>& index,
                     ::amrex::MLNodeHelmDualCstVel& lin_op,
                     const BK19PhysicalParameters& phys_param,
                     BK19IntegratorContext& context, int level, Duration dt,
                     DebugSnapshotProxy& dbg_sn) {
  MultiFab& scratch = context.GetScratch(level);
  const ::amrex::Geometry& geom = context.GetGeometry(level);
  const ::amrex::Periodicity periodicity = geom.periodicity();

  ::amrex::DistributionMapping distribution_map = scratch.DistributionMap();
  ::amrex::BoxArray on_cells = scratch.boxArray();

  // vector field needs one ghost cell width to compute divergence
  MultiFab UV(on_cells, distribution_map, index.momentum.size(),
              one_ghost_cell_width);
  ComputePvFromScratch_(index, UV, scratch, periodicity);

  // construct sigma as in DoEulerBackward_, but without potential temperature
  // factor, since we correct the momentum right away
  MultiFab sigma(on_cells, distribution_map, one_component, no_ghosts);
  sigma.setVal(phys_param.c_p * dt.count());
  MultiFab::Multiply(sigma, scratch, index.PTdensity, 0, one_component,
                     sigma.nGrow());

  // we only need to set sigma in the linear operator, since we only need it
  // for the flux correction
  dbg_sn.SaveData(sigma, "sigma", geom);
  lin_op.setSigma(level, sigma);

  // To compute the fluxes from the old pi we need one ghost cell width
  // Thus, we use periodic boundaries for now and redistribute the pi_n onto the
  // boundary
  // TODO: What happens to pi otherwise?
  MultiFab& pi = context.GetPi(level);
  MultiFab momentum_correction(on_cells, distribution_map,
                               index.momentum.size(), no_ghosts);
  momentum_correction.setVal(0.0);
  // this computes: -sigma Grad(pi)
  lin_op.getFluxes({&momentum_correction}, {&pi});

  MultiFab::Add(scratch, momentum_correction, 0, index.momentum[0],
                index.momentum.size(), no_ghosts);
  dbg_sn.SaveData(
      momentum_correction,
      DebugSnapshot::ComponentNames{"Momentum_corr0", "Momentum_corr1"}, geom);

  RecoverVelocityFromMomentum_(scratch, index);

  // compute update for pi (compressible case), (equation (16) in [BK19])
  ::amrex::BoxArray on_nodes = on_cells;
  on_nodes.surroundingNodes();
  MultiFab div(on_nodes, distribution_map, one_component, no_ghosts);
  div.setVal(0.0);
  lin_op.compDivergence({&div}, {&UV});

  MultiFab dpidP_cells(on_cells, distribution_map, one_component,
                       scratch.nGrow());
  MultiFab dpidP_nodes(on_nodes, distribution_map, one_component, no_ghosts);
  ForEachFab(dpidP_cells, [&](const MFIter& mfi) {
    span<double> dpidP = MakePatchDataView(dpidP_cells[mfi], 0).Span();
    span<double> PTdensity =
        MakePatchDataView(scratch[mfi], index.PTdensity).Span();
    FUB_ASSERT(dpidP.size() == PTdensity.size());
    for (std::ptrdiff_t i = 0; i < dpidP.size(); ++i) {
      dpidP[i] = -dt.count() * (phys_param.gamma - 1.0) / phys_param.Msq *
                 std::pow(PTdensity[i], phys_param.gamma - 2.0);
    }
  });
  AverageCellToNode_(dpidP_nodes, 0, dpidP_cells, 0);
  MultiFab::Multiply(div, dpidP_nodes, 0, 0, one_component, no_ghosts);
  // Note: It would be nice to just multiply by alpha_p=0 in the pseudo
  // incompressible limit. But that does not work, since before we divide by
  // Msq=0.
  if (phys_param.alpha_p > 0.0) {
    MultiFab::Add(pi, div, 0, 0, one_component, no_ghosts);
  }
}

} // namespace

BK19LevelIntegratorOptions::BK19LevelIntegratorOptions(
    const ProgramOptions& options) {
  mlmg_tolerance_rel =
      GetOptionOr(options, "mlmg_tolerance_rel", mlmg_tolerance_rel);
  mlmg_tolerance_abs =
      GetOptionOr(options, "mlmg_tolerance_abs", mlmg_tolerance_abs);
  mlmg_max_iter = GetOptionOr(options, "mlmg_max_iter", mlmg_max_iter);
  mlmg_verbose = GetOptionOr(options, "mlmg_verbose", mlmg_verbose);
  bottom_tolerance_rel =
      GetOptionOr(options, "bottom_tolerance_rel", bottom_tolerance_rel);
  bottom_tolerance_abs =
      GetOptionOr(options, "bottom_tolerance_abs", bottom_tolerance_abs);
  bottom_max_iter = GetOptionOr(options, "bottom_max_iter", bottom_max_iter);
  bottom_verbose = GetOptionOr(options, "bottom_verbose", bottom_verbose);
  always_use_bnorm = GetOptionOr(options, "always_use_bnorm", always_use_bnorm);
  prefix = GetOptionOr(options, "prefix", prefix);
  output_between_steps =
      GetOptionOr(options, "output_between_steps", output_between_steps);
}

BK19LevelIntegrator::BK19LevelIntegrator(
    const CompressibleAdvection<Rank>& equation, AdvectionSolver advection,
    std::shared_ptr<::amrex::MLNodeHelmDualCstVel> linop,
    const BK19PhysicalParameters& physical_parameters,
    const BK19LevelIntegratorOptions& options)
    : AdvectionSolver(std::move(advection)), phys_param_(physical_parameters),
      options_(options), equation_(equation), index_(equation_),
      lin_op_(std::move(linop)) {}

Result<void, TimeStepTooLarge>
BK19LevelIntegrator::AdvanceLevelNonRecursively(int level, Duration dt,
                                                std::pair<int, int> subcycle) {
  AdvectionSolver& advection = GetAdvection();
  BK19IntegratorContext& context = advection.GetContext();
  MultiFab& scratch = context.GetScratch(level);
  MultiFab& pi = context.GetPi(level);
  const ::amrex::Geometry& geom = context.GetGeometry(level);
  const ::amrex::Periodicity periodicity = geom.periodicity();

  DebugStorage& debug = *context.GetPatchHierarchy().GetDebugStorage();
  DebugSnapshotProxy dbgPreStep = debug.AddSnapshot("BK19_pre-step");
  DebugSnapshotProxy dbgAdvFlux = debug.AddSnapshot("BK19_advective_fluxes");
  DebugSnapshotProxy dbgAdv = debug.AddSnapshot("BK19_advect");
  DebugSnapshotProxy dbgAdvB = debug.AddSnapshot("BK19_advect-backward");
  DebugSnapshotProxy dbgAdvBF =
      debug.AddSnapshot("BK19_advect-backward-forward");
  DebugSnapshotProxy dbgAdvBFA =
      debug.AddSnapshot("BK19_advect-backward-forward-advect");
  DebugSnapshotProxy dbgAdvBFAB =
      debug.AddSnapshot("BK19_advect-backward-forward-advect-backward");

  dbgPreStep.SaveData(scratch, GetCompleteVariableNames(), geom);
  dbgPreStep.SaveData(pi, "pi", geom);

  // Save data on current time level for later use
  MultiFab scratch_aux(scratch.boxArray(), scratch.DistributionMap(),
                       scratch.nComp(), no_ghosts);
  scratch_aux.copy(scratch);
  BK19AdvectiveFluxes& Pv = context.GetAdvectiveFluxes(level);

  const Duration half_dt = 0.5 * dt;

  // 1) Compute current Pv and interpolate to face centered quantity
  //    Current Pv is given by: Pv = PTdensity * velocity
  RecomputeAdvectiveFluxes(index_, Pv.on_faces, Pv.on_cells, scratch,
                           periodicity);
  dbgAdvFlux.SaveData(scratch, GetCompleteVariableNames(), geom);
  dbgAdvFlux.SaveData(pi, "pi", geom);

  // 2) Do the advection with the face-centered Pv
  //    Open Question: Coarse Fine Boundary?
  //      - Need an option to do nothing there
  Result<void, TimeStepTooLarge> result =
      Advect_(advection, level, half_dt, subcycle);
  if (!result) {
    return result;
  }

  dbgAdv.SaveData(scratch, GetCompleteVariableNames(), geom);
  dbgAdv.SaveData(pi, "pi", geom);

  // 3) Do the first euler backward integration step for the source term
  context.FillGhostLayerSingleLevel(level);
  MultiFab pi_tmp = DoEulerBackward_(index_, *lin_op_, phys_param_, options_,
                                     context, level, half_dt, dbgAdvB);

  // 4) Recompute Pv at half time
  RecomputeAdvectiveFluxes(index_, Pv.on_faces, Pv.on_cells, scratch,
                           periodicity);

  dbgAdvB.SaveData(Pv.on_cells, DebugSnapshot::ComponentNames{"Pu", "Pv"},
                   geom);
  dbgAdvB.SaveData(Pv.on_faces[0], "Pu_faces", geom, ::amrex::SrcComp(0));
  dbgAdvB.SaveData(Pv.on_faces[1], "Pv_faces", geom, ::amrex::SrcComp(0));
  dbgAdvB.SaveData(scratch, GetCompleteVariableNames(), geom);

  // Copy data from old time level back to scratch
  scratch.copy(scratch_aux);
  context.FillGhostLayerSingleLevel(level);

  // 5) Explicit Euler with old scratch data
  //   - We need a current pi_n here. What is the initial one?
  DoEulerForward_(index_, *lin_op_, phys_param_, context, level, half_dt,
                  dbgAdvBF);

  dbgAdvBF.SaveData(scratch, GetCompleteVariableNames(), geom);
  dbgAdvBF.SaveData(pi, "pi", geom);

  // 6) Do the second advection step with half-time Pv and full time step
  //   - Currently, scratch contains the result of euler forward step,
  //     which started at the old time level.
  context.FillGhostLayerSingleLevel(level);
  result = Advect_(advection, level, dt, subcycle);
  if (!result) {
    return result;
  }

  dbgAdvBFA.SaveData(scratch, GetCompleteVariableNames(), geom);
  dbgAdvBFA.SaveData(pi, "pi", geom);

  // 6) Do the second euler backward integration step for the source term
  MultiFab pi_new = DoEulerBackward_(index_, *lin_op_, phys_param_, options_,
                                     context, level, half_dt, dbgAdvBFAB);

  // Copy pi_n+1 to pi_n
  context.GetPi(level).copy(pi_new);

  dbgAdvBFAB.SaveData(scratch, GetCompleteVariableNames(), geom);

  return boost::outcome_v2::success();
}

} // namespace fub::amrex
