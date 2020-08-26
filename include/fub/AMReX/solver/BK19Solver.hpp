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

#ifndef FUB_AMREX_BK19_SOLVER_HPP
#define FUB_AMREX_BK19_SOLVER_HPP

#include "fub/AMReX/CompressibleAdvectionIntegratorContext.hpp"
#include "fub/equations/CompressibleAdvection.hpp"

#include "fub/AMReX/MLMG/MLNodeHelmholtz.hpp"
#include "fub/ext/ProgramOptions.hpp"
#include "fub/solver/DimensionalSplitLevelIntegrator.hpp"
#include "fub/solver/NoSubcycleSolver.hpp"
#include "fub/split_method/SplittingMethod.hpp"

namespace fub::amrex {

struct BK19SolverOptions {
  BK19SolverOptions() = default;
  BK19SolverOptions(const ProgramOptions& map);

  void Print(SeverityLogger& log);

  bool do_initial_projection = true;
  double mlmg_tolerance_rel = 1.0e-4;
  double mlmg_tolerance_abs = -1.0;
  int mlmg_max_iter = 100;
  int mlmg_verbose = 0;
  double bottom_tolerance_rel = 1.0e-4;
  double bottom_tolerance_abs = -1.0;
  int bottom_max_iter = 20;
  int bottom_verbose = 0;
  int always_use_bnorm = 0;
  std::string prefix = "BK19LevelIntegrator";
  bool output_between_steps = false;
};

struct BK19PhysicalParameters {
  BK19PhysicalParameters() = default;
  BK19PhysicalParameters(const ProgramOptions& map);

  void Print(SeverityLogger& log);

  /// Specific gas constant
  double R_gas{287.4};

  /// Heat capacity ratio
  double gamma{1.4};

  /// Heat capacity at constant pressure
  double c_p{1006.0};

  /// Gravitational acceleration
  double g{10.0}; //  [m / s^2]

  /// Coriolis parameter in beta plane
  double f{0.0};
  std::array<double, 3> k_vect{0.0, 0.0, 1.0};

  double alpha_p{1.0};
  double Msq{0.0};
};

template <int Rank, int VelocityRank = Rank> class BK19Solver {
public:
  ////////////////////////////////////////////////////////////////////////////
  // Typedefs

  using Equation = CompressibleAdvection<Rank, VelocityRank>;
  using IndexMapping = ::fub::IndexMapping<Equation>;
  using Coordinates = Eigen::Matrix<double, Rank, 1>;
  using Complete = ::fub::Complete<Equation>;
  using Conservative = ::fub::Conservative<Equation>;
  using SplittingMethod = ::fub::AnySplitMethod;
  using IntegratorContext = CompressibleAdvectionIntegratorContext;
  using AdvectionLevelIntegrator =
      DimensionalSplitLevelIntegrator<Rank, IntegratorContext, SplittingMethod>;
  using AdvectionSolver = NoSubcycleSolver<AdvectionLevelIntegrator>;
  using LinearOperator = ::amrex::MLNodeHelmholtz;

  ////////////////////////////////////////////////////////////////////////////
  // Constructor

  BK19Solver(const Equation& equation, AdvectionSolver advection,
             std::shared_ptr<::amrex::MLNodeHelmholtz> linop,
             const BK19PhysicalParameters& physical_parameters,
             const BK19SolverOptions& options = BK19SolverOptions());

  ////////////////////////////////////////////////////////////////////////////
  // Accessors

  AdvectionSolver& GetAdvection() noexcept { return advection_; }
  const AdvectionSolver& GetAdvection() const noexcept { return advection_; }

  const std::shared_ptr<GriddingAlgorithm>&
  GetGriddingAlgorithm() const noexcept {
    return advection_.GetGriddingAlgorithm();
  }

  ////////////////////////////////////////////////////////////////////////////
  // Modifiers

  void ResetPatchHierarchy(std::shared_ptr<GriddingAlgorithm> grid) {
    advection_.ResetPatchHierarchy(std::move(grid));
  }

  Duration ComputeStableDt() const { return advection_.ComputeStableDt(); }

  Result<void, TimeStepTooLarge> AdvanceHierarchy(Duration dt);

  void DoInitialProjection();

  void RecomputeAdvectiveFluxes();

  void FillAllGhostLayers();

private:
  ////////////////////////////////////////////////////////////////////////////////
  // Euler Backward Helpers

  void ApplyExplicitCoriolisSourceTerm(Duration dt, double f,
                                       span<const double, 3> k);

  std::vector<::amrex::MultiFab>
  ComputeEulerBackwardRHSAndSetAlphaForLinearOperator(
      Duration dt, const BK19PhysicalParameters& physical_parameters);

  void
  SetSigmaForLinearOperator(Duration dt,
                            const BK19PhysicalParameters& physical_parameters);

  void ApplyDivergenceCorrectionToScratch(
      Duration dt, span<const ::amrex::MultiFab> UV_correction,
      const BK19PhysicalParameters& physical_parameters);

  std::vector<::amrex::MultiFab> DoEulerBackward(Duration dt);
  std::vector<::amrex::MultiFab>
  DoEulerBackward(Duration dt,
                  const BK19PhysicalParameters& physical_parameters);

  ////////////////////////////////////////////////////////////////////////////////
  // Euler Forward Helpers

  void DoEulerForward(Duration dt);

private:
  static constexpr int no_ghosts = 0;
  static constexpr int one_ghost_cell_width = 1;
  static constexpr int one_component = 1;

  BK19PhysicalParameters physical_parameters_;
  BK19SolverOptions options_;
  Equation equation_;
  AdvectionSolver advection_;
  std::shared_ptr<LinearOperator> lin_op_;
};

namespace bk19 {

/// \brief Apply the cell to face average to a specified cell_component on
/// mf_cells and write its result into face_component of mf_faces.
///
/// \tparam Rank  The spatial Rank for the BK19 scheme.
///
/// \param mf_faces  The face centered multifab which shall be written to.
///
/// \param face_component  The component within `mf_faces` which will be written
/// to.
///
/// \param mf_cells  The cell centered multifab which will be read from.
///
/// \param cell_component  The component within mf_cells which will be averaged.
///
/// \param dir  The face direction in which we are averaging to.
///
template <int Rank>
void AverageCellToFace(::amrex::MultiFab& mf_faces, int face_component,
                       const ::amrex::MultiFab& mf_cells, int cell_component,
                       Direction dir);

/// \brief Apply the cell to node average to a specified cell_component on
/// mf_cells and write its result into node_component of mf_nodes.
template <int Rank>
void AverageCellToNode(::amrex::MultiFab& mf_nodes, int node_component,
                       const ::amrex::MultiFab& mf_cells, int cell_component);

template <int Rank, int VelocityRank>
void ComputePvFromScratch(
    ::amrex::MultiFab& dest, const ::amrex::MultiFab& scratch,
    const IndexMapping<CompressibleAdvection<Rank, VelocityRank>>& index);

template <int Rank, int VelocityRank>
void ComputeKCrossMomentum(
    ::amrex::MultiFab& result, const ::amrex::MultiFab& scratch,
    const std::array<double, 3>& k,
    const IndexMapping<CompressibleAdvection<Rank, VelocityRank>>& index);

void ApplyExplicitCoriolisSourceTerm(std::array<span<double>, 2> rhou,
                                     double fac1, double fac2,
                                     span<const double, 3> k);
void ApplyExplicitCoriolisSourceTerm(std::array<span<double>, 3> rhou,
                                     double fac1, double fac2,
                                     span<const double, 3> k);

void ComputeEulerBackwardRHSAndSetAlphaForLinearOperator(
    span<double> dest, span<const double> PTdensity, double mach_squared,
    double alpha_p, Duration dt);

std::vector<::amrex::MultiFab>
CopyScratchFromContext(const IntegratorContext& context);

void CopyScratchToContext(IntegratorContext& context,
                          span<const ::amrex::MultiFab> scratch);

::amrex::Vector<::amrex::MultiFab*> ToPointers(span<::amrex::MultiFab> array);

} // namespace bk19

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
      GetGriddingAlgorithm().GetPatchHierarchy().GetNumberOfLevels();
  for (int level = 0; level < nlevel; ++level) {
    IntegratorContext& context = advection_.GetContext();
    const ::amrex::MultiFab& scratch = context.GetScratch(level);
    ::amrex::MultiFab& Pv_cells = context.GetAdvectiveFluxes(level).on_cells;
    std::array<::amrex::MultiFab, AMREX_SPACEDIM>& Pv_faces =
        context.GetAdvectiveFluxes(level).on_faces;

    // Here we assume that boundary conditions are satisified on scratch
    bk19::ComputePvFromScratch(Pv_cells, scratch, index);

    // Average Pv_i for each velocity direction
    // TODO Something is off here if Rank != Velocity Rank
    constexpr int face_component = 0;
    for (std::size_t dir = 0; dir < std::size_t(Rank); ++dir) {
      const int cell_component = static_cast<int>(dir);
      bk19::AverageCellToFace<Rank>(Pv_faces[dir], face_component, Pv_cells,
                                    cell_component, Direction(dir));
    }
  }
}

template <int Rank, int VelocityRank>
void BK19Solver<Rank, VelocityRank>::DoInitialProjection() {
  std::shared_ptr counters = advection_.GetContext().GetCounterRegistry();
  Timer _ = counters->get_timer("BK19Solver::DoInitialProjection");
  FillAllGhostLayers();
  DoEulerBackward(Duration(1.0));
}

template <int Rank, int VelocityRank>
void BK19Solver<Rank, VelocityRank>::FillAllGhostLayers() {
  const int nlevel =
      GetGriddingAlgorithm().GetPatchHierarchy().GetNumberOfLevels();
  IntegratorContext& context = advection_.GetContext();
  context.FillGhostLayerSingleLevel(0);
  for (int level = 1; level < nlevel; ++level) {
    context.FillGhostLayerTwoLevels(level, level - 1);
  }
}

template <int Rank, int VelocityRank>
std::vector<::amrex::MultiFab>
BK19Solver<Rank, VelocityRank>::DoEulerBackward(Duration dt) {
  return DoEulerBackward(dt, physical_parameters_);
}

template <int Rank, int VelocityRank>
void BK19Solver<Rank, VelocityRank>::ApplyExplicitCoriolisSourceTerm(
    Duration dt, double f, span<const double, 3> k) {
  IntegratorContext& context = advection_.GetContext();
  const int nlevel = context.GetPatchHierarchy().GetNumberOfLevels();
  IndexMapping index = equation_.GetIndexMapping();
  const double factor1 = -dt.count() * f;
  const double factor2 = 1.0 / (1.0 + factor1 * factor1);
  for (int level = 0; level < nlevel; ++level) {
    ::amrex::MultiFab& scratch = context.GetScratch(level);
    ForEachFab(execution::openmp, scratch, [&](const ::amrex::MFIter& mfi) {
      const ::amrex::Box box = mfi.growntilebox();
      if constexpr (VelocityRank == 2) {
        auto momentum_x =
            MakePatchDataView(scratch[mfi], index.momentum[0], box);
        auto momentum_y =
            MakePatchDataView(scratch[mfi], index.momentum[1], box);
        ForEachRow(std::tuple{momentum_x, momentum_y},
                   [factor1, factor2, k](span<double> rhou, span<double> rhov) {
                     bk19::ApplyExplicitCoriolisSourceTerm(
                         std::array<span<double>, 2>{rhou, rhov}, factor1,
                         factor2, k);
                   });
      } else if constexpr (VelocityRank == 3) {
        auto momentum_x =
            MakePatchDataView(scratch[mfi], index.momentum[0], box);
        auto momentum_y =
            MakePatchDataView(scratch[mfi], index.momentum[1], box);
        auto momentum_z =
            MakePatchDataView(scratch[mfi], index.momentum[2], box);
        ForEachRow(std::tuple{momentum_x, momentum_y, momentum_z},
                   [factor1, factor2, k](span<double> rhou, span<double> rhov,
                                         span<double> rhow) {
                     bk19::ApplyExplicitCoriolisSourceTerm(
                         std::array<span<double>, 3>{rhou, rhov, rhow}, factor1,
                         factor2, k);
                   });
      }
    });
  }
}

template <int Rank, int VelocityRank>
std::vector<::amrex::MultiFab> BK19Solver<Rank, VelocityRank>::
    ComputeEulerBackwardRHSAndSetAlphaForLinearOperator(
        Duration dt, const BK19PhysicalParameters& physical_parameters) {
  IntegratorContext& context = advection_.GetContext();
  const int nlevel = context.GetPatchHierarchy().GetNumberOfLevels();
  // first compute the divergence term
  // vector field needs one ghost cell width to compute divergence!
  std::vector<::amrex::MultiFab> rhs_hierarchy =
      ZerosOnNodes(one_component, no_ghosts);
  std::vector<::amrex::MultiFab> UV = ComputePvFromScratch(one_ghost_cells);
  lin_op_->compDivergence(ToPointers(rhs_hierarchy), ToPointers(UV));
  // Now add the diagonal term to rhs if alpha_p > 0.0
  // We also set the alpha coefficients for the linear operator in this case
  if (pyhsical_parameters.alpha_p > 0.0) {
    const double factor1 = -pyhsical_parameters.alpha_p *
                           pyhsical_parameters.Msq / (phys_param.gamma - 1.0) /
                           dt.count();
    FUB_ASSERT(factor1 < 0.0);
    for (int ilvl = 0; ilvl < nlevel; ++ilvl) {
      const std::size_t level = static_cast<std::size_t>(ilvl);
      ::amrex::MultiFab diagonal_term_on_cells =
          AllocateCells(ilvl, one_component, scratch.nGrow());
      ForEachFab(
          execution::openmp, diagonal_term_on_cells, [&](const MFIter& mfi) {
            ::amrex::Box tilebox = mfi.growntilebox();
            auto diagfac = MakePatchDataView(diagfac_cells[mfi], 0, tilebox);
            auto PTdensity =
                MakePatchDataView(scratch[mfi], index.PTdensity, tilebox);
            ForEachRow(
                std::tuple{diagfac, PTdensity},
                [phys_param, dt](span<double> dfac, span<double> PTdens) {
                  for (std::ptrdiff_t i = 0; i < dfac.size(); ++i) {
                    dfac[i] =
                        factor1 * std::pow(PTdens[i], 2.0 - phys_param.gamma);
                  }
                });
          });
      ::amrex::MultiFab diagonal_term_on_nodes =
          AllocateNodes(ilvl, one_component, no_ghosts);
      bk19::AverageCellToNode(diagonal_term_on_nodes, 0, diagonal_term_on_cells,
                              0);
      // This copies the termns into a local array within the linear operator
      lin_op_->setAlpha(ilvl, diagonal_term_on_nodes);

      // Now we weigth the diagonal term with our old solution and add this to
      // the RHS.
      const ::amrex::MultiFab& pi = context.GetPi(ilvl);
      ::amrex::MultiFab::Multiply(diagonal_term_on_nodes, pi, 0, 0,
                                  one_component, no_ghosts);
      ::amrex::MultiFab::Add(rhs_hierarchy[level], diagonal_term_on_nodes, 0, 0,
                             one_component, no_ghosts);
    }
  }
  return rhs;
}

template <int Rank, int VelocityRank>
std::vector<::amrex::MultiFab> BK19Solver<Rank, VelocityRank>::DoEulerBackward(
    Duration dt, const BK19PhysicalParameters& physical_parameters) {

  // If we play with a non-trivial coriolis term we need an explicit source
  // term integration step here
  if (physical_parameters.f > 0) {
    ApplyExplicitCoriolisSourceTerm(dt, physical_parameters.f,
                                    physical_parameters.k);
    FillAllGhostLayers();
  }

  std::vector<::amrex::MultiFab> rhs =
      ComputeEulerBackwardRHSAndSetAlphaForLinearOperator(dt,
                                                          physical_parameters);

  SetSigmaForLinearOperator(dt, physical_parameters);

  // solve elliptic equation for pi
  std::vector<::amrex::MultiFab> pi = ZerosOnNodes(context, one_component);

  {
    Timer _ = counters.get_timer("BK19LevelIntegrator::EulerBackward::solve");
    ::amrex::MLMG nodal_solver(lin_op);
    nodal_solver.setMaxIter(options_.mlmg_max_iter);
    nodal_solver.setVerbose(options_.mlmg_verbose);
    nodal_solver.setBottomVerbose(options_.bottom_verbose);
    nodal_solver.setBottomMaxIter(options_.bottom_max_iter);
    nodal_solver.setBottomToleranceAbs(options_.bottom_tolerance_abs);
    nodal_solver.setBottomTolerance(options_.bottom_tolerance_rel);
    nodal_solver.setAlwaysUseBNorm(options_.always_use_bnorm);
    nodal_solver.solve(AsPointers(pi), AsPointers(rhs),
                       options_.mlmg_tolerance_rel,
                       options_.mlmg_tolerance_abs);
  }

  // compute momentum correction
  std::vector<::amrex::MultiFab> UV_correction =
      ZerosOnCells(context, index.momentum.size());

  // this computes: -sigma Grad(pi)
  lin_op.getFluxes(AsPointers(UV_correction), AsPointers(pi));

  // Given this solution we apply the correction for the momentum on the
  // current scratch.
  // This correction also takes the coriolis force into consideration.
  ApplyDivergenceCorrectionToScratch(dt, UV_correction, physical_parameters);

  return pi;
}

template <int Rank, int VelocityRank>
void BK19Solver<Rank, VelocityRank>::DoEulerForward(Duration dt) {
  IndexMapping index = equation_.GetIndexMapping();
  IntegratorContext& context = advection_.GetContext();
  ::amrex::MLNodeHelmholtz& lin_op = *lin_op_.get();

  // vector field needs one ghost cell width to compute divergence
  //
  std::vector<::amrex::MultiFab> UV;
  for (int level = 0; level < nlevel; ++level) {
    const ::amrex::MultiFab& scratch = context.GetScratch(level);
    const ::amrex::BoxArray& on_cells = scratch.boxArray();

    UV.emplace_back(on_cells, distribution_map, index.momentum.size(),
                    scratch.nGrow());
    ComputePvFromScratch(UV, scratch, index);

    const ::amrex::Geometry& geom = context.GetGeometry(level);
    ::amrex::Box grown_box = geom.growNonPeriodicDomain(scratch.nGrow());
    ::amrex::BoxList boundaries =
        ::amrex::complementIn(grown_box, ::amrex::BoxList{geom.Domain()});
    for (::amrex::Box& box : boundaries) {
      ::amrex::Box ghostbox_with_corners =
          GrowInPeriodicDirection_(box, geom, scratch.nGrow());
      UV.setVal(0.0, ghostbox_with_corners, 0, index.momentum.size(),
                scratch.nGrow());
    }
  }

  // construct sigma as in DoEulerBackward_, but without potential temperature
  // factor, since we correct the momentum right away
  MultiFab sigma(on_cells, distribution_map, one_component, no_ghosts);
  sigma.setVal(phys_param.c_p * dt.count());
  MultiFab::Multiply(sigma, scratch, index.PTdensity, 0, one_component,
                     sigma.nGrow());

  // we only need to set sigma in the linear operator, since we only need it
  // for the flux correction
  dbg_sn.SaveData(sigma, "sigma", geom);
  lin_op_->setSigma(level, sigma);

  MultiFab momentum_correction(on_cells, distribution_map,
                               index.momentum.size(), no_ghosts);
  momentum_correction.setVal(0.0);
  // this computes: -sigma Grad(pi)
  lin_op_->getFluxes({&momentum_correction}, {&pi});

  MultiFab k_cross_rhou(on_cells, distribution_map, index.momentum.size(),
                        no_ghosts);
  ComputeKCrossM_(index.momentum, k_cross_rhou, phys_param.k_vect, scratch);
  const double fac = -dt.count() * phys_param.f;
  MultiFab::Saxpy(momentum_correction, fac, k_cross_rhou, 0, 0,
                  index.momentum.size(), no_ghosts);

  MultiFab::Add(scratch, momentum_correction, 0, index.momentum[0],
                index.momentum.size(), no_ghosts);
  dbg_sn.SaveData(
      momentum_correction,
      DebugSnapshot::ComponentNames{"Momentum_corr0", "Momentum_corr1"}, geom);

  // compute update for pi (compressible case), (equation (16) in [BK19])
  ::amrex::BoxArray on_nodes = on_cells;
  on_nodes.surroundingNodes();
  MultiFab div(on_nodes, distribution_map, one_component, no_ghosts);
  div.setVal(0.0);
  lin_op.compDivergence({&div}, {&UV});

  MultiFab dpidP_cells(on_cells, distribution_map, one_component,
                       scratch.nGrow());
  MultiFab dpidP_nodes(on_nodes, distribution_map, one_component, no_ghosts);
  ForEachFab(execution::openmp, dpidP_cells, [&](const MFIter& mfi) {
    ::amrex::Box tilebox = mfi.growntilebox();
    StridedDataView<double, AMREX_SPACEDIM> dpidP =
        MakePatchDataView(dpidP_cells[mfi], 0, tilebox);
    StridedDataView<double, AMREX_SPACEDIM> PTdensity =
        MakePatchDataView(scratch[mfi], index.PTdensity, tilebox);
    ForEachRow(std::tuple{dpidP, PTdensity},
               [phys_param, dt](span<double> dpi, span<double> PTdens) {
                 for (std::ptrdiff_t i = 0; i < dpi.size(); ++i) {
                   dpi[i] = -dt.count() * (phys_param.gamma - 1.0) /
                            phys_param.Msq *
                            std::pow(PTdens[i], phys_param.gamma - 2.0);
                 }
               });
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

Result<void, TimeStepTooLarge>
BK19Solver<Rank, VelocityRank>::AdvanceHierarchy(Duration dt) {
  std::shared_ptr<CounterRegistry> counters = advection_.GetCounterRegistry();
  Timer measure_everything =
      counters->get_timer("BK19LevelIntegrator::AdvanceLevelNonRecursively");
  AdvectionSolver& advection = GetAdvection();
  CompressibleAdvectionIntegratorContext& context = advection.GetContext();
  const fub::amrex::PatchHierarchy& hier = context.GetPatchHierarchy();
  ::amrex::MultiFab& scratch = context.GetScratch(level);
  ::amrex::MultiFab& pi = *hier.GetPatchLevel(level).nodes;
  const ::amrex::Geometry& geom = context.GetGeometry(level);
  const ::amrex::Periodicity periodicity = geom.periodicity();
  CompressibleAdvectionAdvectiveFluxes& Pv = context.GetAdvectiveFluxes(level);
  const Duration half_dt = 0.5 * dt;

  // Save data on current time level for later use
  std::vector<::amrex::MultiFab> scratch_hierarchy =
      CopyScratchFromContext(context);

  // 1) Compute current Pv and interpolate to face centered quantity
  //    Current Pv is given by: Pv = PTdensity * momentum / density
  RecomputeAdvectiveFluxes();

  // 2) Do the advection with the face-centered Pv
  Result<void, TimeStepTooLarge> result = boost::outcome_v2::success();
  {
    Timer _ = counters->get_timer("BK19Solver::Advection_1");
    result = advection.AdvanceHierarchy(half_dt);
  }
  if (!result) {
    return result;
  }

  // 3) Do the first euler backward integration step for the source term
  {
    Timer _ = counters->get_timer("BK19Solver::EulerBackward_1");
    FillAllGhostLayers();

    std::vector<::amrex::MultiFab> pi_aux = DoEulerBackward(dt);

    // NOTE: the following update of pi in the pseudo-incompressible case is
    // not present in BK19, but a further development in the work of Ray Chow
    if (phys_param_.alpha_p == 0) {
      for (int level = 0; level < nlevel; ++level) {
        hier.GetPatchLevel(level).nodes->copy(pi_aux);
      }
    }
  }

  // 4) Recompute Pv at half time
  FillAllGhostLayers();
  RecomputeAdvectiveFluxes();

  // 5) Explicit Euler with old scratch data
  //   - We need a current pi_n here. What is the initial one?

  {
    Timer _ = counters->get_timer("BK19Solver::RestoreInitialScratchData");
    for (int level = 0; level < nlevel; ++level) {
      ::amrex::MultiFab& scratch = context.GetScratch(level);
      const ::amrex::MultiFab& scratch_aux = scratch_hierarchy[level];
      scratch.copy(scratch_aux);
    }
    FillAllGhostLayers();
  }

  {
    Timer _ = counters->get_timer("BK19Solver::EulerForward");
    // Copy data from old time level back to scratch
    DoEulerForward(half_dt);
  }

  // 6) Do the second advection step with half-time Pv and full time step
  //   - Currently, scratch contains the result of euler forward step,
  //     which started at the old time level.
  {
    Timer _ = counters->get_timer("BK19Solver::Advection_2");
    FillAllGhostLayers();
    result = advection.AdvanceHierarchy(half_dt);
  }
  if (!result) {
    return result;
  }

  // 7) Do the second euler backward integration step for the source term
  {
    Timer _ = counters->get_timer("BK19Solver::EulerBackward_2");
    std::vector<::amrex::MultiFab> pi_new = DoEulerBackward(half_dt);

    // Copy pi_n+1 to pi_n
    for (int ilvl = 0; ilvl < nlevel; ++ilvl) {
      const std::size_t level = static_cast<std::size_t>(ilvl);
      hier.GetPatchLevel(level).nodes->copy(pi_new[level]);
    }
  }

  return boost::outcome_v2::success();
}

} // namespace fub::amrex

#endif