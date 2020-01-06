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
#include "fub/AMReX/bk19/BK19IntegratorContext.hpp"

namespace fub::amrex {

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
void AverageCellToFace_(MultiFab& mf_faces, const MultiFab& mf_cells,
                        int src_component, int dest_component, Direction dir) {
  if constexpr (Rank == 2) {
    FUB_ASSERT(dir == Direction::X || dir == Direction::Y);
    if (dir == Direction::X) {
      ForEachFab(execution::openmp, mf_cells, [&](const MFIter& mfi) {
        auto cells = SliceLast(MakePatchDataView(mf_cells[mfi]), src_component);
        auto faces =
            SliceLast(MakePatchDataView(mf_faces[mfi]), dest_component);
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
      ForEachFab(execution::openmp, mf_cells, [&](const MFIter& mfi) {
        auto cells = SliceLast(MakePatchDataView(mf_cells[mfi]), src_component);
        auto faces =
            SliceLast(MakePatchDataView(mf_faces[mfi]), dest_component);
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

void ComputePvFromScratch_(const IndexMapping<Equation>& index, MultiFab& dest,
                           const MultiFab& scratch) {
  // Shall be: Pv[i] = PTdensity * v[i]
  for (std::size_t i = 0; i < index.momentum.size(); ++i) {
    MultiFab::Copy(dest, scratch, index.PTdensity, i, one_component, dest.nGrow());
    MultiFab::Multiply(dest, scratch, index.velocity[i], i, one_component, dest.nGrow());
  }
}

void RecomputeAdvectiveFluxes_(const IndexMapping<Equation>& index,
                               std::array<MultiFab, Rank>& Pv_faces,
                               MultiFab& Pv_cells, const MultiFab& scratch) {
  ComputePvFromScratch_(index, Pv_cells, scratch);
  for (std::size_t dir = 0; dir < index.velocity.size(); ++dir) {
    AverageCellToFace_(Pv_faces[dir], Pv_cells, 0, dir, Direction(dir));
  }
}

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

::amrex::MultiFab DoEulerBackward_(const Equation& equation,
                                   const IndexMapping<Equation>& index,
                                   ::amrex::MLMG& nodal_solver,
                                   ::amrex::MLNodeHelmDualCstVel& lin_op,
                                   BK19IntegratorContext& context, int level,
                                   Duration dt) {
  MultiFab& scratch = context.GetScratch(level);

  // Construct right hand side
  ::amrex::DistributionMapping distribution_map = scratch.DistributionMap();
  ::amrex::BoxArray on_cells = scratch.boxArray();
  ::amrex::BoxArray on_nodes = on_cells;
  on_nodes.surroundingNodes();
  MultiFab rhs(on_nodes, distribution_map, one_component, no_ghosts);

  // Copy momentum into seperate MultiFab to use compDivergence
  // This assumes f = 0 and N = 0
  // Construct right hand side by: -dt div(momentum)
  // Equation (28) in [BK19]
  //
  // Divergence needs on ghost cell width.
  MultiFab momentum(on_cells, distribution_map, index.momentum.size(),
                    one_ghost_cell_width);
  for (std::size_t i = 0; index.momentum.size(); ++i) {
    MultiFab::Copy(momentum, scratch, index.momentum[i], i, one_component,
                   one_ghost_cell_width);
  }
  momentum.mult(-dt.count());
  lin_op.compDivergence({&rhs}, {&momentum});

  // Construct sigma by: -cp dt^2 (P Theta)^o (Equation (27) in [BK19])
  // MultiFab::Divide(dest, src, src_comp, dest_comp, n_comp, n_grow);
  MultiFab sigma(on_cells, distribution_map, one_component, no_ghosts);
  sigma.setVal(-equation.c_p * dt.count() * dt.count());
  MultiFab::Multiply(sigma, scratch, index.PTdensity, 0, one_component,
                     sigma.nGrow());
  MultiFab::Divide(sigma, scratch, index.PTinverse, 0, one_component,
                   sigma.nGrow());
  lin_op.setSigma(level, sigma);
  MultiFab pi(on_nodes, distribution_map, one_component, no_ghosts);
  nodal_solver.solve({&pi}, {&rhs}, 1e-10, 1e-10);

  MultiFab momentum_correction(on_cells, distribution_map,
                               index.momentum.size(), no_ghosts);
  momentum_correction.setVal(0.0);
  lin_op.getFluxes({&momentum_correction}, {&pi});

  // TODO check sign output of AMReX
  momentum_correction.mult(-1.0 / dt.count(), 0);

  // Fix Momentum
  for (std::size_t i = 0; index.momentum.size(); ++i) {
    MultiFab::Add(scratch, momentum_correction, i, index.momentum[i],
                  one_component, no_ghosts);
  }

  RecoverVelocityFromMomentum_(scratch, index);

  return pi;
}

void DoEulerForward_(const Equation& equation,
                     const IndexMapping<Equation>& index,
                     ::amrex::MLNodeHelmDualCstVel& lin_op,
                     BK19IntegratorContext& context, int level, Duration dt) {
  MultiFab& scratch = context.GetScratch(level);
  ::amrex::BoxArray on_cells = scratch.boxArray();
  ::amrex::DistributionMapping distribution_map = scratch.DistributionMap();

  // Construct sigma as in EulerBackward, but with -cp dt instead of -cp dt^2
  // Maikel: This needs a comment for me, why it is so
  MultiFab sigma(on_cells, distribution_map, one_component, no_ghosts);
  sigma.setVal(-equation.c_p * dt.count());
  MultiFab::Multiply(sigma, scratch, index.PTdensity, 0, one_component,
                     sigma.nGrow());
  MultiFab::Divide(sigma, scratch, index.PTinverse, 0, one_component,
                   sigma.nGrow());
  lin_op.setSigma(level, sigma);

  // To compute the fluxes from the old pi we need one ghost cell width
  // Thus, we use periodic boundaries for now and redistribute the pi_n onto the
  // boundary
  // TODO: What happens to pi otherwise?
  const MultiFab& pi_old = context.GetPi(level);
  MultiFab pi(pi_old.boxArray(), pi_old.DistributionMap(), one_component,
              one_ghost_cell_width);
  // This assumes periodicity in each direction
  pi.ParallelCopy(pi_old, context.GetGeometry(level).periodicity());

  MultiFab momentum_correction(on_cells, distribution_map,
                               index.momentum.size(), no_ghosts);
  momentum_correction.setVal(0.0);
  lin_op.getFluxes({&momentum_correction}, {&pi});
  for (std::size_t i = 0; index.momentum.size(); ++i) {
    MultiFab::Add(scratch, momentum_correction, i, index.momentum[i],
                  one_component, no_ghosts);
  }

  RecoverVelocityFromMomentum_(scratch, index);
}

} // namespace

BK19LevelIntegrator::BK19LevelIntegrator(
    const CompressibleAdvection<Rank>& equation, AdvectionSolver advection,
    std::shared_ptr<::amrex::MLNodeHelmDualCstVel> linop)
    : AdvectionSolver(std::move(advection)), equation_(equation),
      index_(equation_), lin_op_(std::move(linop)),
      nodal_solver_(std::make_shared<::amrex::MLMG>(*lin_op_)) {}

Result<void, TimeStepTooLarge>
BK19LevelIntegrator::AdvanceLevelNonRecursively(int level, Duration dt,
                                                std::pair<int, int> subcycle) {
  AdvectionSolver advection_ = GetAdvection();
  BK19IntegratorContext& context = advection_.GetContext();
  MultiFab& scratch = context.GetScratch(level);

  // Save data on current time level for later use
  MultiFab scratch_aux(scratch.boxArray(), scratch.DistributionMap(),
                       scratch.nComp(), scratch.nGrow());
  scratch_aux.copy(scratch);
  BK19AdvectiveFluxes& Pv = context.GetAdvectiveFluxes(level);

  const Duration half_dt = 0.5 * dt;

  // 1) Compute current Pv and interpolate to face centered quantity
  //    Current Pv is given by: Pv = PTdensity * velocity
  RecomputeAdvectiveFluxes_(index_, Pv.on_faces, Pv.on_cells, scratch);

  // 2) Do the advection with the face-centered Pv
  //    Open Question: Coarse Fine Boundary?
  //      - Need an option to do nothing there
  Result<void, TimeStepTooLarge> result =
      Advect_(advection_, level, half_dt, subcycle);
  if (!result) {
    return result;
  }

  // 3) Do the first euler backward integration step for the source term
  DoEulerBackward_(equation_, index_, *nodal_solver_, *lin_op_,
                   advection_.GetContext(), level, half_dt);

  // 4) Recompute Pv at half time
  RecomputeAdvectiveFluxes_(index_, Pv.on_faces, Pv.on_cells, scratch);
  std::swap(scratch_aux, scratch);
  context.FillGhostLayerSingleLevel(level);

  // 5) Explicit Euler with old scratch data
  //   - We need a current pi_n here. What is the initial one?
  DoEulerForward_(equation_, index_, *lin_op_, advection_.GetContext(), level,
                  half_dt);

  // 6) Do the second advection step with half-time Pv and full time step
  //   - Currently, scratch contains the result of euler forward step,
  //     which started at the old time level.
  context.FillGhostLayerSingleLevel(level);
  result = Advect_(advection_, level, dt, subcycle);
  if (!result) {
    return result;
  }

  // 6) Do the second euler backward integration step for the source term
  MultiFab pi_new =
      DoEulerBackward_(equation_, index_, *nodal_solver_, *lin_op_,
                       advection_.GetContext(), level, half_dt);

  // Copy pi_n+1 to pi_n
  context.GetPi(level).copy(pi_new);

  return boost::outcome_v2::success();
}

} // namespace fub::amrex