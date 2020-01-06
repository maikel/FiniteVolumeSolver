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

inline constexpr int Rank = AMREX_SPACEDIM;

// This works for now but we need something better.
// TODO Introduce some kind of map between variables and indices 
//   IndexMapping<Equation> index = MapIndices(equation);
// with
//   index.density, index.momentum[d], ...
inline constexpr int index_density = 0;
inline constexpr int index_momentum_x = 1;
inline constexpr int index_momentum_y = 2;
inline constexpr int index_momentum_z = 3;
inline constexpr int index_PTdensity = index_momentum_x + AMREX_SPACEDIM;
inline constexpr int index_velocity_x = index_PTdensity + 1;
inline constexpr int index_velocity_y = index_velocity_x + 1;
inline constexpr int index_velocity_z = index_velocity_x + 2;
inline constexpr int index_PTinverse = index_velocity_x + AMREX_SPACEDIM;
using Equation = CompressibleAdvection<Rank>;

// For now the implementation assumes Rank == 2
static_assert(Rank == 2);

namespace {
void AverageCellToFace_(MultiFab& mf_faces, const MultiFab& mf_cells,
                        int src_component, int dest_component, Direction dir) {
  if constexpr (AMREX_SPACEDIM == 2) {
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

void ComputePvFromScratch_(const Equation&, MultiFab& dest,
                           const MultiFab& scratch) {
  // Shall be: Pv = PTdensity * v
  MultiFab::Copy(dest, scratch, 0, index.PTdensity, 1, dest.nGrow());
  MultiFab::Copy(dest, scratch, 1, index.PTdensity, 1, dest.nGrow());
  MultiFab::Multiply(dest, scratch, 0, index.velocity[0], 1, dest.nGrow());
  MultiFab::Multiply(dest, scratch, 1, index.velocity[1], 1, dest.nGrow());
}

void RecomputeAdvectiveFluxes_(const Equation& equation,
                               std::array<MultiFab, Rank>& Pv_faces,
                               MultiFab& Pv_cells, const MultiFab& scratch) {
  ComputePvFromScratch_(equation, Pv_cells, scratch);
  for (int dir = 0; dir < Rank; ++dir) {
    AverageCellToFace_(Pv_faces[dir], Pv_cells, 0, dir, Direction(dir));
  }
}

Result<void, TimeStepTooLarge>
Advect_(BK19LevelIntegrator::AdvectionSolver& advection, int level, Duration dt,
        std::pair<int, int> subcycle) {
  return advection.AdvanceLevelNonRecursively(level, dt, subcycle);
}

::amrex::MultiFab DoEulerBackward_(const Equation& equation,
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
  MultiFab rhs(on_nodes, distribution_map, 1, 0);

  // Copy momentum into seperate MultiFab to use compDivergence
  // This assumes f = 0 and N = 0
  // Construct right hand side by: -dt div(momentum)
  // Equation (28) in [BK19]
  MultiFab momentum(on_cells, distribution_map, AMREX_SPACEDIM, 1);
  MultiFab::Copy(momentum, scratch, index_momentum_x, 0, AMREX_SPACEDIM, 1);
  momentum.mult(-dt.count());
  lin_op.compDivergence({&rhs}, {&momentum});

  // Construct sigma by: -cp dt^2 (P Theta)^o (Equation (27) in [BK19]) 
  MultiFab sigma(on_cells, distribution_map, 1, 0);
  sigma.setVal(-equation.c_p * dt.count() * dt.count());
  // MultiFab::Divide(dest, src, src_comp, dest_comp, n_comp, n_grow);
  MultiFab::Divide(sigma, scratch, index_PTinverse, 0, 1, 0);
  // Maikel: Why is this -cp dt^2 / PTinverse and not -cp dt^2 PTdensity / PTinverse
  // MultiFab::Multiply(sigma, scratch, index_PTdensity, 0, 1, 0);

  lin_op.setSigma(level, sigma);

  MultiFab pi(on_nodes, distribution_map, 1, 0);
  nodal_solver.solve({&pi}, {&rhs}, 1e-10, 1e-10);

  MultiFab momentum_correction(on_cells, distribution_map, AMREX_SPACEDIM, 0);
  momentum_correction.setVal(0.0);
  lin_op.getFluxes({&momentum_correction}, {&pi});

  // TODO check sign output of AMReX
  momentum_correction.mult(-1.0 / dt.count(), 0);

  // Fix Momentum
  MultiFab::Add(scratch, momentum_correction, 0, index_momentum_x, AMREX_SPACEDIM, 0);

  // Recover Velocity from Momentum
  // MultiFab::Copy(dest, src, src_comp, dest_comp, n_comp, n_grow);
  // MultiFab::Divide(dest, src, src_comp, dest_comp, n_comp, n_grow);
  MultiFab::Copy(scratch, scratch, index_momentum_x, index_velocity_x, AMREX_SPACEDIM, 0);
  MultiFab::Divide(scratch, scratch, index_density, index_velocity_x, 1, 0);
  MultiFab::Divide(scratch, scratch, index_density, index_velocity_y, 1, 0);

  return pi;
}

void DoEulerForward_(const Equation& equation,
                     ::amrex::MLNodeHelmDualCstVel& lin_op,
                     BK19IntegratorContext& context, int level, Duration dt) {
  MultiFab& scratch = context.GetScratch(level);
  ::amrex::BoxArray cells = scratch.boxArray();
  ::amrex::DistributionMapping dm = scratch.DistributionMap();

  // Construct sigma as in EulerBackward, but with -cp dt instead of -cp dt^2
  // Maikel: This needs a comment for me, why it is so
  MultiFab sigma(cells, dm, 1, 0);
  sigma.setVal(-equation.c_p * dt.count());
  // sigma = - c_p dt / PTinverse   (Maikel: ????)
  MultiFab::Divide(sigma, scratch, index_PTinverse, 0, 1, 0);

  lin_op.setSigma(level, sigma);

  // To compute the fluxes from the old pi we need one ghost cell width
  // Thus, we use periodic boundaries for now and redistribute the pi_n onto the boundary
  // TODO: What happens to pi otherwise?
  const MultiFab& pi_old = context.GetPi(level);
  MultiFab pi(pi_old.boxArray(), pi_old.DistributionMap(), 1, 1);
  // This assumes periodicity in each direction
  pi.ParallelCopy(pi_old, context.GetGeometry(level).periodicity());
  MultiFab momentum_correction(cells, dm, AMREX_SPACEDIM, 0);
  momentum_correction.setVal(0.0);
  lin_op.getFluxes({&momentum_correction}, {&pi});

  // Fix Momentum
  MultiFab::Add(scratch, momentum_correction, 0, index_momentum_x, AMREX_SPACEDIM, 0);

  // Recover Velocity from Momentum
  // MultiFab::Copy(dest, src, src_comp, dest_comp, n_comp, n_grow);
  MultiFab::Copy(scratch, scratch, index_momentum_x, index_velocity_x, AMREX_SPACEDIM, 0);
  MultiFab::Divide(scratch, scratch, index_density, index_velocity_x, 1, 0);
  MultiFab::Divide(scratch, scratch, index_density, index_velocity_y, 1, 0);
}

} // namespace

BK19LevelIntegrator::BK19LevelIntegrator(
    const CompressibleAdvection<Rank>& equation, AdvectionSolver advection,
    std::shared_ptr<::amrex::MLNodeHelmDualCstVel> linop)
    : AdvectionSolver(std::move(advection)), equation_(equation),
      lin_op_(std::move(linop)),
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
  RecomputeAdvectiveFluxes_(equation_, Pv.on_faces, Pv.on_cells, scratch);

  // 2) Do the advection with the face-centered Pv
  //    Open Question: Coarse Fine Boundary? 
  //      - Need an option to do nothing there
  Result<void, TimeStepTooLarge> result =
      Advect_(advection_, level, half_dt, subcycle);
  if (!result) {
    return result;
  }

  // 3) Do the first euler backward integration step for the source term
  DoEulerBackward_(equation_, *nodal_solver_, *lin_op_, advection_.GetContext(),
                   level, half_dt);

  // 4) Recompute Pv at half time
  RecomputeAdvectiveFluxes_(equation_, Pv.on_faces, Pv.on_cells, scratch);
  std::swap(scratch_aux, scratch);
  context.FillGhostLayerSingleLevel(level);

  // 5) Explicit Euler with old scratch data
  //   - We need a current pi_n here. What is the initial one?
  DoEulerForward_(equation_, *lin_op_, advection_.GetContext(), level, half_dt);

  // 6) Do the second advection step with half-time Pv and full time step
  //   - Currently, scratch contains the result of euler forward step, 
  //     which started at the old time level.
  context.FillGhostLayerSingleLevel(level);
  result = Advect_(advection_, level, dt, subcycle);
  if (!result) {
    return result;
  }

  // 6) Do the second euler backward integration step for the source term
  MultiFab pi_new = DoEulerBackward_(equation_, *nodal_solver_, *lin_op_,
                                     advection_.GetContext(), level, half_dt);

  // Copy pi_n+1 to pi_n
  context.GetPi(level).copy(pi_new);

  return boost::outcome_v2::success();
}

} // namespace fub::amrex