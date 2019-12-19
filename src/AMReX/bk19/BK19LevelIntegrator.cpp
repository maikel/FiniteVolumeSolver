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
inline constexpr int Index_of_P = 7;
using Equation = CompressibleAdvection<Rank>;

// For now the implementation assumes Rank == 2
static_assert(Rank == 2);

namespace {
double EvaluateRhsU_(const Equation& equation, const Complete<Equation>& state,
                     Duration time_step_size) {
  const double f = equation.f;
  const double dt = time_step_size.count();
  const double U = state.velocity[0] * state.PTdensity;
  const double V = state.velocity[1] * state.PTdensity;
  const double dt_times_f = dt * f;
  const double dt_times_f_square = dt_times_f * dt_times_f;
  return (U + dt * f * V) / (1 + dt_times_f_square);
}

double EvaluateRhsV_(const Equation& equation, const Complete<Equation>& state,
                     Duration time_step_size) {
  const double f = equation.f;
  const double dt = time_step_size.count();
  const double U = state.velocity[0] * state.PTdensity;
  const double V = state.velocity[1] * state.PTdensity;
  const double dt_times_f = dt * f;
  const double dt_times_f_square = dt_times_f * dt_times_f;
  return (V - dt * f * U) / (1 + dt_times_f_square);
}

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

void ComputePvFromScratch_(const Equation& equation, MultiFab& dest,
                           const MultiFab& scratch) {
  ForEachFab(execution::openmp, dest, [&](const MFIter& mfi) {
    auto states = MakeView<const Complete<Equation>>(scratch[mfi], equation);
    auto Pv = MakePatchDataView(dest[mfi], 0);
    ForEachIndex(Pv.Box(), [&](int i, int j) {
      const double P = states.PTdensity(i, j);
      const double v = states.velocity(i, j);
      Pv(i, j) = P * v;
    });
  });
}

void RecomputeAdvectiveFluxes_(const Equation& equation,
                               std::array<MultiFab, Rank>& Pv_faces,
                               MultiFab& Pv_cells, const MultiFab& scratch) {
  ComputePvFromScratch_(equation, Pv_cells, scratch);
  for (int dir = 0; dir < Rank; ++dir) {
    AverageCellToFace_(Pv_faces[dir], Pv_cells, 0, 0, Direction(dir));
  }
}

Result<void, TimeStepTooLarge>
Advect_(BK19LevelIntegrator::AdvectionSolver& advection, int level, Duration dt,
        std::pair<int, int> subcycle) {
  return advection.AdvanceLevelNonRecursively(level, dt, subcycle);
}

::amrex::MultiFab DoEulerBackward_(const Equation& equation, ::amrex::MLMG& nodal_solver,
                      ::amrex::MLNodeHelmDualCstVel& lin_op,
                      BK19IntegratorContext& context, int level, Duration dt) {
  MultiFab& scratch = context.GetScratch(level);

  // Construct right hand side
  ::amrex::DistributionMapping dm = scratch.DistributionMap();
  ::amrex::BoxArray cells = scratch.boxArray();
  ::amrex::BoxArray nodes = cells;
  nodes.surroundingNodes();
  MultiFab rhs(nodes, dm, 1, 0);

  // Copy velocity into seperate MultiFab to use compDivergence
  // This assumes f = 0 and N = 0, otherwise you should use EvaluateRhsU_
  MultiFab momentum(cells, dm, AMREX_SPACEDIM, 1);
  MultiFab::Copy(momentum, scratch, 1, 0, AMREX_SPACEDIM, 1);
  momentum.mult(-dt.count());

  lin_op.compDivergence({&rhs}, {&momentum});

  // Construct sigma
  MultiFab sigma(cells, dm, 1, 0);
  sigma.setVal(-equation.c_p * dt.count() * dt.count());

  // TODO try to avoid magic constants by
  // requires some tepmlate magic
  // component = FabComponents(equation);
  // component.PTinerse
  constexpr int index_of_PTinverse = 6;
  MultiFab::Divide(sigma, scratch, index_of_PTinverse, 0, 1, 0);

  lin_op.setSigma(level, sigma);

  MultiFab pi(nodes, dm, 1, 0);
  nodal_solver.solve({&pi}, {&rhs}, 1e-10, 1e-10);

  MultiFab momentum_correction(cells, dm, AMREX_SPACEDIM, 0);
  momentum_correction.setVal(0.0);
  nodal_solver.getFluxes({&momentum_correction});

  // TODO check sign
  momentum_correction.mult(-1.0 / dt.count(), 0);

  MultiFab::Add(scratch, momentum_correction, 0, 1, AMREX_SPACEDIM, 0);
  MultiFab::Copy(scratch, scratch, 1, 4, AMREX_SPACEDIM, 0);
  MultiFab::Divide(scratch, scratch, 4, 0, 1, 0);
  MultiFab::Divide(scratch, scratch, 5, 0, 2, 0);

  context.FillGhostLayerSingleLevel(level);

  return pi;
}

void DoEulerForward_(const Equation& equation, ::amrex::MLNodeHelmDualCstVel& lin_op, BK19IntegratorContext& context,
                     int level, Duration dt) {
  MultiFab& scratch = context.GetScratch(level);
  ::amrex::BoxArray cells = scratch.boxArray();
  ::amrex::DistributionMapping dm = scratch.DistributionMap();

  // Construct sigma
  MultiFab sigma(cells, dm, 1, 0);
  sigma.setVal(-equation.c_p * dt.count());

  // TODO try to avoid magic constants by
  // requires some tepmlate magic
  // component = FabComponents(equation);
  // component.PTinerse
  constexpr int index_of_PTinverse = 6;
  MultiFab::Divide(sigma, scratch, index_of_PTinverse, 0, 1, 0);

  lin_op.setSigma(level, sigma);

  const MultiFab& pi_old = context.GetPi(level);
  MultiFab pi(pi_old.boxArray(), pi_old.DistributionMap(), 1, 1);
  pi.copy(pi_old, context.GetGeometry(level).periodicity());
  MultiFab momentum_correction(cells, dm, AMREX_SPACEDIM, 0);
  momentum_correction.setVal(0.0);
  lin_op.getFluxes({&momentum_correction}, {&pi});

  MultiFab::Add(scratch, momentum_correction, 0, 1, AMREX_SPACEDIM, 0);
  MultiFab::Copy(scratch, scratch, 1, 4, AMREX_SPACEDIM, 0);
  MultiFab::Divide(scratch, scratch, 4, 0, 1, 0);
  MultiFab::Divide(scratch, scratch, 5, 0, 2, 0);

  context.FillGhostLayerSingleLevel(level);
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
  MultiFab& scratch = advection_.GetScratch(level);
  MultiFab scratch_aux(scratch.boxArray(), scratch.DistributionMap(), scratch.nComp(), scratch.nGrow());
  scratch_aux.copy(scratch);
  BK19AdvectiveFluxes& Pv = advection_.GetContext().GetAdvectiveFluxes(level);

  [[maybe_unused]] const Duration half_dt = 0.5 * dt;

  // 1) Compute current Pv and interpolate to face centered quantity
  RecomputeAdvectiveFluxes_(equation_, Pv.on_faces, Pv.on_cells, scratch);

  // 2) Do the advection with the face-centered Pv
  // No accumulation
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

  // 5) Explicit Euler with old scratch data
  DoEulerForward_(equation_, *lin_op_,advection_.GetContext(), level, half_dt);

  // 6) Do the ssecond advection step with half-time Pv
  result = Advect_(advection_, level, dt, subcycle);
  if (!result) {
    return result;
  }

  // 6) Do the second euler backward integration step for the source term
  MultiFab pi_new = DoEulerBackward_(equation_, *nodal_solver_, *lin_op_, advection_.GetContext(),
                   level, half_dt);

  advection_.GetContext().GetPi(level).copy(pi_new);

  return boost::outcome_v2::success();
}

} // namespace fub::amrex