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

#include "fub/AMReX/BK19LevelIntegrator.hpp"

namespace fub::amrex {

using ::amrex::MFIter;
using ::amrex::MultiFab;

inline constexpr int Rank = AMREX_SPACEDIM;
inline constexpr int Index_of_P = 7;
using Equation = CompressibleAdvection<Rank>;

// For now the implementation assumes Rank == 2
static_assert(Rank == 2);

namespace {
double EvaluateRhsU_(const Equation& equation, const Complete& state,
                     Duration time_step_size) {
  const double f = equation.f;
  const double dt = time_step_size.count();
  const double U = state.velocity[0];
  const double V = state.velocity[1];
  const double dt_times_f = dt * f;
  const double dt_times_f_square = dt_times_f * dt_times_f;
  return (U + dt * f * V) / (1 + dt_times_f_square);
}

double EvaluateRhsV_(const Equation& equation, const Complete& state,
                     Duration time_step_size) {
  const double f = equation.f(x);
  const double dt = time_step_size.count();
  const double U = state.velocity[0];
  const double V = state.velocity[1];
  const double dt_times_f = dt * f;
  const double dt_times_f_square = dt_times_f * dt_times_f;
  return (V - dt * f * U) / (1 + dt_times_f_square);
}

double EvaluateRhsW_(const Equation& equation, const Complete& state,
                     Duration time_step_size) {
  const double dt = time_step_size.count();
  const double U = state.velocity[0];
  const double V = state.velocity[1];
  const double dt_times_f = dt * f;
  const double dt_times_f_square = dt_times_f * dt_times_f;
  return (V - dt * f * U) / (1 + dt_times_f_square);
}

void AverageCellToFace_(MultiFab& mf_faces, const MultiFab& mf_cells, int src_component, int dest_component,
                        Direction dir) {
  if constexpr (AMREX_SPACEDIM == 2) {
    FUB_ASSERT(dir == Direction::X || dir == Direction::Y);
    if (dir == Direction::X) {
      ForEachFab(execution::openmp, mf_cells, [&](const MFIter& mfi) {
        auto cells = SliceLast(MakePatchDataView(mf_cells[mfi]), src_component);
        auto faces =
            SliceLast(MakePatchDataView(mf_faces[mfi]), dest_component);
        ForEachIndex(faces.Box(), [](int i, int j) {
          // clang-format off
          faces(i, j) = 1.0 * cells(i - 1, j - 1) + 1.0 * cells(i, j - 1) +
                        2.0 * cells(i - 1,     j) + 2.0 * cells(i,     j) +
                        1.0 * cells(i - 1, j + 1) + 1.0 * cells(i, j + 1);
          faces(i, j) *= 0.125;
          // clang-format on
        });
      });
    } else {
      ForEachFab(exection::openmp, mf_cells, [&](const MFIter& mfi) {
        auto cells = SliceLast(MakePatchDataView(mf_cells[mfi]), src_component);
        auto faces =
            SliceLast(MakePatchDataView(mf_faces[mfi]), dest_component);
        ForEachIndex(faces.Box(), [](int i, int j) {
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

void ComputePvFromScratch_(const Equation& equation, MultiFab& Pv, const MultiFab& scratch) {
  ForEachFab(execution::openmp, Pv, [&](const MFIter& mfi) {
    auto states = MakeView<Complete>(equation, scratch[mfi]);
    auto dest = MakePatchDataView(Pv[mfi], 0);
    ForEachIndex(dest.Box(), [&](int i, int j) {
      const double P = states.PTdensity(i, j);
      const double v = states.velocity(i, j);
      Pv(i, j) = P * v;
    });
  });
}

void RecomputeAdvectiveFluxes_(std::array<MultiFab, Rank>& Pv_faces, MultiFab& Pv_cells, const MultiFab& scratch) {
  ComputePvFromScratch_(Pv_cells, scratch);
  for (int dir = 0; dir < Rank; ++dir) {
    AverageCellToFace_(Pv_faces[dir], Pv_cells, Direction(dir));
  }
}
} // namespace

Result<void, TimeStepTooLargeError>
BK19LevelIntegrator::AdvanceLevel(int level, Duration dt) {
  MultiFab& scratch = advection_.GetScratch(level);
  AdvectiveFluxes& Pv = advection_.GetPv(level);

  const Duration half_dt = 0.5 * dt;

  // 1) Compute current Pv and interpolate to face centered quantity
  RecomputeAdvectiveFluxes_(Pv.on_faces, Pv.on_cells, scratch);

  // 2) Do the advection with these Pv on its faces
  // No accumulation
  Advect_(advection_, half_dt);

  // 3) Do the first euler backward integration step for the source term
  DoEulerBackward_(nodal_solver_, half_dt);

  // 4) Recompute Pv at half time
  RecomputeAdvectiveFluxes_(Pv_faces, Pv_cells, scratch);

  // 5) Do the euler forward step for the source term
  DoEulerForward_(scratch, half_dt);

  // 6) Do the ssecond advection step with half-time Pv
  Advect_(advection_, dt);

  // 6) Do the second euler backward integration step for the source term
  DoEulerBackward_(nodal_solver_, half_dt);
}

} // namespace fub::amrex