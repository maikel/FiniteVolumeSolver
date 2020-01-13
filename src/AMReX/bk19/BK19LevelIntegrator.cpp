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
      ForEachFab(execution::openmp, mf_cells, [&](const MFIter& mfi) {
        auto cells =
            SliceLast(MakePatchDataView(mf_cells[mfi]), cell_component);
        auto all_faces =
            SliceLast(MakePatchDataView(mf_faces[mfi]), face_component);
        // We are cautious and only compute the average for faces which exists
        // and are in range for the cells which exists on our grid
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
      ForEachFab(execution::openmp, mf_cells, [&](const MFIter& mfi) {
        auto cells =
            SliceLast(MakePatchDataView(mf_cells[mfi]), cell_component);
        auto all_faces =
            SliceLast(MakePatchDataView(mf_faces[mfi]), face_component);
        // We are cautious and only compute the average for faces which exists
        // and are in range for the cells which exists on our grid
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

void ComputePvFromScratch_(const IndexMapping<Equation>& index, MultiFab& dest,
                           const MultiFab& scratch) {
  // Shall be: Pv[i] = PTdensity * v[i]
  // Compute Pv_i for each velocity direction
  for (std::size_t i = 0; i < index.momentum.size(); ++i) {
    const int dest_component = static_cast<int>(i);
    MultiFab::Copy(dest, scratch, index.PTdensity, dest_component,
                   one_component, dest.nGrow());
    MultiFab::Multiply(dest, scratch, index.velocity[i], dest_component,
                       one_component, dest.nGrow());
  }
}

void RecomputeAdvectiveFluxes_(const IndexMapping<Equation>& index,
                               std::array<MultiFab, Rank>& Pv_faces,
                               MultiFab& Pv_cells, const MultiFab& scratch) {
  ComputePvFromScratch_(index, Pv_cells, scratch);
  // Average Pv_i for each velocity direction
  constexpr int face_component = 0;
  for (std::size_t dir = 0; dir < index.velocity.size(); ++dir) {
    const int cell_component = static_cast<int>(dir);
    AverageCellToFace_(Pv_faces[dir], face_component, Pv_cells, cell_component,
                       Direction(dir));
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

  ::amrex::DistributionMapping distribution_map = scratch.DistributionMap();
  ::amrex::BoxArray on_cells = scratch.boxArray();
  ::amrex::BoxArray on_nodes = on_cells;
  on_nodes.surroundingNodes();

  // compute RHS for elliptic solve
  MultiFab rhs(on_nodes, distribution_map, one_component, no_ghosts);

  // Copy UV into seperate MultiFab to use compDivergence
  // This assumes f = 0 and N = 0
  // Construct right hand side by: -dt div(UV)
  // Equation (28) in [BK19]
  //
  // Divergence needs on ghost cell width.
  MultiFab UV(on_cells, distribution_map, index.momentum.size(),
              one_ghost_cell_width);
  ComputePvFromScratch_(index, UV, scratch);

  UV.mult(-dt.count());
  lin_op.compDivergence({&rhs}, {&UV});

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

  MultiFab UV_correction(on_cells, distribution_map, index.momentum.size(),
                         no_ghosts);
  UV_correction.setVal(0.0);
  lin_op.getFluxes({&UV_correction}, {&pi});
  // Rupert sagt:
  // Nicht klar, wo die Funktionalität lin_op.getFluxes()  definiert ist.
  // Wird da eine Funktionalität aufgerufen, die schon in AMReX verdrahtet ist?
  // Falls ja, ist das denn auch diejenige, die wir brauchen?

  // Rupert sagt:
  // Achtung: Die Divergenz von  Pv   wird kontrolliert, aber es wird
  // \rho v   am Ende auf das neue Zeitniveau gehoben. Ist sichergestellt,
  // dass die hier verwendeten Routinen das alles richtig machen?

  // TODO check sign output of AMReX
  UV_correction.mult(-1.0 / dt.count(), 0);
  for (std::size_t i = 0; i < index.momentum.size(); ++i) {
    const int UV_component = static_cast<int>(i);
    MultiFab::Multiply(UV_correction, scratch, index.PTinverse, UV_component,
                       one_component, no_ghosts);
    // UV_correction is now a momentum correction. Thus add it.
    MultiFab::Add(scratch, UV_correction, UV_component, index.momentum[i],
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

  ::amrex::DistributionMapping distribution_map = scratch.DistributionMap();
  ::amrex::BoxArray on_cells = scratch.boxArray();

  // Construct sigma as in EulerBackward, but with -cp dt instead of -cp dt^2
  // Maikel: This needs a comment for me, why it is so
  // Rupert sagt:
  // Wenn ich das richtig sehe, liegt das daran, dass die folgenden zwei Zeilen:
  // // TODO check sign output of AMReX
  // momentum_correction.mult(-1.0 / dt.count(), 0);
  // in EulerForward nicht existieren, in EulerBackward aber schon.
  // Siehe auch Kommentar (**) weiter unten.
  MultiFab sigma(on_cells, distribution_map, one_component, no_ghosts);
  sigma.setVal(-equation.c_p * dt.count());
  MultiFab::Multiply(sigma, scratch, index.PTdensity, 0, one_component,
                     sigma.nGrow());
  //   MultiFab::Divide(sigma, scratch, index.PTinverse, 0, one_component,
  //                    sigma.nGrow());
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

  // Rupert sagt:  (**)
  // Hier müssten besagte zwei Zeilen stehen, wenn oben mit dt^2
  // multipliziert worden wäre.

  for (std::size_t i = 0; i < index.momentum.size(); ++i) {
    const int src_component = static_cast<int>(i);
    MultiFab::Add(scratch, momentum_correction, src_component,
                  index.momentum[i], one_component, no_ghosts);
  }

  RecoverVelocityFromMomentum_(scratch, index);

  // Rupert sagt:
  // Achtung, auch hier muss beachtet werden, dass  Pv  nicht dasselbe
  // ist wie  \rho v !   Ich bin mir eben nicht sicher, ob   lin_op.getFluxes()
  // das weiss.
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
  AdvectionSolver& advection_ = GetAdvection();
  BK19IntegratorContext& context = advection_.GetContext();
  MultiFab& scratch = context.GetScratch(level);

  WriteBK19Plotfile writeoutput{};
  writeoutput.plotfilename = "BK19_Pseudo_Incompressible-a/";

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
