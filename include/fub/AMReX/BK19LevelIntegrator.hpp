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

#include "fub/amrex/IntegratorContext.hpp"
#include "fub/solver/DimensionalSplitLevelIntegrator.hpp"
#include "fub/amrex/MLMG/MLNodeHelmDualCstVel.hpp"

namespace fub::amrex {

class BK19LevelIntegrator {
public:
  static constexpr int Rank = AMREX_SPACEDIM;

  using Coordinates = Eigen::Vector<double, Rank>;
  using Equation = CompressibleAdvection<Rank>;
  using Complete = ::fub::Complete<Equation>;
  using Conservative = ::fub::Conservative<Equation>;

  using NodalEllipticSolver = ::amrex::MLNodeLinOp;
  using AdvectionSolver = DimensionalSplitLevelIntegrator<Rank, IntegratorContext>;

  BK19LevelIntegrator(const CompressibleAdvection<Rank>& equation,
      DimensionalSplitLevelIntegrator advection,
      std::shared_ptr<NodalEllipticSolver> nodal_elliptic_solver);

  void ResetPatchHierarchy(std::shared_ptr<GriddingAlgorithm> grid);

  Duration ComputeStableDt(int level);

  double EvaluateRhs(const Complete& state, const Coordinates& x);

  static void AverageCellToFace(::amrex::MultiFab& faces,
                                const ::amrex::MultiFab& cells, int src_comp,
                                int dest_comp, int n_comp, int dir);

  Result<void, TimeStepTooLargeError> AdvanceLevel(int level, Duration dt);

private:
  CompressibleAdvection<Rank> equation_;
  std::shared_ptr<GriddingAlgorithm> grid_;
  AdvectionSolver advection_;
  std::shared_ptr<::amrex::MLNodeLinOp> nodal_solver_;
};

} // namespace fub::amrex