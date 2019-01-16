// Copyright (c) 2018 Maikel Nadolski
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

#ifndef FUB_IDEAL_GAS_KINETIC_DRIVER_HPP
#define FUB_IDEAL_GAS_KINETIC_DRIVER_HPP

#include "fub/SAMRAI/DimensionalSplitBoundaryCondition.hpp"
#include "fub/SAMRAI/DimensionalSplitSystemSolver.hpp"
#include "fub/SAMRAI/DimensionalSplitSystemSourceSolver.hpp"
#include "fub/SAMRAI/GodunovSplitting.hpp"
#include "fub/SAMRAI/ideal_gas/DimensionalSplitTimeIntegrator.hpp"
#include "fub/SAMRAI/ideal_gas/IdealGasKinetics.hpp"
#include "fub/SAMRAI/ideal_gas/KineticSourceTerm.hpp"
#include "fub/core/mdspan.hpp"
#include "fub/core/span.hpp"

#include <SAMRAI/geom/CartesianGridGeometry.h>

#include <boost/optional.hpp>

namespace fub {
namespace ideal_gas {

class KineticDriverGrid {
public:
  using Variable = IdealGasEquation::Variable;

  KineticDriverGrid(DynamicMdSpan<double, 2> buffer, double dx)
      : mdspan_(buffer), dx_{dx} {}

  double& operator()(Variable var, int cell) {
    return mdspan_(static_cast<int>(var), cell);
  }
  const double& operator()(Variable var, int cell) const {
    return mdspan_(static_cast<int>(var), cell);
  }

  DynamicMdSpan<double, 2> Mdspan() noexcept { return mdspan_; }
  DynamicMdSpan<const double, 2> Mdspan() const noexcept { return mdspan_; }

  DynamicExtents<2> Extents() const noexcept { return mdspan_.extents(); }
  int Extent(int dim) const { return mdspan_.extent(dim); }

  double Dx() const noexcept { return dx_; }

private:
  DynamicMdSpan<double, 2> mdspan_;
  double dx_;
};

class KineticDriver {
public:
  KineticDriver(const std::shared_ptr<IdealGasKinetics>& equation,
                const CoordinateRange& coordinates, std::ptrdiff_t n_cells);

  KineticDriver(const std::shared_ptr<IdealGasKinetics>& equation,
                const CoordinateRange& coordinates, std::ptrdiff_t n_cells,
                const std::string& name);

  void InitializeHierarchy(const InitialCondition& init);

  void SetLeftBoundaryCondition(
      const std::shared_ptr<DimensionalSplitBoundaryCondition>& condition);

  void SetRightBoundaryCondition(
      const std::shared_ptr<DimensionalSplitBoundaryCondition>& condition);

  double ComputeStableDt() const;

  void Step(double dt);

  void Advance(double dt);
  void Advance(double dt,
               function_ref<void(const KineticDriver&, double)> feedback);

  double TimePoint() const noexcept { return time_point_; }

  const DimensionalSplitTimeIntegrator&
  GetDimensionalSplitTimeIntegrator() const noexcept {
    return time_integrator_;
  }

  MPI_Comm GetCommunicator() const noexcept {
    return GetPatchHierarchy()->getMPI().getCommunicator();
  }

  std::ptrdiff_t GetNCells() const;

  double GetDx() const noexcept {
    const double* dx = static_cast<const SAMRAI::geom::CartesianGridGeometry&>(
                           *GetPatchHierarchy()->getGridGeometry())
                           .getDx();
    return dx[0];
  }

  const std::shared_ptr<const IdealGasKinetics>& GetEquation() const noexcept {
    return source_term_.GetEquation();
  }

  const std::shared_ptr<SAMRAI::hier::PatchHierarchy>& GetPatchHierarchy() const
      noexcept {
    return hierarchy_;
  }

private:
  DimensionalSplitTimeIntegrator time_integrator_;
  KineticSourceTerm source_term_;
  std::shared_ptr<DimensionalSplitBoundaryConditions> boundary_condition_;
  std::shared_ptr<SAMRAI::hier::PatchHierarchy> hierarchy_;
  double time_point_{0.0};
  GodunovSplitting splitting_{};
  DimensionalSplitSystemSolver hyperbolic_{time_integrator_, splitting_};
  DimensionalSplitSystemSourceSolver solver_{hyperbolic_, source_term_,
                                             splitting_};
};

} // namespace ideal_gas
} // namespace fub

#endif