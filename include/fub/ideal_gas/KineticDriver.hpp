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

#include "fub/ideal_gas/HyperbolicTimeIntegrator.hpp"
#include "fub/ideal_gas/IdealGasKinetics.hpp"
#include "fub/ideal_gas/KineticSourceTerm.hpp"
#include "fub/solver/DimensionalSplitSystemSolver.hpp"
#include "fub/solver/GodunovSplitting.hpp"
#include "fub/solver/HyperbolicSystemSourceSolver.hpp"
#include "fub/solver/SplitBoundaryCondition.hpp"

namespace fub {
namespace ideal_gas {

class KineticDriver {
public:
  KineticDriver(const std::shared_ptr<IdealGasKinetics>& equation,
                const CoordinateRange& coordinates, std::ptrdiff_t n_cells);

  void InitializeHierarchy(const InitialCondition& init);

  void SetLeftBoundaryCondition(
      const std::shared_ptr<const SplitBoundaryCondition>& condition);

  void SetRightBoundaryCondition(
      const std::shared_ptr<const SplitBoundaryCondition>& condition);

  double ComputeStableDt() const;

  void Advance(double dt);

  std::vector<SAMRAI::hier::Patch> GatherAll() const;

private:
  HyperbolicTimeIntegrator time_integrator_;
  KineticSourceTerm source_term_;
  std::shared_ptr<SplitBoundaryConditions> boundary_condition_;
  std::shared_ptr<SAMRAI::hier::PatchHierarchy> hierarchy_;
  double time_point_{0.0};
  GodunovSplitting splitting_{};
  DimensionalSplitSystemSolver hyperbolic_{time_integrator_, splitting_};
  HyperbolicSystemSourceSolver solver_{hyperbolic_, source_term_, splitting_};
};

} // namespace ideal_gas
} // namespace fub

#endif