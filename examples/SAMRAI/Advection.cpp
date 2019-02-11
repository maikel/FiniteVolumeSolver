// Copyright (c) 2019 Maikel Nadolski
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

#include "fub/equations/Advection.hpp"
#include "fub/CartesianCoordinates.hpp"
#include "fub/Eigen.hpp"
#include "fub/SAMRAI/CartesianPatchHierarchy.hpp"
#include "fub/SAMRAI/DimensionalSplitBoundaryCondition.hpp"
#include "fub/SAMRAI/DimensionalSplitFluxMethod.hpp"
#include "fub/SAMRAI/GriddingAlgorithm.hpp"
#include "fub/SAMRAI/HyperbolicSplitLevelIntegrator.hpp"
#include "fub/SAMRAI/HyperbolicSplitPatchIntegrator.hpp"
#include "fub/SAMRAI/RegisterVariables.hpp"
#include "fub/SAMRAI/ScopeGuard.hpp"
#include "fub/tagging/GradientDetector.hpp"

#include <SAMRAI/hier/CoarseFineBoundary.h>
#include <SAMRAI/mesh/GriddingAlgorithm.h>
#include <SAMRAI/pdat/FaceDoubleConstantRefine.h>

#include <SAMRAI/appu/VisitDataWriter.h>
#include <SAMRAI/tbox/PIO.h>

#include <iostream>

namespace fub {
struct WaveData {
  static double wave_package(double x) {
    return std::exp(-20 * x * x) * std::sin(4.0 * M_PI * x);
  }

  void InitializeData(StateView<double, Advection<2>> states,
                      const CartesianCoordinates& coords) const {
    ForEachIndex(Mapping(states), [&](auto... is) {
      Advection<2>::State q{{.mass = wave_package(coords(is...).norm())}};
      states(is...) = q;
    });
  }
};

struct PatchIntegrator : fub::samrai::HyperbolicSplitPatchIntegrator {
  void
  ConservativeUpdateOnPatch(span<SAMRAI::pdat::CellData<double>*> next,
                            span<const SAMRAI::pdat::SideData<double>*> fluxes,
                            span<const SAMRAI::pdat::CellData<double>*> prev,
                            const SAMRAI::hier::Patch& patch, double dt,
                            int dir) override {}

  void
  ReconstructStatesFromCons(span<SAMRAI::pdat::CellData<double>*> states,
                            span<const SAMRAI::pdat::CellData<double>*> cons,
                            const SAMRAI::hier::Patch& patch,
                            int dir) override {}
};

struct FluxMethod : public fub::DimensionalSplitFluxMethod {
  double
  ComputeStableDtOnPatch(span<const SAMRAI::pdat::CellData<double>*> states,
                         const SAMRAI::hier::Patch& patch,
                         int dir) const override {
    return 0;
  }

  void ComputeFluxesOnPatch(span<SAMRAI::pdat::SideData<double>*> fluxes,
                            span<const SAMRAI::pdat::CellData<double>*> states,
                            const SAMRAI::hier::Patch& patch, double dt,
                            int dir) const override {}

  SAMRAI::hier::IntVector
  GetStencilWidth(const SAMRAI::tbox::Dimension& dim) const override {
    return SAMRAI::hier::IntVector::getOne(dim);
  }
};

struct BoundaryCondition_ : public fub::DimensionalSplitBoundaryCondition {
  void SetPhysicalBoundaryCondition(const SAMRAI::hier::Patch& patch,
                                    const SAMRAI::hier::Box& fill_box,
                                    double fill_time, Direction dir,
                                    int side) const override {}

  int GetStencilWidth() const override { return 1; }
};
} // namespace fub

int main(int argc, char** argv) {
  const fub::samrai::ScopeGuard guard(argc, argv);

  const SAMRAI::tbox::Dimension dim(2);
  const SAMRAI::hier::BlockId block_id(0);
  const SAMRAI::hier::Index lower(0, 0);
  const SAMRAI::hier::Index upper(127, 127);

  // Create a Hierarchy of specified extents
  const SAMRAI::hier::Box coarse_index_range(lower, upper, block_id);
  const fub::samrai::CoordinatesRange coordinates{{-1.0, -1.0}, {1.0, 1.0}};

  std::shared_ptr<SAMRAI::hier::PatchHierarchy> hierarchy =
      fub::samrai::CartesianPatchHierarchy(coarse_index_range, coordinates);
  hierarchy->setMaxNumberOfLevels(2);
  hierarchy->setRatioToCoarserLevel(SAMRAI::hier::IntVector(dim, 2), 1);

  const std::array<double, 2> velocity{1.0, 1.0};
  const fub::Advection advection{velocity};
  const fub::samrai::PatchDataIdSet data_ids =
      fub::samrai::RegisterVariables(advection, "Advection");

  using State = fub::Advection<2>::State;
  fub::GradientDetector tagging(std::pair{&State::mass, 1e-2});

  {
    fub::WaveData initial_data{};
    std::shared_ptr<SAMRAI::mesh::GriddingAlgorithm> initial_gridding =
        fub::samrai::GriddingAlgorithm(hierarchy, advection, data_ids.data_ids,
                                       initial_data, tagging);
    initial_gridding->makeCoarsestLevel(0.0);
    initial_gridding->makeFinerLevel(0, 0, 0, 0.0);
  }

  hierarchy->recursivePrint(SAMRAI::tbox::pout, {}, 2);

  fub::PatchIntegrator patch_integrator{};
  fub::FluxMethod flux_method{};
  fub::BoundaryCondition_ boundary_condition{};

  fub::samrai::HyperbolicSplitLevelIntegrator level_integrator(
      data_ids, hierarchy, patch_integrator, flux_method, boundary_condition);

  SAMRAI::hier::VariableDatabase::getDatabase()->printClassData(
      SAMRAI::tbox::pout);

  SAMRAI::appu::VisItDataWriter writer(dim, "VisItWriter", "output");
  writer.registerPlotQuantity("mass", "SCALAR", data_ids.data_ids[0]);
  writer.writePlotData(hierarchy, 0, 0);
}