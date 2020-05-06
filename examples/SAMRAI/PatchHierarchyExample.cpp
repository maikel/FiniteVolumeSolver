// Copyright (c) 2019 Maikel Nadolski
// Copyright (c) 2019 Patrick Denzler
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

#include "fub/SAMRAI.hpp"
#include "fub/Solver.hpp"

#include <SAMRAI/appu/VisItDataWriter.h>
#include <SAMRAI/geom/CartesianGridGeometry.h>
#include <SAMRAI/hier/Patch.h>
#include <SAMRAI/pdat/SideDataFactory.h>

#include "fub/ext/uuid.hpp"

#include "fub/SAMRAI/GriddingAlgorithm.hpp"
#include "fub/SAMRAI/tagging_method/ConstantBox.hpp"
#include "fub/SAMRAI/tagging_method/GradientDetector.hpp"
#include <boost/container/static_vector.hpp>

namespace fub {
void ComputeCoordinates(span<double> x, const SAMRAI::pdat::CellIndex& cell,
                        const SAMRAI::geom::CartesianGridGeometry& geom) {
  const int dim = geom.getDim().getValue();
  span<const double> xlower(geom.getXLower(), dim);
  span<const double> dx(geom.getDx(), dim);
  for (int i = 0; i < dim; ++i) {
    x[i] = xlower[i] + cell[unsigned(i)] * dx[i] + 0.5 * dx[i];
  }
}

std::shared_ptr<const SAMRAI::geom::CartesianGridGeometry>
GetCartesianGridGeometry(const SAMRAI::hier::PatchLevel& level) {
  const SAMRAI::geom::CartesianGridGeometry& base_geom =
      static_cast<const SAMRAI::geom::CartesianGridGeometry&>(
          *level.getGridGeometry());
  return std::static_pointer_cast<SAMRAI::geom::CartesianGridGeometry>(
      base_geom.makeRefinedGridGeometry(MakeUniqueName(),
                                        level.getRatioToLevelZero()));
}

} // namespace fub

struct CircleData {
  using Complete = fub::Complete<fub::Advection2d>;
  fub::samrai::DataDescription data_description_;
  fub::Advection2d equation_;

  //  fub::View<Complete> MakeView(SAMRAI::hier::Patch& patch, const
  //  fub::samrai::DataDescription& desc) {
  //    fub::View<Complete> view;
  //    view.mass =
  //    fub::samrai::MakePatchDataView(patch.getPatchData(desc.data_ids[0]));m
  //  }

  void InitializeData(std::shared_ptr<SAMRAI::hier::PatchLevel> level, const fub::samrai::GriddingAlgorithm&, int,
                      fub::Duration) const {
    std::shared_ptr<const SAMRAI::geom::CartesianGridGeometry> geom =
        fub::GetCartesianGridGeometry(*level);
    for (const std::shared_ptr<SAMRAI::hier::Patch>& patch : *level) {
      SAMRAI::pdat::CellData<double>& data =
          static_cast<SAMRAI::pdat::CellData<double>&>(*patch->getPatchData(0));
      int lower_x = data.getArrayData().getBox().lower(0);
      int lower_y = data.getArrayData().getBox().lower(1);
      int upper_x = data.getArrayData().getBox().upper(0);
      int upper_y = data.getArrayData().getBox().upper(1);
      for (int j = lower_y; j <= upper_y; ++j) {
        for (int i = lower_x; i <= upper_x; ++i) {
          SAMRAI::pdat::CellIndex cell(SAMRAI::hier::Index(i, j));
          double x[2];
          fub::ComputeCoordinates(x, cell, *geom);
          constexpr double r2 = 0.25 * 0.25;
          if (x[0] * x[0] + x[1] * x[1] < r2) {
            data(cell) = 3.0;
          } else {
            data(cell) = 1.0;
          }
        }
      }
    }
  }
};

int main(int argc, char** argv) {
  fub::samrai::ScopeGuard guard(argc, argv);

  // fub::IdealGasMix<1> equation(fub::Burke2012{});
  fub::Advection2d equation{{}};
  using Eq = std::decay_t<decltype(equation)>;

  SAMRAI::tbox::Dimension dim(Eq::Rank());
  fub::samrai::PatchHierarchyOptions hier_opts{SAMRAI::hier::Index(dim, 2), 4};

  const std::array<std::ptrdiff_t, 2> n_cells{128, 128};
  const fub::samrai::CoordinateRange<2> x_range{
      std::array<double, 2>{-1.0, -1.0}, std::array<double, 2>{+1.0, +1.0}};

  std::shared_ptr<SAMRAI::geom::CartesianGridGeometry> geom =
      MakeCartesianGridGeometry(n_cells, x_range);

  fub::samrai::PatchHierarchy ph(equation, geom, hier_opts);
  fub::samrai::DataDescription data_desc = ph.GetDataDescription();

  SAMRAI::hier::Box tagbox(SAMRAI::hier::Index(40, 40),
                           SAMRAI::hier::Index(80, 80),
                           SAMRAI::hier::BlockId(0));
  fub::samrai::ConstantBox constbox{tagbox};

  using State = fub::Advection2d::Complete;
  fub::samrai::GradientDetector gradient{equation,
                                         std::pair{&State::mass, 1e-3}};

  std::vector<int> tb(hier_opts.max_number_of_levels - 1, 4);
  fub::samrai::GriddingAlgorithm ga(std::move(ph),
                                    CircleData{data_desc, equation},
                                    gradient, tb);

  ga.InitializeHierarchy();

  SAMRAI::hier::VariableDatabase* vardb =
      SAMRAI::hier::VariableDatabase::getDatabase();
  vardb->printClassData(std::cout, false);

  ga.GetPatchHierarchy().GetNative()->recursivePrint(std::cout, "", 2);
  // std::cout << std::endl << std::endl <<
  // "---------------------------------------------------------------------" <<
  // std::endl << std::endl; ph1.GetNative()->recursivePrint(std::cout, "", 2);

  fub::samrai::GriddingAlgorithm ga2(ga);
  ga2.GetPatchHierarchy().GetNative()->removePatchLevel(
      ga2.GetPatchHierarchy().GetNative()->getNumberOfLevels() - 1);

  ph = ga.GetPatchHierarchy();
  ph.GetNative()->removePatchLevel(ph.GetNative()->getNumberOfLevels() - 1);

  SAMRAI::appu::VisItDataWriter writer(dim, "VisItWriter", "SAMRAI/PHE");
  writer.registerPlotQuantity("mass", "SCALAR",
                              ph.GetDataDescription().data_ids[0]);

  writer.writePlotData(ga.GetPatchHierarchy().GetNative(), 0);
  writer.writePlotData(ph.GetNative(), 1);
  writer.writePlotData(ga2.GetPatchHierarchy().GetNative(), 2);

  const std::vector<int>& data_ids =
      ga.GetPatchHierarchy().GetDataDescription().data_ids;
  std::vector<SAMRAI::pdat::CellData<double>*> datas(data_ids.size());
  // Generate View
  for (const std::shared_ptr<SAMRAI::hier::Patch>& patch :
       *ga.GetPatchHierarchy().GetNative()->getPatchLevel(0)) {
    std::transform(data_ids.begin(), data_ids.end(), datas.begin(),
                   [&](int id) -> SAMRAI::pdat::CellData<double>* {
                     return static_cast<SAMRAI::pdat::CellData<double>*>(
                         patch->getPatchData(id).get());
                   });
  }

  const SAMRAI::hier::IntVector ghosts(dim, 2);
  std::shared_ptr<SAMRAI::hier::PatchDescriptor> scratch_descriptor =
      std::make_shared<SAMRAI::hier::PatchDescriptor>();
  for (int id : data_desc.data_ids) {
    std::shared_ptr<SAMRAI::hier::Variable> variable{};
    if (vardb->mapIndexToVariable(id, variable)) {
      const std::string& name = variable->getName();
      const int depth =
          static_cast<SAMRAI::pdat::CellVariable<double>*>(variable.get())
              ->getDepth();
      scratch_descriptor->definePatchDataComponent(
          name, std::make_shared<SAMRAI::pdat::CellDataFactory<double>>(
                    depth, ghosts));
    }
  }
  scratch_descriptor->printClassData(SAMRAI::tbox::pout);

  const SAMRAI::hier::IntVector face_ghosts(dim, 1);
  std::shared_ptr<SAMRAI::hier::PatchDescriptor> flux_descriptor =
      std::make_shared<SAMRAI::hier::PatchDescriptor>();
  for (int id :
       fub::span(data_desc.data_ids.data(), data_desc.n_cons_variables)) {
    std::shared_ptr<SAMRAI::hier::Variable> variable{};
    if (vardb->mapIndexToVariable(id, variable)) {
      const std::string& name = variable->getName();
      const int depth =
          static_cast<SAMRAI::pdat::CellVariable<double>*>(variable.get())
              ->getDepth();
      flux_descriptor->definePatchDataComponent(
          name, std::make_shared<SAMRAI::pdat::SideDataFactory<double>>(
                    depth, face_ghosts, true));
    }
  }
  flux_descriptor->printClassData(SAMRAI::tbox::pout);
}