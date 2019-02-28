#include "SAMRAI/hier/Box.h"
#include "SAMRAI/hier/BoxLevel.h"
#include "SAMRAI/hier/VariableDatabase.h"

#include "SAMRAI/geom/CartesianGridGeometry.h"
#include <SAMRAI/geom/CartesianOutersideDoubleWeightedAverage.h>

#include "SAMRAI/pdat/OutersideData.h"
#include "SAMRAI/pdat/OutersideVariable.h"
#include "SAMRAI/pdat/SideData.h"
#include "SAMRAI/pdat/SideVariable.h"

#include <SAMRAI/xfer/CoarsenAlgorithm.h>

#include <memory>

int main(int argc, char** argv) {
  using namespace SAMRAI;
  using namespace tbox;
  using namespace hier;
  using namespace pdat;
  using namespace geom;
  using namespace xfer;
  using namespace std;

  SAMRAI_MPI::init(&argc, &argv);
  SAMRAIManager::initialize();
  SAMRAIManager::startup();

  {
    Dimension dim(2);
    VariableDatabase* vars = VariableDatabase::getDatabase();

    shared_ptr<VariableContext> scratch = vars->getContext("scratch");
    auto flux =
        make_shared<SideVariable<double>>(dim, "flux", std::vector<int>{1, 0});
    auto fluxsum = make_shared<OutersideVariable<double>>(dim, "fluxsum");
    int fid = vars->registerVariableAndContext(flux, scratch,
                                               IntVector::getZero(dim));
    int sid = vars->registerVariableAndContext(fluxsum, scratch,
                                               IntVector::getZero(dim));
    shared_ptr<PatchDescriptor> patch_descriptor = vars->getPatchDescriptor();
    vars->printClassData(std::cout);

    BoxContainer coarse_domain(Box(Index(0, 0), Index(3, 3), BlockId(0)));
    BoxContainer fine_domain(Box(Index(2, 2), Index(5, 5), BlockId(0)));

    double xlo[] = {-1, -1};
    double xup[] = {+1, +1};
    auto geom =
        make_shared<CartesianGridGeometry>("Geometry", xlo, xup, coarse_domain);

    IntVector ratio(dim, 2);
    shared_ptr<PatchHierarchy> hierarchy =
        make_shared<PatchHierarchy>("Hierarchy", geom);
    hierarchy->setMaxNumberOfLevels(2);
    hierarchy->setRatioToCoarserLevel(ratio, 1);

    {
      const int owner_rank = 0;
      BoxLevel layer0(IntVector(dim, 1), geom);
      layer0.addBox(Box(*coarse_domain.begin(), LocalId(0), owner_rank));

      BoxLevel layer1(ratio, geom);
      layer1.addBox(Box(*fine_domain.begin(), LocalId(0), owner_rank));

      hierarchy->makeNewPatchLevel(0, layer0);
      hierarchy->makeNewPatchLevel(1, layer1);
    }

    shared_ptr<PatchLevel> coarse = hierarchy->getPatchLevel(0);
    coarse->allocatePatchData(fid);
    for (shared_ptr<Patch> patch : *coarse) {
      shared_ptr<SideData<double>> data =
          static_pointer_cast<SideData<double>>(patch->getPatchData(fid));
      data->fill(42);
      data->print(data->getBox(), std::cout);
    }

    shared_ptr<PatchLevel> fine = hierarchy->getPatchLevel(1);
    fine->allocatePatchData(sid);
    for (shared_ptr<Patch> patch : *fine) {
      shared_ptr<OutersideData<double>> data =
          static_pointer_cast<OutersideData<double>>(patch->getPatchData(sid));
      data->fill(24);
      data->print(data->getBox(), std::cout);
    }

    CoarsenAlgorithm fluxsum_coarsen{dim};
    fluxsum_coarsen.registerCoarsen(
        fid, sid, make_shared<CartesianOutersideDoubleWeightedAverage>());
    shared_ptr<CoarsenSchedule> sched =
        fluxsum_coarsen.createSchedule(coarse, fine);
    sched->coarsenData();

    for (shared_ptr<Patch> patch : *coarse) {
      shared_ptr<SideData<double>> data =
          static_pointer_cast<SideData<double>>(patch->getPatchData(fid));
      data->print(data->getBox(), std::cout);
    }
  }

  SAMRAIManager::shutdown();
  SAMRAIManager::finalize();
  SAMRAI_MPI::finalize();
}