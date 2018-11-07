#include "fub/initial_data/RiemannProblem.hpp"

namespace fub {

void RiemannProblem::InitializeDataOnPatch(
    const SAMRAI::hier::Patch& patch) const {
  const SAMRAI::hier::Box& box = patch.getBox();
  const SAMRAI::geom::CartesianPatchGeometry& geom =
      *GetCartesianPatchGeometry(patch);
  for (const SAMRAI::hier::Index& index : box) {
    Coordinates x = ComputeCellCoordinates(geom, box, index);
    if (geometry_.ComputeDistanceTo(x) <= 0) {
      FillLeftState(patch, index);
    } else {
      FillRightState(patch, index);
    }
  }
  PostInitialize(patch);
}

} // namespace fub