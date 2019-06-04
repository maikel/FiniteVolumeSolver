#include "fub/equations/ShallowWater.hpp"

namespace fub {

template class FluxMethod<MusclHancock<ShallowWater>>;
template class FluxMethod<MusclHancock<
    ShallowWater, HllMethod<ShallowWater, ShallowWaterSignalVelocities>>>;

} // namespace fub