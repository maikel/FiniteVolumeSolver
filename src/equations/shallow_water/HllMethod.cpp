#include "fub/equations/ShallowWater.hpp"

namespace fub {

template class FluxMethod<Hll<ShallowWater, ShallowWaterSignalVelocities>>;

}