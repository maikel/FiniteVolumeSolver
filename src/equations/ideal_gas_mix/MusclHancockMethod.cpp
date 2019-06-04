#include "fub/equations/IdealGasMix.hpp"

namespace fub {

template class FluxMethod<MusclHancock<
    IdealGasMix<1>,
    Hll<IdealGasMix<1>, EinfeldtSignalVelocities<IdealGasMix<1>>>>>;

template class FluxMethod<MusclHancock<
    IdealGasMix<2>,
    Hll<IdealGasMix<2>, EinfeldtSignalVelocities<IdealGasMix<2>>>>>;

template class FluxMethod<MusclHancock<
    IdealGasMix<3>,
    Hll<IdealGasMix<3>, EinfeldtSignalVelocities<IdealGasMix<3>>>>>;

} // namespace fub