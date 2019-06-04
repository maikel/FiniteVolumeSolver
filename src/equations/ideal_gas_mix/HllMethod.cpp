#include "fub/equations/IdealGasMix.hpp"

namespace fub {

template class FluxMethod<
    Hll<IdealGasMix<1>, EinfeldtSignalVelocities<IdealGasMix<1>>>>;

template class FluxMethod<
    Hll<IdealGasMix<2>, EinfeldtSignalVelocities<IdealGasMix<2>>>>;

template class FluxMethod<
    Hll<IdealGasMix<3>, EinfeldtSignalVelocities<IdealGasMix<3>>>>;
} // namespace fub