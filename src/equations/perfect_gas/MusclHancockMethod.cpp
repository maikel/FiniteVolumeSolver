#include "fub/equations/PerfectGas.hpp"

namespace fub {

template class FluxMethod<MusclHancock<PerfectGas<1>>>;
template class FluxMethod<MusclHancock<PerfectGas<2>>>;
template class FluxMethod<MusclHancock<PerfectGas<3>>>;

} // namespace fub
