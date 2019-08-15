#include "fub/equations/PerfectGas.hpp"

namespace fub {

template class FluxMethod<Godunov<PerfectGas<1>>>;
template class FluxMethod<Godunov<PerfectGas<2>>>;
template class FluxMethod<Godunov<PerfectGas<3>>>;

} // namespace fub
