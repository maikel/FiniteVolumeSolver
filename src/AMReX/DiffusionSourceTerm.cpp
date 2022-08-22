#include "fub/AMReX/DiffusionSourceTerm.hpp"

namespace fub::amrex {

DiffusionSourceTermOptions::DiffusionSourceTermOptions(
    const ProgramOptions& opts) {
  mul = GetOptionOr(opts, "mul", mul);
  reynolds = GetOptionOr(opts, "reynolds", reynolds);
  schmidt = GetOptionOr(opts, "schmidt", schmidt);
  prandtl = GetOptionOr(opts, "prandtl", prandtl);
}

void DiffusionSourceTermOptions::Print(SeverityLogger& log) const {
  BOOST_LOG(log) << fmt::format("  - mul = {}", mul);
  BOOST_LOG(log) << fmt::format("  - reynolds = {}", reynolds);
  BOOST_LOG(log) << fmt::format("  - schmidt = {}", schmidt);
  BOOST_LOG(log) << fmt::format("  - prandtl = {}", prandtl);
}

} // namespace fub::amrex