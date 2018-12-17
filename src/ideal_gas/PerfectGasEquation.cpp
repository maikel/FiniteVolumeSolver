#include "fub/ideal_gas/PerfectGasEquation.hpp"

namespace fub {
namespace ideal_gas {
void PerfectGasEquation::FillFromCons(const CompleteStatePatchData& q,
                                      const ConsStatePatchData& u) const {
  q.density.copy(u.density);
  q.momentum.copy(u.momentum);
  q.energy.copy(u.energy);
  const SAMRAI::hier::Box& box = q.density.getBox();
  for (const SAMRAI::hier::Index& index : box) {
    SAMRAI::pdat::CellIndex i(index);
    const double rho_eps =
        (u.energy(i) - 0.5 * q.momentum(i) * q.momentum(i) / q.density(i));
    constexpr double gamma = 1.4;
    constexpr double gamma_minus_1 = gamma - 1.0;
    q.species(i) = q.density(i);
    q.pressure(i) = rho_eps * gamma_minus_1;
    q.temperature(i) = q.pressure(i) / q.density(i);
    q.speed_of_sound(i) = std::sqrt(gamma * q.temperature(i));
  }
}

void PerfectGasEquation::FillFromPrim(const CompleteStatePatchData& q,
                                      const PrimStatePatchData& w) const {
  q.temperature.copy(w.temperature);
  q.momentum.copy(w.momentum);
  q.pressure.copy(w.pressure);
  const SAMRAI::hier::Box& box = q.density.getBox();
  for (const SAMRAI::hier::Index& index : box) {
    SAMRAI::pdat::CellIndex i(index);
    constexpr double gamma = 1.4;
    constexpr double gamma_minus_1 = gamma - 1.0;
    q.density(i) = q.temperature(i) / q.pressure(i);
    q.momentum(i) *= q.density(i);
    q.species(i) = q.density(i);
    q.energy(i) = q.pressure(i) / gamma_minus_1 +
                  0.5 * q.momentum(i) * q.momentum(i) / q.density(i);
    q.speed_of_sound(i) = std::sqrt(gamma * q.temperature(i));
  }
}

} // namespace ideal_gas
} // namespace fub