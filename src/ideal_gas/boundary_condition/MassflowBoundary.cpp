// Copyright (c) 2018 Maikel Nadolski
//
// Permission is hereby granted, free of charge, to any person obtaining a copy
// of this software and associated documentation files (the "Software"), to deal
// in the Software without restriction, including without limitation the rights
// to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
// copies of the Software, and to permit persons to whom the Software is
// furnished to do so, subject to the following conditions:
//
// The above copyright notice and this permission notice shall be included in
// all copies or substantial portions of the Software.
//
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
// IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
// AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
// OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
// SOFTWARE.

#include "fub/ideal_gas/boundary_condition/MassflowBoundary.hpp"

namespace fub {
namespace ideal_gas {

void MassflowBoundary::SetPhysicalBoundaryCondition(
    const SAMRAI::hier::Patch& patch, const SAMRAI::hier::Box& fill_box,
    [[maybe_unused]] double fill_time, Direction dir, int side) const {
  using Scratch = HyperbolicTimeIntegrator::Scratch;
  IdealGasEquation::CompleteStatePatchData complete =
      integrator_->getCompleteState(patch, Scratch(dir));
  const int direction = static_cast<int>(dir);
  const int sign = side == 0 ? +1 : -1;
  FlameMasterReactor& reactor = kinetics_->GetReactor();
  std::vector<double> fractions(complete.species.getDepth());
  for (const SAMRAI::hier::Index& index : fill_box) {
    SAMRAI::pdat::CellIndex to(index);
    SAMRAI::pdat::CellIndex from(shift(index, Direction(direction), sign));
    CopyMassFractions(fractions, complete.species, from);
    reactor.setMassFractions(fractions);
    reactor.setTemperature(complete.temperature(from));
    reactor.setPressure(complete.pressure(from));
    const double rho = complete.density(from);
    const double gamma = reactor.getCp() / reactor.getCv();
    const double c = complete.speed_of_sound(from);
    const double u = complete.momentum(from) / rho;
    const double M = u / c;
    const double gammaMinus = gamma - 1.0;
    const double gammaPlus = gamma + 1.0;
    const double rGammaPlus = 1.0 / gammaPlus;
    const double c_critical =
        std::sqrt(c * c + 0.5 * gammaMinus * u * u) * std::sqrt(2 * rGammaPlus);
    const double u_n = required_massflow_ / rho / surface_area_;
    const double lambda = u / c_critical;
    const double lambda_n = u_n / c_critical;
    const double gammaQuot = gammaMinus * rGammaPlus;
    const double p = complete.pressure(from);
    const double p0_n =
        p * std::pow(1. - gammaQuot * lambda_n * lambda_n, -gamma / gammaMinus);
    const double p_n =
        p0_n * std::pow(1. - gammaQuot * lambda * lambda, gamma / gammaMinus);
    reactor.setPressureIsentropic(p_n);
    complete.density(to) = reactor.getDensity();
    complete.momentum(to) = reactor.getDensity() * u_n;
    complete.energy(to) = reactor.getDensity() * reactor.getInternalEnergy() +
                          0.5 * reactor.getDensity() * u_n * u_n;
    complete.pressure(to) = reactor.getPressure();
    complete.temperature(to) = reactor.getTemperature();
    complete.speed_of_sound(to) =
        std::sqrt(gamma * reactor.getPressure() / reactor.getDensity());
    span<const double> Y = reactor.getMassFractions();
    for (int s = 0; s < Y.size(); ++s) {
      complete.species(to, s) = reactor.getDensity() * Y[s];
    }
  }
}

} // namespace ideal_gas
} // namespace fub