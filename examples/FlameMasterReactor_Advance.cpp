#include "fub/ideal_gas/FlameMasterReactor.hpp"
#include "fub/ideal_gas/mechanism/Gri30.hpp"
#include "fub/ideal_gas/mechanism/Zhao2008Dme.hpp"

#include "fub/ode_solver/CVodeSolver.hpp"

extern "C" {
#include "TC_defs.h"
#include "TC_interface.h"
}

#include <cstdio>

int main() {
  using namespace fub::ideal_gas;
  Zhao2008Dme mechanism;

  FlameMasterReactor reactor(mechanism);
  std::vector<double> X(reactor.getNSpecies());
  X[Zhao2008Dme::sH2] = 2.0;
  X[Zhao2008Dme::sO2] = 1.0;
  reactor.setMoleFractions(X);
  reactor.setTemperature(1300.);
  reactor.setPressure(101325.);
  reactor.advance(1e-4, [](fub::span<const double> z, double t,
                           FlameMasterReactor* reactor) {
    std::printf("time: %f, temperature: %f, pressure: %f\n", t,
                reactor->getTemperature(), reactor->getPressure());
  });

#ifdef FUB_WITH_SUNDIALS
  std::fill(X.begin(), X.end(), 0.0);
  X[Zhao2008Dme::sH2] = 2.0;
  X[Zhao2008Dme::sO2] = 1.0;
  reactor.setOdeSolver(
      std::make_unique<fub::CVodeSolver>(reactor.getNSpecies() + 1));
  reactor.setMoleFractions(X);
  reactor.setTemperature(1300.);
  reactor.setPressure(101325.);
  reactor.advance(1e-4, [](fub::span<const double> z, double t,
                           FlameMasterReactor* reactor) {
    std::printf("time: %f, temperature: %f, pressure: %f\n", t,
                reactor->getTemperature(), reactor->getPressure());
  });
#endif

  // try {
  //   X[Zhao2008Dme::sH2] = 0.5;
  //   X[Zhao2008Dme::sO2] = 0.5;
  //   reactor.setMoleFractions(X);
  //   reactor.setTemperature(1100.);
  //   reactor.setPressure(101325.);
  //   reactor.advance_tchem(1e-4, [](fub::span<const double> z, double t,
  //                                  FlameMasterReactor* reactor) {
  //     reactor->setMassFractions(z.subspan(1));
  //     reactor->setTemperature(z[0]);
  //     std::printf("time: %f, temperature: %f, pressure: %f\n", t,
  //                 reactor->getTemperature(), reactor->getPressure());
  //   });
  // } catch (...) {
  //   TC_reset();
  //   throw;
  // }
}