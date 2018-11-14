#include "fub/ideal_gas/FlameMasterReactor.hpp"
#include "fub/ideal_gas/mechanism/Gri30.hpp"
#include "fub/ideal_gas/mechanism/Zhao2008Dme.hpp"

extern "C" {
#include "TC_defs.h"
#include "TC_interface.h"
}

#include <cstdio>

int main() {
  using namespace fub::ideal_gas;
  TC_initChem("Zhao2008DME.ckm.chmech", "Zhao2008DME.ckm.chthermo", 1, 1.0);
  Zhao2008Dme mechanism;

  FlameMasterReactor reactor(mechanism);
  std::vector<double> X(reactor.getNSpecies());
  X[Zhao2008Dme::sH2] = 0.5;
  X[Zhao2008Dme::sO2] = 0.5;
  reactor.setMoleFractions(X);
  reactor.setTemperature(1100.);
  reactor.setPressure(101325.);
  reactor.advance(1e-4, [](fub::span<const double> z, double t,
                           FlameMasterReactor* reactor) {
    std::printf("time: %f, temperature: %f, pressure: %f\n", t,
                reactor->getTemperature(), reactor->getPressure());
    return 0;
  });

  try {
    X[Zhao2008Dme::sH2] = 0.5;
    X[Zhao2008Dme::sO2] = 0.5;
    reactor.setMoleFractions(X);
    reactor.setTemperature(1100.);
    reactor.setPressure(101325.);
    reactor.advance_tchem(1e-4, [](fub::span<const double> z, double t,
                                   FlameMasterReactor* reactor) {
      reactor->setMassFractions(z.subspan(1));
      reactor->setTemperature(z[0]);
      std::printf("time: %f, temperature: %f, pressure: %f\n", t,
                  reactor->getTemperature(), reactor->getPressure());
      return 0;
    });
  } catch (...) {
    TC_reset();
    throw;
  }
}