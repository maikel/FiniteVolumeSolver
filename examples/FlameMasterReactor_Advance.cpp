#include "fub/ideal_gas/FlameMasterReactor.hpp"
#include "fub/ideal_gas/mechanism/Gri30.hpp"

#include <cstdio>

int main() {
  using namespace fub::ideal_gas;
  Gri30 mechanism;
  FlameMasterReactor reactor(mechanism);
  std::vector<double> X(reactor.getNSpecies());
  X[Gri30::sH2] = 0.5;
  X[Gri30::sO2] = 0.5;
  reactor.setMoleFractions(X);
  reactor.setTemperature(1100.);
  reactor.setPressure(101325.);
  reactor.advance(1e-4, [](double t, FlameMasterReactor* reactor) {
    std::printf("time: %f, temperature: %f, pressure: %f\n", t,
                reactor->getTemperature(), reactor->getPressure());
    return 0;
  });
}