#include "fub/SAMRAI/ScopeGuard.hpp"
#include "fub/ideal_gas/TChemKinetics.hpp"

extern "C" {
#include "TC_defs.h"
#include "TC_interface.h"
}

#include <cstdio>

int main(int argc, char** argv) {
  using namespace fub::ideal_gas;
  using namespace fub;
  fub::ScopeGuard guard(argc, argv);
  TChemKineticsOptions options;
  options.chemfile = "chem.inp";
  options.thermofile = "therm.dat";
  TChemKinetics equation("IdealGas", SAMRAI::tbox::Dimension(1), options);
  std::vector<double> X(equation.GetNSpecies());
  std::vector<double> TY(equation.GetNSpecies() + 1);
  X[0] = 0.5;
  X[3] = 0.5;
  span<double> Y = make_span(TY).subspan(1);
  TY[0] = 1100.;
  TC_getMl2Ms(X.data(), X.size(), Y.data());
  TC_setThermoPres(101325.0);
  double rho = 0.0;
  TC_getRhoMixMs(TY.data(), TY.size(), &rho);
  FUB_ASSERT(rho > 0.0);
  equation.AdvanceSourceTerm(
      TY, rho, 1e-4, [](span<const double> TandY, double t) {
        std::printf("time: %f, temperature: %f, pressure: %f\n", t, TandY[0],
                    TC_getThermoPres());
        return 0;
      });
}