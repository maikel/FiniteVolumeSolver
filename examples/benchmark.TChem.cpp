#include "fub/SAMRAI/ScopeGuard.hpp"

#include "fub/ideal_gas/FlameMasterReactor.hpp"
#include "fub/ideal_gas/TChemKinetics.hpp"
#include "fub/ideal_gas/mechanism/Burke2012.hpp"
#include "fub/ideal_gas/mechanism/Gri30.hpp"
#include "fub/ideal_gas/mechanism/Zhao2008Dme.hpp"

#include "fub/ode_solver/CVodeSolver.hpp"

extern "C" {
#include "TC_defs.h"
#include "TC_interface.h"
}

#include <benchmark/benchmark.h>

#include <cstdio>

// static void AdvanceSourceTerm_FlameMaster_Burke2012(benchmark::State& state)
// {
//   using namespace fub::ideal_gas;
//   Burke2012 mechanism;
//   FlameMasterReactor reactor(mechanism);
//   std::vector<double> X(reactor.getNSpecies());
//   X[Burke2012::sH2] = 0.5;
//   X[Burke2012::sO2] = 0.5;
//   for (auto _ : state) {
//     reactor.setMoleFractions(X);
//     reactor.setTemperature(1000.);
//     reactor.setPressure(101325.);
//     reactor.advance(1e-4);
//   }
// }
// BENCHMARK(AdvanceSourceTerm_FlameMaster_Burke2012);

// static void AdvanceSourceTerm_FlameMaster_Gri30(benchmark::State& state) {
//   using namespace fub::ideal_gas;
//   Gri30 mechanism;
//   FlameMasterReactor reactor(mechanism);
//   std::vector<double> X(reactor.getNSpecies());
//   X[Gri30::sH2] = 0.5;
//   X[Gri30::sO2] = 0.5;
//   for (auto _ : state) {
//     reactor.setMoleFractions(X);
//     reactor.setTemperature(1100.);
//     reactor.setPressure(101325.);
//     reactor.advance(1e-4);
//   }
// }
// BENCHMARK(AdvanceSourceTerm_FlameMaster_Gri30);

// static void
// AdvanceSourceTerm_FlameMaster_Gri30_NoReac(benchmark::State& state) {
//   using namespace fub::ideal_gas;
//   Gri30 mechanism;
//   FlameMasterReactor reactor(mechanism);
//   std::vector<double> X(reactor.getNSpecies());
//   X[Gri30::sH2] = 0.5;
//   X[Gri30::sO2] = 0.5;
//   for (auto _ : state) {
//     reactor.setMoleFractions(X);
//     reactor.setTemperature(300.);
//     reactor.setPressure(101325.);
//     reactor.advance(1e-4);
//   }
// }
// BENCHMARK(AdvanceSourceTerm_FlameMaster_Gri30_NoReac);

static void AdvanceSourceTerm_FlameMaster_Zhao2008Dme(benchmark::State& state) {
  using namespace fub::ideal_gas;
  Zhao2008Dme mechanism;
  FlameMasterReactor reactor(mechanism);
  std::vector<double> X(reactor.getNSpecies());
  X[Zhao2008Dme::sH2] = 0.5;
  X[Zhao2008Dme::sO2] = 0.5;
  for (auto _ : state) {
    reactor.setMoleFractions(X);
    reactor.setTemperature(1100.);
    reactor.setPressure(101325.);
    reactor.advance(1e-4);
  }
}
BENCHMARK(AdvanceSourceTerm_FlameMaster_Zhao2008Dme);

static void
AdvanceSourceTerm_FlameMaster_Zhao2008Dme_NoReac(benchmark::State& state) {
  using namespace fub::ideal_gas;
  Zhao2008Dme mechanism;
  FlameMasterReactor reactor(mechanism);
  std::vector<double> X(reactor.getNSpecies());
  X[Zhao2008Dme::sH2] = 0.5;
  X[Zhao2008Dme::sO2] = 0.5;
  for (auto _ : state) {
    reactor.setMoleFractions(X);
    reactor.setTemperature(300.);
    reactor.setPressure(101325.);
    reactor.advance(1e-4);
  }
}
BENCHMARK(AdvanceSourceTerm_FlameMaster_Zhao2008Dme_NoReac);

// static void AdvanceSourceTerm_TChem_Burke2012(benchmark::State& state) {
//   using namespace fub;
//   using namespace fub::ideal_gas;
//   TChemKineticsOptions options;
//   options.chemfile = "Burke.inp";
//   options.thermofile = "BurkeTherm.dat";
//   TChemKinetics equation("IdealGas", SAMRAI::tbox::Dimension(1), options);
//   std::vector<double> X(equation.GetNSpecies());
//   std::vector<double> TY(equation.GetNSpecies() + 1);
//   X[1] = 0.5;
//   X[5] = 0.5;
//   span<double> Y = make_span(TY).subspan(1);
//   for (auto _ : state) {
//     TY[0] = 1000.;
//     TC_getMl2Ms(X.data(), X.size(), Y.data());
//     double rho = 0.0;
//     TC_setThermoPres(101325.0);
//     TC_getRhoMixMs(TY.data(), TY.size(), &rho);
//     equation.AdvanceSourceTerm(TY, rho, 1e-4);
//   }
// }
// BENCHMARK(AdvanceSourceTerm_TChem_Burke2012);

// static void AdvanceSourceTerm_TChem_Gri30(benchmark::State& state) {
//   using namespace fub;
//   using namespace fub::ideal_gas;
//   TChemKineticsOptions options;
//   options.chemfile = "chem.inp";
//   options.thermofile = "therm.dat";
//   TChemKinetics equation("IdealGas", SAMRAI::tbox::Dimension(1), options);
//   std::vector<double> X(equation.GetNSpecies());
//   std::vector<double> TY(equation.GetNSpecies() + 1);
//   X[0] = 0.5;
//   X[3] = 0.5;
//   span<double> Y = make_span(TY).subspan(1);
//   for (auto _ : state) {
//     TY[0] = 1100.;
//     TC_getMl2Ms(X.data(), X.size(), Y.data());
//     double rho = 0.0;
//     TC_setThermoPres(101325.0);
//     TC_getRhoMixMs(TY.data(), TY.size(), &rho);
//     equation.AdvanceSourceTerm(TY, rho, 1e-4);
//   }
// }
// BENCHMARK(AdvanceSourceTerm_TChem_Gri30);

// static void AdvanceSourceTerm_TChem_Gri30_NoReac(benchmark::State& state) {
//   using namespace fub;
//   using namespace fub::ideal_gas;
//   TChemKineticsOptions options;
//   options.chemfile = "chem.inp";
//   options.thermofile = "therm.dat";
//   TChemKinetics equation("IdealGas", SAMRAI::tbox::Dimension(1), options);
//   std::vector<double> X(equation.GetNSpecies());
//   std::vector<double> TY(equation.GetNSpecies() + 1);
//   X[0] = 0.5;
//   X[3] = 0.5;
//   span<double> Y = make_span(TY).subspan(1);
//   for (auto _ : state) {
//     TY[0] = 300.;
//     TC_getMl2Ms(X.data(), X.size(), Y.data());
//     double rho = 0.0;
//     TC_setThermoPres(101325.0);
//     TC_getRhoMixMs(TY.data(), TY.size(), &rho);
//     equation.AdvanceSourceTerm(TY, rho, 1e-4);
//   }
// }
// BENCHMARK(AdvanceSourceTerm_TChem_Gri30_NoReac);

using namespace fub::ideal_gas;
Zhao2008Dme mechanism;
TChemReactor reactor(mechanism);

static void AdvanceSourceTerm_TChem_Zhao2008Dme(benchmark::State& state) {
  std::vector<double> X(reactor.GetNSpecies());
  X[Zhao2008Dme::sH2] = 0.5;
  X[Zhao2008Dme::sO2] = 0.5;
  for (auto _ : state) {
    reactor.SetMoleFractions(X);
    reactor.SetPressure(101325.);
    reactor.SetTemperature(1100.);
    reactor.Advance(1e-4);
    FUB_ASSERT(reactor.GetTemperature() > 3000.);
  }
}
BENCHMARK(AdvanceSourceTerm_TChem_Zhao2008Dme);

static void
AdvanceSourceTerm_TChem_Zhao2008Dme_NoReac(benchmark::State& state) {
  std::vector<double> X(reactor.GetNSpecies());
  X[Zhao2008Dme::sH2] = 0.5;
  X[Zhao2008Dme::sO2] = 0.5;
  for (auto _ : state) {
    reactor.SetMoleFractions(X);
    reactor.SetPressure(101325.);
    reactor.SetTemperature(300.);
    reactor.Advance(1e-4);
    FUB_ASSERT(reactor.GetTemperature() < 1000.);
  }
}
BENCHMARK(AdvanceSourceTerm_TChem_Zhao2008Dme_NoReac);

int main(int argc, char** argv) {
  reactor.SetOdeSolver(std::make_unique<fub::CVodeSolver>(reactor.GetNSpecies() + 1));
  fub::ScopeGuard guard(argc, argv);
  benchmark::Initialize(&argc, argv);
  benchmark::RunSpecifiedBenchmarks();
}