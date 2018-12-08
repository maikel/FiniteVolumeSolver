#include "fub/SAMRAI/ScopeGuard.hpp"

#include "fub/ideal_gas/FlameMasterReactor.hpp"
#include "fub/ideal_gas/TChemKinetics.hpp"
#include "fub/ideal_gas/mechanism/Zhao2008Dme.hpp"

#include "fub/ode_solver/CVodeSolver.hpp"

extern "C" {
#include "TC_defs.h"
#include "TC_interface.h"
}

#include <benchmark/benchmark.h>

#include <cstdio>

using namespace fub::ideal_gas;
Zhao2008Dme mechanism;
FlameMasterReactor fm_reactor(mechanism);
TChemReactor tc_reactor(mechanism);

static void AdvanceSourceTerm_FlameMaster_Zhao2008Dme(benchmark::State& state) {
  std::vector<double> X(fm_reactor.getNSpecies());
  X[Zhao2008Dme::sH2] = 2;
  X[Zhao2008Dme::sO2] = 1;
  for (auto _ : state) {
    fm_reactor.setMoleFractions(X);
    fm_reactor.setTemperature(1300.);
    fm_reactor.setPressure(101325.);
    fm_reactor.advance(1e-4);
  }
}
BENCHMARK(AdvanceSourceTerm_FlameMaster_Zhao2008Dme);

static void
AdvanceSourceTerm_FlameMaster_Zhao2008Dme_NoReac(benchmark::State& state) {
  std::vector<double> X(fm_reactor.getNSpecies());
  X[Zhao2008Dme::sH2] = 2;
  X[Zhao2008Dme::sO2] = 1;
  for (auto _ : state) {
    fm_reactor.setMoleFractions(X);
    fm_reactor.setTemperature(300.);
    fm_reactor.setPressure(101325.);
    fm_reactor.advance(1e-4);
  }
}
BENCHMARK(AdvanceSourceTerm_FlameMaster_Zhao2008Dme_NoReac);

static void AdvanceSourceTerm_TChem_Zhao2008Dme(benchmark::State& state) {
  std::vector<double> X(tc_reactor.GetNSpecies());
  X[Zhao2008Dme::sH2] = 2;
  X[Zhao2008Dme::sO2] = 1;
  for (auto _ : state) {
    tc_reactor.SetMoleFractions(X);
    tc_reactor.SetPressure(101325.);
    tc_reactor.SetTemperature(1100.);
    tc_reactor.Advance(1e-4);
    FUB_ASSERT(tc_reactor.GetTemperature() > 3000.);
  }
}
BENCHMARK(AdvanceSourceTerm_TChem_Zhao2008Dme);

static void
AdvanceSourceTerm_TChem_Zhao2008Dme_NoReac(benchmark::State& state) {
  std::vector<double> X(tc_reactor.GetNSpecies());
  X[Zhao2008Dme::sH2] = 2;
  X[Zhao2008Dme::sO2] = 1;
  for (auto _ : state) {
    tc_reactor.SetMoleFractions(X);
    tc_reactor.SetPressure(101325.);
    tc_reactor.SetTemperature(300.);
    tc_reactor.Advance(1e-4);
    FUB_ASSERT(tc_reactor.GetTemperature() < 1000.);
  }
}
BENCHMARK(AdvanceSourceTerm_TChem_Zhao2008Dme_NoReac);

int main(int argc, char** argv) {
  tc_reactor.SetOdeSolver(std::make_unique<fub::CVodeSolver>(tc_reactor.GetNSpecies() + 1));
  fub::ScopeGuard guard(argc, argv);
  benchmark::Initialize(&argc, argv);
  benchmark::RunSpecifiedBenchmarks();
}