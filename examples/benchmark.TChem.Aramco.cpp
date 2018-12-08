#include "fub/SAMRAI/ScopeGuard.hpp"

#include "fub/ideal_gas/FlameMasterReactor.hpp"
#include "fub/ideal_gas/TChemKinetics.hpp"
#include "fub/ideal_gas/mechanism/AramcoMech_DMEonly_74spec.hpp"

#include "fub/ode_solver/CVodeSolver.hpp"

extern "C" {
#include "TC_defs.h"
#include "TC_interface.h"
}

#include <benchmark/benchmark.h>

#include <cstdio>

using Aramco = fub::ideal_gas::AramcoMech_DMEonly_74spec;

using namespace fub::ideal_gas;
Aramco mechanism;
FlameMasterReactor fm_reactor(mechanism);
TChemReactor tc_reactor(mechanism);

static void AdvanceSourceTerm_FlameMaster_Aramco(benchmark::State& state) {
  std::vector<double> X(fm_reactor.getNSpecies());
  X[Aramco::sH2] = 2;
  X[Aramco::sO2] = 1;
  for (auto _ : state) {
    fm_reactor.setMoleFractions(X);
    fm_reactor.setTemperature(1300.);
    fm_reactor.setPressure(101325.);
    fm_reactor.advance(1e-4);
  }
}
BENCHMARK(AdvanceSourceTerm_FlameMaster_Aramco);

static void
AdvanceSourceTerm_FlameMaster_Aramco_NoReac(benchmark::State& state) {
  std::vector<double> X(fm_reactor.getNSpecies());
  X[Aramco::sH2] = 2;
  X[Aramco::sO2] = 1;
  for (auto _ : state) {
    fm_reactor.setMoleFractions(X);
    fm_reactor.setTemperature(300.);
    fm_reactor.setPressure(101325.);
    fm_reactor.advance(1e-4);
  }
}
BENCHMARK(AdvanceSourceTerm_FlameMaster_Aramco_NoReac);

static void AdvanceSourceTerm_TChem_Aramco(benchmark::State& state) {
  std::vector<double> X(tc_reactor.GetNSpecies());
  X[Aramco::sH2] = 2;
  X[Aramco::sO2] = 1;
  for (auto _ : state) {
    tc_reactor.SetMoleFractions(X);
    tc_reactor.SetPressure(101325.);
    tc_reactor.SetTemperature(1100.);
    tc_reactor.Advance(1e-4);
    FUB_ASSERT(tc_reactor.GetTemperature() > 3000.);
  }
}
BENCHMARK(AdvanceSourceTerm_TChem_Aramco);

static void
AdvanceSourceTerm_TChem_Aramco_NoReac(benchmark::State& state) {
  std::vector<double> X(tc_reactor.GetNSpecies());
  X[Aramco::sH2] = 2;
  X[Aramco::sO2] = 1;
  for (auto _ : state) {
    tc_reactor.SetMoleFractions(X);
    tc_reactor.SetPressure(101325.);
    tc_reactor.SetTemperature(300.);
    tc_reactor.Advance(1e-4);
    FUB_ASSERT(tc_reactor.GetTemperature() < 1000.);
  }
}
BENCHMARK(AdvanceSourceTerm_TChem_Aramco_NoReac);

int main(int argc, char** argv) {
  tc_reactor.SetOdeSolver(std::make_unique<fub::CVodeSolver>(tc_reactor.GetNSpecies() + 1));
  fub::ScopeGuard guard(argc, argv);
  benchmark::Initialize(&argc, argv);
  benchmark::RunSpecifiedBenchmarks();
}