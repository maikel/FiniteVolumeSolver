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

#include "fub/ode_solver/CVodeSolver.hpp"
#include "fub/ode_solver/RadauSolver.hpp"

#include <chrono>
#include <cstdio>
#include <vector>

extern "C" {
#include "TC_defs.h"
#include "TC_interface.h"
}

#include "fub/ideal_gas/FlameMasterReactor.hpp"
#include "fub/ideal_gas/mechanism/Zhao2008Dme.hpp"

void IntegrateTChemSystem(const fub::OdeSolver& solver) noexcept {
  using namespace fub;
  TC_setThermoPres(101325.0);
  auto source_term = [](span<double> dydt, span<const double> y, double t) {
    TC_getSrcCV(const_cast<double*>(y.data()), y.size(), dydt.data());
  };
  auto jacobian = [](span<double> jac, span<const double> TandYs, double t) {
    int Ns = TandYs.size() - 1;
    TC_getJacCVTYNanl(const_cast<double*>(TandYs.data()), Ns, jac.data());
  };
  std::vector<double> TY(1 + TC_getNspec());
  // Set initial temperature to 1100K
  TY[0] = 1100.;
  // Set mixture to be a mole ratio of 1:1 of H2 and O2
  TY[2] = 0.5;  // H2
  TY[22] = 0.5; // O2
  TC_getMl2Ms(TY.data() + 1, TC_getNspec(), TY.data() + 1);
  double rhomix = 0.0;
  TC_getRhoMixMs(TY.data(), TY.size(), &rhomix);
  TC_setDens(rhomix);
  auto start = std::chrono::steady_clock::now();
  solver.Integrate(source_term, TY, 0.0, 1e-4, jacobian);
  // solver.Integrate(source_term, TY, 0.0, 1e-4);
  auto end = std::chrono::steady_clock::now();
  std::chrono::duration<double> elapsed = end - start;
  std::printf("elapsed time: %e, T = %f\n", elapsed.count(), TY[0]);
}

void IntegrateFMSystem() noexcept {
  using namespace fub;
  ideal_gas::Zhao2008Dme mechanism;
  ideal_gas::FlameMasterReactor reactor = mechanism;
  std::vector<double> X(reactor.getNSpecies());
  X[ideal_gas::Zhao2008Dme::sH2] = 0.5;
  X[ideal_gas::Zhao2008Dme::sO2] = 0.5;

  auto start = std::chrono::steady_clock::now();
  reactor.setMoleFractions(X);
  reactor.setPressure(101325.0);
  reactor.setTemperature(1100.0);
  reactor.advance(1e-4);
  auto end = std::chrono::steady_clock::now();
  std::chrono::duration<double> elapsed = end - start;
  std::printf("elapsed time: %e, T = %f\n", elapsed.count(),
              reactor.getTemperature());
}

int temperature_root_function(fub::span<double> T,
                              fub::span<const double> TandY, double timepoint) {
  T[0] = TandY[0];
  return 0;
}

int main() {
  TC_initChem("Zhao2008DME.ckm.chmech", "Zhao2008DME.ckm.chthermo", 1, 1.0);
  fub::RadauSolver radau;
#ifdef FUB_WITH_SUNDIALS
  fub::CVodeSolver cvode(1 + TC_getNspec());
  cvode.SetRootFunction(temperature_root_function, 1);
#endif
  for (int i = 0; i < 10; ++i) {
    std::printf("=================\nRadau: ");
    IntegrateTChemSystem(radau);
#ifdef FUB_WITH_SUNDIALS
    std::printf("=================\nCVode: ");
    IntegrateTChemSystem(cvode);
#endif
    std::printf("=================\nFlamemaster: ");
    IntegrateFMSystem();
  }
  TC_reset();
}