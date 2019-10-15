// Copyright (c) 2019 Maikel Nadolski
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

#include "fub/equations/IdealGasMix.hpp"
#include "fub/equations/ideal_gas_mix/mechanism/Burke2012.hpp"

#define CATCH_CONFIG_MAIN
#include <catch2/catch.hpp>

void test_TPX(double T) {
  using namespace fub;
  FlameMasterReactor reactor{Burke2012{}};
  reactor.SetMassFractions("N2:8,O2:2");
  reactor.SetTemperature(T);
  reactor.SetPressure(2.0e5);

  double cp = reactor.GetCp();
  double h = reactor.GetEnthalpy();
  double e = reactor.GetInternalEnergy();

  ArrayXd species(reactor.GetNSpecies(), kDefaultChunkSize);
  for (int i = 0; i < species.rows(); ++i) {
    species.row(i) = Array1d::Constant(reactor.GetMassFractions()[i]);
  }

  reactor.SetMassFractionsArray(species);
  reactor.SetTemperatureArray(Array1d::Constant(T));
  reactor.SetPressureArray(Array1d::Constant(2.0e5));

  Array1d cp_array = reactor.GetCpArray();
  Array1d h_array = reactor.GetEnthalpyArray();
  Array1d e_array = reactor.GetInternalEnergyArray();

  REQUIRE((cp == cp_array).all());
  REQUIRE((h == h_array).all());
  REQUIRE((e == e_array).all());
}

TEST_CASE("TPX") {
  double T0 = 200;
  double TEnd = 3000;
  double dT = 25;

  for (double T = T0; T < TEnd; T += dT) {
    test_TPX(T);
  }
}

void test_SetInternalEnergy(double T) {
  using namespace fub;
  FlameMasterReactor reactor{Burke2012{}};
  reactor.SetMoleFractions("N2:8,O2:2,H2:4");
  reactor.SetTemperature(T);
  reactor.SetPressure(101325.0);

  double e = reactor.GetInternalEnergy();

  ArrayXd species(reactor.GetNSpecies(), kDefaultChunkSize);
  for (int i = 0; i < species.rows(); ++i) {
    species.row(i) = Array1d::Constant(reactor.GetMassFractions()[i]);
  }

  reactor.SetMassFractionsArray(species);
  reactor.SetTemperatureArray(Array1d::Constant(300.0));
  reactor.SetPressureArray(Array1d::Constant(101325.0));

  Array1d ev = Array1d::Constant(e);
//  ev[7] = -406228.18662658607;
  reactor.SetInternalEnergyArray(ev);

  Array1d Tv = Array1d::Constant(T);
//  Tv[7] = 300.0;

  REQUIRE(((reactor.GetTemperatureArray() - Tv).abs() < 1e-6).all());
}

TEST_CASE("SetInternalEnergy") {
  double T0 = 200;
  double TEnd = 3000;
  double dT = 25;

  for (double T = T0; T < TEnd; T += dT) {
    test_SetInternalEnergy(T);
  }
}

