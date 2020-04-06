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
#include "fub/equations/ideal_gas_mix/MusclHancockPrimMethod.hpp"

#define CATCH_CONFIG_MAIN
#include <catch2/catch.hpp>

void test_TPX(double T) {
  using namespace fub;
  Burke2012 burke{};
  FlameMasterReactor reactor{burke};
  reactor.SetMassFractions("N2:8,O2:2");
  reactor.SetTemperature(T);
  reactor.SetPressure(2.0e5);

  double cp = reactor.GetCp();
  double h = reactor.GetEnthalpy();
  double e = reactor.GetInternalEnergy();

  ArrayXd species(reactor.GetNSpecies(), kDefaultChunkSize);
  for (int s = 0; s < species.rows(); ++s) {
    species.row(s) = Array1d::Constant(reactor.GetMassFractions()[s]);
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
  Burke2012 burke{};
  FlameMasterReactor reactor{burke};
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

using fub::IdealGasMix;
using fub::MaskArray;
using fub::Array;
using fub::ArrayXd;
using fub::Array1d;
using fub::Burke2012;
using fub::FlameMasterReactor;

template <int Dim>
using CompleteArray = ::fub::CompleteArray<IdealGasMix<Dim>>;

template <int Dim>
using PrimitiveArray = ::fub::ideal_gas::PrimitiveArray<Dim>;

template <int Rank>
void CompleteFromPrim(IdealGasMix<Rank>& eq,
                      CompleteArray<Rank>& complete,
                      const PrimitiveArray<Rank>& prim, MaskArray mask) {
  FlameMasterReactor& reactor = eq.GetReactor();
  reactor.SetDensityArray(Array1d::Constant(1.0));
  reactor.SetMassFractionsArray(prim.mass_fractions);
  reactor.SetTemperatureArray(prim.temperature);
  reactor.SetPressureArray(prim.pressure);
  eq.CompleteFromReactor(complete, prim.velocity, mask);
}

template <int Rank>
void ToPrim(PrimitiveArray<Rank>& p, const CompleteArray<Rank>& q, MaskArray mask) {
  const Array1d zero = Array1d::Zero();
  p.pressure = mask.select(q.pressure, zero);
  Array1d rho_s = (q.density > 0.0).select(q.density, Array1d::Constant(1.0));
  for (int i = 0; i < Rank; ++i) {
    p.velocity.row(i) = mask.select(q.momentum.row(i) / rho_s, zero);
  }
  p.temperature = mask.select(q.temperature, zero);
  for (int i = 0; i < p.mass_fractions.rows(); ++i) {
    p.mass_fractions.row(i) = mask.select(q.species.row(i) / rho_s, zero);
  }
}

void test_ToPrim_CompleteFromReactor(IdealGasMix<3>& equation, CompleteArray<3>& result, PrimitiveArray<3>& primitive, const CompleteArray<3>& source, MaskArray mask) {
  ToPrim(primitive, source, mask);
  CompleteFromPrim(equation, result, primitive, mask);
  REQUIRE(result.density.isApprox(source.density));
  REQUIRE(result.momentum.isApprox(source.momentum));
  REQUIRE(result.species.isApprox(source.species));
  REQUIRE(result.pressure.isApprox(source.pressure));
  REQUIRE(result.energy.isApprox(source.energy));
  REQUIRE(result.temperature.isApprox(source.temperature));
  REQUIRE(result.c_p.isApprox(source.c_p));
  REQUIRE(result.gamma.isApprox(source.gamma));
}

TEST_CASE("CompleteFromPrim") {
  Burke2012 burke{};
  IdealGasMix<3> ideal_gas{burke};
  CompleteArray<3> source(ideal_gas);
  CompleteArray<3> result(ideal_gas);
  PrimitiveArray<3> primitive{};
  primitive.mass_fractions.resize(ideal_gas.GetReactor().GetNSpecies(), fub::kDefaultChunkSize);

  FlameMasterReactor& reactor = ideal_gas.GetReactor();
  ArrayXd masses(reactor.GetNSpecies(), fub::kDefaultChunkSize);

  MaskArray mask = MaskArray::Constant(true);

  double h2_first = 0.0;
  double h2_last = 1e-20;
  double dh2 = 1e-24;
  Array<double, 3> velocity;
  velocity.row(0) = Array1d::Constant(20.0);
  velocity.row(1) = Array1d::Constant(-21.0);
  velocity.row(2) = Array1d::Constant(10.0);
  for (double h2 = h2_first; h2 < h2_last; h2 += fub::kDefaultChunkSize * dh2) {
    Array1d h2a;
    for (int i = 0; i < fub::kDefaultChunkSize; ++i) {
      h2a[i] = h2 + i * dh2;
    }
    masses.setZero();
    masses.row(Burke2012::sN2) = Array1d::Constant(0.76711787578231128);
    masses.row(Burke2012::sO2) = Array1d::Constant(0.23288212421768872);
    masses.row(Burke2012::sH2) = h2a;
    reactor.SetMassFractionsArray(masses);
    reactor.SetTemperatureArray(Array1d::Constant(320.0));
    reactor.SetPressureArray(Array1d::Constant(101325.0));
    ideal_gas.CompleteFromReactor(source, velocity, mask);
    ideal_gas.CompleteFromCons(source, source, mask);

    test_ToPrim_CompleteFromReactor(ideal_gas, result, primitive, source, mask);
  }
}

