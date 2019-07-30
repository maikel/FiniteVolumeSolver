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

#include "fub/AMReX/PatchHierarchy.hpp"
#include "fub/AMReX/ScopeGuard.hpp"
#include "fub/AMReX/ViewFArrayBox.hpp"
#include "fub/equations/PerfectGas.hpp"

#include <AMReX_FArrayBox.H>

#define CATCH_CONFIG_RUNNER
#include <catch2/catch.hpp>

using Complete = fub::PerfectGas<2>::Complete;
using Conservative = fub::PerfectGas<2>::Conservative;
using CompleteArray = fub::PerfectGas<2>::CompleteArray;
using ConservativeArray = fub::PerfectGas<2>::ConservativeArray;

void InitializeStates(const fub::View<Complete>& states,
                      const fub::PerfectGas<2>& equation, double dx) {
  ForEachIndex(fub::Box<0>(states),
               [&states, dx](std::ptrdiff_t i, std::ptrdiff_t j) {
                 const double x = double(i) * dx;
                 const double y = double(j) * dx;
                 states.density(i, j) = 1.0;
                 states.momentum(i, j, 0) = 0.0;
                 states.momentum(i, j, 1) = 0.0;
                 states.energy(i, j) = 1.0 + 0.5 * std::sin(x * x + y * y);
               });
  CompleteFromCons(equation, states, AsConst(AsCons(states)));
}

TEST_CASE("simd single state") {
  fub::PerfectGas<2> equation{};
  std::array<Complete, 4> states{};

  states[0].density = 1.0;
  states[0].momentum.fill(0);
  states[0].energy = 2.0;
  equation.CompleteFromCons(states[0], AsCons(states[0]));

  states[1] = states[0];

  states[2].density = 1.0;
  states[2].momentum.fill(0);
  states[2].energy = 1.0;
  equation.CompleteFromCons(states[2], AsCons(states[2]));

  states[3] = states[2];

  fub::EinfeldtSignalVelocities<fub::PerfectGas<2>> signals{};
  fub::HllMethod hll(equation, signals);
  fub::MusclHancockMethod muscl(equation, hll);

  const fub::Duration dt(1.0);
  const double dx = 1.0;
  const fub::Direction dir = fub::Direction::X;
  Conservative flux{};
  muscl.ComputeNumericFlux(flux, states, dt, dx, dir);

  std::array<CompleteArray, 4> array_states;

  array_states[0].density.fill(1.0);
  array_states[0].momentum.fill(0.0);
  array_states[0].energy.fill(2.0);
  equation.CompleteFromCons(array_states[0], AsCons(array_states[0]));

  array_states[1] =  array_states[0];

  array_states[2].density.fill(1.0);
  array_states[2].momentum.fill(0);
  array_states[2].energy.fill(1.0);
  equation.CompleteFromCons(array_states[2], AsCons(array_states[2]));

  array_states[3] =  array_states[2];

  ConservativeArray array_flux{};
  muscl.ComputeNumericFlux(array_flux, array_states, dt, dx, dir);

  REQUIRE(flux.density == array_flux.density(0));
  REQUIRE(flux.momentum(0) == array_flux.momentum(0, 0));
  REQUIRE(flux.momentum(1) == array_flux.momentum(1, 0));
  REQUIRE(flux.energy == array_flux.energy(0));
}

TEST_CASE("simd impl produces the exact same results as the naive impl") {
  fub::PerfectGas<2> equation{};
  fub::amrex::DataDescription data_desc =
      fub::amrex::MakeDataDescription(equation);

  SECTION("Direction X") {
    amrex::Box cell_box{{AMREX_D_DECL(0, 0, 0)}, {AMREX_D_DECL(15, 15, 0)}};
    amrex::Box face_box = amrex::convert(cell_box, {AMREX_D_DECL(1, 0, 0)});
    cell_box.grow(0, 1);

    amrex::FArrayBox states_data(cell_box, data_desc.n_state_components);
    const double dx = 0.0125;
    fub::View<Complete> states = fub::amrex::MakeView<Complete>(
        states_data, equation, fub::amrex::AsIndexBox<2>(cell_box));
    InitializeStates(states, equation, dx);

    amrex::FArrayBox fluxes_data(face_box, data_desc.n_cons_components);
    amrex::FArrayBox simd_fluxes_data(face_box, data_desc.n_cons_components);

    fub::HllMethod hll(equation,
                       fub::EinfeldtSignalVelocities<fub::PerfectGas<2>>{});

    const fub::Duration dt(
        hll.ComputeStableDt(AsConst(states), dx, fub::Direction::X));

    fub::View<Conservative> fluxes = fub::amrex::MakeView<Conservative>(
        fluxes_data, equation, fub::amrex::AsIndexBox<2>(face_box));

    hll.ComputeNumericFluxes(fluxes, AsConst(states), dt, dx, fub::Direction::X);

    fluxes = fub::amrex::MakeView<Conservative>(
        simd_fluxes_data, equation, fub::amrex::AsIndexBox<2>(face_box));

    hll.ComputeNumericFluxes(fub::execution::simd, fluxes, AsConst(states),
                            dt, dx, fub::Direction::X);

    const double* flux_first = fluxes_data.dataPtr();
    const double* flux_last = flux_first + fluxes_data.size();
    const double* simd_flux_first = simd_fluxes_data.dataPtr();
    REQUIRE(std::equal(flux_first, flux_last, simd_flux_first));
  }

  SECTION("Direction Y") {
    amrex::Box cell_box{{AMREX_D_DECL(0, 0, 0)}, {AMREX_D_DECL(15, 15, 0)}};
    amrex::Box face_box = amrex::convert(cell_box, {AMREX_D_DECL(0, 1, 0)});
    cell_box.grow(1, 1);

    amrex::FArrayBox states_data(cell_box, data_desc.n_state_components);
    const double dx = 0.0125;
    fub::View<Complete> states = fub::amrex::MakeView<Complete>(
        states_data, equation, fub::amrex::AsIndexBox<2>(cell_box));
    InitializeStates(states, equation, dx);

    amrex::FArrayBox fluxes_data(face_box, data_desc.n_cons_components);
    amrex::FArrayBox simd_fluxes_data(face_box, data_desc.n_cons_components);

    fub::HllMethod hll(equation,
                       fub::EinfeldtSignalVelocities<fub::PerfectGas<2>>{});

    const fub::Duration dt(
        hll.ComputeStableDt(AsConst(states), dx, fub::Direction::Y));

    fub::View<Conservative> fluxes = fub::amrex::MakeView<Conservative>(
        fluxes_data, equation, fub::amrex::AsIndexBox<2>(face_box));

    hll.ComputeNumericFluxes(fluxes, AsConst(states), dt,
                             dx, fub::Direction::Y);

    fluxes = fub::amrex::MakeView<Conservative>(
        simd_fluxes_data, equation, fub::amrex::AsIndexBox<2>(face_box));

    hll.ComputeNumericFluxes(fub::execution::simd, fluxes, AsConst(states),
                             dt, dx, fub::Direction::Y);

    const double* flux_first = fluxes_data.dataPtr();
    const double* flux_last = flux_first + fluxes_data.size();
    const double* simd_flux_first = simd_fluxes_data.dataPtr();
    REQUIRE(std::equal(flux_first, flux_last, simd_flux_first));
  }
}

int main(int argc, char* argv[]) {
  fub::amrex::ScopeGuard _(argc, argv);
  return Catch::Session().run(argc, argv);
}
