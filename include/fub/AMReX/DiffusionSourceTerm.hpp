// Copyright (c) 2020 Maikel Nadolski
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

#ifndef FUB_AMREX_DIFFUSION_SOURCE_TERM_HPP
#define FUB_AMREX_DIFFUSION_SOURCE_TERM_HPP

#include "fub/AMReX/IntegratorContext.hpp"
#include "fub/AMReX/TimeIntegrator.hpp"
#include "fub/AMReX/ViewFArrayBox.hpp"

namespace fub::amrex {

struct DiffusionSourceTermOptions {
  double mul{3.0};
  double reynolds{1.0 / 0.003};
  double schmidt{1.0};
  double prandtl{1.0};
};

template <typename EulerEquation> class DiffusionSourceTerm {
public:
  using Conservative = ::fub::Conservative<EulerEquation>;
  using Complete = ::fub::Complete<EulerEquation>;
  using ConservativeArray = ::fub::ConservativeArray<EulerEquation>;
  using CompleteArray = ::fub::CompleteArray<EulerEquation>;

  DiffusionSourceTerm(
      const EulerEquation& equation,
      const DiffusionSourceTermOptions& options = DiffusionSourceTermOptions());

  Duration ComputeStableDt(int level) const;

  Result<void, TimeStepTooLarge>
  AdvanceLevel(IntegratorContext& simulation_data, int level, Duration dt,
               const ::amrex::IntVect& ngrow = ::amrex::IntVect()) const;

  // ComputeDiffusionFluxes(::amrex::MultiFab& fluxes, const ::amrex::MultiFab&
  // scratch, double dxinv);

  void ComputeDiffusionFlux(ConservativeArray& fluxes,
                            span<const ConservativeArray, 2> states,
                            double dx) const;

private:
  EulerEquation equation_;
  DiffusionSourceTermOptions options_;

  struct Constants {
    Constants(const DiffusionSourceTermOptions& o) {
      mul = o.mul;
      schmidt_invers = 1.0 / o.schmidt;
      reynolds_invers = 1.0 / o.reynolds;
      prandtl_invers = 1.0 / o.prandtl;
      mul_schmidt_inv = schmidt_invers * mul;
      mul_prandtl_inv = prandtl_invers * mul;
      pippo_four_three_mu_effective = 4.0 / 3.0 * mul * reynolds_invers;
      mu_Pr_effective = mul * reynolds_invers * prandtl_invers;
      mu_Sc_effective = mul * reynolds_invers * schmidt_invers;
    }

    double mul;
    double schmidt_invers;
    double reynolds_invers;
    double prandtl_invers;
    double mul_schmidt_inv;
    double mul_prandtl_inv;
    double pippo_four_three_mu_effective;
    double mu_Pr_effective;
    double mu_Sc_effective;
  } constants_{options_};

  Array1d Enthalpy(const ConservativeArray& q, const Array1d& rhoinvers) const
      noexcept;
};

template <typename EulerEquation>
DiffusionSourceTerm<EulerEquation>::DiffusionSourceTerm(
    const EulerEquation& equation, const DiffusionSourceTermOptions& options)
    : equation_(equation), options_(options) {}

template <typename EulerEquation>
Result<void, TimeStepTooLarge> DiffusionSourceTerm<EulerEquation>::AdvanceLevel(
    IntegratorContext& simulation_data, int level, Duration dt,
    const ::amrex::IntVect&) const {
  simulation_data.ApplyBoundaryCondition(level, Direction::X);
  ::amrex::MultiFab& scratch = simulation_data.GetScratch(level);
  const ::amrex::MultiFab& fluxes =
      simulation_data.GetFluxes(level, Direction::X);
  ::amrex::MultiFab fluxes_diffusion(fluxes.boxArray(),
                                     fluxes.DistributionMap(), fluxes.nComp(),
                                     fluxes.nGrowVect());
  std::array<ConservativeArray, 2> stencil{ConservativeArray(equation_),
                                           ConservativeArray(equation_)};
  ConservativeArray flux(equation_);
  const ::amrex::Geometry& geom = simulation_data.GetGeometry(level);
  const double dx = geom.CellSize(0);
  const double dxinv = 1.0 / dx;
  ForEachFab(fluxes_diffusion, [&](const ::amrex::MFIter& mfi) {
    const ::amrex::Box face_tilebox = mfi.growntilebox();
    const ::amrex::Box cell_validbox = scratch[mfi].box();
    const auto [cell_box, face_box] = GetCellsAndFacesInStencilRange(
        cell_validbox, face_tilebox, 1, Direction::X);
    View<const Conservative> patch_data =
        MakeView<const Conservative>(scratch[mfi], equation_, cell_box);
    View<const Conservative> lefts = Shrink(patch_data, Direction::X, {0, 1});
    View<const Conservative> rights = Shrink(patch_data, Direction::X, {1, 0});
    View<Conservative> fluxes_diff =
        MakeView<Conservative>(fluxes_diffusion[mfi], equation_, face_box);
    ForEachRow(std::tuple{fluxes_diff, lefts, rights},
               [&](const Row<Conservative>& row_flux,
                   const Row<const Conservative>& row_left,
                   const Row<const Conservative>& row_right) {
                 ViewPointer inL = Begin(row_left);
                 ViewPointer inR = Begin(row_right);
                 ViewPointer out = Begin(row_flux);
                 ViewPointer last = End(row_left);
                 int n = static_cast<int>(get<0>(last) - get<0>(inL));
                 while (n >= kDefaultChunkSize) {
                   Load(stencil[0], inL);
                   Load(stencil[1], inR);
                   ComputeDiffusionFlux(flux, stencil, dxinv);
                   Store(out, flux);
                   Advance(inL, kDefaultChunkSize);
                   Advance(inR, kDefaultChunkSize);
                   Advance(out, kDefaultChunkSize);
                   n = static_cast<int>(get<0>(last) - get<0>(inL));
                 }
                 LoadN(stencil[0], inL, n);
                 LoadN(stencil[1], inR, n);
                 ComputeDiffusionFlux(flux, stencil, dxinv);
                 StoreN(out, flux, n);
               });
  });
  ForwardIntegrator(execution::simd)
      .UpdateConservatively(scratch, scratch, fluxes_diffusion, geom, dt,
                            Direction::X);
  Reconstruction(execution::simd, equation_).CompleteFromCons(scratch, scratch);
  simulation_data.ApplyBoundaryCondition(level, Direction::X);
  return boost::outcome_v2::success();
}

template <typename EulerEquation>
Duration DiffusionSourceTerm<EulerEquation>::ComputeStableDt([
    [maybe_unused]] int level) const {
  return Duration(std::numeric_limits<double>::max());
}

template <typename EulerEquation>
Array1d DiffusionSourceTerm<EulerEquation>::Enthalpy(
    const ConservativeArray& q, const Array1d& rhoinvers) const noexcept {
  const Array1d rhou = q.momentum.row(0);
  const Array1d rhoe = q.energy;
  const double gamm = equation_.gamma;
  return ((gamm * rhoe - gamm * 0.5 * rhou * rhou * rhoinvers) * rhoinvers);
}

template <typename EulerEquation>
void DiffusionSourceTerm<EulerEquation>::ComputeDiffusionFlux(
    ConservativeArray& flux, span<const ConservativeArray, 2> states,
    double dxinv) const {
  FUB_ASSERT(!states[0].density.isNaN().any());
  FUB_ASSERT(!states[1].density.isNaN().any());
  const Array1d rhoinvL = 1.0 / states[0].density;
  const Array1d rhoinvR = 1.0 / states[1].density;
  const Array1d uL = states[0].momentum.row(0) * rhoinvL;
  const Array1d uR = states[1].momentum.row(0) * rhoinvR;
  const Array1d u = 0.5 * (uL + uR);
  const Array1d du_dx = (uR - uL) * dxinv;
  const Array1d hL = Enthalpy(states[0], rhoinvL);
  const Array1d hR = Enthalpy(states[1], rhoinvR);
  const Array1d dh_dx = (hR - hL) * dxinv;

  flux.momentum.row(0) = -constants_.pippo_four_three_mu_effective * du_dx;
  flux.energy = -(constants_.pippo_four_three_mu_effective * u * du_dx +
                  constants_.mu_Pr_effective * dh_dx);
  for (int s = 0; s < states[0].species.rows(); ++s) {
    const Array1d YL = states[0].species.row(s) * rhoinvL;
    const Array1d YR = states[1].species.row(s) * rhoinvR;
    const Array1d dY_dx = (YR - YL) * dxinv;
    flux.species.row(s) = -constants_.mu_Sc_effective * dY_dx;
  }
}

} // namespace fub::amrex

#endif