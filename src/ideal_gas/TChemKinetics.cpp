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

#include "fub/ideal_gas/TChemKinetics.hpp"
#include "src/solver/ode_solver/Radau.hpp"

extern "C" {
#include "TC_defs.h"
#include "TC_interface.h"
}

#include <algorithm>
#include <numeric>

namespace fub {
namespace ideal_gas {
namespace {
int reference_counter_ = 0;

int InitTCAndGetNSpecies_(const TChemKineticsOptions& options) {
  if (reference_counter_ == 0) {
    TC_initChem(options.chemfile.c_str(), options.thermofile.c_str(), 1, 0.5);
  }
  return TC_getNspec();
}

double Dot_(span<const double> x, span<const double> y) {
  FUB_ASSERT(x.size() == y.size());
  int size = x.size();
  double total = 0.0;
  for (int i = 0; i < size; ++i) {
    total += x[i] * y[i];
  }
  return total;
}

double GetUSpecMs_(span<const double> TandY, span<double> ui) {
  TC_getUspecMs(TandY[0], ui.size(), ui.data());
  span<const double> Y = TandY.subspan(1);
  return Dot_(Y, ui);
}

void Normalize_(span<double> x) {
  const double total = std::accumulate(x.begin(), x.end(), 0.0);
  std::transform(x.begin(), x.end(), x.begin(),
                 [total](double xi) { return xi / total; });
}

double GetTemperatureFromInternalEnergy_(double target, span<double> TandY,
                                         span<double> uis, double dTtol) {
  double dT = 0.0;
  double Tnew = 300.0;
  double Cvnew = 0.0;
  TandY[0] = Tnew;
  double Unew = GetUSpecMs_(TandY, uis);
  TC_getMs2CvMixMs(TandY.data(), TandY.size(), &Cvnew);

  double Utop = Unew;
  double Ubot = Unew;
  double Uold = Unew;
  double Ttop = Tnew;
  double Tbot = Tnew;
  double Told = Tnew;

  bool unstablePhase = false;
  double Tunstable = -1;

  double UConvErr;

  // Newton iteration
  // This is exactly like Cantera's setState_UV implementation,
  // except that we do not check for upper/lower temperature limits
  for (int i = 0; i < 1000; i++) {
    Told = Tnew;
    Uold = Unew;

    double cvd = Cvnew;
    if (cvd < 0) {
      unstablePhase = true;
      Tunstable = Tnew;
    }

    dT = std::max(-100., std::min(100., (target - Uold) / cvd));
    Tnew = Told + dT;

    // This limits the step size to make the algorithm convergent
    // See Cantera for details
    if ((dT > 0 && unstablePhase) || (dT <= 0 && !unstablePhase)) {
      if (Ubot < target && Tnew < (0.75 * Tbot + 0.25 * Told)) {
        dT = .75 * (Tbot - Told);
      }
    } else if (Utop > target && Tnew > (.75 * Ttop + .25 * Told)) {
      dT = .75 * (Ttop - Told);
      Tnew = Told + dT;
    }

    // Set the new temperature, but try to stay in the stable region
    // with cv > 0
    for (int its = 0; its < 10; its++) {
      Tnew = Told + dT;
      TandY[0] = Tnew;
      Unew = GetUSpecMs_(TandY, uis);
      TC_getMs2CvMixMs(TandY.data(), TandY.size(), &Cvnew);

      if (Cvnew < 0) {
        Tunstable = Tnew;
        dT *= .25;
      } else {
        break;
      }
    }

    if (Unew == target) {
      return Tnew;
    } else if (Unew > target && (Utop < target || Unew < Utop)) {
      Utop = Unew;
      Ttop = Tnew;
    } else if (Unew < target && (Ubot > target || Unew > Ubot)) {
      Ubot = Unew;
      Tbot = Tnew;
    }

    // Check for convergence
    double Uerr = target - Unew;
    double acvd = std::max(std::abs(cvd), 1e-5);
    double denom = std::max(std::abs(target), acvd * dTtol);
    UConvErr = std::abs(Uerr / denom);
    if (UConvErr < 1e-5 * dTtol || std::abs(dT) < dTtol) {
      return Tnew;
    }
  }
  throw std::runtime_error(
      "TChemKinetics: No Convergence in GetTemperatureFromInternalEnergy");
}
} // namespace

TChemKinetics::TChemKinetics(std::string prefix, SAMRAI::tbox::Dimension dim,
                             TChemKineticsOptions options) noexcept
    : IdealGasKinetics(prefix, dim, InitTCAndGetNSpecies_(options)) {}

TChemKinetics::~TChemKinetics() noexcept {
  --reference_counter_;
  if (reference_counter_ == 0) {
    TC_reset();
  }
}

void TChemKinetics::FillFromCons(const CompleteStatePatchData& q,
                                 const ConsStatePatchData& cons) const {
  q.density.copy(cons.density);
  q.momentum.copy(cons.momentum);
  q.energy.copy(cons.energy);
  q.species.copy(cons.species);
  SAMRAI::hier::Box intersection(q.density.getDim());
  q.density.getGhostBox().intersect(cons.density.getGhostBox(), intersection);
  std::vector<double> TandY(q.species.getDepth() + 1);
  std::vector<double> buffer(q.species.getDepth());
  span<double> Y = make_span(TandY).subspan(1);
  for (const SAMRAI::hier::Index& index : intersection) {
    SAMRAI::pdat::CellIndex cell(index);
    const double rho = q.density(cell);
    CopyMassFractions(Y, q.species, cell);
    Normalize_(Y);
    const double rhou = q.momentum(cell);
    const double rhoE = q.energy(cell);
    TC_setDens(rho);
    const double U = (rhoE - 0.5 * rhou * rhou / rho) / rho;
    const double T = GetTemperatureFromInternalEnergy_(U, TandY, buffer, 1e-5);
    q.temperature(cell) = T;
    double cp = 0.0;
    TandY[0] = T;
    TC_getMs2CpMixMs(TandY.data(), TandY.size(), &cp);
    double cv = 0.0;
    TC_getMs2CvMixMs(TandY.data(), TandY.size(), &cv);
    const double R = cp - cv;
    const double p = rho * R * T;
    q.pressure(cell) = p;
    const double gamma = cp / cv;
    q.speed_of_sound(cell) = std::sqrt(gamma * p / rho);
  }
}

void TChemKinetics::FillFromPrim(const CompleteStatePatchData& q,
                                 const PrimStatePatchData& prim) const {
  q.temperature.copy(prim.temperature);
  q.momentum.copy(prim.momentum);
  q.pressure.copy(prim.pressure);
  q.species.copy(prim.species);
  SAMRAI::hier::Box intersection(q.density.getDim());
  q.density.getGhostBox().intersect(prim.pressure.getGhostBox(), intersection);
  std::vector<double> TandY(q.species.getDepth() + 1);
  std::vector<double> buffer(q.species.getDepth());
  span<double> Y = make_span(TandY).subspan(1);
  for (const SAMRAI::hier::Index& index : intersection) {
    SAMRAI::pdat::CellIndex cell(index);
    const double p = q.pressure(cell);
    TC_setThermoPres(p);
    CopyMassFractions(Y, q.species, cell);
    Normalize_(Y);
    TC_getMl2Ms(Y.data(), Y.size(), Y.data());
    TandY[0] = q.temperature(cell);
    double rho;
    TC_getRhoMixMs(TandY.data(), TandY.size(), &rho);
    q.density(cell) = rho;
    const double rhou = q.momentum(cell);
    // internal energy = (total energy - kinetic energy) / density
    const double eps = GetUSpecMs_(TandY, buffer);
    const double rhoE = rho * eps + 0.5 * rhou * rhou / rho;
    q.energy(cell) = rhoE;
    double cp = 0.0;
    double cv = 0.0;
    TC_getMs2CpMixMs(TandY.data(), TandY.size(), &cp);
    TC_getMs2CvMixMs(TandY.data(), TandY.size(), &cv);
    const double gamma = cp / cv;
    q.speed_of_sound(cell) = std::sqrt(gamma * p / rho);
    for (int s = 0; s < Y.size(); ++s) {
      q.species(cell, s) = rho * Y[s];
    }
  }
}

void TChemKinetics::AdvanceSourceTerm(span<double> TandY, double rho,
                                      double time_step_size) const {
  auto odeRhs = [](span<double> dZdt, span<const double> Z, double time_point) {
    TC_getSrcCV(const_cast<double*>(Z.data()), Z.size(), dZdt.data());
  };
  auto odeJac = [](span<double> dfdZ, span<const double> Z, double time_point) {
    TC_getJacCVTYNanl(const_cast<double*>(Z.data()), Z.size() - 1, dfdZ.data());
  };
  TC_setDens(rho);
  Radau::integrate_jacobian(odeRhs, TandY, 0.0, time_step_size, odeJac);
  // Radau::integrate(odeRhs, TandY, 0.0, time_step_size);
}

void TChemKinetics::AdvanceSourceTerm(
    span<double> TandY, double rho, double time_step_size,
    function_ref<int(span<const double>, double)> feedback) const {
  auto odeRhs = [](span<double> dZdt, span<const double> Z, double time_point) {
    TC_getSrcCV(const_cast<double*>(Z.data()), Z.size(), dZdt.data());
  };
  auto odeJac = [](span<double> dfdZ, span<const double> Z, double time_point) {
    TC_getJacCVTYNanl(const_cast<double*>(Z.data()), Z.size() - 1, dfdZ.data());
  };
  TC_setDens(rho);
  Radau::integrate_jacobian_feedback(odeRhs, TandY, 0.0, time_step_size, odeJac,
                                     feedback);
}

void TChemKinetics::AdvanceSourceTerm(const CompleteStatePatchData& q,
                                      double time_step_size) const {
  const SAMRAI::hier::Box& box = q.density.getGhostBox();
  std::vector<double> TandY(q.species.getDepth() + 1);
  std::vector<double> ui(q.species.getDepth());
  span<double> Y = make_span(TandY).subspan(1);
  for (const SAMRAI::hier::Index& index : box) {
    SAMRAI::pdat::CellIndex cell(index);
    TC_setThermoPres(q.pressure(cell));
    CopyMassFractions(Y, q.species, cell);
    Normalize_(Y);
    TandY[0] = q.temperature(cell);
    AdvanceSourceTerm(TandY, q.density(cell), time_step_size);
    const double u = q.momentum(cell) / q.density(cell);
    double rho;
    TC_getRhoMixMs(TandY.data(), TandY.size(), &rho);
    q.density(cell) = rho;
    q.temperature(cell) = TandY[0];
    const double p = TC_getThermoPres();
    q.pressure(cell) = p;
    q.energy(cell) = rho * GetUSpecMs_(TandY, ui) + 0.5 * rho * u * u;
    q.momentum(cell) = rho * u;
    double cp = 0.0;
    double cv = 0.0;
    TC_getMs2CpMixMs(TandY.data(), TandY.size(), &cp);
    TC_getMs2CvMixMs(TandY.data(), TandY.size(), &cv);
    const double gamma = cp / cv;
    q.speed_of_sound(cell) = std::sqrt(gamma * p / rho);
    for (int s = 0; s < q.species.getDepth(); ++s) {
      q.species(cell, s) = rho * Y[s];
    }
  }
}

} // namespace ideal_gas
} // namespace fub