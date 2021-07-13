// Copyright (c) 2021 Christian Zenker
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

#ifndef FUB_EQUATIONS_PERFECT_GAS_INITIALIZE_SHOCK_HPP
#define FUB_EQUATIONS_PERFECT_GAS_INITIALIZE_SHOCK_HPP

/// \file
///
/// This file contains a Feedback function that can be applied to the cutcell
/// Integratorcontext. It applies a Riemann shock state from a given Mach Number
/// and an average state. So far it is only implemented for 2D Perfect Gas
/// equations.

#include "fub/equations/PerfectGas.hpp"

#include "fub/AMReX/cutcell/AverageState.hpp"
#include "fub/AMReX/cutcell/IntegratorContext.hpp"

#include "fub/ext/Log.hpp"
#include "fub/ext/ProgramOptions.hpp"

namespace fub::amrex::cutcell {

namespace feedback_functions {
/// \namespace
/// \brief This namespace includes all classes and functions that can be passed
/// as feedback functions.

struct ShockOptions {
  ShockOptions() = default;
  ShockOptions(const ProgramOptions& options);

  void Print(SeverityLogger& log);

  double shock_mach_number{1.0}; ///< Shock Machnumber
  // double shock_x_location{0.0};  ///< Shock location
  ::amrex::Box
      average_post_shock_box{}; ///< average box for the post shock state
  Duration shock_time{10.0};    ///< time at which the shock should be applied
  // Direction dir = Direction::X;

  // problem when simulation is started from checkpoint after shock was applied!
  bool shock_was_applied{false};
};

template <int Rank> class ShockFeedback {
public:
  ShockFeedback(const PerfectGas<Rank>& equation, const ShockOptions& opts)
      : equation_(equation), options_(opts) {}

  void ComputePreShockState(PerfectGas<2>& equation,
                            Complete<PerfectGas<2>>& pre_shock,
                            const Complete<PerfectGas<2>>& post_shock,
                            double M_S,
                            const fub::Array<double, 2, 1>& normal) {
    const double g = equation.gamma;
    const double gp1 = g + 1.0;
    const double gm1 = g - 1.0;
    const double M_post = equation.Machnumber(post_shock);
    const double M2 = (M_post - M_S) * (M_post - M_S);

    const double rho_c = gp1 * M2 / (gm1 * M2 + 2.0);
    const double rho_post = post_shock.density;
    const double rho_pre = rho_post * rho_c;

    const double p_c = (2.0 * g * M2 - gm1) / gp1;
    const double p_post = post_shock.pressure;
    const double p_pre = p_post * p_c;

    const double t = rho_post / rho_pre;
    const double a_post = post_shock.speed_of_sound;
    const double u_S = M_S * a_post;
    const double u_post = equation.Velocity(post_shock).matrix().norm();
    const double u_pre = u_S * (1.0 - t) + u_post * t;

    pre_shock = equation.CompleteFromPrim(rho_pre, u_pre * normal, p_pre);
  }

  void operator()(amrex::cutcell::IntegratorContext& context, int, Duration,
                  std::pair<int, int>) {

    constexpr int Dim = 2;
    FUB_ASSERT(Rank == Dim);
    FUB_ASSERT(AMREX_SPACEDIM == 2);

    int coarsest_level = 0;
    fub::amrex::cutcell::PatchHierarchy hier = context.GetPatchHierarchy();
    const Duration timepoint = hier.GetTimePoint(coarsest_level);

    if ((timepoint < options_.shock_time) || options_.shock_was_applied) {
      return;
    }

    // get all relevant data
    ::amrex::MultiFab& scratch = context.GetScratch(coarsest_level);
    const ::amrex::Geometry& geom = context.GetGeometry(coarsest_level);

    // get average state right from the shock location == post shock state
    Complete<PerfectGas<Dim>> post_shock_state(equation_);
    AverageState(post_shock_state, hier, coarsest_level,
                 options_.average_post_shock_box);

    // compute pre shock state from post shock state
    Complete<PerfectGas<Dim>> pre_shock_state(equation_);
    const fub::Array<double, AMREX_SPACEDIM, 1> normal{
        AMREX_D_DECL(1.0, 0.0, 0.0)};

    ComputePreShockState(equation_, pre_shock_state, post_shock_state,
                         options_.shock_mach_number, normal);

    fub::SeverityLogger log = fub::GetInfoLogger();
    BOOST_LOG(log) << "Shock was applied at time " << timepoint.count() << "\n";
    BOOST_LOG(log) << "Post-Shock-State:\n"
                   << "\tdensity: " << post_shock_state.density << " [kg/m^3]\n"
                   << "\tvelocity: "
                   << equation_.Velocity(post_shock_state).transpose()
                   << " [m/s]\n"
                   << "\tpressure: " << post_shock_state.pressure << " [Pa]";

    BOOST_LOG(log) << "Calculated Pre-Shock-State:\n"
                   << "\tdensity: " << pre_shock_state.density << " [kg/m^3]\n"
                   << "\tvelocity: "
                   << equation_.Velocity(pre_shock_state).transpose()
                   << " [m/s]\n"
                   << "\tpressure: " << pre_shock_state.pressure << " [Pa]";

    // write pre shock state in all cells left from the shock location
    const ::amrex::IntVect bigEnd = options_.average_post_shock_box.bigEnd();
    const ::amrex::IntVect smallEnd = geom.Domain().smallEnd();
    const ::amrex::Box pre_shock_box{smallEnd, bigEnd};
    const ::amrex::MultiFab& alphas =
        hier.GetEmbeddedBoundary(coarsest_level)->getVolFrac();

    Complete<PerfectGas<Dim>> state{equation_};
    ForEachFab(execution::seq, scratch, [&](const ::amrex::MFIter& mfi) {
      ::amrex::FArrayBox& fab = scratch[mfi];
      const ::amrex::FArrayBox& alpha = alphas[mfi];
      ::amrex::Box box_to_fill = mfi.growntilebox() & pre_shock_box;
      if (!box_to_fill.isEmpty()) {
        auto states = MakeView<Complete<PerfectGas<Dim>>>(fab, equation_,
                                                          mfi.growntilebox());
        ForEachIndex(AsIndexBox<AMREX_SPACEDIM>(box_to_fill), [&](auto... is) {
          Index<AMREX_SPACEDIM> dest{is...};
          ::amrex::IntVect iv{
              AMREX_D_DECL(int(dest[0]), int(dest[1]), int(dest[2]))};
          if (alpha(iv) > 0.0) {
            Load(state, states, dest);
            Store(states, pre_shock_state, dest);
          }
        });
      }
    });
    options_.shock_was_applied = true;
  }

private:
  PerfectGas<Rank> equation_;
  ShockOptions options_;
};

} // namespace feedback_functions

} // namespace fub::amrex::cutcell

#endif // FUB_EQUATIONS_PERFECT_GAS_INITIALIZE_SHOCK_HPP