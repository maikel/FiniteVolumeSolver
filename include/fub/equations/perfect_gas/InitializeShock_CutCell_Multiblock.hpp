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

#ifndef FUB_EQUATIONS_PERFECT_GAS_INITIALIZE_SHOCK_CUTCELL_MULTIBLOCK_HPP
#define FUB_EQUATIONS_PERFECT_GAS_INITIALIZE_SHOCK_CUTCELL_MULTIBLOCK_HPP

/// \file
///
/// This file contains a Feedback function that can be applied to the cutcell
/// Integratorcontext. It applies a Riemann shock state from a given Mach Number
/// and an average state. So far it is only implemented for 2D Perfect Gas
/// equations.

#include "fub/equations/PerfectGas.hpp"

#include "fub/AMReX/cutcell/AverageState.hpp"
#include "fub/AMReX/multi_block/MultiBlockIntegratorContext2.hpp"

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
  ShockFeedback(const PerfectGas<1>& tube_equation, const PerfectGas<Rank>& plena_equation, const ShockOptions& opts)
      : tube_equation_(tube_equation), plena_equation_(plena_equation), options_(opts) {}

  void ComputePreShockState(PerfectGas<Rank>& equation,
                            Complete<PerfectGas<Rank>>& pre_shock,
                            const Complete<PerfectGas<Rank>>& post_shock,
                            double M_S,
                            const fub::Array<double, Rank, 1>& normal) {
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

  void operator()(amrex::MultiBlockIntegratorContext2& context, Duration) {
    
    static constexpr int Tube_Rank = 1;
    static constexpr int Plenum_Rank = 2;
    FUB_ASSERT(AMREX_SPACEDIM == 2);
    
    int coarsest_level = 0;

    const Duration timepoint = context.GetTimePoint(coarsest_level);

    if ((timepoint < options_.shock_time) || options_.shock_was_applied) {
      return;
    }

    auto tubes = context.Tubes();
    auto plena = context.Plena();

    // all the following stuff is only implemented for one tube and plena !!
    FUB_ASSERT(tubes.size() == 1);
    FUB_ASSERT(plena.size() == 1);

    // get all relevant data
    fub::amrex::cutcell::PatchHierarchy plena_hier = plena[0].GetPatchHierarchy();
    ::amrex::MultiFab& plena_scratch =  plena[0].GetScratch(coarsest_level);
    ::amrex::MultiFab& tube_scratch = tubes[0].GetScratch(coarsest_level);

    // get average state right from the shock location == post shock state
    Complete<PerfectGas<Plenum_Rank>> post_shock_state(plena_equation_);
    AverageState(post_shock_state, plena_hier, coarsest_level,
                 options_.average_post_shock_box);

    // compute pre shock state from post shock state
    Complete<PerfectGas<Plenum_Rank>> pre_shock_state(plena_equation_);
    const fub::Array<double, Plenum_Rank, 1> normal{
        AMREX_D_DECL(1.0, 0.0, 0.0)};

    ComputePreShockState(plena_equation_, pre_shock_state, post_shock_state,
                         options_.shock_mach_number, normal);

    fub::SeverityLogger log = fub::GetInfoLogger();
    BOOST_LOG(log) << "Shock was applied at time " << timepoint.count() << "\n";
    BOOST_LOG(log) << "Post-Shock-State:\n"
                   << "\tdensity: " << post_shock_state.density << " [kg/m^3]\n"
                   << "\tvelocity: "
                   << plena_equation_.Velocity(post_shock_state).transpose()
                   << " [m/s]\n"
                   << "\tpressure: " << post_shock_state.pressure << " [Pa]";

    BOOST_LOG(log) << "Calculated Pre-Shock-State:\n"
                   << "\tdensity: " << pre_shock_state.density << " [kg/m^3]\n"
                   << "\tvelocity: "
                   << plena_equation_.Velocity(pre_shock_state).transpose()
                   << " [m/s]\n"
                   << "\tpressure: " << pre_shock_state.pressure << " [Pa]";

    Complete<PerfectGas<Tube_Rank>> tube_pre_shock_state(tube_equation_);
    Eigen::Array<double, Tube_Rank, 1> tube_velocity = Eigen::Array<double, Tube_Rank, 1>::Zero();
    tube_velocity[0] = plena_equation_.Velocity(pre_shock_state)[0];
    tube_pre_shock_state = tube_equation_.CompleteFromPrim(pre_shock_state.density, tube_velocity, pre_shock_state.pressure);


    // write pre shock state in all cells left from the shock location in the plenum
    const ::amrex::IntVect bigEnd = options_.average_post_shock_box.bigEnd();
    const int scratch_gcw = plena[0].GetOptions().scratch_gcw;
    const ::amrex::IntVect smallEnd{-scratch_gcw, -scratch_gcw};
    const ::amrex::Box pre_shock_box{smallEnd, bigEnd};
    const ::amrex::MultiFab& alphas =
        plena_hier.GetEmbeddedBoundary(coarsest_level)->getVolFrac();

    // Complete<PerfectGas<Plenum_Rank>> state{plena_equation_};
    ForEachFab(execution::seq, plena_scratch, [&](const ::amrex::MFIter& mfi) {
      ::amrex::FArrayBox& fab = plena_scratch[mfi];
      const ::amrex::FArrayBox& alpha = alphas[mfi];
      ::amrex::Box box_to_fill = mfi.growntilebox() & pre_shock_box;
      if (!box_to_fill.isEmpty()) {
        auto states = MakeView<Complete<PerfectGas<Plenum_Rank>>>(fab, plena_equation_,
                                                          mfi.growntilebox());
        ForEachIndex(AsIndexBox<AMREX_SPACEDIM>(box_to_fill), [&](auto... is) {
          Index<AMREX_SPACEDIM> dest{is...};
          ::amrex::IntVect iv{
              AMREX_D_DECL(int(dest[0]), int(dest[1]), int(dest[2]))};
          if (alpha(iv) > 0.0) {
            // Load(state, states, dest);
            // old state doesn't matter, only write new state
            Store(states, pre_shock_state, dest);
          }
        });
      }
    });

    // write pre shock state in all cells left from the tube
    // const int tube_scratch_gcw = tubes[0].GetOptions().scratch_gcw;
    const ::amrex::Box& tube_pre_shock_box = tubes[0].GetGeometry(coarsest_level).Domain();
    // const ::amrex::IntVect tube_bigEnd{0};
    // const ::amrex::IntVect tube_smallEnd{-tube_scratch_gcw};
    // const ::amrex::Box tube_pre_shock_box{tube_smallEnd, tube_bigEnd};

    ForEachFab(execution::seq, tube_scratch, [&](const ::amrex::MFIter& mfi) {
      ::amrex::FArrayBox& tube_fab = tube_scratch[mfi];
      ::amrex::Box tube_box_to_fill = mfi.growntilebox() & tube_pre_shock_box;
      if (!tube_box_to_fill.isEmpty()) {
        auto tube_states = MakeView<Complete<PerfectGas<Tube_Rank>>>(tube_fab, tube_equation_,
                                                          mfi.growntilebox());
        ForEachIndex(AsIndexBox<Tube_Rank>(tube_box_to_fill), [&](auto... is) {
          Index<Tube_Rank> tube_dest{is...};
          // Load(state, states, dest);
          // old state doesn't matter, only write new state
          Store(tube_states, tube_pre_shock_state, tube_dest);
        });
      }
    });

    options_.shock_was_applied = true;
  }

private:
  PerfectGas<1> tube_equation_;
  PerfectGas<Rank> plena_equation_;
  ShockOptions options_;
};

} // namespace feedback_functions

} // namespace fub::amrex::cutcell

#endif // FUB_EQUATIONS_PERFECT_GAS_INITIALIZE_SHOCK_CUTCELL_MULTIBLOCK_HPP