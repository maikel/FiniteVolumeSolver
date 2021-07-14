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

#ifndef FUB_AMREX_BOUNDARY_CONDITION_CUTCELL_SHOCK_VALVE_HPP
#define FUB_AMREX_BOUNDARY_CONDITION_CUTCELL_SHOCK_VALVE_HPP

/// \file
///
/// This boundary condition simulates a valve that admits a specific mass flow
/// when open and models a transmissive boundary when closed.


#include "fub/AMReX/ForEachFab.hpp"
#include "fub/AMReX/GriddingAlgorithm.hpp"

#include "fub/AMReX/boundary_condition/TransmissiveBoundary.hpp"
#include "fub/AMReX/cutcell/boundary_condition/MassflowBoundary_PerfectGas.hpp"
#include "fub/equations/perfect_gas/InitializeShock.hpp"

#include "fub/Duration.hpp"

#include <boost/serialization/access.hpp>

#include "fub/ext/Log.hpp"
#include "fub/ext/ProgramOptions.hpp"

#include <limits>
#include <string>

namespace fub::amrex::cutcell {
/// \ingroup BoundaryCondition
///
struct ShockValveOptions {
  ShockValveOptions() = default;
  ShockValveOptions(const ProgramOptions& vm);

  void Print(SeverityLogger& log);

  std::string channel{"ShockValve"};
  MassflowBoundary_PerfectGasOptions massflow_boundary{};
  feedback_functions::ShockOptions shock_options{};
};

enum class ShockValveState { open, closed };
/// \ingroup BoundaryCondition
///
struct ShockValve {
  ShockValveState state{ShockValveState::open};
  Duration last_shock{std::numeric_limits<double>::lowest()};
};
} // namespace fub::amrex::cutcell

namespace boost::serialization {
template <typename Archive>
void serialize(Archive& ar, ::fub::amrex::cutcell::ShockValve& valve,
               unsigned int /* version */) {
  // clang-format off
  int state = static_cast<int>(valve.state);
  ar & state;
  valve.state = static_cast<::fub::amrex::cutcell::ShockValveState>(state);
  ar & valve.last_shock;
  // clang-format on
}
} // namespace boost::serialization

namespace fub::amrex::cutcell {
/// \ingroup BoundaryCondition
///
class ShockValveBoundary {
public:
  ShockValveBoundary(const PerfectGas<2>& equation, ShockValveOptions options);

  [[nodiscard]] const ShockValveOptions& GetOptions() const noexcept;

  [[nodiscard]] const ShockValve& GetValve() const noexcept;

  void FillBoundary(::amrex::MultiFab& mf, const GriddingAlgorithm& gridding,
                    int level);

  void FillBoundary(::amrex::MultiFab& mf, const GriddingAlgorithm& gridding,
                    int level, Direction dir);

private:
  template <typename Archive>
  friend void serialize(Archive& ar, ShockValveBoundary& boundary,
                        unsigned int /* version */) {
    // clang-format off
    ar & boundary.valve_;
    // clang-format on
  }

  ShockValveOptions options_;
  PerfectGas<2> equation_;
  ShockValve valve_{};
};
} // namespace fub::amrex::cutcell

#endif // FINITEVOLUMESOLVER_CUTCELL_SHOCK_VALVE_HPP
