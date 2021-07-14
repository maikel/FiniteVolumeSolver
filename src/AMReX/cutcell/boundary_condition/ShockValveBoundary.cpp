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

#include "fub/AMReX/cutcell/boundary_condition/ShockValveBoundary.hpp"
#include "fub/AMReX/ForEachIndex.hpp"

#include "fub/ext/ProgramOptions.hpp"

#include <utility>

#include <boost/log/common.hpp>
#include <boost/log/sources/severity_logger.hpp>
#include <boost/log/trivial.hpp>

namespace fub::amrex::cutcell {

void ShockValveOptions::Print(SeverityLogger& log) {
  BOOST_LOG(log) << fmt::format("Shock Valve '{}' Options:", channel);
  BOOST_LOG(log) << fmt::format(" --- Massflow part:");
  massflow_boundary.Print(log);
  BOOST_LOG(log) << fmt::format(" --- ShockFeedback part:");
  shock_options.Print(log);
}

ShockValveOptions::ShockValveOptions(const ProgramOptions& options) {
  FUB_GET_OPTION_VAR(options, channel);
  massflow_boundary = GetOptions(options, "massflow_boundary");
  shock_options = GetOptions(options, "schock_feedback");
}

ShockValveBoundary::ShockValveBoundary(const PerfectGas<2>& equation,
                                       ShockValveOptions options)
    : options_{std::move(options)}, equation_{equation} {}

const ShockValveOptions& ShockValveBoundary::GetOptions() const noexcept {
  return options_;
}

const ShockValve& ShockValveBoundary::GetValve() const noexcept {
  return valve_;
}

namespace {

void ChangeState_(ShockValveState& state, const ::amrex::Geometry&,
                  const GriddingAlgorithm& grid, int level,
                  Duration& /* last_shock_change */,
                  const ShockValveOptions& options, PerfectGas<2>& /* eq */) {
  const Duration current_time = grid.GetPatchHierarchy().GetTimePoint(0);
  boost::log::sources::severity_logger<boost::log::trivial::severity_level> log(
      boost::log::keywords::severity = boost::log::trivial::info);
  BOOST_LOG_SCOPED_LOGGER_TAG(log, "Channel", options.channel);
  BOOST_LOG_SCOPED_LOGGER_TAG(log, "Time", current_time.count());
  BOOST_LOG_SCOPED_LOGGER_TAG(log, "Level", level);
  switch (state) {
  case ShockValveState::closed:
    BOOST_LOG(log) << "Shock valve is transmissive!";
    break;
  case ShockValveState::open:
    if (options.shock_options.shock_time < current_time) {
      state = ShockValveState::closed;
      BOOST_LOG(log) << "Shock valve has changed to transmissive!";
    }
  default:
    break;
  }
}

} // namespace

void ShockValveBoundary::FillBoundary(::amrex::MultiFab& mf,
                                      const GriddingAlgorithm& grid, int level,
                                      Direction dir) {
  if (dir == Direction::X) {
    FillBoundary(mf, grid, level);
  }
}

void ShockValveBoundary::FillBoundary(::amrex::MultiFab& mf,
                                      const GriddingAlgorithm& grid,
                                      int level) {
  const ::amrex::Geometry& geom = grid.GetPatchHierarchy().GetGeometry(level);
  const int ngrow = mf.nGrow(0);
  ::amrex::Box grown_box = geom.growNonPeriodicDomain(ngrow);
  ::amrex::BoxList boundaries =
      ::amrex::complementIn(grown_box, ::amrex::BoxList{geom.Domain()});
  if (boundaries.isEmpty()) {
    return;
  }

  // Change State Machine if neccessary
  ChangeState_(valve_.state, geom, grid, level, valve_.last_shock, options_,
               equation_);

  // ReflectiveBoundary closed(execution::seq, equation_, Direction::X, 0);
  TransmissiveBoundary closed{
      Direction::X, 0}; // TODO Test if cutcell implementation is needed!!!
  MassflowBoundary_PerfectGas inflow_boundary(equation_,
                                              options_.massflow_boundary);
  Complete<PerfectGas<2>> state{equation_};
  switch (valve_.state) {
  case ShockValveState::closed:
    closed.FillBoundary(mf, geom);
    break;
  case ShockValveState::open:
    inflow_boundary.FillBoundary(mf, grid, level);
  }
}
} // namespace fub::amrex::cutcell
