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

#ifndef FUB_AMREX_BOUNDARY_CONDITION_GENERIC_PRESSURE_VALVE_HPP
#define FUB_AMREX_BOUNDARY_CONDITION_GENERIC_PRESSURE_VALVE_HPP

#include "fub/AMReX/GriddingAlgorithm.hpp"
#include "fub/AMReX/MultiFabUtilities.hpp"
#include "fub/AMReX/boundary_condition/ConstantBoundary.hpp"
#include "fub/AMReX/boundary_condition/ReflectiveBoundary.hpp"

#include "fub/ext/Log.hpp"

namespace fub::amrex {

struct GenericPressureValveBoundaryOptions {
  std::string prefix{"GenericPressureValve"};
  double forward_efficiency{1.0};
  double backward_efficiency{0.0};
  double open_below_pressure{1.0};
  Direction dir{Direction::X};
  int side{0};
};

template <typename EulerEquation, typename InflowFunction>
class GenericPressureValveBoundary {
public:
  GenericPressureValveBoundary(
      const EulerEquation& equation, InflowFunction fn,
      const GenericPressureValveBoundaryOptions& options);

  void FillBoundary(::amrex::MultiFab& mf, const GriddingAlgorithm& gridding,
                    int level);

  void FillBoundary(::amrex::MultiFab& mf, const GriddingAlgorithm& gridding,
                    int level, Direction dir);

  std::optional<Duration> GetTimePointWhenOpened() const noexcept {
    return t_opened_;
  }

private:
  EulerEquation equation_;
  InflowFunction inflow_function_;
  GenericPressureValveBoundaryOptions options_;
  ConstantBoundary<EulerEquation> constant_boundary_;
  ReflectiveBoundary<execution::SequentialTag, EulerEquation>
      reflective_boundary_;

  IndexMapping<EulerEquation> comps_{equation_};
  std::optional<Duration> t_opened_{};
};

template <typename EulerEquation, typename InflowFunction>
GenericPressureValveBoundary<EulerEquation, InflowFunction>::
    GenericPressureValveBoundary(
        const EulerEquation& equation, InflowFunction fn,
        const GenericPressureValveBoundaryOptions& options)
    : equation_(equation), inflow_function_(std::move(fn)),
      options_(options), constant_boundary_{options_.dir, options_.side,
                                            equation_,
                                            Complete<EulerEquation>(equation_)},
      reflective_boundary_(equation_, options_.dir, options_.side) {}

template <typename EulerEquation, typename InflowFunction>
void GenericPressureValveBoundary<EulerEquation, InflowFunction>::FillBoundary(
    ::amrex::MultiFab& mf, const GriddingAlgorithm& gridding, int level,
    Direction dir) {
  if (dir == options_.dir) {
    FillBoundary(mf, gridding, level);
  }
}

template <typename EulerEquation, typename InflowFunction>
void GenericPressureValveBoundary<EulerEquation, InflowFunction>::FillBoundary(
    ::amrex::MultiFab& mf, const GriddingAlgorithm& gridding, int level) {
  const ::amrex::Geometry& geom =
      gridding.GetPatchHierarchy().GetGeometry(level);
  const ::amrex::Box domain_box = geom.Domain();
  const ::amrex::Box inner_box =
      GetInnerBox(domain_box, options_.side, options_.dir, 1);
  const double inner_pressure =
      GetMeanValueInBox(mf, inner_box, comps_.pressure);
  const Duration t = gridding.GetTimePoint();
  if (inner_pressure <= options_.open_below_pressure) {
    if (!t_opened_) {
      t_opened_ = t;
      SeverityLogger log = GetInfoLogger();
      BOOST_LOG_SCOPED_LOGGER_TAG(log, "Channel", options_.prefix);
      BOOST_LOG_SCOPED_LOGGER_TAG(log, "Time", t);
      BOOST_LOG_SCOPED_LOGGER_TAG(log, "Level", level);
      BOOST_LOG(log) << fmt::format("Pressure valve has opened because inner pressure is {}.", inner_pressure);
    }
    const Duration t_diff = t - *t_opened_;
    KineticState<EulerEquation> kinetic_state(equation_);
    std::invoke(inflow_function_, equation_, kinetic_state, t_diff, mf,
                gridding, level);
    constexpr int N = EulerEquation::Rank();
    Array<double, N, 1> zero = Array<double, N, 1>::Zero();
    euler::CompleteFromKineticState(equation_, constant_boundary_.state,
                                    kinetic_state, zero);
    euler::IsentropicExpansionWithoutDissipation(
        equation_, constant_boundary_.state, constant_boundary_.state,
        inner_pressure, options_.forward_efficiency);
    constant_boundary_.FillBoundary(mf, gridding, level);
  } else {
    if (t_opened_ && inner_pressure > 1.1 * options_.open_below_pressure) {
      t_opened_.reset();
      SeverityLogger log = GetInfoLogger();
      BOOST_LOG_SCOPED_LOGGER_TAG(log, "Channel", options_.prefix);
      BOOST_LOG_SCOPED_LOGGER_TAG(log, "Time", t);
      BOOST_LOG_SCOPED_LOGGER_TAG(log, "Level", level);
      BOOST_LOG(log) << fmt::format("Pressure valve has closed because inner pressure is {}.", inner_pressure);
    }
    reflective_boundary_.FillBoundary(mf, gridding, level);
  }
}

} // namespace fub::amrex

#endif