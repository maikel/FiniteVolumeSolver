// Copyright (c) 2021 Maikel Nadolski
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

#ifndef FUB_AMREX_BOUNDARY_CONDITION_GENERIC_PRESSURE_VALVE_HPP
#define FUB_AMREX_BOUNDARY_CONDITION_GENERIC_PRESSURE_VALVE_HPP

#include "fub/AMReX/GriddingAlgorithm.hpp"
#include "fub/AMReX/MultiFabUtilities.hpp"
#include "fub/AMReX/boundary_condition/ConstantBoundary.hpp"
#include "fub/AMReX/boundary_condition/ReflectiveBoundary.hpp"
#include "fub/equations/perfect_gas_mix/PlenaControl.hpp"

#include "fub/ext/Log.hpp"

namespace fub::amrex {

struct GenericPressureValveBoundaryOptions {
  std::string prefix{"GenericPressureValve"};
  Direction dir{Direction::X};
  int side{0};
};

struct ChangeTOpened_ReducedModelDemo {
  template <typename EulerEquation>
  [[nodiscard]] std::optional<Duration>
  operator()(EulerEquation&, std::optional<Duration> t_opened,
             double inner_pressure,
             const perfect_gas_mix::gt::PlenumState& compressor_state,
             const GriddingAlgorithm& gridding, int) const noexcept {
    const double p_ref = compressor_state.pressure;
    if (!t_opened && inner_pressure <= p_ref) {
      t_opened = gridding.GetTimePoint();
    } else if (inner_pressure > 1.1 * p_ref) {
      t_opened.reset();
    }
    return t_opened;
  }
};

struct ChangeTOpened_Klein {
  template <typename EulerEquation>
  [[nodiscard]] std::optional<Duration>
  operator()(EulerEquation&, std::optional<Duration> t_opened,
             double inner_pressure,
             const perfect_gas_mix::gt::PlenumState& compressor_state,
             const GriddingAlgorithm& gridding, int) const noexcept {
    if (!t_opened && inner_pressure <= compressor_state.pressure) {
      t_opened = gridding.GetTimePoint();
    } else if (inner_pressure > compressor_state.pressure) {
      t_opened.reset();
    }
    return t_opened;
  }
};

struct ChangeTOpened_DeflagrationValve {
  template <typename EulerEquation>
  [[nodiscard]] std::optional<fub::Duration>
  operator()(EulerEquation&, std::optional<fub::Duration>, double,
             const fub::perfect_gas_mix::gt::PlenumState&,
             const fub::amrex::GriddingAlgorithm& gridding,
             int) const noexcept {
    return gridding.GetTimePoint();
  }
};

struct ChangeTOpened_SwitchDeflagrationToSECValve {
  template <typename EulerEquation>
  [[nodiscard]] std::optional<Duration>
  operator()(EulerEquation&, std::optional<Duration> t_opened,
             double inner_pressure,
             const perfect_gas_mix::gt::PlenumState& compressor_state,
             const GriddingAlgorithm& gridding, int) const noexcept {
    if (!compressor_state.SEC_Mode) {
      return gridding.GetTimePoint();
    } else {
      if (!t_opened && inner_pressure <= compressor_state.pressure) {
        t_opened = gridding.GetTimePoint();
      } else if (inner_pressure > compressor_state.pressure) {
        t_opened.reset();
      }
      return t_opened;
    }
  }
};

struct IsBlockedIfLargePressure {
  template <typename EulerEquation>
  [[nodiscard]] bool
  operator()(EulerEquation&, std::optional<Duration> /* t_opened */,
             double inner_pressure,
             const perfect_gas_mix::gt::PlenumState& compressor_state,
             const GriddingAlgorithm& /* gridding */,
             int /* level */) const noexcept {
    return inner_pressure > compressor_state.pressure;
  }
};

struct IsNeverBlocked {
  template <typename EulerEquation>
  [[nodiscard]] bool
  operator()(EulerEquation&, std::optional<fub::Duration> /* t_opened */,
             double, const fub::perfect_gas_mix::gt::PlenumState&,
             const fub::amrex::GriddingAlgorithm& /* gridding */,
             int /* level */) const noexcept {
    return false;
  }
};

struct IsNeverBlockedTilTargetPressure {
  template <typename EulerEquation>
  [[nodiscard]] bool
  operator()(EulerEquation&, std::optional<Duration> /* t_opened */,
             double inner_pressure,
             const perfect_gas_mix::gt::PlenumState& compressor_state,
             const GriddingAlgorithm& /* gridding */,
             int /* level */) const noexcept {
    if (!compressor_state.SEC_Mode) {
      return false;
    } else {
      return inner_pressure > compressor_state.pressure;
    }
  }
};

template <typename EulerEquation, typename InflowFunction,
          typename ChangeTOpenedT,
          typename IsBlocked = IsBlockedIfLargePressure>
class GenericPressureValveBoundary {
public:
  GenericPressureValveBoundary(
      const EulerEquation& equation,
      std::shared_ptr<const perfect_gas_mix::gt::ControlState> control_state,
      InflowFunction fn,
      const GenericPressureValveBoundaryOptions& options = {});

  void FillBoundary(::amrex::MultiFab& mf, const GriddingAlgorithm& gridding,
                    int level);

  void FillBoundary(::amrex::MultiFab& mf, const GriddingAlgorithm& gridding,
                    int level, Direction dir);

  std::optional<Duration> GetTimePointWhenOpened() const noexcept {
    return t_opened_;
  }

  const std::shared_ptr<const perfect_gas_mix::gt::ControlState>&
  GetSharedControlState() const noexcept {
    return control_state_;
  }

private:
  void ChangeTOpened(double inner_pressure, const GriddingAlgorithm& gridding,
                     int level);

  EulerEquation equation_;
  GenericPressureValveBoundaryOptions options_;
  ConstantBoundary<EulerEquation> constant_boundary_;
  ReflectiveBoundary<execution::SequentialTag, EulerEquation>
      reflective_boundary_;

  std::shared_ptr<const perfect_gas_mix::gt::ControlState> control_state_;
  IndexMapping<EulerEquation> comps_{equation_};
  std::optional<Duration> t_opened_{};

#if __has_cpp_attribute(no_unique_address)
  [[no_unique_address]] InflowFunction inflow_function_{};
  [[no_unique_address]] ChangeTOpenedT change_t_opened_{};
  [[no_unique_address]] IsBlocked is_blocked_{};
#else
  InflowFunction inflow_function_{};
  ChangeTOpenedT change_t_opened_{};
  IsBlocked is_blocked_{};
#endif
};

template <typename Equation, typename InflowFunction>
struct PressureValveBoundary_ReducedModelDemo
    : public GenericPressureValveBoundary<Equation, InflowFunction,
                                          ChangeTOpened_ReducedModelDemo> {
  using GenericPressureValveBoundary<
      Equation, InflowFunction,
      ChangeTOpened_ReducedModelDemo>::GenericPressureValveBoundary;
};

template <typename Equation, typename InflowFunction>
PressureValveBoundary_ReducedModelDemo(
    const Equation&, std::shared_ptr<const perfect_gas_mix::gt::ControlState>,
    InflowFunction, const GenericPressureValveBoundaryOptions&)
    -> PressureValveBoundary_ReducedModelDemo<Equation, InflowFunction>;

template <typename Equation, typename InflowFunction>
PressureValveBoundary_ReducedModelDemo(
    const Equation&, std::shared_ptr<const perfect_gas_mix::gt::ControlState>,
    InflowFunction)
    -> PressureValveBoundary_ReducedModelDemo<Equation, InflowFunction>;

template <typename Equation, typename InflowFunction>
struct PressureValveBoundary_Klein
    : public GenericPressureValveBoundary<Equation, InflowFunction,
                                          ChangeTOpened_Klein> {
  using GenericPressureValveBoundary<
      Equation, InflowFunction,
      ChangeTOpened_Klein>::GenericPressureValveBoundary;
};

template <typename Equation, typename InflowFunction>
PressureValveBoundary_Klein(
    const Equation&, std::shared_ptr<const perfect_gas_mix::gt::ControlState>,
    InflowFunction, const GenericPressureValveBoundaryOptions&)
    -> PressureValveBoundary_Klein<Equation, InflowFunction>;

template <typename Equation, typename InflowFunction>
PressureValveBoundary_Klein(
    const Equation&, std::shared_ptr<const perfect_gas_mix::gt::ControlState>,
    InflowFunction) -> PressureValveBoundary_Klein<Equation, InflowFunction>;

template <typename Equation, typename InflowFunction>
struct DeflagrationValve
    : public GenericPressureValveBoundary<Equation, InflowFunction,
                                          ChangeTOpened_DeflagrationValve,
                                          IsNeverBlocked> {
  using GenericPressureValveBoundary<
      Equation, InflowFunction, ChangeTOpened_DeflagrationValve,
      IsNeverBlocked>::GenericPressureValveBoundary;
};

template <typename Equation, typename InflowFunction>
DeflagrationValve(const Equation&,
                  std::shared_ptr<const perfect_gas_mix::gt::ControlState>,
                  InflowFunction, const GenericPressureValveBoundaryOptions&)
    -> DeflagrationValve<Equation, InflowFunction>;

template <typename Equation, typename InflowFunction>
DeflagrationValve(const Equation&,
                  std::shared_ptr<const perfect_gas_mix::gt::ControlState>,
                  InflowFunction)
    -> DeflagrationValve<Equation, InflowFunction>;

template <typename Equation, typename InflowFunction>
struct SwitchDeflagrationToSECValve
    : public GenericPressureValveBoundary<
          Equation, InflowFunction, ChangeTOpened_SwitchDeflagrationToSECValve,
          IsNeverBlockedTilTargetPressure> {
  using GenericPressureValveBoundary<
      Equation, InflowFunction, ChangeTOpened_SwitchDeflagrationToSECValve,
      IsNeverBlockedTilTargetPressure>::GenericPressureValveBoundary;
};

template <typename Equation, typename InflowFunction>
SwitchDeflagrationToSECValve(
    const Equation&, std::shared_ptr<const perfect_gas_mix::gt::ControlState>,
    InflowFunction, const GenericPressureValveBoundaryOptions&)
    -> SwitchDeflagrationToSECValve<Equation, InflowFunction>;

template <typename Equation, typename InflowFunction>
SwitchDeflagrationToSECValve(
    const Equation&, std::shared_ptr<const perfect_gas_mix::gt::ControlState>,
    InflowFunction) -> SwitchDeflagrationToSECValve<Equation, InflowFunction>;

template <typename EulerEquation, typename InflowFunction,
          typename ChangeTOpened, typename IsBlocked>
GenericPressureValveBoundary<EulerEquation, InflowFunction, ChangeTOpened,
                             IsBlocked>::
    GenericPressureValveBoundary(
        const EulerEquation& equation,
        std::shared_ptr<const perfect_gas_mix::gt::ControlState> control_state,
        InflowFunction fn, const GenericPressureValveBoundaryOptions& options)
    : equation_(equation),
      options_(options), constant_boundary_{options_.dir, options_.side,
                                            equation_,
                                            Complete<EulerEquation>(equation_)},
      reflective_boundary_(equation_, options_.dir, options_.side),
      control_state_(std::move(control_state)),
      inflow_function_(std::move(fn)) {}

template <typename EulerEquation, typename InflowFunction,
          typename ChangeTOpened, typename IsBlocked>
void GenericPressureValveBoundary<
    EulerEquation, InflowFunction, ChangeTOpened,
    IsBlocked>::FillBoundary(::amrex::MultiFab& mf,
                             const GriddingAlgorithm& gridding, int level,
                             Direction dir) {
  if (dir == options_.dir) {
    FillBoundary(mf, gridding, level);
  }
}
template <typename EulerEquation, typename InflowFunction,
          typename ChangeTOpenedT, typename IsBlocked>
void GenericPressureValveBoundary<
    EulerEquation, InflowFunction, ChangeTOpenedT,
    IsBlocked>::ChangeTOpened(double inner_pressure,
                              const GriddingAlgorithm& gridding, int level) {
  std::optional<Duration> t_opened_new =
      std::invoke(change_t_opened_, equation_, t_opened_, inner_pressure,
                  control_state_->compressor, gridding, level);
  // if (inner_pressure <= options_.open_below_pressure) {
  const Duration t = gridding.GetTimePoint();
  if (!t_opened_ && t_opened_new) {
    SeverityLogger log = GetInfoLogger();
    BOOST_LOG_SCOPED_LOGGER_TAG(log, "Channel", options_.prefix);
    BOOST_LOG_SCOPED_LOGGER_TAG(log, "Time", t.count());
    BOOST_LOG_SCOPED_LOGGER_TAG(log, "Level", level);
    BOOST_LOG(log) << fmt::format(
        "Pressure valve has opened because inner pressure is {}.",
        inner_pressure);
  } else if (t_opened_ && !t_opened_new) {
    SeverityLogger log = GetInfoLogger();
    BOOST_LOG_SCOPED_LOGGER_TAG(log, "Channel", options_.prefix);
    BOOST_LOG_SCOPED_LOGGER_TAG(log, "Time", t.count());
    BOOST_LOG_SCOPED_LOGGER_TAG(log, "Level", level);
    BOOST_LOG(log) << fmt::format(
        "Pressure valve has closed because inner pressure is {}.",
        inner_pressure);
  }
  t_opened_ = t_opened_new;
}

template <typename EulerEquation, typename InflowFunction,
          typename ChangeTOpenedT, typename IsBlocked>
void GenericPressureValveBoundary<
    EulerEquation, InflowFunction, ChangeTOpenedT,
    IsBlocked>::FillBoundary(::amrex::MultiFab& mf,
                             const GriddingAlgorithm& gridding, int level) {
  const ::amrex::Geometry& geom =
      gridding.GetPatchHierarchy().GetGeometry(level);
  const ::amrex::Box domain_box = geom.Domain();
  const ::amrex::Box inner_box =
      GetInnerBox(domain_box, options_.side, options_.dir, 1);
  const double inner_pressure =
      GetMeanValueInBox(mf, inner_box, comps_.pressure);

  ChangeTOpened(inner_pressure, gridding, level);

  bool is_blocked =
      std::invoke(is_blocked_, equation_, t_opened_, inner_pressure,
                  control_state_->compressor, gridding, level);
  if (is_blocked) {
    reflective_boundary_.FillBoundary(mf, gridding, level);
  } else {
    FUB_ASSERT(t_opened_);
    const Duration t = gridding.GetTimePoint();
    const Duration t_diff = t - *t_opened_;
    std::invoke(inflow_function_, equation_, constant_boundary_.state,
                control_state_->compressor, inner_pressure, t_diff, mf,
                gridding, level);
    constant_boundary_.FillBoundary(mf, gridding, level);
  }
}

} // namespace fub::amrex

#endif