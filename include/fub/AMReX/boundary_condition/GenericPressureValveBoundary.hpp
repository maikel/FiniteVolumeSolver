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

#include "fub/AMReX/MultiFabUtilities.hpp"

namespace fub::amrex {

struct GenericPressureValveBoundaryOptions {
  double forward_efficiency{1.0};
  double backward_efficiency{0.0};
  double open_below_pressure{1.0};
  Direction dir{Direction::X};
  int side{0};
};

template <typename EulerEquation> class GenericPressureValveBoundary {
public:
  using InflowFunction =
      std::function<KineticState<EulerEquation>(IntegratorContext&, Duration)>;

  enum class State { open, closed };

  GenericPressureValveBoundary(
      const EulerEquation& equation, InflowFunction fn,
      const GenericPressureValveBoundaryOptions& options);

  void FillBoundary(::amrex::MultiFab& mf, const GriddingAlgorithm& gridding,
                    int level);

  void FillBoundary(::amrex::MultiFab& mf, const GriddingAlgorithm& gridding,
                    int level, Direction dir);

  State GetState() const noexcept { return state_; }
  std::optional<Duration> GetTimePointWhenOpened() const noexcept {
    return t_opened_;
  }

private:
  EulerEquation equation_;
  InflowFunction inflow_funtion_;
  GenericPressureValveBoundaryOptions options_;

  IndexMapping<EulerEquation> comps_{equation_};
  State state_{State::closed};
  Duration t_opened_{};
};

template <typename EulerEquation>
void GenericPressureValveBoundary<EulerEquation>::FillBoundary(
    ::amrex::MultiFab& mf, const GriddingAlgorithm& gridding, int level) {
  ::amrex::MultiFab& cell_data =
      gridding.GetPatchHierarchy().GetPatchLevel(level).data;
  const ::amrex::Geometry& geom =
      gridding.GetPatchHierarchy().GetGeometry(level);
  const ::amrex::Box domain_box = geom.Domain();
  const ::amrex::Box inner_box =
      GetInnerBox(domain_box, options_.side, options_.dir, 1);
  const double inner_pressure =
      GetMeanValueInBox(cell_data, inner_box, comps_.pressure);
  if (inner_pressure <= options_.open_below_pressure) {

  }
  
}

} // namespace fub::amrex

#endif