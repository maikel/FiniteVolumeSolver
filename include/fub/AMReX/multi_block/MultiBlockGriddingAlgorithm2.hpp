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

#ifndef FUB_AMREX_COUPLED_GRIDDING_ALGORITHM2_HPP
#define FUB_AMREX_COUPLED_GRIDDING_ALGORITHM2_HPP

#include "fub/AMReX/boundary_condition/BoundarySet.hpp"
#include "fub/AMReX/cutcell/boundary_condition/BoundarySet.hpp"
#include "fub/AMReX/multi_block/MultiBlockBoundary2.hpp"

#include "fub/AMReX/GriddingAlgorithm.hpp"
#include "fub/AMReX/cutcell/GriddingAlgorithm.hpp"

#include "fub/output/CounterOutput.hpp"

#include "fub/core/span.hpp"
#include <range/v3/view/enumerate.hpp>
#include <range/v3/view/zip.hpp>

namespace fub {
namespace amrex {

/// \ingroup GriddingAlgorithm
class MultiBlockGriddingAlgorithm2 {
public:
  template <typename TubeEquation, typename PlenumEquation>
  MultiBlockGriddingAlgorithm2(
      const TubeEquation& tube_equation, const PlenumEquation& plenum_equation,
      std::vector<std::shared_ptr<GriddingAlgorithm>> tubes,
      std::vector<std::shared_ptr<cutcell::GriddingAlgorithm>> plena,
      std::vector<BlockConnection> connectivity);

  MultiBlockGriddingAlgorithm2(const MultiBlockGriddingAlgorithm2& other);
  MultiBlockGriddingAlgorithm2&
  operator=(const MultiBlockGriddingAlgorithm2& other);

  MultiBlockGriddingAlgorithm2(MultiBlockGriddingAlgorithm2&& other) noexcept =
      default;

  MultiBlockGriddingAlgorithm2&
  operator=(MultiBlockGriddingAlgorithm2&& other) noexcept = default;

  void PreAdvanceHierarchy();

  std::ptrdiff_t GetCycles() const noexcept { return plena_[0]->GetCycles(); }

  Duration GetTimePoint() const noexcept { return plena_[0]->GetTimePoint(); }

  [[nodiscard]] span<const std::shared_ptr<GriddingAlgorithm>> GetTubes() const
      noexcept;
  [[nodiscard]] span<const std::shared_ptr<cutcell::GriddingAlgorithm>>
  GetPlena() const noexcept;

  [[nodiscard]] span<const BlockConnection> GetConnectivity() const noexcept;
  [[nodiscard]] span<AnyMultiBlockBoundary>
  GetBoundaries(int level = 0) noexcept {
    return boundaries_[static_cast<std::size_t>(level)];
  }

  void RegridAllFinerLevels(int which_level);

private:
  std::vector<std::shared_ptr<GriddingAlgorithm>> tubes_;
  std::vector<std::shared_ptr<cutcell::GriddingAlgorithm>> plena_;
  std::vector<BlockConnection> connectivity_;
  std::vector<std::vector<AnyMultiBlockBoundary>> boundaries_;
  std::vector<AnyBoundaryCondition> tube_bcs_{};
  std::vector<cutcell::AnyBoundaryCondition> plenum_bcs_{};
};

template <typename TubeEquation, typename PlenumEquation>
MultiBlockGriddingAlgorithm2::MultiBlockGriddingAlgorithm2(
    const TubeEquation& tube_equation, const PlenumEquation& plenum_equation,
    std::vector<std::shared_ptr<GriddingAlgorithm>> tubes,
    std::vector<std::shared_ptr<cutcell::GriddingAlgorithm>> plena,
    std::vector<BlockConnection> connectivity)
    : tubes_{std::move(tubes)}, plena_{std::move(plena)},
      connectivity_(std::move(connectivity)) {
  const int nlevel = plena_[0]->GetPatchHierarchy().GetMaxNumberOfLevels();
  boundaries_.resize(static_cast<std::size_t>(nlevel));
  for (auto&& [level, boundaries] : ranges::view::enumerate(boundaries_)) {
    for (const BlockConnection& conn : connectivity_) {
      boundaries.emplace_back(
          MultiBlockBoundary2(tube_equation, plenum_equation), *this, conn,
          conn.ghost_cell_width, level);
    }
  }
  // Add multi block boundaries to the tubes and plena
  for (auto&& [tube_id, tube] : ranges::view::enumerate(tubes_)) {
    amrex::BoundarySet boundary;
    boundary.conditions.push_back(tube->GetBoundaryCondition());
    // Add all multi block boundaries to boundary.conditions
    for (auto& boundaries : boundaries_) {
      for (auto&& [conn, bc] : ranges::view::zip(connectivity_, boundaries)) {
        if (conn.tube.id == tube_id) {
          boundary.conditions.emplace_back(bc);
        }
      }
    }
    tube->GetBoundaryCondition() = boundary;
  }
  for (auto&& [plenum_id, plenum] : ranges::view::enumerate(plena_)) {
    amrex::cutcell::BoundarySet boundary;
    boundary.conditions.push_back(plenum->GetBoundaryCondition());
    // Add all multi block boundaries to boundary.conditions
    for (auto& boundaries : boundaries_) {
      for (auto&& [conn, bc] : ranges::view::zip(connectivity_, boundaries)) {
        if (conn.plenum.id == plenum_id) {
          boundary.conditions.emplace_back(bc);
        }
      }
    }
    plenum->GetBoundaryCondition() = boundary;
  }
}

} // namespace amrex

template <typename PrintDuration>
class CounterOutput<amrex::MultiBlockGriddingAlgorithm2, PrintDuration>
    : public OutputAtFrequencyOrInterval<amrex::MultiBlockGriddingAlgorithm2> {
public:
  CounterOutput(const ProgramOptions& po,
                std::chrono::steady_clock::time_point ref)
      : OutputAtFrequencyOrInterval<amrex::MultiBlockGriddingAlgorithm2>(po),
        reference_{ref} {}

  CounterOutput(const ProgramOptions& po)
      : OutputAtFrequencyOrInterval<amrex::MultiBlockGriddingAlgorithm2>(po),
        reference_{std::chrono::steady_clock::now()} {}

  CounterOutput(std::chrono::steady_clock::time_point reference,
                std::vector<std::ptrdiff_t> frequencies,
                std::vector<Duration> intervals)
      : OutputAtFrequencyOrInterval<amrex::MultiBlockGriddingAlgorithm2>(
            std::move(frequencies), std::move(intervals)),
        reference_(reference) {}

  void operator()(const amrex::MultiBlockGriddingAlgorithm2& grid) override {
    std::chrono::steady_clock::time_point now =
        std::chrono::steady_clock::now();
    auto diff =
        std::chrono::duration_cast<std::chrono::nanoseconds>(now - reference_);
    std::vector<CounterResult> statistics = grid.GetPlena()[0]
                                                ->GetPatchHierarchy()
                                                .GetCounterRegistry()
                                                ->gather_statistics();
    if (statistics.size()) {
      SeverityLogger log = GetInfoLogger();
      BOOST_LOG_SCOPED_LOGGER_TAG(log, "Channel", "CounterOutput");
      BOOST_LOG_SCOPED_LOGGER_TAG(log, "Time", grid.GetTimePoint().count());
      BOOST_LOG(log) << print_statistics<PrintDuration>(statistics, diff.count());
    }
  }

private:
  std::chrono::steady_clock::time_point reference_;
};

} // namespace fub

#endif // FUB_AMREX_COUPLED_GRIDDING_ALGORITHM_HPP
