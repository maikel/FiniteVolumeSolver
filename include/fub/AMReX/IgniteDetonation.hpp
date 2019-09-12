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

#ifndef FUB_AMREX_IGNITE_DETONATION_HPP
#define FUB_AMREX_IGNITE_DETONATION_HPP

#include "fub/AMReX/GriddingAlgorithm.hpp"
#include "fub/TimeStepError.hpp"
#include "fub/equations/IdealGasMix.hpp"
#include "fub/ext/outcome.hpp"

#include <boost/program_options.hpp>

namespace fub::amrex {

struct IgniteDetonationOptions {
  IgniteDetonationOptions() = default;

  explicit IgniteDetonationOptions(const boost::program_options::variables_map& vm);

  static boost::program_options::options_description
  GetCommandLineOptions();

  double measurement_position{1.0};
  double equivalence_ratio_criterium{0.95};
  double temperature_low{300.0};
  double temperature_high{2000.0};
  double ramp_width{0.05};
  double ignite_position{0.0};
  Duration ignite_interval{0.0};
};

class IgniteDetonation {
public:
  static constexpr int Rank = 1;

  IgniteDetonation(const IdealGasMix<1>& eq,
                   std::shared_ptr<GriddingAlgorithm> grid,
                   const IgniteDetonationOptions& opts = {});

  /////////////////////////////////////////////////////////////////////////
  // member functions needed for being a source term

  void
  ResetHierarchyConfiguration(std::shared_ptr<amrex::GriddingAlgorithm> grid) {
    gridding_ = std::move(grid);
  }

  [[nodiscard]] Duration ComputeStableDt() const noexcept;

  [[nodiscard]] Result<void, TimeStepTooLarge> AdvanceLevel(int level,
                                                            Duration dt);

private:
  IdealGasMix<1> equation_;
  std::shared_ptr<GriddingAlgorithm> gridding_;
  IgniteDetonationOptions options_;
  Duration last_ignition_{-std::numeric_limits<double>::infinity()};
};

} // namespace fub::amrex

#endif
