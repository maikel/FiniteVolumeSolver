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

#ifndef FUB_AMREX_BOUNDARY_CONDITION_ISENTROPIC_HPP
#define FUB_AMREX_BOUNDARY_CONDITION_ISENTROPIC_HPP

#include "fub/AMReX/GriddingAlgorithm.hpp"
#include "fub/equations/IdealGasMix.hpp"

#include "fub/ext/Log.hpp"
#include "fub/ext/ProgramOptions.hpp"

#include <fmt/format.h>

#include <optional>

namespace fub::amrex {

/// \ingroup BoundaryCondition
///
struct IsentropicPressureBoundaryOptions {
  IsentropicPressureBoundaryOptions() = default;
  IsentropicPressureBoundaryOptions(const ProgramOptions& options);

  template <typename Log> void Print(Log& log) {
    BOOST_LOG(log) << fmt::format(" - channel_name = {}", channel_name);
    BOOST_LOG(log) << fmt::format(" - outer_pressure = {} [Pa]",
                                  outer_pressure);
    BOOST_LOG(log) << fmt::format(" - efficiency = {} [-]", efficiency);
    if (boundary_section) {
      std::array<int, AMREX_SPACEDIM> lower, upper;
      std::copy_n(boundary_section->smallEnd().getVect(), AMREX_SPACEDIM,
                  lower.data());
      std::copy_n(boundary_section->bigEnd().getVect(), AMREX_SPACEDIM,
                  upper.data());
      BOOST_LOG(log) << fmt::format(
          " - boundary_section = {{{{{}}}, {{{}}}}} [-]", lower, upper);
    } else {
      BOOST_LOG(log) << " - boundary_section = {} [-]";
    }
    BOOST_LOG(log) << fmt::format(" - direction = {} [-]", int(direction));
    BOOST_LOG(log) << fmt::format(" - side = {} [-]", side);
  }

  std::string channel_name{"IsentropicPressureBoundary"};
  double outer_pressure = 101325.0;
  double efficiency = 1.0;
  std::optional<::amrex::Box> boundary_section{};
  Direction direction = Direction::X;
  int side = 0;
};

/// \ingroup BoundaryCondition
///
/// \brief This boundary models an isentropic pressure expansion for the
/// one-dimensional ideal gas equations for mixtures.
class IsentropicPressureBoundary {
public:
  IsentropicPressureBoundary(const IdealGasMix<1>& eq,
                             const IsentropicPressureBoundaryOptions& options);

  IsentropicPressureBoundary(const IdealGasMix<1>& eq, double outer_pressure,
                             Direction dir, int side);

  void FillBoundary(::amrex::MultiFab& mf, const GriddingAlgorithm& gridding,
                    int level);

  void FillBoundary(::amrex::MultiFab& mf, const GriddingAlgorithm& gridding,
                    int level, Direction dir);

  void FillBoundary(::amrex::MultiFab& mf, const ::amrex::Geometry& geom);

private:
  IdealGasMix<1> equation_;
  double outer_pressure_;
  Direction dir_;
  int side_;
};

} // namespace fub::amrex

#endif
