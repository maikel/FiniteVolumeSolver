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

#include "fub/AMReX/cutcell/boundary_condition/TurbineMassflowBoundary.hpp"

namespace fub::amrex::cutcell {
TurbineMassflowBoundaryOptions::TurbineMassflowBoundaryOptions(
    const ProgramOptions& options) {
  channel_name = GetOptionOr(options, "channel_name", channel_name);
  boundary_section = GetOptionOr(options, "boundary_section", boundary_section);
  relative_surface_area =
      GetOptionOr(options, "relative_surface_area", relative_surface_area);
  massflow_correlation =
      GetOptionOr(options, "massflow_correlation", massflow_correlation);
  dir = GetOptionOr(options, "dir", dir);
  side = GetOptionOr(options, "side", side);
  int mode_v = static_cast<int>(mode);
  mode_v = GetOptionOr(options, "mode", mode_v);
  if (mode_v != static_cast<int>(TurbineMassflowMode::cellwise) &&
      mode_v != static_cast<int>(TurbineMassflowMode::average_inner_state) &&
      mode_v != static_cast<int>(TurbineMassflowMode::average_outer_state) &&
      mode_v != static_cast<int>(TurbineMassflowMode::average_massflow)) {
    throw std::runtime_error(fmt::format(
        "TurbineMassflowBoundary: Invalid mode (= {}) from options.", mode_v));
  } else {
    mode = static_cast<TurbineMassflowMode>(mode_v);
  }
  coarse_average_mirror_box = GetOptionOr(options, "coarse_average_mirror_box",
                                          coarse_average_mirror_box);
}

void TurbineMassflowBoundaryOptions::Print(SeverityLogger& log) const {
  BOOST_LOG(log) << fmt::format(" - channel_name = '{}'", channel_name);
  std::array<int, AMREX_SPACEDIM> lower, upper;
  std::copy_n(boundary_section.smallEnd().getVect(), AMREX_SPACEDIM,
              lower.data());
  std::copy_n(boundary_section.bigEnd().getVect(), AMREX_SPACEDIM,
              upper.data());
  BOOST_LOG(log) << fmt::format(" - boundary_section = {{{{{}}}, {{{}}}}} [-]",
                                lower, upper);
  BOOST_LOG(log) << fmt::format(" - relative_surface_area = {} [-]",
                                relative_surface_area);
  BOOST_LOG(log) << fmt::format(" - massflow_correlation = {} [-]",
                                massflow_correlation);
  static constexpr std::string_view mode_names[] = {
      "cellwise", "average_inner_state", "average_outer_state",
      "average_massflow"};
  const int mode_v = static_cast<int>(mode);
  FUB_ASSERT(0 <= mode_v && mode_v < 4);
  BOOST_LOG(log) << fmt::format(" - mode = {} ({}) [0, 1, 2, 3]", mode_v,
                                mode_names[mode_v]);
  std::copy_n(coarse_average_mirror_box.smallEnd().getVect(), AMREX_SPACEDIM,
              lower.data());
  std::copy_n(coarse_average_mirror_box.bigEnd().getVect(), AMREX_SPACEDIM,
              upper.data());
  BOOST_LOG(log) << fmt::format(
      " - coarse_average_mirror_box = {{{{{}}}, {{{}}}}} [-]", lower, upper);
  BOOST_LOG(log) << fmt::format(" - dir = {} [-]", static_cast<int>(dir));
  BOOST_LOG(log) << fmt::format(" - side = {} [-]", side);
}
} // namespace fub::amrex::cutcell