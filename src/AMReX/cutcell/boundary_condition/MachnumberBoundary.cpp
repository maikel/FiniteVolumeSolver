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

#include "fub/AMReX/cutcell/boundary_condition/MachnumberBoundary.hpp"

namespace fub::amrex::cutcell {
MachnumberBoundaryOptions::MachnumberBoundaryOptions(
    const ProgramOptions& options) {
  channel_name = GetOptionOr(options, "channel_name", channel_name);
  boundary_section = GetOptionOr(options, "boundary_section", boundary_section);
  required_mach_number =
      GetOptionOr(options, "required_mach_number", required_mach_number);
  dir = GetOptionOr(options, "dir", dir);
  side = GetOptionOr(options, "side", side);
}

void MachnumberBoundaryOptions::Print(SeverityLogger& log) const {
  BOOST_LOG(log) << fmt::format(" - channel_name = '{}'", channel_name);
  std::array<int, AMREX_SPACEDIM> lower, upper;
  std::copy_n(boundary_section.smallEnd().getVect(), AMREX_SPACEDIM,
              lower.data());
  std::copy_n(boundary_section.bigEnd().getVect(), AMREX_SPACEDIM,
              upper.data());
  BOOST_LOG(log) << fmt::format(" - boundary_section = {{{{{}}}, {{{}}}}} [-]",
                                lower, upper);
  BOOST_LOG(log) << fmt::format(" - required_mach_number = {} [-]",
                                required_mach_number);
  BOOST_LOG(log) << fmt::format(" - dir = {} [-]", static_cast<int>(dir));
  BOOST_LOG(log) << fmt::format(" - side = {} [-]", side);
}
} // namespace fub::amrex::cutcell