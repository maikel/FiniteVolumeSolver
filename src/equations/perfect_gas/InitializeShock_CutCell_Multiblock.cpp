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

#include "fub/equations/perfect_gas/InitializeShock_CutCell_Multiblock.hpp"

namespace fub::amrex::cutcell::feedback_functions {

ShockOptions::ShockOptions(const ProgramOptions& options) {
  FUB_GET_OPTION_VAR(options, shock_mach_number);
  // FUB_GET_OPTION_VAR(options, shock_x_location);
  FUB_GET_OPTION_VAR(options, average_post_shock_box);
  FUB_GET_OPTION_VAR(options, shock_time);
  // FUB_GET_OPTION_VAR(options, dir);
}

void ShockOptions::Print(SeverityLogger& log) {
  FUB_PRINT_OPTION_VAR(log, shock_mach_number);
  // FUB_PRINT_OPTION_VAR(log, shock_x_location);
  FUB_PRINT_OPTION_VAR(log, shock_time);
  // FUB_PRINT_OPTION_VAR(log, dir);
  std::array<int, AMREX_SPACEDIM> lower, upper;
  std::copy_n(average_post_shock_box.smallEnd().getVect(), AMREX_SPACEDIM,
              lower.data());
  std::copy_n(average_post_shock_box.bigEnd().getVect(), AMREX_SPACEDIM,
              upper.data());
  BOOST_LOG(log) << fmt::format(" - average_post_shock_box = {{{{{}}}, {{{}}}}} [-]",
                                lower, upper);
}

} // namespace fub::amrex::cutcell::feedback_functions
