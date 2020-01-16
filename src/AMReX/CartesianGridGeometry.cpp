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

#include "fub/AMReX/CartesianGridGeometry.hpp"

namespace fub::amrex {

CartesianGridGeometry::CartesianGridGeometry(const ProgramOptions& options) {
  cell_dimensions = GetOptionOr(options, "cell_dimensions", cell_dimensions);
  std::array<double, AMREX_SPACEDIM> xlo{};
  xlo = GetOptionOr(options, "x_lower", xlo);
  std::array<double, AMREX_SPACEDIM> xup{AMREX_D_DECL(1.0, 1.0, 1.0)};
  xup = GetOptionOr(options, "x_upper", xup);
  coordinates = ::amrex::RealBox(xlo, xup);
  periodicity = GetOptionOr(options, "periodicity", periodicity);
}

} // namespace fub::amrex
