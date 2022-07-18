// Copyright (c) 2021 Maikel Nadolski
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

#ifndef FUB_AMREX_CUTCELL_HGRID_FLUX_METHOD_HPP
#define FUB_AMREX_CUTCELL_HGRID_FLUX_METHOD_HPP

#include "fub/Duration.hpp"
#include "fub/Direction.hpp"

#include "fub/AMReX/cutcell/IntegratorContext.hpp"


namespace fub::amrex::cutcell {

template <typename Equation> class HGridFluxMethod {
public:
  explicit HGridFluxMethod(const Equation& eq);

  // This function performs a half time step with a non-conservative first order
  // accurate dimensionally split h-grid method. Then it stores the resulting
  // tangential mass flows on the embedded boundaries as reference values.
  void PreAdvanceLevel(IntegratorContext& simulation_data, int level,
                       Duration dt);

  // This function uses the one-dimensional second order accurate h-grid method.
  // On the embedded boundary it enforces a mass flows that was computed in a
  // prior call to PreAdvanceLevel.
  void ComputeNumericFluxes(IntegratorContext& simulation_data, int level,
                            Duration dt, Direction dir);

  // Returns the a stable time step size as if for regular grids.
  Duration ComputeStableDt(const IntegratorContext& simulation_data, int level,
                           Direction dir);

  const Equation& GetEquation() const noexcept;

private:
  Equation equation_;
  std::vector<::amrex::MultiFab> reference_mass_flow_per_level_;
};

extern template <> HGridFluxMethod<fub::PerfectGas<2>>;

} // namespace fub::amrex::cutcell

#endif