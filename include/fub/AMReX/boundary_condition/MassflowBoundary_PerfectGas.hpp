// Copyright (c) 2021 Christian Zenker
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

#ifndef FUB_AMREX_MASSFLOW_BOUNDARY_PERFECTGAS_HPP
#define FUB_AMREX_MASSFLOW_BOUNDARY_PERFECTGAS_HPP

#include "fub/AMReX/GriddingAlgorithm.hpp"

#include "fub/Direction.hpp"
#include "fub/equations/PerfectGas.hpp"

#include <AMReX.H>

namespace fub::amrex {
/// \ingroup BoundaryCondition
///
struct MassflowBoundary_PerfectGasOptions {
  MassflowBoundary_PerfectGasOptions() = default;
  MassflowBoundary_PerfectGasOptions(const ProgramOptions& options);

  void Print(SeverityLogger& log) const;

  std::string channel_name{"MassflowBoundary_PerfectGas"};
  ::amrex::Box coarse_inner_box{};
  double required_massflow = 0.0;
  double surface_area = 1.0;
  Direction dir = Direction::X;
  int side = 0;
};

/// \ingroup BoundaryCondition
///
/// \brief This boundary models an inflow boundary with constant mean mass flow.
template <int Dim> class MassflowBoundary_PerfectGas {
public:
  MassflowBoundary_PerfectGas(
      const PerfectGas<Dim>& eq,
      const MassflowBoundary_PerfectGasOptions& options);

  void FillBoundary(::amrex::MultiFab& mf, const GriddingAlgorithm& grid,
                    int level);

  void FillBoundary(::amrex::MultiFab& mf, const GriddingAlgorithm& grid,
                    int level, Direction dir);

  void FillBoundary(::amrex::MultiFab& mf, const ::amrex::Geometry& geom,
                    const Complete<PerfectGas<Dim>>& state);

  void ComputeBoundaryState(Complete<PerfectGas<Dim>>& boundary,
                            const Complete<PerfectGas<Dim>>& inner);

private:
  PerfectGas<Dim> equation_;
  MassflowBoundary_PerfectGasOptions options_;
};

extern template class MassflowBoundary_PerfectGas<1>;
extern template class MassflowBoundary_PerfectGas<2>;
extern template class MassflowBoundary_PerfectGas<3>;

} // namespace fub::amrex

#endif // FUB_AMREX_MASSFLOW_BOUNDARY_PERFECTGAS_HPP
