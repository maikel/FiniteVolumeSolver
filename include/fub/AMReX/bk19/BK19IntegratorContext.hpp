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

#ifndef FUB_AMREX_BK19_INTEGRATOR_CONTEXT_HPP
#define FUB_AMREX_BK19_INTEGRATOR_CONTEXT_HPP

#include "fub/AMReX/IntegratorContext.hpp"

namespace fub::amrex {

struct BK19AdvectiveFluxes {
  ::amrex::MultiFab on_cells;
  std::array<::amrex::MultiFab, AMREX_SPACEDIM> on_faces;
};

class BK19IntegratorContext : public IntegratorContext {
public:
  BK19IntegratorContext(std::shared_ptr<GriddingAlgorithm> gridding,
                    HyperbolicMethod method);

  BK19IntegratorContext(std::shared_ptr<GriddingAlgorithm> gridding,
                    HyperbolicMethod method, int cell_gcw, int face_gcw);

  /// \brief Deeply copies a context and all its distributed data for all MPI
  /// ranks.
  BK19IntegratorContext(const BK19IntegratorContext&);
  BK19IntegratorContext operator=(const BK19IntegratorContext&);

  BK19IntegratorContext(BK19IntegratorContext&&) = default;
  BK19IntegratorContext& operator=(BK19IntegratorContext&&) = default;

  BK19AdvectiveFluxes& GetAdvectiveFluxes(int level);
  const BK19AdvectiveFluxes& GetAdvectiveFluxes(int level) const;

  /// \brief Replaces the underlying gridding algorithm with the specified one.
  void ResetHierarchyConfiguration(
      std::shared_ptr<GriddingAlgorithm> gridding) override;

  /// \brief Whenever the gridding algorithm changes the data hierarchy this
  /// function will regrid all distributed helper variables managed by the
  /// context.
  ///
  /// \param[in] level  The level number of the coarsest level which changed its
  /// shape. Regrid all levels finer than level.
  void ResetHierarchyConfiguration(int level = 0) override;

private:
  std::vector<BK19AdvectiveFluxes> Pv_;
};

} // namespace fub::amrex

#endif
