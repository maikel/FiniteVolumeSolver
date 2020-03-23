// Copyright (c) 2020 Stefan Vater
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

#ifndef FUB_AMREX_VISCOSITY_SOURCE_TERM_HPP
#define FUB_AMREX_VISCOSITY_SOURCE_TERM_HPP

#include "fub/AMReX/IntegratorContext.hpp"
#include "fub/TimeStepError.hpp"
#include "fub/equations/PerfectGas.hpp"
#include "fub/ext/outcome.hpp"
#include <AMReX_MLTensorOp.H>

namespace fub::amrex {

template <int R> class ViscositySourceTerm {
public:
  static constexpr int Rank = R;
  static constexpr std::size_t sRank = static_cast<std::size_t>(Rank);

  ViscositySourceTerm(const PerfectGas<Rank>& eq,
                  const std::shared_ptr<GriddingAlgorithm>& grid);

//   ViscositySourceTerm(const ViscositySourceTerm& other);
//   ViscositySourceTerm& operator=(const ViscositySourceTerm& other);
//
//   ViscositySourceTerm(ViscositySourceTerm&& other) noexcept = default;
//   ViscositySourceTerm& operator=(ViscositySourceTerm&& other) noexcept = default;

  /////////////////////////////////////////////////////////////////////////
  // member functions needed for being a source term

  void
  ResetHierarchyConfiguration(const std::shared_ptr<GriddingAlgorithm>& grid);

  Duration ComputeStableDt(int level);

  Result<void, TimeStepTooLarge>
  AdvanceLevel(IntegratorContext& simulation_data, int level, Duration dt);

  //////////////////////////////////////////////////////////////////////////
  // additional member functions to get class data

private:
  PerfectGas<Rank> equation_;

  std::unique_ptr<::amrex::MLTensorOp> m_reg_solve_op;
  std::unique_ptr<::amrex::MLTensorOp> m_reg_apply_op;

  // DiffusionOp verbosity
  int m_verbose = 0;

  // Options to control MLMG behavior
  int m_mg_verbose = 0;
  int m_mg_cg_verbose = 0;
  int m_mg_max_iter = 100;
  int m_mg_cg_maxiter = 100;
  int m_mg_max_fmg_iter = 0;
  int m_mg_max_coarsening_level = 100;
  int m_mg_maxorder = 2;
  ::amrex::Real m_mg_rtol = 1.0e-11;
  ::amrex::Real m_mg_atol = 1.0e-14;
  std::string m_bottom_solver_type = "bicgstab";

};

extern template class ViscositySourceTerm<1>;
extern template class ViscositySourceTerm<2>;
extern template class ViscositySourceTerm<3>;

} // namespace fub::amrex

#endif
