// Copyright (c) 2020 Maikel Nadolski
// Copyright (c) 2019 Stefan Vater
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

#include "fub/AMReX/solver/BK19Solver.hpp"
#include "fub/ext/Vc.hpp"

namespace fub::amrex::bk19 {

std::vector<::amrex::MultiFab>
CopyScratchFromContext(const IntegratorContext& context) {
  const int nlevel = context.GetPatchHierarchy().GetNumberOfLevels();
  std::vector<::amrex::MultiFab> scratch_vector;
  scratch_vector.reserve(nlevel);
  for (int level = 0; level < nlevel; ++level) {
    const ::amrex::MultiFab& scratch = context.GetScratch(level);
    ::amrex::MultiFab& scratch_aux = scratch_vector.emplace_back(
        scratch.boxArray(), scratch.DistributionMap(), scratch.nComp(),
        no_ghosts);
    scratch_aux.back().copy(scratch);
  }
  return scratch_vector;
}

void ApplyExplicitCoriolisSourceTerm(std::array<span<double>, 2> momentum,
                                     double factor1, double factor2,
                                     span<const double, 3> k) {
  Vc::Vector<double> fac1(factor1);
  Vc::Vector<double> fac2(factor2);
  Vc::Vector<double> k2(k[2]);
  double* out_x = momentum[0].begin();
  double* out_y = momentum[1].begin();
  const double* end = momentum[0].end();
  std::ptrdiff_t n = end - out_x;
  while (n >= size) {
    Vc::Vector<double> rhou(out_x, Vc::Unaligned);
    Vc::Vector<double> rhov(out_y, Vc::Unaligned);
    Vc::Vector rhou_next = fac2 * (rhou - fac1 * k2 * rhov);
    Vc::Vector rhov_next = fac2 * (rhov + fac1 * k2 * rhou);
    rhou_next.store(out_x, Vc::Unaligned);
    rhov_next.store(out_y, Vc::Unaligned);
    AdvanceBy(size, out_x, out_y);
    n = end - out_x;
  }
  const auto mask = Vc::Vector<double>([](int i) { return i; }) < int(n);
  const Vc::Vector<double> rhou = mask_load(out_x, mask);
  const Vc::Vector<double> rhov = mask_load(out_y, mask);
  Vc::Vector rhou_next = fac2 * (rhou - fac1 * k2 * rhov);
  Vc::Vector rhov_next = fac2 * (rhov + fac1 * k2 * rhou);
  rhou_next.store(out_x, mask, Vc::Unaligned);
  rhov_next.store(out_y, mask, Vc::Unaligned);
}

void ApplyExplicitCoriolisSourceTerm(std::array<span<double>, 3> rhou,
                                     double fac1, double fac2,
                                     span<const double, 3> k) {
  Vc::Vector<double> fac1(factor1);
  Vc::Vector<double> fac2(factor2);
  Vc::Vector<double> k0(k[0]);
  Vc::Vector<double> k1(k[1]);
  Vc::Vector<double> k2(k[2]);
  double* out_x = momentum[0].begin();
  double* out_y = momentum[1].begin();
  double* out_z = momentum[2].begin();
  const double* end = momentum[0].end();
  std::ptrdiff_t n = end - out_x;
  while (n >= size) {
    Vc::Vector<double> rhou(out_x, Vc::Unaligned);
    Vc::Vector<double> rhov(out_y, Vc::Unaligned);
    Vc::Vector<double> rhow(out_y, Vc::Unaligned);
    Vc::Vector rhou_next = fac2 * (rhou + fac1 * (k1 * rhow - k2 * rhov));
    Vc::Vector rhov_next = fac2 * (rhov + fac1 * (k2 * rhou - k0 * rhow));
    Vc::Vector rhow_next = fac2 * (rhow + fac1 * (k0 * rhov - k1 * rhou));
    rhou_next.store(out_x, Vc::Unaligned);
    rhov_next.store(out_y, Vc::Unaligned);
    AdvanceBy(size, out_x, out_y, out_z);
    n = end - out_x;
  }
  const auto mask = Vc::Vector<double>([](int i) { return i; }) < int(n);
  const Vc::Vector<double> rhou = mask_load(out_x, mask);
  const Vc::Vector<double> rhov = mask_load(out_y, mask);
  const Vc::Vector<double> rhow = mask_load(out_z, mask);
  Vc::Vector rhou_next = fac2 * (rhou + fac1 * (k1 * rhow - k2 * rhov));
  Vc::Vector rhov_next = fac2 * (rhov + fac1 * (k2 * rhou - k0 * rhow));
  Vc::Vector rhow_next = fac2 * (rhow + fac1 * (k0 * rhov - k1 * rhou));
  rhou_next.store(out_x, mask, Vc::Unaligned);
  rhov_next.store(out_y, mask, Vc::Unaligned);
  rhow_next.store(out_z, mask, Vc::Unaligned);
}

} // namespace fub::amrex::bk19