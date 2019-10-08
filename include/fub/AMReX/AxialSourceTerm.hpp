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

#ifndef FUB_AMREX_AXIAL_SOURCE_TERM_HPP
#define FUB_AMREX_AXIAL_SOURCE_TERM_HPP

#include "fub/AMReX/GriddingAlgorithm.hpp"
#include "fub/TimeStepError.hpp"
#include "fub/equations/IdealGasMix.hpp"
#include "fub/ext/outcome.hpp"

namespace fub::amrex {

class AxialSourceTerm {
public:
  static constexpr int Rank = 1;

  AxialSourceTerm(const IdealGasMix<1>& eq,
                  std::function<double(double)> diameter,
                  std::shared_ptr<GriddingAlgorithm> gridding);

  AxialSourceTerm(const AxialSourceTerm& other);
  AxialSourceTerm& operator=(const AxialSourceTerm& other);

  AxialSourceTerm(AxialSourceTerm&& other) noexcept = default;
  AxialSourceTerm& operator=(AxialSourceTerm&& other) noexcept = default;

  /////////////////////////////////////////////////////////////////////////
  // member functions needed for being a source term

  void PreAdvanceLevel(int level, Duration dt, int subcycle);
  //  void PostAdvanceLevel(int level, Duration dt, int subcycle);

  void ResetHierarchyConfiguration(
      std::shared_ptr<amrex::GriddingAlgorithm>&& gridding);

  void ResetHierarchyConfiguration(
      const std::shared_ptr<amrex::GriddingAlgorithm>& gridding);

  Duration ComputeStableDt();

  Result<void, TimeStepTooLarge> AdvanceLevel(int level, Duration dt);

  //////////////////////////////////////////////////////////////////////////
  // additional member functions to get class data

  const IdealGasMix<1>& GetEquation() const;

  span<const ::amrex::MultiFab> GetAxialVariations() const;

  Duration GetTimePoint() const;

  std::ptrdiff_t GetCycles() const;

  amrex::PatchHierarchy& GetPatchHierarchy();
  const amrex::PatchHierarchy& GetPatchHierarchy() const;

private:
  std::function<double(double)> diameter_;
  IdealGasMix<1> equation_;
  std::shared_ptr<GriddingAlgorithm> gridding_;
  std::vector<::amrex::MultiFab> Ax_;
};

} // namespace fub::amrex

#endif