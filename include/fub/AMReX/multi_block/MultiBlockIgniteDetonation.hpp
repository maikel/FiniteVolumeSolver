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

#ifndef FUB_AMREX_MULTI_BLOCK_IGNITE_DETONATION_HPP
#define FUB_AMREX_MULTI_BLOCK_IGNITE_DETONATION_HPP

#include "fub/AMReX/IgniteDetonation.hpp"
#include "fub/AMReX/multi_block/MultiBlockIntegratorContext2.hpp"

#include <boost/serialization/access.hpp>
#include <boost/serialization/vector.hpp>

#include <vector>

namespace fub::amrex {

class MultiBlockIgniteDetonation {
public:
  static constexpr int Rank = 1;

  MultiBlockIgniteDetonation(const IdealGasMix<1>& equation,
                             std::size_t n_tubes, int max_refinement_level,
                             const IgniteDetonationOptions& opts = {});

  MultiBlockIgniteDetonation(const IdealGasMix<1>& equation,
                             std::size_t n_tubes, int max_refinement_level,
                             const std::vector<IgniteDetonationOptions>& opts);

  [[nodiscard]] static Duration ComputeStableDt(const MultiBlockIntegratorContext2&, int level) noexcept;

  [[nodiscard]] Result<void, TimeStepTooLarge>
  AdvanceLevel(MultiBlockIntegratorContext2& context, int level, Duration dt,
               const ::amrex::IntVect& ngrow = ::amrex::IntVect(0));

private:
  std::vector<IgniteDetonation> source_terms_;
  int max_number_levels_{1};

  friend class boost::serialization::access;
  template <class Archive>
  void serialize(Archive& ar, unsigned int /* version */) {
    for (auto&& source_term : source_terms_) {
      ar & source_term;
    }
  }
};

} // namespace fub::amrex

#endif
